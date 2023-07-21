#' Generate taxa-level individual change boxplot pairs
#'
#' This function generates a boxplot pair that shows the change in abundance of taxa at different levels (Phylum, Family, or Genus) between two time points.
#' It allows users to compare the changes in abundance within and between different groups and strata.
#'
#' @name generate_taxa_change_scatterplot_pair
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param taxa.level A character vector specifying the taxa level(s) to include in the analysis. Default is c('Phylum', 'Family', 'Genus').
#' @param prev.filter A numeric value specifying the minimum prevalence threshold for filtering taxa.
#' @param abund.filter A numeric value specifying the minimum abundance threshold for filtering taxa.
#' @param pdf A logical value indicating whether to save the boxplot to a PDF file. Default is TRUE.
#' @param file.ann A character string specifying the file annotation. Default is NULL.
#' @param ... Additional arguments passed to \code{ggplot2::ggsave()}.
#'
#' @return A ggplot object representing the taxa-level individual change boxplot pair.
#'
#' @examples
#' # Load required libraries and data
#' library(microbiome)
#' library(vegan)
#' library(tidyverse)
#' data(peerj32)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#' peerj32.obj$meta.dat <- peerj32.obj$meta.dat %>% select(all_of("subject")) %>% distinct()
#' %>% mutate(cons = runif(n(),0,5)) %>% left_join(peerj32.obj$meta.dat,by = "subject") %>%
#' column_to_rownames("sample")
#' # Generate the boxplot pair
#' plot_list_all <- generate_taxa_change_scatterplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "cons",
#'   strata.var = "sex",
#'   change.base = "1",
#'   change.func = "lfc",
#'   feature.level = "Genus",
#'   feature.dat.type = "other",
#'   features.plot = NULL,
#'   top.k.plot = 8,
#'   top.k.func = "mean",
#'   prev.filter = 0.1,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_taxa_change_scatterplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           change.func = "difference",
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter = 0.1,
           abund.filter = 0.1,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    feature.dat.type <- match.arg(feature.dat.type)

    mStat_validate_data(data.obj)

    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      otu_tab <-
        load_data_obj_count(mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm)
    } else{
      otu_tab <- load_data_obj_count(data.obj)
    }

    tax_tab <- load_data_obj_taxonomy(data.obj) %>%
      as.data.frame() %>%
      {if("original" %in% feature.level) mutate(., original = rownames(.)) else .} %>%
      select(all_of(feature.level))

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
        time.var, group.var, strata.var, subject.var
      )))

    if (is.null(group.var)) {
      group.var = "ALL"
      meta_tab$ALL <- ""
    }

    if (is.null(strata.var)) {
      strata.var = "ALL2"
      meta_tab$ALL2 <- ""
    }

    # Define the colors
    if (is.null(palette)) {
      colors <- c("#92c5de", "#0571b0", "#f4a582", "#ca0020")
    } else {
      colors <- palette
    }

    theme_function <- switch(
      theme.choice,
      prism = ggprism::theme_prism(),
      classic = theme_classic(),
      gray = theme_gray(),
      bw = theme_bw(),
      ggprism::theme_prism()
    ) # 根据用户选择设置主题

    # 使用用户自定义主题（如果提供），否则使用默认主题
    theme_to_use <-
      if (!is.null(custom.theme))
        custom.theme else
          theme_function

    aes_function <- if (!is.null(strata.var)){
      aes(shape = !!sym(strata.var), color = !!sym(strata.var))
    } else {
      aes(color = !!sym(time.var))
    }

    # 将 OTU 表与分类表合并
    otu_tax <-
      cbind(otu_tab, tax_tab)

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list <- lapply(feature.level, function(feature.level) {
      # Filter taxa based on prevalence and abundance
      otu_tax_filtered <- otu_tax %>%
        gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        group_by_at(vars(!!sym(feature.level))) %>%
        summarise(total_count = mean(count),
                  prevalence = sum(count > 0) / n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count, -prevalence) %>%
        left_join(otu_tax, by = feature.level)

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax_filtered %>%
        gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        group_by_at(vars(sample, !!sym(feature.level))) %>%
        summarise(count = sum(count)) %>%
        spread(key = "sample", value = "count")

      compute_function <- function(top.k.func) {
        if (is.function(top.k.func)) {
          results <-
            top.k.func(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
        } else {
          switch(top.k.func,
                 "mean" = {
                   results <-
                     rowMeans(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                              na.rm = TRUE)
                 },
                 "sd" = {
                   results <-
                     rowSds(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                            na.rm = TRUE)
                   names(results) <- rownames(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
                 },
                 stop("Invalid function specified"))
        }

        return(results)
      }

      if (is.null(features.plot) &&
          !is.null(top.k.plot) && !is.null(top.k.func)) {
        features.plot <- names(sort(compute_function(top.k.func), decreasing = TRUE)[1:top.k.plot])
      }

      # 转换计数为数值类型
      otu_tax_agg_numeric <-
        mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      otu_tab_norm <- otu_tax_agg_numeric %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      otu_tab_norm_agg <- otu_tax_agg_numeric %>%
        gather(-!!sym(feature.level), key = "sample", value = "count") %>%
        inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      taxa.levels <-
        otu_tab_norm_agg %>% select(all_of(feature.level)) %>% distinct() %>% pull()

      # 首先，把数据分为两个子集，一个为change.base，一个为change.after
      df_t0 <- otu_tab_norm_agg %>% filter(!!sym(time.var) == change.base)
      df_ts <- otu_tab_norm_agg %>% filter(!!sym(time.var) != change.base)

      # 然后，使用inner_join合并这两个子集，基于Phylum、subject和sex
      df <- inner_join(df_ts, df_t0, by = c(feature.level, subject.var), suffix = c("_ts", "_t0"), relationship = "many-to-many")

      # 最后，计算新的count值
      if (is.function(change.func)) {
        df <- df %>% mutate(new_count = change.func(count_ts, count_t0))
      } else if (change.func == "lfc") {
        # 对于对数折叠变化("lfc")，我们需要插补数据
        # 首先，为每个分类计算非零最小值的一半
        half_nonzero_min_time_2 <- df %>%
          filter(count_ts > 0) %>%
          group_by(!!sym(feature.level)) %>%
          summarize(half_nonzero_min = min(count_ts) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- df %>%
          filter(count_t0 > 0) %>%
          group_by(!!sym(feature.level)) %>%
          summarize(half_nonzero_min = min(count_t0) / 2,
                    .groups = "drop")

        # 然后，用这些值来插补数据
        df <- left_join(df, half_nonzero_min_time_2, by = feature.level, suffix = c("_t0", "_ts"))
        df <- left_join(df, half_nonzero_min_time_1, by = feature.level, suffix = c("_t0", "_ts"))
        df$count_ts[df$count_ts == 0] <- df$half_nonzero_min_ts[df$count_ts == 0]
        df$count_t0[df$count_t0 == 0] <- df$half_nonzero_min_t0[df$count_t0 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message("Imputation was performed using half the minimum nonzero count for each taxa at different time points.")

        df <- df %>% mutate(new_count = log2(count_ts) - log2(count_t0))
      } else if (change.func == "relative difference"){
        df <- df %>%
          mutate(new_count = case_when(
            count_ts == 0 & count_t0 == 0 ~ 0,
            TRUE ~ (count_ts - count_t0) / (count_ts + count_t0)
          ))
      } else {
        df <- df %>% mutate(new_count = count_ts - count_t0)
      }

      df <- df %>% left_join(meta_tab %>% select(-all_of(time.var)) %>% distinct(), by = c(subject.var))

      df <- df %>% setNames(ifelse(names(.) == paste0(time.var,"_ts"), time.var, names(.)))

      # 提前判断change.func的类型，如果是自定义函数则给出特定的标签，否则保持原样
      ylab_label <- if (feature.dat.type != "other") {
        if (is.function(change.func)) {
          paste0("Change in Relative Abundance", " (custom function)")
        } else {
          paste0("Change in Relative Abundance", " (", change.func, ")")
        }
      } else {
        if (is.function(change.func)) {
          paste0("Change in Abundance", " (custom function)")
        } else {
          paste0("Change in Abundance", " (", change.func, ")")
        }
      }

      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 5) {
          features.plot <- taxa.levels[1:4]
        } else {
          features.plot <- taxa.levels
        }
      }

        scatterplot <-
          ggplot(df %>% filter(!!sym(feature.level) %in% features.plot),
                 aes(
                   x = !!sym(group.var),
                   y = new_count,
                   fill = !!sym(strata.var),
                   color = !!sym(strata.var),  # Add this line to separate color by strata.var
                   group = !!sym(strata.var)   # Add this line to separate smoothing line by strata.var
                 )) +
          geom_smooth(se = TRUE, method = 'lm') +
          geom_point(aes_function,data = df %>% filter(!!sym(feature.level) %in% features.plot),
                     size = 4) +
          scale_shape_manual(values = c(21, 22, 24, 25)) +
          scale_fill_manual(values = colors) +
          ylab(ylab_label) +
          #scale_linetype_manual(values = c("solid", "dashed")) + # 设置曲线类型
          scale_color_manual(values = colors, guide = guide_legend(override.aes = list(size = 0))) +
          theme_to_use +
          theme(
            axis.text.x = element_text(
              vjust = 0.5,
              hjust = 1,
              size = base.size
            ),
            axis.text.y = element_text(color = "black", size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.title.x = element_text(size = base.size),
            legend.position = "right",
            legend.direction = "vertical",
            legend.box = "vertical",
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            axis.text = ggplot2::element_text(color = "black",
                                              size = base.size),
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_text(size = 16),
            plot.title = element_text(hjust = 0.5, size = 20)
          ) +
          facet_nested_wrap(as.formula(paste(".~",feature.level)), scales = "fixed")


        if (group.var == "ALL"){
          scatterplot <- scatterplot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
        }

        # Save the stacked dotplot as a PDF file
        if (pdf) {
          pdf_name <- paste0(
            "taxa_indiv_change_scatterplot_v2",
            "_",
            "subject_",
            subject.var,
            "_",
            "time_",
            time.var,
            "_",
            "feature_level_",
            feature.level,
            "_",
            "prev_filter_",
            prev.filter,
            "_",
            "abund_filter_",
            abund.filter
          )
          if (!is.null(group.var)) {
            pdf_name <- paste0(pdf_name, "_", "group_", group.var)
          }
          if (!is.null(strata.var)) {
            pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
          }
          if (!is.null(file.ann)) {
            pdf_name <- paste0(pdf_name, "_", file.ann)
          }
          pdf_name <- paste0(pdf_name, ".pdf")
          pdf(pdf_name, width = pdf.wid, height = pdf.hei)
          # Use lapply to print each ggplot object in the list to a new PDF page
          print(scatterplot)
          # Close the PDF device
          dev.off()
        }

      return(scatterplot)
    })

    return(plot_list)
  }
