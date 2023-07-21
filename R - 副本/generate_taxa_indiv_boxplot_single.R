#' Generate boxplot pair of individual taxa abundance
#'
#' This function creates a boxplot pair of individual taxa abundance at a specified
#' taxonomic level between two or more time points. It takes a data object containing
#' OTU, taxonomy, and metadata tables as input, along with metadata variables, and
#' outputs a ggplot2 boxplot.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param group.var A character string specifying the group variable in the metadata. Default is NULL.
#' @param strata.var A character string specifying the strata variable in the metadata. Default is NULL.
#' @param feature.level A character string specifying the taxonomic level for the boxplot. Default is 'Phylum'.
#' @param prev.filter A numeric value specifying the prevalence threshold for filtering taxa.
#' @param abund.filter A numeric value specifying the abundance threshold for filtering taxa.
#' @param pdf A logical value indicating whether to save the plot as a PDF file. Default is TRUE.
#' @param file.ann An optional character string to be used as a file annotation when saving the PDF. Default is NULL.
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A ggplot2 boxplot of individual taxa abundance at the specified taxonomic level.
#'
#' @examples
#' # Load required libraries and example data
#' library(microbiome)
#' library(tidyverse)
#' data(peerj32)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#'
#' # Generate the boxplot pair
#' generate_taxa_indiv_boxplot_single(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = NULL,
#'   strata.var = NULL,
#'   feature.level = c("Phylum"),
#'   feature.dat.type = "proportion",
#'   Transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.1,
#'   base.size = 20,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = c(ggsci::pal_npg()(9),ggsci::pal_jama()(7),ggsci::pal_lancet()(9)),
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_indiv_boxplot_single(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   Transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_taxa_indiv_boxplot_single <-
  function(data.obj,
           subject.var,
           time.var,
           t.level,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           top.k.plot = NULL,
           top.k.func = NULL,
           Transform = c("identity", "sqrt", "log"),
           prev.filter = 0.05,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "prism",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    feature.dat.type <- match.arg(feature.dat.type)

    mStat_validate_data(data.obj)

    if (!is.character(subject.var))
      stop("`subject.var` should be a character string.")
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    if (!is.null(time.var) &!is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    meta_tab <- load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

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

    aes_function <-  aes(
      x = !!sym(group.var),
      y = value,
      fill = !!sym(group.var)
    )

    # 设置颜色，根据 time.var 的唯一值数量生成颜色列表
    if (is.null(palette)){
      col <-
        c(
          "#E31A1C",
          "#1F78B4",
          "#FB9A99",
          "#33A02C",
          "#FDBF6F",
          "#B2DF8A",
          "#A6CEE3",
          "#BA7A70",
          "#9D4E3F",
          "#829BAB"
        )
    } else{
      col <- palette
    }

    theme_function <- switch(theme.choice,
                             prism = ggprism::theme_prism(),
                             classic = theme_classic(),
                             gray = theme_gray(),
                             bw = theme_bw(),
                             ggprism::theme_prism()) # 根据用户选择设置主题

    # 使用用户自定义主题（如果提供），否则使用默认主题
    theme_to_use <- if (!is.null(custom.theme)) custom.theme else theme_function

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list_all <- lapply(feature.level, function(feature.level) {
      otu_tax <-
        cbind(otu_tab,
              tax_tab %>% select(all_of(feature.level)))

      otu_tax_filtered <- otu_tax %>%
        gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        group_by_at(vars(!!sym(feature.level))) %>%
        summarise(total_count = mean(value),
                  prevalence = sum(value > 0) / n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count, -prevalence) %>%
        left_join(otu_tax, by = feature.level)

      otu_tax_agg <- otu_tax_filtered %>%
        gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        group_by_at(vars(sample, !!sym(feature.level))) %>%
        summarise(value = sum(value)) %>%
        spread(key = "sample", value = "value")

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

      otu_tax_agg_numeric <- otu_tax_agg %>%
        gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        mutate(value = as.numeric(value))

      otu_tax_agg_merged <-
        left_join(otu_tax_agg_numeric, meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        select(one_of(c("sample",
                        feature.level,
                        subject.var,
                        time.var,
                        group.var,
                        strata.var,
                        "value")))

      # Apply transformation
      if (feature.dat.type %in% c("count","proportion")){
        # Apply transformation
        if (Transform %in% c("identity", "sqrt", "log")) {
          if (Transform == "identity") {
            # No transformation needed
          } else if (Transform == "sqrt") {
            otu_tax_agg_merged$value <- sqrt(otu_tax_agg_merged$value)
          } else if (Transform == "log") {
            # Find the half of the minimum non-zero proportion for each taxon
            min_half_nonzero <- otu_tax_agg_merged %>%
              group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              summarise(min_half_value = min(value[value > 0]) / 2) %>%
              ungroup()
            # Replace zeros with the log of the half minimum non-zero proportion
            otu_tax_agg_merged <- otu_tax_agg_merged %>%
              group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              ungroup() %>%
              left_join(min_half_nonzero, by = feature.level) %>%
              mutate(value = ifelse(value == 0, log10(min_half_value), log10(value))) %>%
              select(-min_half_value)
          }
        }
      }

      taxa.levels <-
        otu_tax_agg_merged %>% select(feature.level) %>% distinct() %>% pull()

      n_subjects <- length(unique(otu_tax_agg_merged[[subject.var]]))

      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      plot_list <- lapply(taxa.levels, function(tax) {

        sub_otu_tax_agg_merged <- otu_tax_agg_merged %>% filter(!!sym(feature.level) == tax)

        boxplot <-
          ggplot(sub_otu_tax_agg_merged,
                 aes_function) +
          geom_violin(trim = FALSE, alpha = 0.8) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.1
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.1,
            fill = "white"
          ) +
          geom_jitter(width = 0.1, alpha = 0.1, size = 1) +
          scale_fill_manual(values = col) +
          {
            if (feature.dat.type == "other"){
              labs(
                x = time.var,
                y = "Abundance",
                title = tax
              )
            } else {
              labs(
                x = time.var,
                y = paste("Relative Abundance(", Transform, ")"),
                title = tax
              )
            }
          } +
          theme_to_use +
          theme(
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0, "cm"),
            plot.title = element_text(hjust = 0.5, size = 20),
            strip.text.x = element_text(size = 12, color = "black"),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(color = "black", size = base.size),
            axis.text.y = element_text(color = "black", size = (base.size-2)),
            axis.title.x = element_text(size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_text(size = 16)
          )

        if (!is.null(group.var)) {
          if (is.null(strata.var)) {
          } else {
            boxplot <-
              boxplot + ggh4x::facet_nested(as.formula(paste("~", strata.var)), scales = "free", space = "free") + theme(panel.spacing = unit(0,"lines"))
          }
        }

        if (feature.dat.type != "other"){
          # 添加对Y轴刻度的修改
          if (Transform == "sqrt") {
            boxplot <- boxplot + scale_y_continuous(
              labels = function(x) sapply(x, function(i) as.expression(substitute(a^b, list(a = i, b = 2))))
            )
          } else if (Transform == "log") {
            boxplot <- boxplot + scale_y_continuous(
              labels = function(x) sapply(x, function(i) as.expression(substitute(10^a, list(a = i))))
            )
          }
        }

        return(boxplot)
      })

      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_boxplot_single",
          "_",
          "subject_",
          subject.var,
          "_",
          "feature_level_",
          feature.level,
          "_",
          "transform_",
          Transform,
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
        pdf_name <- paste0(pdf_name,"_", feature.level, ".pdf")
        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Use lapply to print each ggplot object in the list to a new PDF page
        lapply(plot_list, print)
        # Close the PDF device
        dev.off()
      }

      return(plot_list)
    })
    return(plot_list_all)
  }
