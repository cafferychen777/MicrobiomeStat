#' Check if a Variable is Continuous Numeric
#'
#' This function checks if a given variable is continuous numeric by checking if it is numeric and has at least 10 unique values.
#'
#' @param x A variable that you want to check.
#'
#' @return Logical. Returns TRUE if the variable is continuous numeric, FALSE otherwise.
#'
#' @examples
#' is_continuous_numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) # TRUE
#' is_continuous_numeric(c(1, 2, 3, 4, 5)) # FALSE
#' is_continuous_numeric(c('a', 'b', 'c')) # FALSE
#'
#' @export
is_continuous_numeric <- function(x) {
  if (is.numeric(x) && length(unique(x)) >= 10) {
    # 10是你选择的阈值，你可以根据需要调整
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Generate Individual Change Scatterplot Pairs for Taxonomic Composition Data
#'
#' This function generates scatterplots to visualize the change in taxonomic composition of samples between two time points in a longitudinal study.
#' It also provides options for grouping and stratifying data, and selecting the top k features based on a user-defined function.
#'
#' @param data.obj A list object containing the input data.
#' @param subject.var A string indicating the variable for subject identifiers.
#' @param time.var A string indicating the variable for time points.
#' @param group.var Optional string specifying the variable for groups.
#' @param strata.var Optional string specifying the variable for strata.
#' @param change.base A string indicating the base time point for change computation.
#' @param change.func A string or function to compute the change in abundance. If a string, it should be one of "difference", "relative difference", or "lfc" (log fold change). If a function, it should take two arguments representing the abundances at two time points and return a numeric value indicating the change.
#' @param feature.level A string indicating the taxonomic level to plot.
#' @param features.plot A character vector of features to include in the plot. If NULL, top features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type A string indicating the type of data in the input object. Options are "count", "proportion", "other".
#' @param top.k.plot An integer indicating the top K features to plot based on the function specified in `top.k.func`.
#' @param top.k.func A function to determine the top K features to plot.
#' @param prev.filter A numeric value indicating the minimum prevalence for a feature to be included in the plot.
#' @param abund.filter A numeric value indicating the minimum abundance for a feature to be included in the plot.
#' @param base.size A numeric value specifying the base font size for the plot.
#' @param theme.choice A string specifying the theme of the plot. Default is "bw".
#' @param custom.theme A custom ggplot theme if the user wants to apply it. Default is NULL.
#' @param palette A character vector specifying the color palette. Default is NULL.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF file. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF file. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of ggplot objects, one for each taxonomic level.
#' @details
#' This function generates a scatterplot of the change in taxa abundances between two time points in a longitudinal study.
#' The scatterplot can be stratified by a group variable and/or other variables.
#' It also allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#'  library(vegan)
#'  data(peerj32.obj)
#'
#'  # Generate the scatterplot pairs
#'  plot_list <- generate_taxa_indiv_change_scatterplot_pair(
#'    data.obj = peerj32.obj,
#'    subject.var = "subject",
#'    time.var = "time",
#'    group.var = "cons",
#'    strata.var = "sex",
#'    change.base = "1",
#'    taxa.level = "Phylum",
#'    top.k.plot = 3,
#'    prev.filter = 0.1,
#'    abund.filter = 0.01
#'  )
#' }
#'
#' @export
generate_taxa_indiv_change_scatterplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           change.func = "difference",
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
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

    mStat_validate_data(data.obj)

    feature.dat.type <- match.arg(feature.dat.type)

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
      {if("original" %in% feature.level) dplyr::mutate(., original = rownames(.)) else .} %>%
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

    plot_list_all <- lapply(feature.level, function(feature.level) {
      # Filter taxa based on prevalence and abundance
      otu_tax_filtered <- otu_tax %>%
        tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        dplyr::group_by_at(vars(!!sym(feature.level))) %>%
        dplyr::summarise(total_count = mean(count),
                  prevalence = sum(count > 0) / dplyr::n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count, -prevalence) %>%
        dplyr::left_join(otu_tax, by = feature.level)

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax_filtered %>%
        tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        dplyr::group_by_at(vars(sample, !!sym(feature.level))) %>%
        dplyr::summarise(count = sum(count)) %>%
        tidyr::spread(key = "sample", value = "count")

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
                     matrixStats::rowSds(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
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
        dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      otu_tab_norm <- otu_tax_agg_numeric %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      otu_tab_norm_agg <- otu_tax_agg_numeric %>%
        tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
        dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      taxa.levels <-
        otu_tab_norm_agg %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # 首先，把数据分为两个子集，一个为change.base，一个为change.after
      df_t0 <- otu_tab_norm_agg %>% filter(!!sym(time.var) == change.base)
      df_ts <- otu_tab_norm_agg %>% filter(!!sym(time.var) != change.after)

      # 然后，使用dplyr::inner_join合并这两个子集，基于Phylum、subject和sex
      df <- dplyr::inner_join(df_ts, df_t0, by = c(feature.level, subject.var), suffix = c("_ts", "_t0"), relationship = "many-to-many")

      # 最后，计算新的count值
      if (is.function(change.func)) {
        df <- df %>% dplyr::mutate(new_count = change.func(count_ts, count_t0))
      } else if (change.func == "lfc") {
        # 对于对数折叠变化("lfc")，我们需要插补数据
        # 首先，为每个分类计算非零最小值的一半
        half_nonzero_min_time_2 <- df %>%
          filter(count_ts > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(count_ts) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- df %>%
          filter(count_t0 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(count_t0) / 2,
                    .groups = "drop")

        # 然后，用这些值来插补数据
        df <- dplyr::left_join(df, half_nonzero_min_time_2, by = feature.level, suffix = c("_t0", "_ts"))
        df <- dplyr::left_join(df, half_nonzero_min_time_1, by = feature.level, suffix = c("_t0", "_ts"))
        df$count_ts[df$count_ts == 0] <- df$half_nonzero_min_ts[df$count_ts == 0]
        df$count_t0[df$count_t0 == 0] <- df$half_nonzero_min_t0[df$count_t0 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message("Imputation was performed using half the minimum nonzero count for each taxa at different time points.")

        df <- df %>% dplyr::mutate(new_count = log2(count_ts) - log2(count_t0))
      } else if (change.func == "relative difference"){
        df <- df %>%
          dplyr::mutate(new_count = case_when(
            count_ts == 0 & count_t0 == 0 ~ 0,
            TRUE ~ (count_ts - count_t0) / (count_ts + count_t0)
          ))
      } else {
        df <- df %>% dplyr::mutate(new_count = count_ts - count_t0)
      }

      df <- df %>% dplyr::left_join(meta_tab %>% select(-all_of(time.var)) %>% dplyr::distinct(), by = c(subject.var))

      df <- df %>% setNames(ifelse(names(.) == paste0(time.var,"_ts"), time.var, names(.)))

      # 提前判断change.func的类型，如果是自定义函数则给出特定的标签，否则保持原样
      ylab_label <- if (feature.dat.type != "other") {
        if (is.function(change.func)) {
          paste0("Change in Relative Abundance", " (custom function)")
        } else {
          paste0("Change in Relative Abundance", " (", change.func, ")")
        }
      }
      else {
        if (is.function(change.func)) {
          paste0("Change in Abundance", " (custom function)")
        } else {
          paste0("Change in Abundance", " (", change.func, ")")
        }
      }

      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      plot_list <- lapply(taxa.levels, function(tax) {
        scatterplot <-
          ggplot(df %>% filter(!!sym(feature.level) == tax),
                 aes(
                   x = !!sym(group.var),
                   y = new_count,
                   fill = !!sym(strata.var),
                   color = !!sym(strata.var),  # Add this line to tidyr::separate color by strata.var
                   group = !!sym(strata.var)   # Add this line to tidyr::separate smoothing line by strata.var
                 )) +
          geom_smooth(se = TRUE, method = 'lm') +
          geom_point(aes_function,data = df %>% filter(!!sym(feature.level) == tax),
                     size = 4) +
          scale_shape_manual(values = c(21, 22, 24, 25)) +
          scale_fill_manual(values = colors) +
          ylab(ylab_label) +
          ggtitle(tax) +
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
          )

        if (group.var == "ALL"){
          scatterplot <- scatterplot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
        }

        return(scatterplot)
      })

      # Save the stacked dotplot as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_change_scatterplot",
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
        lapply(plot_list, print)
        # Close the PDF device
        dev.off()
      }
      return(plot_list)
    })

    return(plot_list_all)
  }
