#' @title Generate Boxplot of Individual Taxa Abundance Over Time
#'
#' @description This function creates a boxplot showing the abundance distribution of individual taxa at a specified taxonomic level over time from longitudinal data. It takes a MicrobiomeStat data object as input.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param t0.level The base level for time points in longitudinal data.
#' @param ts.levels The levels for time points in longitudinal data.
#' @param group.var Optional grouping variable in metadata.
#' @param strata.var Optional stratification variable in metadata.
#' @param feature.level Taxonomic level(s) for boxplots.
#' @param features.plot A character vector specifying the taxa to be plotted. If NULL (default), the top k taxa by mean abundance will be plotted.
#' @param feature.dat.type A character string specifying the type of the data in feature.dat. Options are "count", "proportion", or "other".
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
#' @param Transform Transformation to apply before plotting. Default is "log10".
#' @param prev.filter Prevalence threshold for filtering taxa. Default 0.05.
#' @param abund.filter Abundance threshold for filtering taxa. Default 0.01.
#' @param base.size Base font size for the generated plots.
#' @param theme.choice Plot theme choice (default: "bw").
#' @param custom.theme Custom ggplot2 theme (optional).
#' @param palette Color palette used for the plots.
#' @param pdf Logical, if TRUE save plot as a multi-page PDF file. Default is TRUE.
#' @param file.ann Optional string for file annotation to add to PDF name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A ggplot object showing the abundance distribution of taxa over time.
#'
#' @examples
#' # Load required libraries and example data
#' library(tidyverse)
#'
#' # Generate the boxplot pair
#' data(ecam.obj)
#' generate_taxa_indiv_boxplot_long(
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
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(peerj32.obj)
#'
#' generate_taxa_indiv_boxplot_long(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "other",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   Transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.1,
#'   base.size = 20,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_taxa_indiv_boxplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
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
    if (!is.character(time.var))
      stop("`time.var` should be a character string.")
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    # 提取数据
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

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

    line_aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = value,
        group = !!sym(subject.var),
        color = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = value,
        group = !!sym(subject.var)
      )
    }

    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = value,
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = value,
        fill = !!sym(time.var)
      )
    }

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
              dplyr::group_by(!!sym(feature.level)) %>%
              summarise(min_half_value = min(value[value > 0]) / 2) %>%
              dplyr::ungroup()
            # Replace zeros with the log of the half minimum non-zero proportion
            otu_tax_agg_merged <- otu_tax_agg_merged %>%
              left_join(min_half_nonzero, by = feature.level) %>%
              mutate(value = ifelse(value == 0, log10(min_half_value), log10(value))) %>%
              select(-min_half_value)
          }
        }
      }

      taxa.levels <-
        otu_tax_agg_merged %>% select(feature.level) %>% distinct() %>% pull()

      n_subjects <- length(unique(otu_tax_agg_merged[[subject.var]]))
      n_times <- length(unique(otu_tax_agg_merged[[time.var]]))

      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      plot_list <- lapply(taxa.levels, function(tax) {

        sub_otu_tax_agg_merged <- otu_tax_agg_merged %>% filter(!!sym(feature.level) == tax)

        # 在数据处理部分创建一个新的数据框
        average_sub_otu_tax_agg_merged <- NULL
        if (n_times > 10 || n_subjects > 25) {
          if (!is.null(group.var) && !is.null(strata.var)) {
            average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
              summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
              dplyr::ungroup() %>%
              mutate(!!sym(subject.var) := "ALL")
          } else if (!is.null(group.var)) {
            average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(group.var), !!sym(time.var)) %>%
              summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
              dplyr::ungroup() %>%
              mutate(!!sym(subject.var) := "ALL")
          } else {
            average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(time.var)) %>%
              summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
              dplyr::ungroup() %>%
              mutate(!!sym(subject.var) := "ALL")
          }
        }

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
          geom_line(
            line_aes_function,
            alpha = 1,
            linewidth = 0.6,
            color = "black",
            linetype = "dashed", # 更改线条类型为虚线
            data = if (!is.null(average_sub_otu_tax_agg_merged)) average_sub_otu_tax_agg_merged else sub_otu_tax_agg_merged
          ) +
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
              boxplot <-
                boxplot + ggh4x::facet_nested(as.formula(paste("~", group.var)), scales = "fixed")
            } else {
              boxplot <-
                boxplot + ggh4x::facet_nested(as.formula(paste("~", strata.var, "+", group.var)), scales = "free", space = "free") + theme(panel.spacing = unit(0,"lines"))
            }
          }

        # Add geom_jitter() if the number of unique time points or subjects is greater than 10
        if (n_subjects > 20 || n_times > 10) {
          boxplot <- boxplot + geom_jitter(width = 0.1, alpha = 0.1, size = 1)
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
          "taxa_indiv_boxplot_long",
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
