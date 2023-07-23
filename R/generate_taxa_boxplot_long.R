#' Generate Taxa Boxplots for Longitudinal Data
#'
#' This function generates boxplots to visualize the taxonomic composition of samples for longitudinal data.
#' It also provides options for grouping and stratifying data.
#'
#' @param data.obj A list object containing the input data.
#' @param subject.var A string indicating the variable for subject identifiers.
#' @param time.var A string indicating the variable for time points.
#' @param t0.level A string indicating the reference time point. Default is NULL.
#' @param ts.levels A vector indicating the time points to include in the analysis. Default is NULL.
#' @param group.var A string indicating the variable for group identifiers. Default is NULL.
#' @param strata.var A string indicating the variable for strata identifiers. Default is NULL.
#' @param feature.level A string indicating the taxonomic level to plot.
#' @param feature.dat.type A string indicating the type of data in the input object.
#' @param features.plot A character vector of features to include in the plot. If NULL, top features will be selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot An integer indicating the top K features to plot based on the function specified in `top.k.func`.
#' @param top.k.func A function to determine the top K features to plot.
#' @param Transform A string indicating the transformation to apply to the data before plotting. Options are "identity", "sqrt", "log". Default is "identity".
#' @param prev.filter A numeric value indicating the minimum prevalence for a feature to be included in the plot.
#' @param abund.filter A numeric value indicating the minimum abundance for a feature to be included in the plot.
#' @param base.size A numeric value specifying the base size of the plot. Default is 16.
#' @param theme.choice A string specifying the theme of the plot. Default is "prism".
#' @param custom.theme A ggplot2 theme object for user-defined theme. Default is NULL.
#' @param palette A character vector specifying the color palette. Default is NULL.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of ggplot objects for each taxonomic level.
#' @details
#' This function generates a boxplot of taxa abundances for longitudinal data.
#' The boxplot can be stratified by a group variable and/or other variables.
#' It allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
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
#' generate_taxa_boxplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = NULL,
#'   feature.level = c("Phylum"),
#'   features.plot = sample(unique(ecam.obj$feature.ann[,"Phylum"]),1),
#'   feature.dat.type = "proportion",
#'   Transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = c(ggsci::pal_npg()(9),ggsci::pal_jama()(7),ggsci::pal_lancet()(9)),
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_boxplot_long(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   Transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_taxa_boxplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
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
      {
        if ("original" %in% feature.level)
          mutate(., original = rownames(.))
        else
          .
      } %>%
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
    if (is.null(palette)) {
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

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list <- lapply(feature.level, function(feature.level) {
      otu_tax <-
        cbind(otu_tab,
          tax_tab %>% select(all_of(feature.level)))

      otu_tax_filtered <- otu_tax %>%
        gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        group_by_at(vars(!!sym(feature.level))) %>%
        summarise(total_count = mean(value),
                  prevalence = sum(value > 0) / n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count,-prevalence) %>%
        left_join(otu_tax, by = feature.level)

      otu_tax_agg <- otu_tax_filtered %>%
        gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        group_by_at(vars(sample,!!sym(feature.level))) %>%
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
        gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        mutate(value = as.numeric(value))

      otu_tax_agg_merged <-
        left_join(otu_tax_agg_numeric,
                  meta_tab %>% rownames_to_column("sample"),
                  by = "sample") %>%
        select(one_of(
          c(
            "sample",
            feature.level,
            subject.var,
            time.var,
            group.var,
            strata.var,
            "value"
          )
        ))

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
              summarise(min_half_value = min(value[value > 0]) / 2) %>%
              ungroup()
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

      n_subjects <-
        length(unique(otu_tax_agg_merged[[subject.var]]))
      n_times <- length(unique(otu_tax_agg_merged[[time.var]]))

      sub_otu_tax_agg_merged <- otu_tax_agg_merged

      # 在数据处理部分创建一个新的数据框
      average_sub_otu_tax_agg_merged <- NULL
      if (n_times > 10 || n_subjects > 25) {
        if (!is.null(group.var) & !is.null(strata.var)) {
          average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
            group_by(
              !!sym(strata.var),
              !!sym(group.var),
              !!sym(time.var),
              !!sym(feature.level)
            ) %>%
            summarise(across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
            ungroup() %>%
            mutate(!!sym(subject.var) := "ALL")
        } else if (!is.null(group.var)) {
          average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
            group_by(!!sym(group.var),
                     !!sym(time.var),
                     !!sym(feature.level)) %>%
            summarise(across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
            ungroup() %>%
            mutate(!!sym(subject.var) := "ALL")
        } else {
          average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
            group_by(!!sym(time.var)) %>%
            summarise(across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
            ungroup() %>%
            mutate(!!sym(subject.var) := "ALL")
        }
      }

      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 6) {
          features.plot <- taxa.levels[1:6]
        } else {
          features.plot <- taxa.levels
        }
      }

      boxplot <-
        ggplot(sub_otu_tax_agg_merged %>% filter(!!sym(feature.level) %in% features.plot),
               aes_function) +
        geom_violin(trim = FALSE, alpha = 0.8) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
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
          linetype = "dashed",
          # 更改线条类型为虚线
          data = if (!is.null(average_sub_otu_tax_agg_merged))
            average_sub_otu_tax_agg_merged  %>% filter(!!sym(feature.level) %in% features.plot)
          else
            sub_otu_tax_agg_merged %>% filter(!!sym(feature.level) %in% features.plot)
        ) +
        scale_fill_manual(values = col) +
        {
          if (feature.dat.type == "other"){
            labs(x = time.var,
                 y = "Abundance")
          } else {
            labs(x = time.var,
                 y = paste("Relative Abundance(", Transform, ")"))
          }
        } +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          plot.title = element_text(hjust = 0.5, size = 20),
          strip.text.x = element_text(size = base.size, color = "black"),
          strip.text.y = element_text(size = base.size, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = (base.size -
                                                                2)),
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
            boxplot + ggh4x::facet_nested_wrap(as.formula(paste(
              "~", feature.level, "+", group.var
            )), scales = "fixed",
            nrow = ifelse(
              ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 4, (length(taxa.levels) + 1) / 4) < 1,
              1,
              ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 4, (length(taxa.levels) + 1) / 4)
            ))
        } else {
          boxplot <-
            boxplot + ggh4x::facet_nested_wrap(as.formula(
              paste("~", feature.level, "+", group.var, "+", strata.var)
            ),
            scales = "fixed",
            nrow = ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 2, (length(taxa.levels) + 1) / 2))
        }
      }

      # Add geom_jitter() if the number of unique time points or subjects is greater than 10
      if (n_subjects > 20 || n_times > 10) {
        boxplot <- boxplot + geom_jitter(width = 0.1,
                                         alpha = 0.1,
                                         size = 1)
      }

      if (feature.dat.type != "other"){
        # 添加对Y轴刻度的修改
        if (Transform == "sqrt") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  a ^ b, list(a = i, b = 2)
                )))
          )
        } else if (Transform == "log") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  10 ^ a, list(a = i)
                )))
          )
        }
      }

      return(boxplot)
    })

    # Save the plots as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_indiv_boxplot_long_v2",
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
      pdf_name <- paste0(pdf_name, "_", feature.level, ".pdf")
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Use lapply to print each ggplot object in the list to a new PDF page
      lapply(plot_list, print)
      # Close the PDF device
      dev.off()
    }

    return(plot_list)
  }
