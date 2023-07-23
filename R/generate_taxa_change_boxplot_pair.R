#' Generate Individual Taxa Change Boxplot Pair
#'
#' This function generates boxplot pairs visualizing the changes in abundance of individual taxa over time.
#' The boxplots show the change in abundance for each taxon for different groups and strata.
#' It also allows for a prevalence and abundance filter to be applied to the data, and can optionally save the plots as a PDF.
#'
#' @param data.obj A data object containing the data to be plotted.
#' @param subject.var A character string specifying the subject variable.
#' @param time.var A character string specifying the time variable.
#' @param group.var A character string specifying the group variable. Default is NULL.
#' @param strata.var A character string specifying the strata variable. Default is NULL.
#' @param change.base The time point to be used as the base for calculating change in abundance.
#' @param change.func The function to be used for calculating change in abundance. Default is "relative difference".
#' @param feature.level The feature level at which to plot the data.
#' @param feature.dat.type The type of the feature data. Can be "count", "proportion", or "other". Default is "count".
#' @param features.plot A character vector specifying the features to be plotted. Default is NULL.
#' @param top.k.plot The top k features to be plotted. Default is NULL.
#' @param top.k.func The function to be used for selecting the top k features. Default is NULL.
#' @param prev.filter The prevalence filter to be applied to the data.
#' @param abund.filter The abundance filter to be applied to the data.
#' @param base.size The base size for the plot. Default is 16.
#' @param theme.choice The theme choice for the plot. Default is "prism".
#' @param custom.theme An optional custom theme for the plot. Default is NULL.
#' @param palette A character vector specifying the color palette for the plot. Default is NULL.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann An optional character string to be appended to the file name of the PDF. Default is NULL.
#' @param pdf.wid The width of the PDF. Default is 11.
#' @param pdf.hei The height of the PDF. Default is 8.5.
#' @param ... Additional arguments passed to the underlying functions.
#'
#'
#' @return A list of ggplot objects, each of which is a boxplot visualizing the changes in abundance of individual taxa over time.
#'
#' @examples
#' # Load required libraries and data
#' library(vegan)
#' library(tidyverse)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' # Generate the boxplot pair
#' generate_taxa_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   change.base = "1",
#'   change.func = "lfc",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 6,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.01,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_taxa_change_boxplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           change.func = "relative difference",
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter,
           abund.filter,
           base.size = 16,
           theme.choice = "prism",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    # Data validation
    mStat_validate_data(data.obj)

    feature.dat.type <- match.arg(feature.dat.type)

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

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      ))) %>% rownames_to_column("sample")

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

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list <- lapply(feature.level, function(feature.level) {
      # Merge OTU table with taxonomy table
      otu_tax <-
        cbind(otu_tab, tax_tab %>% select(all_of(c(feature.level))))

      # Filter taxa based on prevalence and abundance
      otu_tax_filtered <- otu_tax %>%
        gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        group_by_at(vars(!!sym(feature.level))) %>%
        summarise(total_value = mean(value),
                  prevalence = sum(value > 0) / n()) %>%
        filter(prevalence >= prev.filter, total_value >= abund.filter) %>%
        select(-all_of(c("total_value", "prevalence"))) %>%
        left_join(otu_tax, by = feature.level)

      # Aggregate OTU table
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

      # Convert values to numeric and add sample ID
      otu_tax_agg_numeric <- otu_tax_agg %>%
        gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        mutate(value = as.numeric(value))

      # Add metadata to the aggregated OTU table
      otu_tax_agg_merged <-
        left_join(otu_tax_agg_numeric, meta_tab, by = "sample") %>%
        select(all_of(
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

      change.after <-
        unique(otu_tax_agg_merged %>% select(all_of(c(time.var))))[unique(otu_tax_agg_merged %>% select(all_of(c(time.var)))) != change.base]

      # Calculate the change in abundance for each taxa
      # 拆分成一个列表，每个time值都有一个独立的tibble
      split_data <-
        split(otu_tax_agg_merged,
              f = otu_tax_agg_merged %>%
                dplyr::group_by(!!sym(time.var)) %>% select(all_of(c(time.var))))

      # 提取split_data中的第一个和第二个表
      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # 将这两个表连接在一起，以便计算差值
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # 计算value的差值
      if (is.function(change.func)) {
        combined_data <-
          combined_data %>% mutate(value_diff = change.func(value_time_2, value_time_1))
      } else if (change.func == "lfc") {
        half_nonzero_min_time_2 <- combined_data %>%
          filter(value_time_2 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          summarize(half_nonzero_min = min(value_time_2) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- combined_data %>%
          filter(value_time_1 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          summarize(half_nonzero_min = min(value_time_1) / 2,
                    .groups = "drop")

        combined_data <-
          left_join(
            combined_data,
            half_nonzero_min_time_2,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data <-
          left_join(
            combined_data,
            half_nonzero_min_time_1,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data$value_time_2[combined_data$value_time_2 == 0] <-
          combined_data$half_nonzero_min_time_2[combined_data$value_time_2 == 0]
        combined_data$value_time_1[combined_data$value_time_1 == 0] <-
          combined_data$half_nonzero_min_time_1[combined_data$value_time_1 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message(
          "Imputation was performed using half the minimum nonzero proportion for each taxon at different time points."
        )

        combined_data <-
          combined_data %>% mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (change.func == "relative difference") {
        combined_data <- combined_data %>%
          mutate(value_diff = case_when(
            value_time_2 == 0 & value_time_1 == 0 ~ 0,
            TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
          ))
      } else {
        combined_data <-
          combined_data %>% mutate(value_diff = value_time_2 - value_time_1)
      }

      combined_data <-
        combined_data %>% left_join(meta_tab %>% filter(!!sym(time.var) == change.after), by = subject.var)

      taxa.levels <-
        combined_data %>% select(all_of(c(feature.level))) %>% distinct() %>% pull()

      if (is.null(group.var)) {
        group.var = "group"
        combined_data$group <- "ALL"
      }

      if (!is.null(features.plot)) {
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      } else {
        if (length(taxa.levels) >= 5) {
          taxa.levels <- taxa.levels[1:4]
        } else {
          taxa.levels <- taxa.levels
        }
      }

      # Create the boxplot
      boxplot <-
        ggplot(
          combined_data %>% filter(!!sym(feature.level) %in% taxa.levels),
          aes(
            x = !!sym(group.var),
            y = value_diff,
            fill = !!sym(group.var)
          )
        ) +
        geom_violin(trim = F, alpha = 0.8) +
        geom_jitter(width = 0.1,
                    alpha = 0.3,
                    size = 1.5) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
          fill = "white"
        ) +
        scale_alpha_manual(values = c(0.5, 0.5)) +
        scale_fill_manual(values = col) +
        labs(x = group.var,
             y = ylab_label) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          plot.title = element_text(hjust = 0.5, size = 20),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16),
          ...
        )

      if (!is.null(strata.var)) {
        boxplot <- boxplot +
          ggh4x::facet_nested_wrap(
            as.formula(paste("~",feature.level,"+",strata.var)),
            scales = "fixed",
            nrow = ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 2, (length(taxa.levels) + 1) / 2)
          )
      } else {
        boxplot <- boxplot +
          ggh4x::facet_nested_wrap(
            as.formula(paste("~",feature.level)),
            scales = "fixed",
            nrow = ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 2, (length(taxa.levels) + 1) / 2)
          )
      }

      if (group.var == "group" &&
          unique(combined_data$group)[1] == "ALL") {
        boxplot <- boxplot +
          theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank()
          )
      }

      return(boxplot)

    })

    # Save the plots as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_indiv_change_boxplot_pair_v2",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "change_base_",
        change.base,
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
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Use lapply to print each ggplot object in the list to a new PDF page
      lapply(plot_list, print)
      # Close the PDF device
      dev.off()
    }

    # Return the boxplot for display
    return(plot_list)
  }
