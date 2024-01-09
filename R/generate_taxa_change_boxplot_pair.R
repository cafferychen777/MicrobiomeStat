#' Generate Individual Taxa Change Boxplot Pair
#'
#' This function generates boxplot pairs visualizing the changes in abundance of individual taxa over time.
#' The boxplots show the change in abundance for each taxon for different groups and strata.
#' It also allows for a prevalence and abundance filter to be applied to the data, and can optionally save the plots as a PDF.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string specifying the subject variable.
#' @param time.var A character string specifying the time variable.
#' @param group.var A character string specifying the group variable. Default is NULL.
#' @param strata.var A character string specifying the strata variable. Default is NULL.
#' @param change.base The time point value to be used as the baseline for calculating change in abundance between time points. Should be one of the values from `time.var`. The change will be calculated as the difference between this baseline time point and the subsequent time point.
#' @param feature.change.func The method or function used to calculate the change in feature abundance between time points.
#' The following options are supported:
#'
#' - "relative change": Computes the relative change as (time_2 - time_1) / (time_2 + time_1). If both values are zero, the result is zero.
#' - "log fold change": Computes the log2 fold change between time points. Zero values are imputed as half the minimum nonzero value of the respective feature at the given time point before taking the logarithm.
#' - "absolute change": Computes the absolute difference between time points.
#' - A custom function: The provided function should take two numeric vectors as input (values at time 1 and time 2) and return a numeric vector of differences. Users should ensure that their function handles zero values appropriately.
#'
#' If an unrecognized value or no value is provided for `feature.change.func`, the default behavior will be to compute the absolute difference between time points.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot, when `features.plot` is NULL.
#' Default is NULL, in which case all features passing filters will be plotted.
#' @param top.k.func Function to use for selecting top k features, when `features.plot` is NULL.
#' Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
#' case features will be selected by abundance.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size The base size for the plot. Default is 16.
#' @param theme.choice Plot theme choice. Can be one of:
#'   - "prism": ggprism::theme_prism()
#'   - "classic": theme_classic()
#'   - "gray": theme_gray()
#'   - "bw": theme_bw()
#' Default is "bw".
#' @param custom.theme A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. To use a custom theme, first create a theme object with ggplot2::theme(), then pass it to this argument. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=16, color="red"),
#'   legend.position = "none"
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. Default is NULL, which will use the default theme based on `theme.choice`.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
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
#' \dontrun{
#' # Load required libraries and data
#' library(vegan)
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
#'   feature.change.func = "relative change",
#'   feature.level = "original",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 6,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#'
#' # Generate the boxplot pair
#' generate_taxa_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "relative change",
#'   feature.level = "original",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 6,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_change_boxplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           feature.change.func = "relative change",
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter,
           abund.filter,
           base.size = 16,
           theme.choice = "bw",
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

    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      ))) %>% rownames_to_column("sample")

    # Assuming mStat_get_theme function is already defined
    # Replace the existing theme selection code with this:
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    col <- mStat_get_palette(palette)

    ylab_label <- if (feature.dat.type != "other") {
      if (is.function(feature.change.func)) {
        paste0("Change in Relative Abundance", " (custom function)")
      } else {
        paste0("Change in Relative Abundance", " (", feature.change.func, ")")
      }
    }
    else {
      if (is.function(feature.change.func)) {
        paste0("Change in Abundance", " (custom function)")
      } else {
        paste0("Change in Abundance", " (", feature.change.func, ")")
      }
    }

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list <- lapply(feature.level, function(feature.level) {

      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm
      }

      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        tibble::rownames_to_column(feature.level)

      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Convert values to numeric and add sample ID
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Add metadata to the aggregated OTU table
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab, by = "sample") %>%
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
      if (is.function(feature.change.func)) {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = feature.change.func(value_time_2, value_time_1))
      } else if (feature.change.func == "log fold change") {
        half_nonzero_min_time_2 <- combined_data %>%
          filter(value_time_2 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(value_time_2) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- combined_data %>%
          filter(value_time_1 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(value_time_1) / 2,
                    .groups = "drop")

        combined_data <-
          dplyr::left_join(
            combined_data,
            half_nonzero_min_time_2,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )

        combined_data <-
          dplyr::left_join(
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
          combined_data %>% dplyr::mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (feature.change.func == "relative change") {
        combined_data <- combined_data %>%
          dplyr::mutate(value_diff = dplyr::case_when(
            value_time_2 == 0 & value_time_1 == 0 ~ 0,
            TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
          ))
      } else if (feature.change.func == "absolute change"){
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      } else {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      }

      combined_data <-
        combined_data %>% dplyr::left_join(meta_tab %>% filter(!!sym(time.var) == change.after), by = subject.var)

      taxa.levels <-
        combined_data %>% select(all_of(c(feature.level))) %>% dplyr::distinct() %>% dplyr::pull()

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

    names(plot_list) <- feature.level

    # Save the plots as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_boxplot_pair",
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
