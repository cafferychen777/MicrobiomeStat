#' Generate Individual Change Boxplot Pairs for Taxonomic Composition Data
#'
#' This function generates boxplots to visualize the change in taxonomic composition of samples between two time points in a longitudinal study.
#' It also provides options for grouping and stratifying data, and selecting the top k features based on a user-defined function.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A string indicating the variable for subject identifiers.
#' @param time.var A string indicating the variable for time points.
#' @param group.var Optional string specifying the variable for groups.
#' @param strata.var Optional string specifying the variable for strata.
#' @param change.base A string indicating the base time point for change computation. This should match one of the time points present in the metadata for the 'time.var' variable.
#' @param feature.change.func Specifies the method or function used to compute the change between two time points.
#' The following options are available:
#'
#' - "absolute change": Computes the difference between the abundance values at the two time points (`value_time_2` and `value_time_1`).
#'
#' - "log fold change": Computes the log2 fold change between the two time points. For zero abundances, imputation is performed using half of the minimum nonzero value for each feature level at the respective time point before taking the logarithm.
#'
#' - "relative change": Computes the relative change as `(value_time_2 - value_time_1) / (value_time_2 + value_time_1)`. If both time points have an abundance of 0, the change is defined as 0.
#'
#' - A custom function: If a user-defined function is provided, it should take two numeric vectors as input corresponding to the abundances at the two time points (`value_time_1` and `value_time_2`) and return a numeric vector of the computed change. This custom function will be applied directly.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
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
#' @param base.size A numeric value specifying the base font size of the plot.
#' @param theme.choice
#' Plot theme choice. Specifies the visual style of the plot. Can be one of the following pre-defined themes:
#'   - "prism": Utilizes the ggprism::theme_prism() function from the ggprism package, offering a polished and visually appealing style.
#'   - "classic": Applies theme_classic() from ggplot2, providing a clean and traditional look with minimal styling.
#'   - "gray": Uses theme_gray() from ggplot2, which offers a simple and modern look with a light gray background.
#'   - "bw": Employs theme_bw() from ggplot2, creating a classic black and white plot, ideal for formal publications and situations where color is best minimized.
#'   - "light": Implements theme_light() from ggplot2, featuring a light theme with subtle grey lines and axes, suitable for a fresh, modern look.
#'   - "dark": Uses theme_dark() from ggplot2, offering a dark background, ideal for presentations or situations where a high-contrast theme is desired.
#'   - "minimal": Applies theme_minimal() from ggplot2, providing a minimalist theme with the least amount of background annotations and colors.
#'   - "void": Employs theme_void() from ggplot2, creating a blank canvas with no axes, gridlines, or background, ideal for custom, creative plots.
#' Each theme option adjusts various elements like background color, grid lines, and font styles to match the specified aesthetic.
#' Default is "bw", offering a universally compatible black and white theme suitable for a wide range of applications.
#' @param custom.theme
#' A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. Custom themes are useful for creating publication-ready figures with specific formatting requirements.
#'
#' To use a custom theme, create a theme object with ggplot2::theme(), including any desired customizations. Common customizations for publication-ready figures might include adjusting text size for readability, altering line sizes for clarity, and repositioning or formatting the legend. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=14, face="bold"),        # Bold axis titles with larger font
#'   axis.text = ggplot2::element_text(size=12),                      # Slightly larger axis text
#'   legend.position = "top",                                         # Move legend to the top
#'   legend.background = ggplot2::element_rect(fill="lightgray"),     # Light gray background for legend
#'   panel.background = ggplot2::element_rect(fill="white", colour="black"), # White panel background with black border
#'   panel.grid.major = ggplot2::element_line(colour = "grey90"),     # Lighter color for major grid lines
#'   panel.grid.minor = ggplot2::element_blank(),                     # Remove minor grid lines
#'   plot.title = ggplot2::element_text(size=16, hjust=0.5)           # Centered plot title with larger font
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. If `custom.theme` is NULL (the default), the theme is determined by `theme.choice`. This flexibility allows for both easy theme selection for general use and detailed customization for specific presentation or publication needs.
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
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of ggplot objects, one for each taxonomic level.
#' @details
#' This function generates a boxplot of the change in taxa abundances between two time points in a longitudinal study.
#' The boxplot can be stratified by a group variable and/or other variables.
#' It allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and data
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' # Generate the boxplot pair
#' generate_taxa_indiv_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 20,
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
#' generate_taxa_indiv_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 20,
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
generate_taxa_indiv_change_boxplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           feature.change.func = "relative change",
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter = 0.01,
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

    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    plot_list_all <- lapply(feature.level, function(feature.level) {

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
        rownames_to_column(feature.level)

      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Convert values to numeric and add sample ID
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Add metadata to the aggregated OTU table
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab, by = "sample") %>%
        select(all_of(c("sample",
               feature.level,
               subject.var,
               time.var,
               group.var,
               strata.var,
               "value")))

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
        combined_data <- combined_data %>% dplyr::mutate(value_diff = feature.change.func(value_time_2, value_time_1))
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

        combined_data <- dplyr::left_join(combined_data, half_nonzero_min_time_2, by = feature.level, suffix = c("_time_1", "_time_2"))
        combined_data <- dplyr::left_join(combined_data, half_nonzero_min_time_1, by = feature.level, suffix = c("_time_1", "_time_2"))
        combined_data$value_time_2[combined_data$value_time_2 == 0] <- combined_data$half_nonzero_min_time_2[combined_data$value_time_2 == 0]
        combined_data$value_time_1[combined_data$value_time_1 == 0] <- combined_data$half_nonzero_min_time_1[combined_data$value_time_1 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message("Imputation was performed using half the minimum nonzero proportion for each taxon at different time points.")

        combined_data <- combined_data %>% dplyr::mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (feature.change.func == "relative change"){
        combined_data <- combined_data %>%
          dplyr::mutate(value_diff = dplyr::case_when(
            value_time_2 == 0 & value_time_1 == 0 ~ 0,
            TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
          ))
      } else if (feature.change.func == "absolute change"){
        combined_data <- combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      } else {
        combined_data <- combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      }

      combined_data <-
        combined_data %>% dplyr::left_join(meta_tab %>% filter(!!sym(time.var) == change.after), by = subject.var)

      taxa.levels <-
        combined_data %>% select(all_of(c(feature.level))) %>% dplyr::distinct() %>% dplyr::pull()

      if (is.null(group.var)) {
        group.var = "group"
        combined_data$group <- "ALL"
      }

      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      plot_list <- lapply(taxa.levels, function(tax) {
        # Create the boxplot
        boxplot <-
          ggplot(
            combined_data %>% filter(!!sym(feature.level) == tax),
            aes(
              x = !!sym(group.var),
              y = value_diff,
              fill = !!sym(group.var)
            )
          ) +
          #geom_violin(trim = F, alpha = 0.8) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.1
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.1,
            #fill = "white"
          ) +
          geom_jitter(width = 0.1,
                      alpha = 0.3,
                      size = 1.7) +
          scale_alpha_manual(values = c(0.5, 0.5)) +
          scale_fill_manual(values = col) +
          labs(
            x = group.var,
            y = ylab_label,
            title = tax
          ) +
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
            legend.title = ggplot2::element_text(size = 16)
          )

        if (!is.null(strata.var)) {
          boxplot <- boxplot +
            ggh4x::facet_nested(cols = vars(!!sym(strata.var)),
                       scales = "fixed",
                       space = "free")
        }

        if (group.var == "group" && unique(combined_data$group)[1] == "ALL") {
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
          "taxa_indiv_change_boxplot_pair",
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

      names(plot_list) <- taxa.levels

      return(plot_list)
    })
    names(plot_list_all) <- feature.level
    # Return the boxplot for display
    return(plot_list_all)
  }
