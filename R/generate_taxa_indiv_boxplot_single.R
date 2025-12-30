#' Generate Taxa Boxplots for Single Time Point
#'
#' This function generates boxplots to visualize the taxonomic composition of samples for a single time point in a longitudinal study.
#' It provides options for grouping and stratifying data, and selecting the top k features based on a user-defined function.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A string indicating the variable for subject identifiers.
#' @param time.var A string indicating the variable for time points. If NULL, the function assumes that data for a single time point is provided.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var Optional string specifying the variable for groups.
#' @param strata.var Optional string specifying the variable for strata.
#' @param feature.level A string specifying the taxonomic level to plot.
#' @param features.plot A character vector or named list specifying which feature IDs (e.g. OTU IDs) to plot.
#' If a character vector, the same features will be plotted for all taxonomic levels.
#' If a named list, each element should correspond to a taxonomic level in `feature.level`,
#' with the value being a character vector of feature IDs for that level.
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
#' @param transform A string indicating the transformation to apply to the data before plotting. Options are:
#' - "identity": No transformation (default)
#' - "sqrt": Square root transformation
#' - "log": Logarithmic transformation. Zeros are replaced with half of the minimum non-zero value for each taxon before log transformation.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size A numeric value specifying the base font size for the plot.
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
#' @param pdf.wid A numeric value specifying the width of the PDF file. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF file. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of ggplot objects for each taxonomic level.
#' @details
#' This function generates a boxplot of taxa abundances for a single time point in a longitudinal study.
#' The boxplot can be stratified by a group variable and/or other variables.
#' It allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#' # Generate the boxplot pair
#' data(ecam.obj)
#' generate_taxa_indiv_boxplot_single(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t.level = NULL,
#'   group.var = "antiexposedall",
#'   strata.var = NULL,
#'   feature.level = c("Phylum","Family"),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 20,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_indiv_boxplot_single(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Family"),
#'   features.plot = list(Family = c("Bacteroidaceae", "Lachnospiraceae")),
#'   feature.dat.type = "count",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
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
#' }
#' @export
generate_taxa_indiv_boxplot_single <-
  function(data.obj,
           subject.var,
           time.var,
           t.level = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           top.k.plot = NULL,
           top.k.func = NULL,
           transform = c("sqrt", "identity", "log"),
           prev.filter = 0.01,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    # Match and validate input arguments
    feature.dat.type <- match.arg(feature.dat.type)
    transform <- match.arg(transform)

    # Validate the input data object
    mStat_validate_data(data.obj)

    # If subject variable is not provided, create a default one
    if (is.null(subject.var)){
      data.obj$meta.dat$subject.id <- rownames(data.obj$meta.dat)
      subject.var <- "subject.id"
    }

    # Check if input variables are properly specified
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    # Subset the data if time variable and level are provided
    if (!is.null(time.var) &!is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    # Select relevant variables from metadata (filter out NULL vars)
    meta_tab <- select_meta_vars(data.obj$meta.dat, subject.var, group.var, time.var, strata.var)

    # Define aesthetic function for plotting
    # This function determines how the data will be mapped to visual properties in the plot
    aes_function <-  aes(
      x = !!sym(group.var),
      y = value,
      fill = !!sym(group.var)
    )

    # Get color palette for the plot
    col <- mStat_get_palette(palette)

    # Get the appropriate theme for the plot
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Adjust filtering parameters based on input conditions
    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    # For "other" data type, check if data contains negative values
    # If so, adjust abundance filter to handle negative values appropriately
    if (feature.dat.type == "other") {
      # Check if any feature table contains negative values
      has_negative <- FALSE
      if (!is.null(data.obj$feature.tab)) {
        has_negative <- any(data.obj$feature.tab < 0, na.rm = TRUE)
      }
      if (!has_negative && !is.null(data.obj$feature.agg.list)) {
        for (agg_table in data.obj$feature.agg.list) {
          if (any(agg_table < 0, na.rm = TRUE)) {
            has_negative <- TRUE
            break
          }
        }
      }

      if (has_negative) {
        message("Note: Negative values detected in 'other' data type. Abundance filtering is disabled to preserve all features.")
        abund.filter <- -Inf  # Set to negative infinity to include all features regardless of abundance
      }
    }

    # Normalize count data if necessary
    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    # Generate plots for each taxonomic level
    plot_list_all <- lapply(feature.level, function(feature.level) {

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Select appropriate feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Filter and prepare the feature table
      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        rownames_to_column(feature.level)

      # Determine features to plot for this specific taxonomic level
      current_features_plot <- NULL

      # Handle features.plot parameter (support both vector and list formats)
      if (!is.null(features.plot)) {
        if (is.list(features.plot)) {
          # If features.plot is a list, extract features for current level
          if (feature.level %in% names(features.plot)) {
            current_features_plot <- features.plot[[feature.level]]
          }
        } else {
          # If features.plot is a vector, use it for all levels (backward compatibility)
          current_features_plot <- features.plot
        }
      }

      # Select top k features if no specific features are provided
      if (is.null(current_features_plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
        computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
        current_features_plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Reshape the data for plotting
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Merge feature data with metadata
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        select(one_of(c("sample",
                        feature.level,
                        subject.var,
                        time.var,
                        group.var,
                        strata.var,
                        "value")))

      # Apply data transformation if specified
      if (feature.dat.type %in% c("count","proportion")){
        if (transform %in% c("identity", "sqrt", "log")) {
          if (transform == "identity") {
            # No transformation needed
          } else if (transform == "sqrt") {
            otu_tax_agg_merged$value <- sqrt(otu_tax_agg_merged$value)
          } else if (transform == "log") {
            # For log transformation, we need to handle zeros
            # We replace zeros with half of the minimum non-zero value for each taxon
            min_half_nonzero <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::summarise(min_half_value = min(value[value > 0]) / 2) %>%
              dplyr::ungroup()
            otu_tax_agg_merged <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::ungroup() %>%
              dplyr::left_join(min_half_nonzero, by = feature.level) %>%
              dplyr::mutate(value = ifelse(value == 0, log10(min_half_value), log10(value))) %>%
              select(-min_half_value)
          }
        }
      }

      # Determine if there are many samples
      manysample <- calc_sample_count(data.obj$meta.dat, time.var, group.var) > 8
      wds <- if(manysample) 7 else 3

      # Get unique taxa levels
      taxa.levels <-
        otu_tax_agg_merged %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()

      # Count number of subjects
      n_subjects <- length(unique(otu_tax_agg_merged[[subject.var]]))

      # Filter taxa levels if specified
      if (!is.null(current_features_plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% current_features_plot]
      }

      # Determine the appropriate Y-axis label for "other" data type
      # Check the overall aggregated data for negative values, not individual taxa
      if (feature.dat.type == "other") {
        overall_data_range <- range(otu_tax_agg_merged$value, na.rm = TRUE)
        overall_has_negative <- any(otu_tax_agg_merged$value < 0, na.rm = TRUE)

        # Intelligent label inference based on overall data characteristics
        if (overall_has_negative && overall_data_range[1] > -15 && overall_data_range[2] < 15) {
          global_y_label <- "Log-transformed Abundance"
        } else if (overall_has_negative) {
          global_y_label <- "Transformed Abundance"
        } else if (all(otu_tax_agg_merged$value >= 0 &
                       otu_tax_agg_merged$value <= 1, na.rm = TRUE)) {
          global_y_label <- "Normalized Abundance (0-1)"
        } else {
          global_y_label <- "Abundance"
        }
      }

      # Generate individual plots for each taxon
      plot_list <- lapply(taxa.levels, function(tax) {

        sub_otu_tax_agg_merged <- otu_tax_agg_merged %>% filter(!!sym(feature.level) == tax)

        # Remove outliers from the data before plotting
        sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
          dplyr::group_by(!!sym(group.var)) %>%
          dplyr::mutate(
            Q1 = quantile(value, 0.25, na.rm = TRUE),
            Q3 = quantile(value, 0.75, na.rm = TRUE),
            IQR = Q3 - Q1,
            lower_fence = Q1 - 1.5 * IQR,
            upper_fence = Q3 + 1.5 * IQR
          ) %>%
          dplyr::filter(value >= lower_fence & value <= upper_fence) %>%
          dplyr::select(-Q1, -Q3, -IQR, -lower_fence, -upper_fence) %>%
          dplyr::ungroup()

        # Create the boxplot
        boxplot <-
          ggplot(sub_otu_tax_agg_merged,
                 aes_function) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.3
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.3
          ) +
          geom_jitter(
            width = 0.35,
            alpha = 0.6,
            size = 2,
            aes(color = !!sym(group.var))
          ) +
          scale_fill_manual(values = col) +
          scale_color_manual(values = col) +
          {
            if (feature.dat.type == "other"){
              labs(
                x = group.var,
                y = global_y_label,  # Use the globally determined label
                title = tax
              )
            } else {
              labs(
                x = group.var,  # Fixed: add X-axis label for consistency
                y = paste("Relative Abundance (", transform, ")", sep = ""),
                title = tax
              )
            }
          } +
          theme_to_use +
          theme(
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0, "cm"),
            plot.title = element_text(hjust = 0.5, size = 20),
            strip.text.x = element_text(size = base.size, color = "black"),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(color = "black", size = base.size * 0.8),
            axis.text.y = element_text(color = "black", size = (base.size-2)),
            axis.title.x = element_text(size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_text(size = 16)
          )

        # Adjust x-axis labels if there are more than 2 groups
        if(length(unique(data.obj$meta.dat[[group.var]])) > 2){
          boxplot <- boxplot +
            scale_x_discrete(guide = guide_axis(n.dodge = 2))
        }

        # Add faceting if strata variable is provided
        if (!is.null(group.var)) {
          if (!is.null(strata.var)) {
            boxplot <-
              boxplot + ggh4x::facet_nested(as.formula(paste("~", strata.var)), scales = "free", space = "free") + theme(panel.spacing = unit(0,"lines"))
          }
        }

        # Modify y-axis scale based on the transformation
        if (feature.dat.type != "other"){
          if (transform == "sqrt") {
            boxplot <- boxplot + scale_y_continuous(
              labels = function(x) sapply(x, function(i) as.expression(substitute(a^b, list(a = i, b = 2))))
            )
          } else if (transform == "log") {
            boxplot <- boxplot + scale_y_continuous(
              labels = function(x) sapply(x, function(i) as.expression(substitute(10^a, list(a = i))))
            )
          }
        }

        return(boxplot)
      })

      # Save the plots as a PDF file if requested
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
          transform,
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

      names(plot_list) <- taxa.levels
      return(plot_list)
    })

    names(plot_list_all) <- feature.level
    return(plot_list_all)
  }