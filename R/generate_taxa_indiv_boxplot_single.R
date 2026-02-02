#' @title Generate Individual Taxa Boxplots for Single Time Point
#'
#' @description Creates boxplots showing the abundance distribution of individual taxa
#' at specified taxonomic levels for a single time point.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param t.level Character string specifying the time level/value to subset data to.
#'   Default NULL does not subset data.
#' @param features.plot A character vector or named list specifying which feature IDs to plot.
#'   If a named list, each element should correspond to a taxonomic level in `feature.level`.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param transform Transformation to apply: "identity" (default), "sqrt", or "log".
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
           time.var = NULL,
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

    # Check if input variables are properly specified
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    # Subset the data if time variable and level are provided
    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    # Select relevant variables from metadata (filter out NULL vars)
    meta_tab <- select_meta_vars(data.obj$meta.dat, group.var, time.var, strata.var)

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
      otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

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
        tidyr::gather(key = "sample", value = "value", -all_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Merge feature data with metadata
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        select(all_of(c("sample",
                        feature.level,
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

      # Count number of samples
      n_subjects <- length(unique(otu_tax_agg_merged[["sample"]]))

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
