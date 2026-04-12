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
           feature.level,
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
    data.obj <- mStat_validate_data(data.obj)

    # Check if input variables are properly specified
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    context <- mStat_prepare_taxa_single_context(
      data.obj = data.obj,
      time.var = time.var,
      t.level = t.level,
      group.var = group.var,
      strata.var = strata.var
    )
    data.obj <- context$data.obj

    # Select relevant variables from metadata (filter out NULL vars)
    meta_tab <- context$meta_tab
    placeholder_group <- mStat_ensure_group_placeholder(
      meta_tab,
      group.var = group.var,
      value = "ALL",
      column_name = "ALL"
    )
    meta_tab <- placeholder_group$df
    resolved_group_var <- placeholder_group$group.var

    # Define aesthetic function for plotting
    # This function determines how the data will be mapped to visual properties in the plot
    aes_function <-  aes(
      x = !!sym(resolved_group_var),
      y = value,
      fill = !!sym(resolved_group_var)
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

    abund.filter <- mStat_adjust_other_abundance_filter(
      data.obj = data.obj,
      feature.dat.type = feature.dat.type,
      abund.filter = abund.filter
    )

    # Normalize count data if necessary.
    analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

    # Generate plots for each taxonomic level
    plot_list_all <- lapply(feature.level, function(feature.level) {

      # Aggregate data by taxonomy if necessary
      otu_tax_agg <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter)

      # Determine features to plot for this specific taxonomic level
      current_features_plot <- mStat_resolve_selected_features(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        features.plot = if (is.list(features.plot)) features.plot[[feature.level]] else features.plot,
        top.k.plot = top.k.plot,
        top.k.func = top.k.func
      )

      otu_tax_agg_merged <-
        mStat_prepare_taxa_long_data(
          feature.dat = otu_tax_agg,
          feature.level = feature.level,
          value_col = "value",
          meta.dat = meta_tab
        ) %>%
        select(all_of(c("sample",
                        feature.level,
                        time.var,
                        resolved_group_var,
                        strata.var,
                        "value")))

      otu_tax_agg_merged <- mStat_transform_taxa_long_values(
        long.df = otu_tax_agg_merged,
        feature.level = feature.level,
        value_col = "value",
        feature.dat.type = feature.dat.type,
        transform = transform
      )

      # Determine if there are many samples
      manysample <- calc_sample_count(data.obj$meta.dat, time.var, group.var) > 8
      wds <- if(manysample) 7 else 3

      # Get unique taxa levels
      taxa.levels <-
        otu_tax_agg_merged %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()
      taxa.levels <- mStat_resolve_selected_features(
        feature.level = feature.level,
        features.plot = current_features_plot,
        taxa.levels = taxa.levels,
        fallback_n = length(taxa.levels)
      )

      # Count number of samples
      n_subjects <- length(unique(otu_tax_agg_merged[["sample"]]))

      global_y_label <- mStat_get_taxa_value_ylabel(feature.dat.type, transform)

      # Generate individual plots for each taxon
      plot_list <- lapply(taxa.levels, function(tax) {

        sub_otu_tax_agg_merged <- otu_tax_agg_merged %>% filter(!!sym(feature.level) == tax)

        # Remove outliers from the data before plotting
        sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
          dplyr::group_by(!!sym(resolved_group_var)) %>%
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
            aes(color = !!sym(resolved_group_var))
          ) +
          scale_fill_manual(values = col) +
          scale_color_manual(values = col) +
          labs(
            x = group.var,
            y = global_y_label,
            title = tax
          ) +
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
        if(!is.null(group.var) && length(unique(stats::na.omit(data.obj$meta.dat[[group.var]]))) > 2){
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
        pdf_name <- mStat_append_pdf_group_suffixes(
          pdf_name = pdf_name,
          group.var = group.var,
          strata.var = strata.var
        )
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }
        pdf_name <- paste0(pdf_name,"_", feature.level, ".pdf")
        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Print each ggplot object to a new PDF page
        for (plot_obj in plot_list) {
          print(plot_obj)
        }
        # Close the PDF device
        dev.off()
      }

      names(plot_list) <- taxa.levels
      return(plot_list)
    })

    names(plot_list_all) <- feature.level
    return(plot_list_all)
  }
