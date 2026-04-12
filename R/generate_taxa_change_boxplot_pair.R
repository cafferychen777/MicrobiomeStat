#' Generate Taxa Change Boxplots for Paired Data
#'
#' Creates boxplots showing changes in taxa abundance between paired time points.
#' Supports multiple change metrics and grouping/stratification.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param features.plot Character vector of specific feature IDs to plot.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#'
#' @return A list of ggplot objects showing abundance changes.
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
           feature.level,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
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
    # Validate input data
    data.obj <- mStat_validate_data(data.obj)

    # Match and validate the feature data type
    feature.dat.type <- match.arg(feature.dat.type)

    # Validate input variables
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

    # Extract metadata
    meta_tab <- mStat_prepare_meta_tab(
      meta.dat = data.obj$meta.dat,
      vars = list(subject.var, time.var, group.var, strata.var)
    )

    # Get plot theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Get color palette
    col <- mStat_get_palette(palette)

    ylab_label <- mStat_get_taxa_change_ylabel(
      feature.dat.type = feature.dat.type,
      feature.change.func = feature.change.func
    )

    # Adjust filtering parameters if necessary
    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    # Normalize count data if necessary.
    analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

    # Generate plots for each taxonomic level
    plot_list <- lapply(feature.level, function(feature.level) {

      # Aggregate data by taxonomy if necessary
      otu_tax_agg <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter)

      current_features_plot <- mStat_resolve_selected_features(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        features.plot = features.plot,
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

      pair_change <- mStat_prepare_taxa_pair_change_data(
        long.df = otu_tax_agg_merged,
        feature.level = feature.level,
        subject.var = subject.var,
        time.var = time.var,
        change.base = change.base,
        feature.change.func = feature.change.func,
        context = "taxa change boxplotting"
      )
      combined_data <- pair_change$combined_data
      change.after <- pair_change$change.after

      # Add metadata for the change time point
      combined_data <- mStat_attach_pair_metadata(
        df = combined_data,
        meta_tab = meta_tab,
        subject.var = subject.var,
        time.var = time.var,
        mode = "followup_time",
        change.after = change.after
      )

      # Get unique taxa levels
      taxa.levels <-
        combined_data %>% select(all_of(c(feature.level))) %>% dplyr::distinct() %>% dplyr::pull()

      placeholder_group <- mStat_ensure_group_placeholder(
        combined_data,
        group.var = group.var,
        value = "ALL",
        column_name = "ALL"
      )
      combined_data <- placeholder_group$df
      resolved_group_var <- placeholder_group$group.var

      # Select features to plot
      taxa.levels <- mStat_resolve_selected_features(
        feature.level = feature.level,
        features.plot = current_features_plot,
        taxa.levels = taxa.levels,
        fallback_n = 4
      )

      # Create the boxplot
      boxplot <-
        ggplot(
          combined_data %>% filter(!!sym(feature.level) %in% taxa.levels),
          aes(
            x = !!sym(resolved_group_var),
            y = value_diff,
            fill = !!sym(resolved_group_var)
          )
        ) +
        geom_jitter(width = 0.1,
                    alpha = 0.3,
                    size = 1.5) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
        ) +
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

      # Add facets if strata variable is provided
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

      # Adjust theme if there's only one group
      if (is.null(group.var)) {
        boxplot <- boxplot +
          theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank()
          )
      }

      return(boxplot)

    })

    # Name the plot list with feature levels
    names(plot_list) <- feature.level

    # Save the plots as a PDF file if specified
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
      pdf_name <- mStat_append_pdf_group_suffixes(
        pdf_name = pdf_name,
        group.var = group.var,
        strata.var = strata.var
      )
      if (!is.null(file.ann)) {
        pdf_name <- paste0(pdf_name, "_", file.ann)
      }
      pdf_name <- paste0(pdf_name, ".pdf")
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Print each ggplot object to a new PDF page
      for (plot_obj in plot_list) {
        print(plot_obj)
      }
      # Close the PDF device
      dev.off()
    }

    # Return the boxplot for display
    return(plot_list)
  }
