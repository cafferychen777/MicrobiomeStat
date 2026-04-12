#' Generate Taxa Boxplots for Single Time Point
#'
#' Creates boxplots showing taxa abundance distributions at a single time point.
#' Supports grouping, stratification, and various transformations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param t.level Character string specifying the time level to subset data to.
#'   Default NULL uses all data.
#' @param features.plot Character vector of specific feature IDs to plot.
#'   If NULL, features are selected based on top.k.plot and top.k.func.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#' @param transform Transformation to apply: "identity", "sqrt", or "log".
#' @param point.alpha Numeric value (0-1) for jitter point transparency. Default 0.6.
#'
#' @return A list of ggplot objects for each taxonomic level.
#'
#' @examples
#' \dontrun{
#' # Generate the boxplot pair
#' data(ecam.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = "diet",
#'   strata.var = NULL,
#'   feature.level = c("Phylum"),
#'   features.plot = sample(unique(ecam.obj$feature.ann[,"Phylum"]),3),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   point.alpha = 0.6,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Phylum"),
#'   features.plot = sample(unique(ecam.obj$feature.ann[,"Phylum"]),3),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   point.alpha = 0.6,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = NULL,
#'   strata.var = NULL,
#'   feature.level = c("Order", "Phylum", "Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   point.alpha = 0.6,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   base.size = 12,
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
generate_taxa_boxplot_single <-
  function(data.obj,
           time.var = NULL,
           t.level = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           transform = c("sqrt", "identity", "log"),
           prev.filter = 0.01,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           point.alpha = 0.6,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    # Match and validate the feature data type
    feature.dat.type <- match.arg(feature.dat.type)

    # Match and validate the transformation method
    transform <- match.arg(transform)

    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Validate input variables
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    # Create a subject identifier from sample names (rownames)
    data.obj$meta.dat$subject.id <- rownames(data.obj$meta.dat)
    subject.var <- "subject.id"

    context <- mStat_prepare_taxa_single_context(
      data.obj = data.obj,
      time.var = time.var,
      t.level = t.level,
      group.var = group.var,
      strata.var = strata.var,
      subject.var = subject.var
    )
    data.obj <- context$data.obj

    # Extract metadata
    meta_tab <-
      select_meta_vars(data.obj$meta.dat, subject.var, group.var, time.var, strata.var)
    placeholder_group <- mStat_ensure_group_placeholder(
      meta_tab,
      group.var = group.var,
      value = "ALL",
      column_name = "ALL"
    )
    meta_tab <- placeholder_group$df
    resolved_group_var <- placeholder_group$group.var

    # Get color palette
    col <- mStat_get_palette(palette)

    # Get plot theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

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

      # Define aesthetic mapping for the plot
      aes_function <- aes(
        x = !!sym(feature.level),
        y = value,
        fill = !!sym(resolved_group_var)
      )

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
            resolved_group_var,
            strata.var,
            "value"
          )
        ))

      otu_tax_agg_merged <- mStat_transform_taxa_long_values(
        long.df = otu_tax_agg_merged,
        feature.level = feature.level,
        value_col = "value",
        feature.dat.type = feature.dat.type,
        transform = transform
      )

      # Get unique taxa levels
      taxa.levels <-
        otu_tax_agg_merged %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # Count number of subjects
      n_subjects <-
        length(unique(stats::na.omit(otu_tax_agg_merged[[subject.var]])))

      sub_otu_tax_agg_merged <- otu_tax_agg_merged

      # Select features to plot
      current_features_plot <- mStat_resolve_selected_features(
        feature.level = feature.level,
        features.plot = current_features_plot,
        taxa.levels = taxa.levels,
        fallback_n = 6
      )

      # Get number of group levels if group variable is specified
      if (!is.null(group.var)){
        group.levels <- sub_otu_tax_agg_merged %>%
          select(!!sym(resolved_group_var)) %>%
          pull() %>%
          stats::na.omit() %>%
          unique() %>%
          length()
      }

      # Create box plot
      boxplot <-
        ggplot(sub_otu_tax_agg_merged %>% filter(!!sym(feature.level) %in% current_features_plot),
               aes_function) +
        # Add error bars
        stat_boxplot(
          geom = "errorbar",
          position = position_dodge(width = 0.8),
          width = 0.2
        ) +
        # Add boxplots without outliers
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.2,
          outlier.shape = NA,  # Hide outliers since we show all points
          alpha = 0.8  # Make boxplot slightly transparent so points are visible
        ) +
        # Draw jitter points last so they appear on top of boxplots
        geom_jitter(
          aes(fill = !!sym(resolved_group_var)),  # Need fill aesthetic for dodging
          position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
          alpha = point.alpha,
          size = 2,
          color = "gray30"
        ) +
        scale_fill_manual(values = col) +
        labs(y = mStat_get_taxa_value_ylabel(feature.dat.type, transform)) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          strip.text.x = element_text(size = base.size, color = "black"),
          strip.text.y = element_text(size = base.size, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = (base.size -
                                                                                   2)),
          axis.text.y = element_text(color = "black", size = (base.size -
                                                                2)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      # Add facets if strata variable is specified.
      if (!is.null(strata.var)) {
        boxplot <- boxplot + ggh4x::facet_nested_wrap(
          as.formula(paste('~', strata.var)),
          scales = "free_x"
        )
      }

      # Modify y-axis scale based on transformation
      if (feature.dat.type != "other") {
        if (transform == "sqrt") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  a ^ b, list(a = i, b = 2)
                )))
          )
        } else if (transform == "log") {
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

    # Save plots as PDF if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_boxplot_single",
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
      pdf_name <- paste0(pdf_name, "_", feature.level, ".pdf")
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Print each ggplot object to a new PDF page
      for (plot_obj in plot_list) {
        print(plot_obj)
      }
      # Close the PDF device
      dev.off()
    }

    # Name the plot list with feature levels
    names(plot_list) <- feature.level

    # Return the list of plots
    return(plot_list)
  }
