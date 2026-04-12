#' @title Generate Longitudinal Spaghetti Plots of Taxonomic Composition
#'
#' @description Generates spaghetti plots showing taxonomic composition trends over time,
#' with individual subject trajectories and group-level mean lines.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param ... Additional parameters to be passed.
#'
#' @return A ggplot object containing the longitudinal line plot.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data("subset_T2D.obj")
#' generate_taxa_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "sample_body_site",
#'   feature.level = c("Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "npg",
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' }
#' @export
generate_taxa_spaghettiplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
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

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Check if input variables are of the correct type
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

    context <- mStat_prepare_taxa_long_context(
      data.obj = data.obj,
      subject.var = subject.var,
      time.var = time.var,
      group.var = group.var,
      strata.var = strata.var,
      t0.level = t0.level,
      ts.levels = ts.levels
    )
    data.obj <- context$data.obj
    meta_tab <- context$meta_tab

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)
    text_sizes <- mStat_get_spaghettiplot_text_sizes(base.size)

    # Adjust filters if necessary based on input parameters
    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    # Normalize count data if necessary.
    analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

    # Generate plots for each feature level
    plot_list <- lapply(feature.level, function(feature.level){

      # Aggregate data by taxonomy if necessary
      otu_tax_agg <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter)

      selected_features <- mStat_resolve_selected_features(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        features.plot = features.plot,
        top.k.plot = top.k.plot,
        top.k.func = top.k.func
      )

      df <- mStat_prepare_taxa_long_data(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        value_col = "count",
        meta.dat = meta_tab
      )

      placeholder_group <- mStat_ensure_group_placeholder(df, group.var = group.var)
      df <- placeholder_group$df
      resolved_group_var <- placeholder_group$group.var

      # Calculate mean counts
      mean_group_vars <- c(time.var, resolved_group_var, strata.var)
      mean_group_vars <- mean_group_vars[!sapply(mean_group_vars, is.null)]
      mean_df <- mStat_summarize_mean_by_groups(
        long.df = df,
        feature.level = feature.level,
        group_vars = mean_group_vars,
        value_col = "count",
        mean_col = "mean_count"
      )

      # Get unique taxa levels
      taxa.levels <-
        df %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # Select features to plot
      selected_features <- mStat_resolve_selected_features(
        feature.level = feature.level,
        features.plot = selected_features,
        taxa.levels = taxa.levels,
        fallback_n = 4
      )

      # Subset data for selected features
      sub_df <- df %>% filter(!!sym(feature.level) %in% selected_features)
      sub_df.mean <- mean_df %>% filter(!!sym(feature.level) %in% selected_features)

      # Get number of strata levels if strata variable is specified
      if (!is.null(strata.var)){
        strata.levels <- df %>% select(!!sym(strata.var)) %>% pull() %>% unique() %>% length()
      }

      # Create the spaghetti plot
      lineplot <- ggplot() +
        geom_line(
          data = sub_df,
          aes(
            x = .data[[time.var]],
            y = .data[["count"]],
            group = .data[[subject.var]],
            color = .data[[resolved_group_var]]
          ),
          alpha = 0.5
        ) +
        geom_line(
          data = sub_df.mean,
          aes(
            x = .data[[time.var]],
            y = .data[["mean_count"]],
            group = .data[[resolved_group_var]],
            color = .data[[resolved_group_var]]
          ),
          size = 2
        ) +
        geom_point(
          data = sub_df.mean,
          aes(
            x = .data[[time.var]],
            y = .data[["mean_count"]],
            group = .data[[resolved_group_var]],
            color = .data[[resolved_group_var]]
          ),
          size = 3
        ) +
        scale_color_manual(values = col) +
        labs(
          x = time.var,
          y = mStat_get_taxa_value_ylabel(feature.dat.type),
          color = if (is.null(group.var)) NULL else group.var
        ) +
        {
          if (!is.null(strata.var)) {
            ggh4x::facet_nested_wrap(as.formula(paste('~', feature.level,"+",strata.var)), scales = "free_y", ncol = strata.levels*2)  # Use facet_wrap with strata.var as the faceting variable
          } else {
            ggh4x::facet_nested_wrap(as.formula(paste('~',feature.level)), scales = "free_y", ncol = 3)
          }
        } +
        theme_to_use +
        theme(
          plot.title = element_text(
            hjust = 0.5,
            size = text_sizes$title
          ),
          panel.spacing.x = unit(0, "cm"),
          axis.title.x = element_text(size = text_sizes$axis.title),
          axis.title.y = element_text(size = text_sizes$axis.title),
          axis.text.x = element_text(size = text_sizes$axis.text),
          axis.text.y = element_text(size = text_sizes$axis.text),
          legend.title = element_text(size = text_sizes$legend.title),
          legend.text = element_text(size = text_sizes$legend.text)
        )

      # Remove legend for single group case
      if (is.null(group.var)) {
        lineplot <- lineplot + theme(legend.position = "none")
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_spaghettiplot_long",
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
        # Print the plot to the PDF
        print(lineplot)
        # Close the PDF device
        dev.off()
      }

      return(lineplot)
    })

    # Name the plots in the list
    names(plot_list) <- feature.level

    # Return the list of plots
    return(plot_list)
  }
