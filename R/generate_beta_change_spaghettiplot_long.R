#' Generate Beta Diversity Change Spaghetti Plot for Longitudinal Data
#'
#' Creates spaghetti plots showing beta diversity changes from baseline over time.
#' Individual subject trajectories and group means are visualized.
#'
#' @note Subjects without a baseline time point (t0.level) will be excluded.
#'
#' @name generate_beta_change_spaghettiplot_long
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A spaghetti plot displaying the beta diversity change over time, stratified by the specified grouping and/or strata variables (if provided). The plot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_beta_change_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' generate_beta_change_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = NULL,
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_beta_change_spaghettiplot_long <-
  function(data.obj,
           dist.obj = NULL,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           dist.name = c("BC"),
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    # Exit the function if no distance metric is specified
    if (is.null(dist.name)){
      return()
    }

    # Input validation to ensure correct variable types
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

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Set the theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)
    text_sizes <- mStat_get_spaghettiplot_text_sizes(base.size, variant = "wide_panel")

    # Inform users about the requirement for baseline time points
    message("Note: This function requires subjects to have a baseline time point (specified by `t0.level`) to calculate the change in beta diversity. Subjects without a baseline time point will be excluded from the visualization.")

    # Generate plots for each specified distance metric
    plot_list <- lapply(dist.name,function(dist.name){
      # Calculate beta diversity if not provided
      if (is.null(dist.obj)) {
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
        dist.obj <-
          mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        # Adjust distances for covariates if specified
        if (!is.null(adj.vars)){
          dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
        }
      } else {
        prepared_context <- mStat_prepare_precomputed_beta_context(
          dist.obj = dist.obj,
          dist.name = dist.name,
          data.obj = data.obj,
          time.var = time.var,
          t0.level = t0.level,
          ts.levels = ts.levels,
          process_time = TRUE
        )
        data.obj <- prepared_context$data.obj
        dist.obj <- prepared_context$dist.obj
        meta_tab <- mStat_extract_dist_metadata(
          dist.obj = dist.obj,
          dist.name = dist.name,
          vars = c(subject.var, time.var, group.var, strata.var),
          data.obj = data.obj
        )
      }

      # Check if the specified distance metric exists in the distance object.
      if (is.null(dist.obj[[dist.name]])) {
        stop("dist.obj does not contain ", dist.name, call. = FALSE)
      }

      meta_tab <- mStat_meta_to_tibble(meta_tab, sample_col = "sample")

      resolved_time <- mStat_resolve_followup_timepoints(
        values = meta_tab[[time.var]],
        time.var = time.var,
        t0.level = t0.level,
        ts.levels = ts.levels,
        context = "beta change spaghetti plot"
      )
      change.base <- resolved_time$t0.level

      # Calculate pairwise distances between baseline and follow-up time points for each subject
      long.df <- mStat_prepare_beta_change_long_data(
        dist.matrix = dist.obj[[dist.name]],
        meta.dat = meta_tab,
        subject.var = subject.var,
        time.var = time.var,
        change.base = change.base
      )

      # Add group and strata information to the long-format data
      long.df <- mStat_attach_change_metadata(
        change.df = long.df,
        meta.dat = meta_tab,
        by = c(subject.var, time.var),
        vars = c(group.var, strata.var)
      )

      placeholder_group <- mStat_ensure_group_placeholder(long.df, group.var = group.var)
      long.df <- placeholder_group$df
      resolved_group_var <- placeholder_group$group.var

      # Calculate mean distances for each time point and group (and strata if applicable)
      if (is.null(strata.var)) {
        long.df.mean <- long.df %>%
          dplyr::group_by(!!sym(time.var), !!sym(resolved_group_var)) %>%
          dplyr::summarize(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop")
      } else {
        long.df.mean <- long.df %>%
          dplyr::group_by(!!sym(time.var), !!sym(resolved_group_var), !!sym(strata.var)) %>%
          dplyr::summarize(mean_distance = mean(distance, na.rm = TRUE), .groups = "drop")
      }

      long.df <- long.df %>% dplyr::arrange(!!sym(subject.var), !!sym(time.var))

      # Set y-axis label based on whether distances are adjusted for covariates
      if (is.null(adj.vars)) {
        y_label <- paste(dist.name, "Distance from Baseline")
      } else {
        y_label <- paste(dist.name, "Distance from Baseline (adjusted by:", paste(adj.vars, collapse = ", "), ")")
      }

      # Create the spaghetti plot
      p <- ggplot() +
        geom_point(
          data = long.df,
          aes(
            x = .data[[time.var]],
            y = .data[["distance"]],
            group = .data[[subject.var]],
            color = .data[[resolved_group_var]]
          ),
          alpha = 0.3,
          size = 3
        ) +
        geom_line(
          data = long.df.mean,
          aes(
            x = .data[[time.var]],
            y = .data[["mean_distance"]],
            group = .data[[resolved_group_var]],
            color = .data[[resolved_group_var]]
          ),
          size = 2
        ) +
        geom_point(
          data = long.df.mean,
          aes(
            x = .data[[time.var]],
            y = .data[["mean_distance"]],
            group = .data[[resolved_group_var]],
            color = .data[[resolved_group_var]]
          ),
          size = 5
        ) +
        labs(x = time.var, y = y_label, color = group.var) +
        scale_color_manual(
          values = col
        ) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = text_sizes$axis.title),
          axis.title.y = element_text(size = text_sizes$axis.title),
          axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = text_sizes$axis.text),
          axis.text.y = element_text(size = text_sizes$axis.text),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = text_sizes$legend.text),
          legend.title = ggplot2::element_text(size = text_sizes$legend.title),
          legend.key.size = unit(10, "mm"),
          legend.key.spacing = unit(2, "mm")
        )

      # Remove legend if there's only one group
      if (is.null(group.var)){
        p <- p + theme(legend.position = "none")
      }

      # Add faceting by strata if a strata variable is specified
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(
          cols = vars(!!sym(strata.var)),
          scales = "free",
          space = "free"
        )
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0("beta_change_spaghettiplot_",
                           dist.name,
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var)

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

        ggsave(
          filename = pdf_name,
          plot = p,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      }

      return(p)
    })

    # Name the plots in the list according to the distance metrics
    names(plot_list) <- dist.name
    return(plot_list)
  }
