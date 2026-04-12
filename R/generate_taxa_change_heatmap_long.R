#' @title Generate Taxonomic Change Heatmap for Longitudinal Data
#'
#' @description Creates heatmaps showing taxa abundance changes from baseline across
#' multiple time points. Uses pheatmap with hierarchical clustering and group annotations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param features.plot Character vector of specific feature IDs to plot.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#' @param cluster.rows Logical, whether to cluster rows. Default TRUE.
#' @param cluster.cols Logical, whether to cluster columns. Default FALSE.
#'
#' @return A list of ggplot heatmap objects.
#'
#' @seealso \code{\link{pheatmap}}, \code{\link{mStat_normalize_data}}
#'
#' @examples
#' \dontrun{
#' library(pheatmap)
#' data(ecam.obj)
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   feature.level = c("Family","Class"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   palette = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   feature.level = c("Family","Class"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   cluster.rows = FALSE,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   palette = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_change_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Phylum"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 50,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Phylum"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   cluster.rows = FALSE,
#'   top.k.plot = 50,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#'
#' @return An object of class pheatmap, the generated heatmap plot
#' @export
#' @importFrom dplyr distinct pull
#' @seealso \code{\link{pheatmap}}
generate_taxa_change_heatmap_long <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              t0.level = NULL,
                                              ts.levels = NULL,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              feature.level,
                                              feature.change.func = "relative change",
                                              feature.dat.type = c("count", "proportion", "other"),
                                              features.plot = NULL,
                                              top.k.plot = NULL,
                                              top.k.func = NULL,
                                              prev.filter = 0.01,
                                              abund.filter = 0.01,
                                              base.size = 10,
                                              palette = NULL,
                                              cluster.rows = NULL,
                                              cluster.cols = NULL,
                                              pdf = TRUE,
                                              file.ann = NULL,
                                              pdf.wid = 11,
                                              pdf.hei = 8.5,
                                              ...) {

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Validate the input data object
  data.obj <- mStat_validate_data(data.obj)

  # Check if the input variables are of the correct type
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

  # Process the time variable in the data object
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  # Capture original factor levels before any data extraction
  fl <- mStat_capture_factor_levels(data.obj, group.var, strata.var)
  data.obj <- fl$data.obj

  # Extract relevant columns from the metadata
  meta_tab <- data.obj$meta.dat %>%
    as.data.frame() %>%
    select(all_of(c(subject.var,group.var,time.var,strata.var)))

  placeholder_group <- mStat_ensure_group_placeholder(
    meta_tab,
    group.var = group.var,
    value = "ALL",
    column_name = "ALL"
  )
  meta_tab <- placeholder_group$df
  resolved_group_var <- placeholder_group$group.var
  has_group <- !is.null(group.var)

  # Set default values for clustering if not provided
  cluster.cols <- mStat_resolve_optional_flag(cluster.cols, FALSE, "cluster.cols")
  cluster.rows <- mStat_resolve_optional_flag(cluster.rows, TRUE, "cluster.rows")

  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(resolved_group_var) := interaction(!!sym(resolved_group_var), !!sym(strata.var), sep = .STRATA_SEP))
  }

  meta_tab_with_sample <- mStat_meta_with_sample(meta_tab)

  # Adjust filtering parameters based on input conditions
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  # Normalize count data if necessary.
  data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

  # Generate plots for each taxonomic level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

    # Select top k features if specified
    selected_features <- mStat_resolve_selected_features(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      features.plot = features.plot,
      top.k.plot = top.k.plot,
      top.k.func = top.k.func
    )

    df_mean_value <- mStat_prepare_taxa_long_data(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      value_col = "value",
      meta.dat = meta_tab_with_sample
    ) %>%
      mStat_summarize_mean_by_groups(
        feature.level = feature.level,
        group_vars = c(resolved_group_var, time.var),
        value_col = "value",
        mean_col = "mean_value"
      )

    resolved_time <- mStat_resolve_followup_timepoints(
      values = meta_tab[[time.var]],
      time.var = time.var,
      t0.level = t0.level,
      ts.levels = ts.levels,
      context = "taxa change heatmap"
    )
    t0.level <- resolved_time$t0.level
    ts.levels <- resolved_time$ts.levels

    wide_data <- mStat_prepare_change_heatmap_long_matrix(
      summary.df = df_mean_value,
      feature.level = feature.level,
      group.var = resolved_group_var,
      time.var = time.var,
      baseline_time = t0.level,
      followup_times = ts.levels,
      feature.change.func = feature.change.func,
      value_col = "mean_value"
    )

    # Save original column names
    original_colnames <- colnames(wide_data)

    # Remove columns where all values are NA
    wide_data <-
      wide_data[, colSums(is.na(wide_data)) != nrow(wide_data), drop = FALSE]

    if (!all(is.finite(wide_data))) {
      wide_data[!is.finite(wide_data)] <- NA_real_
    }

    # Save new column names
    new_colnames <- colnames(wide_data)

    # Find out removed columns
    removed_colnames <- setdiff(original_colnames, new_colnames)

    # Notify user about removed columns
    if (length(removed_colnames) > 0) {
      message(
        "The following combinations were all NA, therefore they have been removed: ",
        paste(removed_colnames, collapse = ", ")
      )
    } else {
      message("No columns have been removed.")
    }

    # Prepare annotation for columns
    annotation_col <- meta_tab %>%
      select(!!sym(time.var), !!sym(resolved_group_var)) %>%
      filter(!!sym(time.var) != t0.level) %>%
      tibble::as_tibble() %>%
      dplyr::distinct() %>%
      dplyr::mutate(group_time = paste(!!sym(resolved_group_var), !!sym(time.var), sep = "_")) %>%
      filter(group_time %in% new_colnames) %>%
      tibble::column_to_rownames("group_time")

    # Sort annotation by group and time
    annotation_col_sorted <-
      annotation_col[order(annotation_col[[resolved_group_var]], annotation_col[[time.var]]),]

    # Handle strata variable if present
    if (!is.null(strata.var)) {
      group_var_for_restore <- if (has_group) group.var else resolved_group_var
      annotation_col_sorted <- annotation_col_sorted %>%
        tidyr::separate(!!sym(resolved_group_var),
                 into = c(group_var_for_restore, strata.var),
                 sep = .STRATA_SEP)
      annotation_col_sorted <- mStat_restore_factor_levels(
        annotation_col_sorted, fl$levels, group_var_for_restore, strata.var)

      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]], annotation_col_sorted[[group_var_for_restore]], annotation_col_sorted[[time.var]]), ]

    }

    # Handle the case when no grouping variable is provided
    if (!has_group) {
      annotation_col_sorted <-
        annotation_col_sorted %>% select(all_of(c(time.var)))
    }

    # Sort wide data according to annotation
    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted), drop = FALSE]

    # Filter features if specified
    if (!is.null(selected_features)) {
      wide_data_sorted <-
        wide_data_sorted[rownames(wide_data_sorted) %in% selected_features,, drop = FALSE]
    }

    all_na_rows <- rowSums(!is.na(wide_data_sorted)) == 0
    if (any(all_na_rows)) {
      wide_data_sorted <- wide_data_sorted[!all_na_rows, , drop = FALSE]
    }

    # Calculate gaps for heatmap
    group_var_for_restore <- if (has_group) group.var else resolved_group_var
    if (has_group) {
      gaps <-
        cumsum(table(annotation_col_sorted[[group_var_for_restore]]))[-length(table(annotation_col_sorted[[group_var_for_restore]]))]
      if (length(gaps) == 0) {
        gaps <- NULL
      }
    } else {
      gaps <- NULL
    }

    if (!is.null(strata.var)){
      gaps <-
        cumsum(table(annotation_col_sorted[[strata.var]]))[-length(table(annotation_col_sorted[[strata.var]]))]
      if (length(gaps) == 0) {
        gaps <- NULL
      }
    }

    if (!is.null(gaps) && ncol(wide_data_sorted) <= 1) {
      gaps <- NULL
    }

    # Set up color scheme for heatmap
    n_colors <- 100
    col <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")

    # Find the maximum absolute value in the data
    matrix_values <- as.numeric(wide_data_sorted)
    finite_values <- matrix_values[is.finite(matrix_values)]
    if (length(finite_values) == 0) {
      max_abs_val <- 1
    } else {
      max_abs_val <- max(abs(finite_values))
      if (!is.finite(max_abs_val) || max_abs_val <= 0) {
        max_abs_val <- 1
      }
    }

    # Keep white centered at zero with a fixed midpoint
    zero_pos <- floor(n_colors / 2)

    # Create color vectors
    my_col <-
      c(
        colorRampPalette(col[1:3])(zero_pos),
        colorRampPalette(col[3:5])(n_colors - zero_pos + 1)
      )

    # Create break points for color scale
    break_points <-
      seq(-max_abs_val, max_abs_val, length.out = length(my_col) + 1)

    # Get color palette
    color_vector <- mStat_get_palette(palette)

    # Set up annotation colors
    if (!is.null(strata.var) && !is.null(group.var)){
      # Get unique levels for group and strata variables
      # Check if variables are already factors to preserve their order
      if (is.factor(annotation_col_sorted[[group.var]])) {
        group_levels <- levels(annotation_col_sorted[[group.var]])
      } else {
        group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      }
      
      if (is.factor(annotation_col_sorted[[strata.var]])) {
        strata_levels <- levels(annotation_col_sorted[[strata.var]])
      } else {
        strata_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      }

      # Assign colors to group and strata levels
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

      # Create a list of annotation colors
      annotation_colors_list <- setNames(
        list(group_colors, strata_colors),
        c(group.var, strata.var)
      )
    } else if (!is.null(group.var)){
      # Get unique levels for group variable
      # Check if variable is already a factor to preserve its order
      if (is.factor(annotation_col_sorted[[group.var]])) {
        group_levels <- levels(annotation_col_sorted[[group.var]])
      } else {
        group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      }
      
      # Assign colors to group levels
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      
      # Create annotation color list
      annotation_colors_list <- setNames(
        list(group_colors),
        c(group.var)
      )
    } else {
      annotation_colors_list <- NULL
    }

    heatmap_gaps <- if (!is.null(gaps) && length(gaps) > 0 && ncol(wide_data_sorted) > 1) {
      gaps
    } else {
      NULL
    }

    # Plot stacked heatmap
    heatmap_plot <- pheatmap::pheatmap(
      mat = wide_data_sorted[order(rowMeans(abs(wide_data_sorted), na.rm = TRUE), decreasing = TRUE), , drop = FALSE],
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = heatmap_gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col,
      breaks = break_points,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_long",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "group_",
        group.var,
        "_",
        "strata_",
        strata.var,
        "_",
        "taxa_",
        feature.level,
        "_",
        "prev_filter_",
        prev.filter,
        "_",
        "abund_filter_",
        abund.filter
      )

      if (!is.null(file.ann)) {
        pdf_name <- paste0(pdf_name, "_", file.ann)
      }

      pdf_name <- paste0(pdf_name, ".pdf")

      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_heatmap_plot
      )
    }
    return(gg_heatmap_plot)
  })

  names(plot_list) <- feature.level

  # Return the heatmap plot for display
  return(plot_list)
}
