#' @title Generate Taxa Heatmap for Longitudinal Data
#'
#' @description Generates hierarchical clustered heatmaps for longitudinal microbiome data
#' using the pheatmap package, with options for group and strata annotations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to compute the top k taxa. Default is NULL (uses mean).
#' @param cluster.rows Logical indicating if rows should be clustered. Default is TRUE.
#' @param cluster.cols Logical indicating if columns should be clustered. Default is FALSE.
#' @param ... Additional parameters to be passed to the pheatmap() function.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[2:18],
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Family","Phylum","Genus", "Class"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[2:18],
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Family","Phylum","Genus", "Class"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   cluster.rows = FALSE,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "sample_body_site",
#'   strata.var = "subject_gender",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 15
#' )
#' }
#' @export
#' @importFrom dplyr distinct pull
#' @seealso \code{\link{pheatmap}}
generate_taxa_heatmap_long <- function(data.obj,
                                       subject.var,
                                       time.var,
                                       t0.level,
                                       ts.levels,
                                       group.var = NULL,
                                       strata.var = NULL,
                                       feature.level,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       features.plot = NULL,
                                       top.k.plot = NULL,
                                       top.k.func = NULL,
                                       prev.filter = 0.01,
                                       abund.filter = 0.01,
                                       base.size = 10,
                                       palette = NULL,
                                       cluster.cols = NULL,
                                       cluster.rows = NULL,
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       ...) {
  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Validate the input data object
  data.obj <- mStat_validate_data(data.obj)

  # Input validation for key variables
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

  # Process time variable in the data object
  data.obj <-
    mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  # Capture original factor levels before any data extraction
  fl <- mStat_capture_factor_levels(data.obj, group.var, strata.var)
  data.obj <- fl$data.obj

  # Extract relevant metadata
  meta_tab <-
    data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
      subject.var, group.var, time.var, strata.var
    )))

  placeholder_group <- mStat_ensure_group_placeholder(
    meta_tab,
    group.var = group.var,
    value = "ALL",
    column_name = "ALL"
  )
  meta_tab <- placeholder_group$df
  resolved_group_var <- placeholder_group$group.var

  # If a strata variable is provided, create an interaction term with the group variable
  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(resolved_group_var) := interaction(!!sym(resolved_group_var), !!sym(strata.var), sep = .STRATA_SEP))
  }

  # Set default clustering options if not specified
  cluster.cols <- mStat_resolve_optional_flag(cluster.cols, FALSE, "cluster.cols")
  cluster.rows <- mStat_resolve_optional_flag(cluster.rows, TRUE, "cluster.rows")

  # Determine whether to apply abundance and prevalence filters
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  # Normalize count data if necessary.
  analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

  # Create a list to store plots for each taxonomic level
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

    # Convert counts to numeric
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    # Prepare the data for heatmap visualization
    otu_tab_norm <-
      otu_tax_agg_numeric %>%
      filter(!is.na(!!sym(feature.level))) %>%
      mStat_as_taxa_feature_matrix(
        feature.level = feature.level,
        feature_in_column = TRUE
      )

    wide_data <- mStat_prepare_taxa_long_data(
      feature.dat = otu_tax_agg_numeric,
      feature.level = feature.level,
      value_col = "value",
      meta.dat = meta_tab,
      sample_col = "sample"
    ) %>%
      dplyr::filter(!is.na(.data[[feature.level]])) %>%
      mStat_summarize_mean_by_groups(
        feature.level = feature.level,
        group_vars = c(resolved_group_var, time.var),
        value_col = "value",
        mean_col = "mean_value"
      ) %>%
      tidyr::unite("group_time", dplyr::all_of(c(resolved_group_var, time.var)), sep = "_") %>%
      tidyr::pivot_wider(names_from = "group_time", values_from = "mean_value") %>%
      mStat_as_taxa_feature_matrix(
        feature.level = feature.level,
        feature_in_column = TRUE
      )

    # Prepare column annotations for the heatmap
    annotation_col <- meta_tab %>%
      select(!!sym(time.var), !!sym(resolved_group_var)) %>%
      tibble::as_tibble() %>%
      dplyr::distinct() %>%
      dplyr::mutate(group_time = paste(!!sym(resolved_group_var), !!sym(time.var), sep = "_")) %>%
      tibble::column_to_rownames("group_time")

    # Sort the annotation columns
    annotation_col_sorted <-
      annotation_col[order(annotation_col[[resolved_group_var]], annotation_col[[time.var]]),]

    # Handle strata variable if present
    if (!is.null(strata.var)) {
      group_var_for_restore <- if (is.null(group.var)) resolved_group_var else group.var
      annotation_col_sorted <- annotation_col_sorted %>%
        tidyr::separate(!!sym(resolved_group_var),
                        into = c(group_var_for_restore, strata.var),
                        sep = .STRATA_SEP) %>%
        dplyr::select(!!sym(time.var), !!sym(group_var_for_restore), !!sym(strata.var))
      annotation_col_sorted <- mStat_restore_factor_levels(
        annotation_col_sorted, fl$levels, group_var_for_restore, strata.var)

      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]],
                                    annotation_col_sorted[[group_var_for_restore]],
                                    annotation_col_sorted[[time.var]]),]
    }

    # Sort the data according to the annotation
    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted), drop = FALSE]

    group_var_for_restore <- if (is.null(group.var)) resolved_group_var else group.var

    # Calculate gaps for the heatmap
    if (!is.null(group.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[group_var_for_restore]]))[-length(table(annotation_col_sorted[[group_var_for_restore]]))]
    } else {
      gaps <- NULL
    }

    if (!is.null(strata.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[strata.var]]))[-length(table(annotation_col_sorted[[strata.var]]))]
    }

    heatmap_gaps <- if (!is.null(gaps) && length(gaps) > 0 && ncol(wide_data_sorted) > 1) {
      gaps
    } else {
      NULL
    }

    # Define color palette for the heatmap
    col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
    my_col <- colorRampPalette(col)
    n_colors <- 100

    # Filter features if specified
    if (!is.null(current_features_plot)) {
      wide_data_sorted <-
        wide_data_sorted[rownames(wide_data_sorted) %in% current_features_plot, , drop = FALSE]
    }

    # Define colors for group and strata variables
    color_vector <- mStat_get_palette(palette)
    
    # Check if group.var is already a factor to preserve its order
    if (is.factor(annotation_col_sorted[[group_var_for_restore]])) {
      group_levels <- levels(annotation_col_sorted[[group_var_for_restore]])
    } else {
      group_levels <-
        annotation_col_sorted %>% dplyr::select(all_of(c(group_var_for_restore))) %>% distinct() %>% pull()
    }
    group_colors <-
      setNames(color_vector[1:length(group_levels)], group_levels)

    if (!is.null(strata.var)){
      # Check if strata.var is already a factor to preserve its order
      if (is.factor(annotation_col_sorted[[strata.var]])) {
        strata_levels <- levels(annotation_col_sorted[[strata.var]])
      } else {
        strata_levels <-
          annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      }
      strata_colors <-
        setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)
    }

    # Create annotation color list
    if (!is.null(strata.var)){
      annotation_colors_list <- setNames(list(group_colors, strata_colors),
                                         c(group_var_for_restore, strata.var))
    } else {
      annotation_colors_list <- setNames(list(group_colors),
                                         c(group_var_for_restore))
    }

    if (is.null(group.var)){
      annotation_col_sorted <- annotation_col_sorted %>% select(-all_of("ALL"))
      annotation_colors_list <- NULL
    }

    # Generate heatmap for average values
    heatmap_plot_average <- pheatmap::pheatmap(
      mat = wide_data_sorted[order(rowMeans(wide_data_sorted, na.rm = TRUE), decreasing = TRUE), , drop = FALSE],
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = heatmap_gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col(n_colors),
      ...
    )

    gg_heatmap_plot_average <- as.ggplot(heatmap_plot_average)

    # Filter features for individual heatmap if specified
    if (!is.null(current_features_plot)) {
      otu_tab_norm <-
        otu_tab_norm[rownames(otu_tab_norm) %in% current_features_plot, , drop = FALSE]
    }

    # Prepare metadata for individual heatmap
    if (is.null(group.var)){
      meta_tab <- meta_tab %>% select(-all_of("ALL"))
      annotation_colors_list <- NULL
    } else {
      meta_tab <-
        data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
          subject.var, group.var, time.var, strata.var
        )))
    }

    # Generate heatmap for individual samples
    heatmap_plot_indiv <- pheatmap::pheatmap(
      mat = otu_tab_norm[order(rowMeans(otu_tab_norm, na.rm = TRUE), decreasing = TRUE), , drop = FALSE],
      annotation_col = meta_tab,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      gaps_col = heatmap_gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col(n_colors),
      ...
    )

    gg_heatmap_plot_indiv <- as.ggplot(heatmap_plot_indiv)

    # Save the average heatmap as a PDF file if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_long_average",
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

      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_heatmap_plot_average
      )
    }

    # Save the individual heatmap as a PDF file if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_long_indiv",
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

      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_heatmap_plot_indiv
      )
    }

    # Combine individual and average heatmaps into a list
    sub.plot_list <- list(gg_heatmap_plot_indiv, gg_heatmap_plot_average)
    names(sub.plot_list) <- c("indiv","average")

    return(sub.plot_list)
  })

  # Name the plots in the list by feature level
  names(plot_list) <- feature.level

  # Return the list of heatmap plots
  return(plot_list)
}
