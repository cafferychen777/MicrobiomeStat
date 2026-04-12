#' @title Generate Taxa Change Heatmap for Paired Data
#'
#' @description Creates heatmaps showing taxa abundance changes between paired time points.
#' Generates both individual and group-averaged views using pheatmap.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param features.plot Character vector of specific feature IDs to plot.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#' @param cluster.rows Logical, whether to cluster rows. Default TRUE.
#' @param cluster.cols Logical, whether to cluster columns. Default NULL.
#'
#' @return A list containing individual and average heatmap objects.
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(pheatmap)
#' data(peerj32.obj)
#' generate_taxa_change_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = NULL,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_change_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = FALSE,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_pairs.obj)
#' generate_taxa_change_heatmap_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = NULL,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_change_heatmap_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = FALSE,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_change_heatmap_pair <- function(data.obj,
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
                                              base.size = 10,
                                              palette = NULL,
                                              cluster.rows = NULL,
                                              cluster.cols = NULL,
                                              pdf = TRUE,
                                              file.ann = NULL,
                                              pdf.wid = 11,
                                              pdf.hei = 8.5,
                                              ...) {
  # Validate the input data object
  data.obj <- mStat_validate_data(data.obj)

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Extract relevant metadata
  meta_tab <- mStat_prepare_meta_tab(
    meta.dat = data.obj$meta.dat,
    vars = list(time.var, group.var, strata.var, subject.var)
  )

  default_cluster_cols <- if (!is.null(group.var)) {
    !(is.numeric(meta_tab %>% select(all_of(group.var)) %>% dplyr::pull()) &&
        !is.integer(meta_tab %>% select(all_of(group.var)) %>% dplyr::pull()))
  } else {
    TRUE
  }
  cluster.cols <- mStat_resolve_optional_flag(cluster.cols, default_cluster_cols, "cluster.cols")
  cluster.rows <- mStat_resolve_optional_flag(cluster.rows, TRUE, "cluster.rows")

  # Adjust filtering parameters based on input conditions
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  # Normalize count data if necessary.
  analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

  # Generate plots for each feature level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    otu_tax_agg <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter)

    # Select top k features if specified
    selected_features <- mStat_resolve_selected_features(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      features.plot = features.plot,
      top.k.plot = top.k.plot,
      top.k.func = top.k.func
    )

    merged_data <- mStat_prepare_taxa_long_data(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      value_col = "value",
      meta.dat = meta_tab,
      join = "inner"
    )

    pair_change <- mStat_prepare_taxa_pair_change_data(
      long.df = merged_data,
      feature.level = feature.level,
      subject.var = subject.var,
      time.var = time.var,
      change.base = change.base,
      feature.change.func = feature.change.func,
      context = "taxa change heatmap plotting"
    )
    combined_data <- pair_change$combined_data

    # Create a matrix of value differences
    value_diff_matrix <- combined_data %>%
      select(all_of(feature.level), !!sym(subject.var), value_diff) %>%
      tidyr::pivot_wider(names_from = !!sym(subject.var), values_from = value_diff) %>%
      mStat_as_taxa_feature_matrix(
        feature.level = feature.level,
        feature_in_column = TRUE
      )

    # Prepare metadata for annotation
    subject_meta <- data.frame(
      colnames(value_diff_matrix),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    colnames(subject_meta) <- subject.var
    subject_meta <- mStat_attach_subject_level_metadata(
      df = subject_meta,
      meta.dat = meta_tab,
      subject.var = subject.var,
      vars = c(group.var, strata.var)
    )

    aligned_subjects <- mStat_align_subject_metadata_to_matrix(
      value_matrix = value_diff_matrix,
      meta_tab = subject_meta,
      subject.var = subject.var,
      keep_vars = c(group.var, strata.var)
    )
    value_diff_matrix <- aligned_subjects$value_matrix
    sorted_meta_tab <- aligned_subjects$sorted_meta

    # Sort metadata by group variable if present
    if (!is.null(group.var)) {
      sorted_meta_tab <-
        sorted_meta_tab[order(sorted_meta_tab %>% select(!!sym(group.var)) %>% as.matrix()),]
    }

    # Rearrange the columns of the value difference matrix
    value_diff_matrix <- value_diff_matrix[, sorted_meta_tab[[subject.var]], drop = FALSE]

    # Prepare annotation columns for the heatmap
    if (!is.null(strata.var) && !is.null(group.var)){
      annotation_cols <-
        sorted_meta_tab %>%
        select(all_of(c(group.var, strata.var, subject.var))) %>%
        tibble::column_to_rownames(var = subject.var)
      annotation_col_sorted <-
        annotation_cols[order(annotation_cols[[strata.var]], annotation_cols[[group.var]]),]
    } else if (!is.null(group.var)){
      annotation_cols <-
        sorted_meta_tab %>%
        select(all_of(c(subject.var, group.var))) %>%
        tibble::column_to_rownames(var = subject.var)
      annotation_col_sorted <- annotation_cols
    } else {
      annotation_col_sorted <- NULL
    }

   # Rearrange value difference matrix columns if annotation is present
   if (!is.null(annotation_col_sorted)){
     value_diff_matrix <-
       value_diff_matrix[, rownames(annotation_col_sorted), drop = FALSE]
   }

    # Determine gaps for the heatmap based on strata or group variables
    if (!is.null(strata.var)) {
      strata_counts <- table(sorted_meta_tab[[strata.var]])
      gaps <- cumsum(strata_counts)[-length(strata_counts)]
    } else if (!is.null(group.var)) {
      group_counts <- table(sorted_meta_tab[[group.var]])
      gaps <- cumsum(group_counts)[-length(group_counts)]
    } else {
      gaps <- NULL
    }

    # Filter features to plot if specified
    if (!is.null(selected_features)) {
      value_diff_matrix <-
        value_diff_matrix[rownames(value_diff_matrix) %in% selected_features,, drop = FALSE]
    }

    # Set up color scheme for the heatmap
    n_colors <- 100
    col <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
    matrix_values <- as.numeric(value_diff_matrix)
    finite_values <- matrix_values[is.finite(matrix_values)]
    if (length(finite_values) == 0) {
      max_abs_val <- 1
    } else {
      max_abs_val <- max(abs(finite_values))
      if (!is.finite(max_abs_val) || max_abs_val <= 0) {
        max_abs_val <- 1
      }
    }

    zero_pos <- floor(n_colors / 2)
    my_col <-
      c(
        colorRampPalette(col[1:3])(zero_pos),
        colorRampPalette(col[3:5])(n_colors - zero_pos + 1)
      )
    break_points <-
      seq(-max_abs_val, max_abs_val, length.out = length(my_col) + 1)

    # Set up color palette for annotations
    color_vector <- mStat_get_palette(palette)

    if (!is.null(strata.var) && !is.null(group.var)){
      group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      strata_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

      annotation_colors_list <- setNames(
        list(group_colors, strata_colors),
        c(group.var, strata.var)
      )
    } else if(!is.null(group.var)){
      group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)

      annotation_colors_list <- setNames(
        list(group_colors),
        c(group.var)
      )
    } else {
      annotation_colors_list <- NULL
    }

    heatmap_gaps <- if (!is.null(gaps) && length(gaps) > 0 && ncol(value_diff_matrix) > 1) {
      gaps
    } else {
      NULL
    }

    # Generate the heatmap for individual data
    heatmap_plot <- pheatmap::pheatmap(
      mat = value_diff_matrix[order(rowMeans(abs(value_diff_matrix), na.rm = TRUE), decreasing = TRUE), , drop = FALSE],
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      annotation_legend = TRUE,
      show_colnames = TRUE,
      show_rownames = TRUE,
      border_color = NA,
      silent = TRUE,
      gaps_col = heatmap_gaps,
      fontsize = base.size,
      color = my_col,
      breaks = break_points,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Calculate average value differences for group-level heatmap
    combined_with_meta <- dplyr::left_join(
      combined_data %>% dplyr::select(all_of(feature.level), !!sym(subject.var), value_diff),
      sorted_meta_tab,
      by = subject.var,
      relationship = "many-to-one"
    )

    if (!is.null(strata.var)){
      average_value_diff_matrix <- combined_with_meta %>%
        dplyr::group_by(!!sym(feature.level), !!sym(group.var), !!sym(strata.var)) %>%
        dplyr::summarize(mean_value_diff = mean(value_diff), .groups = 'drop') %>%
        tidyr::pivot_wider(names_from = all_of(c(group.var, strata.var)), values_from = mean_value_diff) %>%
        mStat_as_taxa_feature_matrix(
          feature.level = feature.level,
          feature_in_column = TRUE
        )

      create_annotation_df <- function(colnames_vec, delimiter = "_") {
        split_names <- strsplit(colnames_vec, delimiter)
        if (length(unique(sapply(split_names, length))) != 1) {
          stop("All column names must have the same number of delimiters.")
        }

        annotation_df <- do.call(rbind.data.frame, split_names)
        colnames(annotation_df) <- c(group.var, strata.var)
        annotation_df[] <- lapply(annotation_df, function(x) factor(x, levels = unique(x)))
        annotation_matrix <- as.data.frame(annotation_df)
        rownames(annotation_matrix) <- colnames_vec
        annotation_matrix
      }

      annotation_df <- create_annotation_df(colnames(average_value_diff_matrix))
      annotation_df <- annotation_df[order(annotation_df[,strata.var]),]
    } else if (!is.null(group.var)){
      average_value_diff_matrix <- combined_with_meta %>%
        dplyr::group_by(!!sym(feature.level), !!sym(group.var)) %>%
        dplyr::summarize(mean_value_diff = mean(value_diff), .groups = 'drop') %>%
        tidyr::pivot_wider(names_from = all_of(c(group.var)), values_from = mean_value_diff) %>%
        mStat_as_taxa_feature_matrix(
          feature.level = feature.level,
          feature_in_column = TRUE
        )

      annotation_df <- data.frame(Var1 = colnames(average_value_diff_matrix)) %>%
        dplyr::mutate(Var2 = Var1) %>%
        tibble::column_to_rownames("Var2") %>%
        dplyr::rename(!!sym(group.var) := Var1)
    } else {
      average_value_diff_matrix <- combined_with_meta %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarize(mean_value_diff = mean(value_diff), .groups = 'drop') %>%
        mStat_as_taxa_feature_matrix(
          feature.level = feature.level,
          feature_in_column = TRUE
        )

      annotation_df <- NULL
    }

    # Filter features for plotting if specified
    if (!is.null(selected_features)) {
      average_value_diff_matrix <-
        average_value_diff_matrix[rownames(average_value_diff_matrix) %in% selected_features,, drop = FALSE]
    }

    # Reorder columns of average_value_diff_matrix to match annotation_df if applicable
    if (!is.null(annotation_df)){
      average_value_diff_matrix <- average_value_diff_matrix[, rownames(annotation_df), drop = FALSE]
    }

    # Generate the heatmap for average data
    average_heatmap_plot <- pheatmap::pheatmap(
      mat = average_value_diff_matrix[order(rowMeans(abs(average_value_diff_matrix), na.rm = TRUE), decreasing = TRUE), , drop = FALSE],
      annotation_col = annotation_df,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      annotation_legend = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      border_color = NA,
      silent = TRUE,
      gaps_col = NULL,
      fontsize = base.size,
      color = my_col,
      breaks = break_points
    )

    gg_average_heatmap_plot <- as.ggplot(average_heatmap_plot)

    # Convert feature.change.func to string if it's a custom function
    if (is.function(feature.change.func)) {
      feature.change.func = "custom function"
    }

    # Save individual heatmap as PDF if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_pair_indiv",
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
      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_heatmap_plot
      )
    }

    # Save average heatmap as PDF if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_pair_average",
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
      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_average_heatmap_plot
      )
    }
    
    # Create a list of both heatmap plots
    sub.plot_list <- list(gg_average_heatmap_plot, gg_heatmap_plot)
    names(sub.plot_list) <- c("average", "indiv")

    return(sub.plot_list)
  })

  # Name the plot list with feature levels and return
  names(plot_list) <- feature.level
  return(plot_list)
}
