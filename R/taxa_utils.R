#' @title Taxa Data Utilities
#' @description Helper functions for retrieving, aggregating, and filtering
#'   feature (OTU/ASV/taxa) tables. Centralizes the aggregation + extraction +
#'   filtering pipeline that was previously duplicated across 30+ files.
#' @name taxa_utils
#' @keywords internal
NULL


#' Retrieve a filtered feature table for a given taxonomic level
#'
#' Performs the three-step pipeline that appears in virtually every
#' \code{generate_taxa_*} function:
#' \enumerate{
#'   \item Aggregate features to \code{feature.level} if not already cached.
#'   \item Extract the aggregated table (or the raw table for \code{"original"}).
#'   \item Optionally filter by prevalence / abundance, then convert rownames
#'         to a column named \code{feature.level}.
#' }
#'
#' @param data.obj A MicrobiomeStat data object containing at least
#'   \code{feature.tab} and optionally \code{feature.agg.list}.
#' @param feature.level Character string specifying the taxonomic level to
#'   extract. Use \code{"original"} to retrieve \code{data.obj$feature.tab}
#'   without aggregation.
#' @param prev.filter Numeric prevalence threshold passed to
#'   \code{\link{mStat_filter}}. Set to 0 to skip filtering. Default 0.
#' @param abund.filter Numeric abundance threshold passed to
#'   \code{\link{mStat_filter}}. Set to 0 to skip filtering. Default 0.
#' @param feature.col Logical; if \code{TRUE} (default), rownames are converted
#'   to a leading column named \code{feature.level} via
#'   \code{tibble::rownames_to_column}. Set to \code{FALSE} when the downstream
#'   consumer expects a matrix-like object with rownames (e.g. \code{linda()}).
#'
#' @return When \code{feature.col = TRUE}, a data.frame with a leading column
#'   named \code{feature.level}. When \code{FALSE}, a data.frame with rownames.
#'
#' @keywords internal
get_taxa_data <- function(data.obj,
                          feature.level,
                          prev.filter = 0,
                          abund.filter = 0,
                          feature.col = TRUE) {

  # Step 1: Aggregate if needed (idempotent, cached in data.obj$feature.agg.list)
  if (is.null(data.obj$feature.agg.list[[feature.level]]) &
      feature.level != "original") {
    data.obj <-
      mStat_aggregate_by_taxonomy(data.obj = data.obj,
                                  feature.level = feature.level)
  }

  # Step 2: Extract the appropriate table
  if (feature.level != "original") {
    otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
  } else {
    otu_tax_agg <- data.obj$feature.tab
  }

  # Step 3: Filter
  result <- otu_tax_agg %>%
    as.data.frame() %>%
    mStat_filter(prev.filter = prev.filter,
                 abund.filter = abund.filter)

  if (feature.col) {
    result <- tibble::rownames_to_column(result, feature.level)
  }

  result
}


#' @keywords internal
mStat_as_taxa_feature_matrix <- function(feature.dat,
                                         feature.level,
                                         feature_in_column = FALSE) {
  feature.df <- as.data.frame(
    feature.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (feature_in_column) {
    if (!feature.level %in% colnames(feature.df)) {
      stop(
        "The feature column `", feature.level, "` is missing from the taxa data.",
        call. = FALSE
      )
    }
    feature.df <- tibble::column_to_rownames(feature.df, var = feature.level)
  }

  feature_ids <- rownames(feature.df)
  if (is.null(feature_ids) || anyNA(feature_ids) || anyDuplicated(feature_ids)) {
    stop(
      "Taxa data must have unique, non-missing feature row names.",
      call. = FALSE
    )
  }

  data.matrix(feature.df)
}


#' @keywords internal
mStat_normalize_feature_matrix_by_sample_total <- function(feature.mat) {
  feature.mat <- data.matrix(feature.mat)
  sample_totals <- colSums(feature.mat, na.rm = TRUE)
  normalized <- matrix(
    0,
    nrow = nrow(feature.mat),
    ncol = ncol(feature.mat),
    dimnames = dimnames(feature.mat)
  )

  nonzero_samples <- is.finite(sample_totals) & sample_totals != 0
  if (any(nonzero_samples)) {
    normalized[, nonzero_samples] <- sweep(
      feature.mat[, nonzero_samples, drop = FALSE],
      2,
      sample_totals[nonzero_samples],
      "/"
    )
  }

  normalized
}


#' @keywords internal
mStat_as_taxa_composition_matrix <- function(feature.dat,
                                             feature.level,
                                             feature_in_column = TRUE) {
  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.dat,
    feature.level = feature.level,
    feature_in_column = feature_in_column
  )

  mStat_normalize_feature_matrix_by_sample_total(feature.mat)
}


#' @keywords internal
mStat_meta_with_sample <- function(meta.dat,
                                   sample_col = "sample") {
  meta.df <- as.data.frame(
    meta.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (sample_col %in% colnames(meta.df)) {
    return(meta.df)
  }

  tibble::rownames_to_column(meta.df, var = sample_col)
}


#' @keywords internal
mStat_resolve_optional_flag <- function(flag,
                                        default,
                                        name) {
  if (is.null(flag)) {
    return(default)
  }

  if (!is.logical(flag) || length(flag) != 1 || is.na(flag)) {
    stop("`", name, "` should be a single TRUE/FALSE value or NULL.", call. = FALSE)
  }

  flag
}


#' @keywords internal
mStat_normalize_count_data_if_needed <- function(data.obj,
                                                 feature.dat.type,
                                                 message_text = paste(
                                                   "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
                                                 )) {
  if (!identical(feature.dat.type, "count")) {
    return(data.obj)
  }

  message(message_text)
  mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
}


#' @keywords internal
mStat_adjust_other_abundance_filter <- function(data.obj,
                                                feature.dat.type,
                                                abund.filter) {
  if (!identical(feature.dat.type, "other")) {
    return(abund.filter)
  }

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
    return(-Inf)
  }

  abund.filter
}


#' @keywords internal
mStat_get_taxa_change_ylabel <- function(feature.dat.type,
                                         feature.change.func) {
  base_label <- if (identical(feature.dat.type, "other")) {
    "Change in Abundance"
  } else {
    "Change in Relative Abundance"
  }

  suffix <- if (is.function(feature.change.func)) {
    "custom function"
  } else {
    as.character(feature.change.func)
  }

  paste0(base_label, " (", suffix, ")")
}


#' @keywords internal
mStat_get_taxa_value_ylabel <- function(feature.dat.type,
                                        transform = NULL) {
  if (identical(feature.dat.type, "other")) {
    return("Abundance")
  }

  if (is.null(transform)) {
    return("Relative Abundance")
  }

  paste0("Relative Abundance (", transform, ")")
}


#' @keywords internal
mStat_get_taxa_dotplot_size_range <- function(taxa.levels) {
  if (taxa.levels <= 2) {
    return(c(40, 57))
  } else if (taxa.levels <= 4) {
    return(c(35, 42))
  } else if (taxa.levels <= 6) {
    return(c(30, 37))
  } else if (taxa.levels <= 8) {
    return(c(25, 33))
  } else if (taxa.levels < 10) {
    return(c(17, 20))
  } else if (taxa.levels < 20) {
    return(c(10, 15))
  } else if (taxa.levels < 30) {
    return(c(8, 13))
  } else if (taxa.levels < 40) {
    return(c(6, 10))
  } else if (taxa.levels < 50) {
    return(c(4, 8))
  }

  c(1, 4)
}


#' @keywords internal
mStat_get_spaghettiplot_text_sizes <- function(base.size,
                                              variant = c("default", "wide_panel")) {
  variant <- match.arg(variant)

  if (identical(variant, "wide_panel")) {
    return(list(
      title = base.size * 1.25,
      axis.title = base.size * 1.5,
      axis.text = base.size * 2,
      legend.title = base.size * 2,
      legend.text = base.size * 2
    ))
  }

  list(
    title = base.size * 1.25,
    axis.title = base.size * 0.75,
    axis.text = base.size * 0.5,
    legend.title = base.size,
    legend.text = base.size * 0.75
  )
}


#' @keywords internal
mStat_prepare_meta_tab <- function(meta.dat,
                                   vars,
                                   sample_col = "sample") {
  vars <- unlist(vars, recursive = TRUE, use.names = FALSE)
  vars <- vars[!vapply(vars, is.null, logical(1))]
  vars <- unique(as.character(vars))

  if (length(vars) == 0) {
    return(mStat_meta_with_sample(meta.dat, sample_col = sample_col))
  }

  selected_meta <- as.data.frame(
    meta.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>% dplyr::select(dplyr::all_of(vars))

  mStat_meta_with_sample(selected_meta, sample_col = sample_col)
}


#' @keywords internal
mStat_prepare_taxa_long_data <- function(feature.dat,
                                         feature.level,
                                         value_col = "value",
                                         meta.dat = NULL,
                                         sample_col = "sample",
                                         feature_in_column = TRUE,
                                         join = c("left", "inner")) {
  join <- match.arg(join)

  feature.df <- as.data.frame(
    feature.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!feature_in_column) {
    feature.df <- tibble::rownames_to_column(feature.df, var = feature.level)
  }

  long.df <- feature.df %>%
    tidyr::pivot_longer(
      cols = -all_of(feature.level),
      names_to = sample_col,
      values_to = value_col
    ) %>%
    dplyr::mutate(!!value_col := as.numeric(.data[[value_col]]))

  if (is.null(meta.dat)) {
    return(long.df)
  }

  meta.df <- mStat_meta_with_sample(meta.dat = meta.dat, sample_col = sample_col)

  if (join == "inner") {
    dplyr::inner_join(long.df, meta.df, by = sample_col)
  } else {
    dplyr::left_join(long.df, meta.df, by = sample_col)
  }
}


#' @keywords internal
mStat_prepare_taxa_single_context <- function(data.obj,
                                              time.var = NULL,
                                              t.level = NULL,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              subject.var = NULL,
                                              add_subject_id = FALSE) {
  if (add_subject_id) {
    data.obj$meta.dat$subject.id <- rownames(data.obj$meta.dat)
    subject.var <- "subject.id"
  }

  if (!is.null(time.var) && !is.null(t.level)) {
    data.obj <- mStat_subset_by_meta_values(data.obj, time.var, t.level)
  }

  list(
    data.obj = data.obj,
    meta_tab = select_meta_vars(data.obj$meta.dat, subject.var, group.var, time.var, strata.var),
    subject.var = subject.var
  )
}


#' @keywords internal
mStat_prepare_taxa_long_context <- function(data.obj,
                                            subject.var,
                                            time.var,
                                            group.var = NULL,
                                            strata.var = NULL,
                                            t0.level = NULL,
                                            ts.levels = NULL) {
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  list(
    data.obj = data.obj,
    meta_tab = select_meta_vars(data.obj$meta.dat, subject.var, group.var, time.var, strata.var)
  )
}


#' @keywords internal
mStat_prepare_taxa_pair_change_data <- function(long.df,
                                                feature.level,
                                                subject.var,
                                                time.var,
                                                change.base = NULL,
                                                feature.change.func,
                                                context) {
  pair_slices <- mStat_prepare_pair_time_slices(
    df = long.df,
    time.var = time.var,
    change.base = change.base,
    context = context
  )

  combined_data <- pair_slices$data_time_1 %>%
    dplyr::inner_join(
      pair_slices$data_time_2,
      by = c(feature.level, subject.var),
      suffix = c("_time_1", "_time_2")
    ) %>%
    dplyr::mutate(value_diff = compute_taxa_change(
      value_after = value_time_2,
      value_before = value_time_1,
      method = feature.change.func,
      feature_id = .data[[feature.level]]
    ))

  list(
    combined_data = combined_data,
    change.base = pair_slices$change.base,
    change.after = pair_slices$change.after
  )
}


#' @keywords internal
mStat_attach_pair_metadata <- function(df,
                                       meta_tab,
                                       subject.var,
                                       time.var,
                                       mode = c("followup_time", "subject"),
                                       change.after = NULL) {
  mode <- match.arg(mode)

  meta_source <- if (mode == "followup_time") {
    if (is.null(change.after)) {
      stop("`change.after` is required when mode = 'followup_time'.", call. = FALSE)
    }
    mStat_prepare_subject_metadata(
      meta.dat = meta_tab,
      subject.var = subject.var,
      vars = setdiff(colnames(meta_tab), c("sample", subject.var, time.var)),
      time.var = time.var,
      time.value = change.after
    )
  } else {
    mStat_prepare_subject_metadata(
      meta.dat = meta_tab,
      subject.var = subject.var,
      vars = setdiff(colnames(meta_tab), c("sample", subject.var, time.var))
    )
  }

  if (is.null(meta_source)) {
    return(df)
  }

  dplyr::left_join(df, meta_source, by = subject.var, relationship = "many-to-one")
}


#' @keywords internal
mStat_align_subject_metadata_to_matrix <- function(value_matrix,
                                                   meta_tab,
                                                   subject.var,
                                                   keep_vars) {
  aligned_meta <- meta_tab %>%
    dplyr::filter(.data[[subject.var]] %in% colnames(value_matrix)) %>%
    dplyr::select(dplyr::all_of(c(subject.var, keep_vars))) %>%
    dplyr::distinct(.data[[subject.var]], .keep_all = TRUE) %>%
    tibble::as_tibble()

  order_index <- match(colnames(value_matrix), aligned_meta[[subject.var]])

  suppressWarnings({
    sorted_meta <- aligned_meta[order_index, , drop = FALSE]
    rownames(sorted_meta) <- sorted_meta[[subject.var]]
  })

  list(
    value_matrix = value_matrix[, rownames(sorted_meta), drop = FALSE],
    sorted_meta = sorted_meta
  )
}


#' @keywords internal
mStat_resolve_selected_features <- function(feature.dat = NULL,
                                            feature.level,
                                            features.plot = NULL,
                                            top.k.plot = NULL,
                                            top.k.func = NULL,
                                            taxa.levels = NULL,
                                            fallback_n = NULL) {
  selected_features <- features.plot

  if (is.null(selected_features) && !is.null(top.k.plot) && !is.null(top.k.func)) {
    if (is.null(feature.dat)) {
      stop("`feature.dat` is required when selecting top-k features.", call. = FALSE)
    }
    computed_values <- compute_function(top.k.func, feature.dat, feature.level)
    top_k_n <- min(top.k.plot, length(computed_values))
    selected_features <- names(sort(computed_values, decreasing = TRUE)[seq_len(top_k_n)])
  }

  if (is.null(selected_features) && !is.null(taxa.levels) && !is.null(fallback_n)) {
    selected_features <- taxa.levels[seq_len(min(fallback_n, length(taxa.levels)))]
  }

  selected_features
}


#' @keywords internal
mStat_ensure_group_placeholder <- function(df,
                                           group.var,
                                           value = "ALL",
                                           column_name = ".mstat_group") {
  if (!is.null(group.var)) {
    return(list(df = df, group.var = group.var))
  }

  resolved_column_name <- column_name
  suffix <- 1
  while (resolved_column_name %in% colnames(df)) {
    resolved_column_name <- paste0(column_name, "_", suffix)
    suffix <- suffix + 1
  }

  df[[resolved_column_name]] <- value
  list(df = df, group.var = resolved_column_name)
}


#' @keywords internal
mStat_summarize_mean_by_groups <- function(long.df,
                                           feature.level,
                                           group_vars,
                                           value_col = "value",
                                           mean_col = "mean_value",
                                           mean_transform = NULL) {
  summary.df <- long.df %>%
    dplyr::group_by(dplyr::across(all_of(c(feature.level, group_vars)))) %>%
    dplyr::summarise(.mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  if (!is.null(mean_transform)) {
    summary.df$.mean_value <- mean_transform(summary.df$.mean_value)
  }

  summary.df %>% dplyr::rename(!!mean_col := .mean_value)
}


#' @keywords internal
mStat_filter_test_result_features <- function(result.df,
                                              features.plot = NULL) {
  if (is.null(features.plot)) {
    return(result.df)
  }

  result.df %>% dplyr::filter(.data[["Variable"]] %in% features.plot)
}


#' @keywords internal
mStat_filter_taxa_features <- function(feature.dat,
                                       feature.level,
                                       features.plot = NULL) {
  if (is.null(features.plot)) {
    return(feature.dat)
  }

  feature.dat[feature.dat[[feature.level]] %in% features.plot, , drop = FALSE]
}


#' @keywords internal
mStat_prepare_stack_levels <- function(composition.mat,
                                       feature.level,
                                       feature.number,
                                       other_first = TRUE,
                                       other_inclusive = FALSE) {
  avg_abund <- rowMeans(composition.mat)
  other.abund.cutoff <- sort(avg_abund, decreasing = TRUE)[feature.number]

  feature.df <- composition.mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(feature.level)

  if (!is.na(other.abund.cutoff)) {
    if (other_inclusive) {
      feature.df[[feature.level]][avg_abund <= other.abund.cutoff] <- "Other"
    } else {
      feature.df[[feature.level]][avg_abund < other.abund.cutoff] <- "Other"
    }
  }

  long.df <- feature.df %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), sum), .groups = "drop") %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(feature.level),
      names_to = "sample",
      values_to = "value"
    )

  overall_summary <- long.df %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::summarise(overall_mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(is_other = .data[[feature.level]] == "Other")

  if (other_first) {
    overall_summary <- overall_summary %>% dplyr::arrange(!is_other, overall_mean)
  } else {
    overall_summary <- overall_summary %>% dplyr::arrange(is_other, overall_mean)
  }

  new_levels <- overall_summary[[feature.level]]
  if (other_first && !is.na(other.abund.cutoff) && "Other" %in% new_levels) {
    new_levels <- c("Other", setdiff(new_levels, "Other"))
  }

  list(
    long.df = long.df,
    new_levels = new_levels,
    other.abund.cutoff = other.abund.cutoff
  )
}


#' @keywords internal
mStat_prepare_stacked_positions <- function(long.df,
                                            feature.level,
                                            id_var,
                                            ordered_levels,
                                            terminal_ids,
                                            cumulative_col = "cumulative_value",
                                            next_cumulative_col = "next_cumulative_value") {
  long.df %>%
    dplyr::group_by(!!rlang::sym(id_var)) %>%
    dplyr::mutate(!!rlang::sym(feature.level) := factor(.data[[feature.level]], levels = ordered_levels)) %>%
    dplyr::arrange(match(.data[[feature.level]], ordered_levels), .by_group = TRUE) %>%
    dplyr::mutate(!!cumulative_col := 1 - cumsum(.data[["value"]])) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::mutate(
      !!next_cumulative_col := dplyr::if_else(
        .data[[id_var]] %in% terminal_ids,
        NA_real_,
        dplyr::lead(.data[[cumulative_col]])
      )
    ) %>%
    dplyr::ungroup()
}


#' @keywords internal
mStat_prepare_average_stack_data <- function(long.df,
                                             feature.level,
                                             group.var,
                                             time.var,
                                             ordered_levels,
                                             terminal_time_values,
                                             bar_width = 0.6,
                                             bar_spacing = bar_width / 2) {
  df_average <- mStat_summarize_mean_by_groups(
    long.df = long.df,
    feature.level = feature.level,
    group_vars = c(group.var, time.var),
    value_col = "value",
    mean_col = "mean_value"
  ) %>%
    dplyr::mutate(!!rlang::sym(feature.level) := factor(.data[[feature.level]], levels = ordered_levels)) %>%
    dplyr::arrange(match(.data[[feature.level]], ordered_levels), .data[[group.var]], .data[[time.var]]) %>%
    dplyr::group_by(.data[[group.var]], .data[[time.var]]) %>%
    dplyr::mutate(cumulative_mean_value = 1 - cumsum(.data[["mean_value"]])) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data[[feature.level]]) %>%
    dplyr::mutate(
      next_cumulative_mean_value = dplyr::if_else(
        .data[[time.var]] %in% terminal_time_values,
        NA_real_,
        dplyr::lead(.data[["cumulative_mean_value"]])
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      joint_factor = interaction(.data[[time.var]], .data[[group.var]]),
      x_offset = ifelse(.data[["cumulative_mean_value"]] == 0,
                        (bar_width + bar_spacing) / 2,
                        -(bar_width + bar_spacing) / 2)
    )

  df_average$joint_factor <- droplevels(df_average$joint_factor)
  df_average$joint_factor_numeric <- match(df_average$joint_factor, levels(df_average$joint_factor))

  list(
    df = df_average,
    labels = sub("\\..*", "", levels(df_average$joint_factor))
  )
}


#' @keywords internal
mStat_summarize_grouped_taxa_long <- function(long.df,
                                              feature.level,
                                              group_vars = NULL,
                                              value_col = "value",
                                              mean_col = "mean_abundance",
                                              prevalence_col = "prevalence",
                                              mean_transform = NULL) {
  summary.df <- long.df %>%
    dplyr::group_by(dplyr::across(all_of(c(group_vars, feature.level)))) %>%
    dplyr::summarise(
      .mean_value = mean(.data[[value_col]], na.rm = TRUE),
      .prevalence = mean(.data[[value_col]] != 0, na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(mean_transform)) {
    summary.df$.mean_value <- mean_transform(summary.df$.mean_value)
  }

  summary.df %>%
    dplyr::rename(
      !!mean_col := .mean_value,
      !!prevalence_col := .prevalence
    )
}


#' @keywords internal
mStat_prepare_change_heatmap_long_matrix <- function(summary.df,
                                                     feature.level,
                                                     group.var,
                                                     time.var,
                                                     baseline_time,
                                                     followup_times,
                                                     feature.change.func,
                                                     value_col = "mean_value") {
  wide.df <- summary.df %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(time.var), values_from = !!rlang::sym(value_col))

  for (ts in followup_times) {
    change_col_name <- paste0("change_", ts)
    wide.df[[change_col_name]] <- compute_taxa_change(
      value_after = wide.df[[as.character(ts)]],
      value_before = wide.df[[as.character(baseline_time)]],
      method = feature.change.func,
      feature_id = wide.df[[feature.level]],
      verbose = identical(ts, followup_times[[1]])
    )
  }

  heatmap.df <- wide.df %>%
    dplyr::select(dplyr::all_of(c(feature.level, group.var, paste0("change_", followup_times))))

  feature_ids <- unique(heatmap.df[[feature.level]])
  group_levels <- unique(heatmap.df[[group.var]])
  heatmap_mat <- matrix(
    NA_real_,
    nrow = length(feature_ids),
    ncol = length(group_levels) * length(followup_times),
    dimnames = list(feature_ids, as.vector(outer(group_levels, followup_times, paste, sep = "_")))
  )

  for (group_value in group_levels) {
    group_df <- heatmap.df[heatmap.df[[group.var]] == group_value, , drop = FALSE]
    row_index <- match(group_df[[feature.level]], feature_ids)
    heatmap_mat[row_index, paste(group_value, followup_times, sep = "_")] <-
      as.matrix(group_df[, paste0("change_", followup_times), drop = FALSE])
  }

  heatmap_mat
}


#' @keywords internal
mStat_transform_taxa_long_values <- function(long.df,
                                             feature.level,
                                             value_col = "value",
                                             feature.dat.type = c("count", "proportion", "other"),
                                             transform = c("identity", "sqrt", "log")) {
  feature.dat.type <- match.arg(feature.dat.type)
  transform <- match.arg(transform)

  if (!feature.dat.type %in% c("count", "proportion") || transform == "identity") {
    return(long.df)
  }

  if (transform == "sqrt") {
    long.df[[value_col]] <- sqrt(long.df[[value_col]])
    return(long.df)
  }

  min_half_nonzero <- long.df %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::filter(sum(.data[[value_col]], na.rm = TRUE) != 0) %>%
    dplyr::summarise(
      min_half_value = min(.data[[value_col]][.data[[value_col]] > 0], na.rm = TRUE) / 2,
      .groups = "drop"
    )

  long.df %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::filter(sum(.data[[value_col]], na.rm = TRUE) != 0) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(min_half_nonzero, by = feature.level) %>%
    dplyr::mutate(
      !!value_col := ifelse(
        .data[[value_col]] == 0,
        log10(min_half_value),
        log10(.data[[value_col]])
      )
    ) %>%
    dplyr::select(-min_half_value)
}


#' @keywords internal
mStat_prepare_taxa_clr_long_data <- function(feature.dat,
                                             feature.level,
                                             prev.filter = 0,
                                             abund.filter = -Inf,
                                             sample_col = "sample",
                                             value_col = "value",
                                             feature_in_column = FALSE) {
  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.dat,
    feature.level = feature.level,
    feature_in_column = feature_in_column
  )

  sample_totals <- colSums(feature.mat, na.rm = TRUE)
  valid_samples <- is.finite(sample_totals) & sample_totals > 0

  if (!any(valid_samples)) {
    return(tibble::tibble(
      !!feature.level := character(),
      !!sample_col := character(),
      !!value_col := numeric()
    ))
  }

  if (any(!valid_samples)) {
    warning(
      "Found ", sum(!valid_samples),
      " samples with zero total abundance during CLR preparation; these samples were excluded.",
      call. = FALSE
    )
  }

  sample_ids <- colnames(feature.mat)[valid_samples]
  transformed_samples <- lapply(sample_ids, function(sample_id) {
    sample_values <- feature.mat[, sample_id]
    positive_values <- sample_values[sample_values > 0 & is.finite(sample_values)]

    pseudocount <- min(positive_values, na.rm = TRUE) / 2
    imputed_values <- ifelse(sample_values == 0, pseudocount, sample_values)
    geometric_mean <- exp(mean(log(imputed_values), na.rm = TRUE))

    log(imputed_values / geometric_mean)
  })

  clr.mat <- do.call(cbind, transformed_samples)
  rownames(clr.mat) <- rownames(feature.mat)
  colnames(clr.mat) <- sample_ids

  clr.filtered <- mStat_filter(
    clr.mat,
    prev.filter = prev.filter,
    abund.filter = abund.filter
  )

  mStat_prepare_taxa_long_data(
    feature.dat = clr.filtered,
    feature.level = feature.level,
    value_col = value_col,
    sample_col = sample_col,
    feature_in_column = FALSE
  )
}


#' @keywords internal
.mStat_build_taxa_linda_formula <- function(group.var = NULL,
                                            adj.vars = NULL,
                                            time.var = NULL,
                                            subject.var = NULL,
                                            random_slopes = TRUE) {
  adj.vars_str <- if (!is.null(adj.vars)) paste(adj.vars, collapse = " + ") else NULL

  if (is.null(time.var)) {
    if (is.null(group.var)) {
      fixed_effects <- if (!is.null(adj.vars_str)) adj.vars_str else "1"
    } else if (!is.null(adj.vars_str)) {
      fixed_effects <- paste(adj.vars_str, "+", group.var)
    } else {
      fixed_effects <- group.var
    }
    random_effects <- paste("(1 |", subject.var, ")")
  } else {
    if (is.null(group.var)) {
      fixed_effects <- if (!is.null(adj.vars_str)) paste(adj.vars_str, "+", time.var) else time.var
    } else if (!is.null(adj.vars_str)) {
      fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
    } else {
      fixed_effects <- paste(group.var, "*", time.var)
    }

    random_effects <- if (random_slopes) {
      paste("(1 +", time.var, "|", subject.var, ")")
    } else {
      paste("(1 |", subject.var, ")")
    }
  }

  paste(fixed_effects, random_effects, sep = " + ")
}


#' @keywords internal
.mStat_prune_zero_total_samples <- function(feature.dat,
                                            meta.dat) {
  keep_samples <- colSums(feature.dat) > 0

  if (!any(keep_samples)) {
    return(NULL)
  }

  list(
    feature.dat = feature.dat[, keep_samples, drop = FALSE],
    meta.dat = meta.dat[keep_samples, , drop = FALSE]
  )
}


#' @keywords internal
.mStat_get_group_reference_level <- function(meta.dat,
                                             group.var) {
  if (is.null(group.var)) {
    return(NULL)
  }

  group_values <- meta.dat[[group.var]]
  if (!(is.factor(group_values) || is.character(group_values))) {
    return(NULL)
  }

  levels(as.factor(group_values))[1]
}


#' @keywords internal
.mStat_run_linda <- function(feature.dat,
                             meta.dat,
                             formula,
                             feature.dat.type,
                             extra_args = list(),
                             prev.filter = NULL,
                             mean.abund.filter = NULL,
                             fallback_formula = NULL,
                             fallback_message = NULL,
                             muffle_all_filtered_warning = FALSE) {
  base_args <- list(
    feature.dat = feature.dat,
    meta.dat = meta.dat,
    feature.dat.type = feature.dat.type
  )

  if (!is.null(prev.filter)) {
    base_args$prev.filter <- prev.filter
  }

  if (!is.null(mean.abund.filter)) {
    base_args$mean.abund.filter <- mean.abund.filter
  }

  call_linda <- function(formula_text) {
    linda_args <- c(
      base_args,
      list(formula = paste("~", formula_text)),
      extra_args
    )

    if (!muffle_all_filtered_warning) {
      return(do.call(linda, linda_args))
    }

    withCallingHandlers(
      do.call(linda, linda_args),
      warning = function(w) {
        if (grepl("All features were filtered out", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  if (is.null(fallback_formula)) {
    return(call_linda(formula))
  }

  tryCatch(
    call_linda(formula),
    error = function(e) {
      message("Error in linda: ", conditionMessage(e))
      if (!is.null(fallback_message)) {
        message(fallback_message)
      }
      call_linda(fallback_formula)
    }
  )
}


#' @keywords internal
.mStat_extract_pair_linda_outputs <- function(linda_output,
                                              group_var = NULL,
                                              time_var = NULL,
                                              reference_level = NULL) {
  result_list <- list()
  output_names <- names(linda_output)

  if (is.null(group_var)) {
    for (df_name in output_names) {
      result_list[[df_name]] <- linda_output[[df_name]]
    }
    return(result_list)
  }

  group_prefix_pattern <- paste0("^", mStat_escape_regex(group_var))
  matching_dfs <- grep(group_prefix_pattern, output_names, value = TRUE)
  for (df_name in matching_dfs) {
    group_value <- strsplit(df_name, split = ":", fixed = TRUE)[[1]][1]
    group_value <- sub(group_prefix_pattern, "", group_value)

    is_interaction <- !is.null(time_var) && grepl(paste0(":", time_var), df_name, fixed = TRUE)

    label <- if (is.null(reference_level)) {
      if (is_interaction) {
        df_name
      } else {
        group_var
      }
    } else {
      base_label <- if (nzchar(group_value)) {
        paste0(group_value, " vs ", reference_level, " (Reference)")
      } else {
        group_var
      }

      if (is_interaction) {
        paste0(base_label, " [Interaction]")
      } else {
        paste0(base_label, " [Main Effect]")
      }
    }

    result_list[[label]] <- linda_output[[df_name]]
  }

  result_list
}


#' @keywords internal
.mStat_extract_trend_linda_outputs <- function(linda_output,
                                               group_var = NULL,
                                               time_var = NULL,
                                               reference_level = NULL) {
  result_list <- list()
  output_names <- names(linda_output)

  if (is.null(group_var)) {
    if (!is.null(time_var) && time_var %in% output_names) {
      result_list[[time_var]] <- linda_output[[time_var]]
    }
    return(result_list)
  }

  if (!is.null(time_var) && time_var %in% output_names) {
    result_list[[time_var]] <- linda_output[[time_var]]
  }

  group_prefix_pattern <- paste0("^", mStat_escape_regex(group_var))
  matching_dfs <- grep(group_prefix_pattern, output_names, value = TRUE)
  for (df_name in matching_dfs) {
    term_head <- strsplit(df_name, split = ":", fixed = TRUE)[[1]][1]
    group_value <- sub(group_prefix_pattern, "", term_head)
    is_interaction <- !is.null(time_var) && grepl(paste0(":", time_var), df_name, fixed = TRUE)

    if (nzchar(group_value)) {
      base_label <- paste0(group_value, " vs ", reference_level, " (Reference)")
      label <- if (is_interaction) {
        paste0(base_label, " [Interaction]")
      } else {
        paste0(base_label, " [Main Effect]")
      }
    } else {
      label <- if (is_interaction) {
        paste0(group_var, ":", time_var)
      } else {
        group_var
      }
    }

    result_list[[label]] <- linda_output[[df_name]]
  }

  result_list
}


#' @keywords internal
mStat_summarize_taxa_features <- function(feature.dat,
                                          feature.level,
                                          feature_in_column = FALSE) {
  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.dat,
    feature.level = feature.level,
    feature_in_column = feature_in_column
  )
  feature_ids <- rownames(feature.mat)

  tibble::tibble(
    !!feature.level := feature_ids,
    avg_abundance = unname(rowMeans(feature.mat, na.rm = TRUE)),
    prevalence = unname(rowMeans(feature.mat != 0, na.rm = TRUE))
  )
}


#' @keywords internal
mStat_format_linda_feature_results <- function(result.df,
                                               feature.level,
                                               feature.stats,
                                               include_significant = FALSE) {
  required_cols <- c("log2FoldChange", "lfcSE", "pvalue", "padj")
  missing_cols <- setdiff(required_cols, colnames(result.df))
  if (length(missing_cols) != 0) {
    stop(
      "LinDA-style results are missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (include_significant && !"reject" %in% colnames(result.df)) {
    stop(
      "LinDA-style results must include `reject` when `include_significant = TRUE`.",
      call. = FALSE
    )
  }

  formatted <- tibble::rownames_to_column(
    as.data.frame(result.df, stringsAsFactors = FALSE),
    var = feature.level
  ) %>%
    dplyr::left_join(feature.stats, by = feature.level)

  select_cols <- c(
    feature.level,
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj",
    if (include_significant) "reject",
    "avg_abundance",
    "prevalence"
  )

  formatted %>%
    dplyr::select(all_of(select_cols)) %>%
    dplyr::rename(
      Variable = all_of(feature.level),
      Coefficient = log2FoldChange,
      SE = lfcSE,
      P.Value = pvalue,
      Adjusted.P.Value = padj,
      Mean.Abundance = avg_abundance,
      Prevalence = prevalence
    ) %>%
    {
      if (include_significant) {
        dplyr::rename(., Significant = reject)
      } else {
        .
      }
    }
}
