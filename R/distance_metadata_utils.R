#' @title Distance Metadata Utilities
#' @description Internal helpers for attaching and recovering metadata from
#'   precomputed distance objects. Distances naturally store sample labels, not
#'   arbitrary metadata, so these helpers make the metadata contract explicit.
#' @name distance_metadata_utils
#' @keywords internal
NULL


#' @keywords internal
mStat_select_metadata_columns <- function(meta.dat, vars = NULL) {
  meta.dat <- as.data.frame(meta.dat, stringsAsFactors = FALSE)

  vars <- unique(unlist(vars, use.names = FALSE))
  vars <- vars[!is.na(vars)]

  if (length(vars) == 0) {
    return(meta.dat[, 0, drop = FALSE])
  }

  missing_vars <- setdiff(vars, colnames(meta.dat))
  if (length(missing_vars) != 0) {
    stop(
      "The following metadata variables are missing: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }

  meta.dat[, vars, drop = FALSE]
}


#' @keywords internal
mStat_get_dist_labels <- function(dist.matrix) {
  labels <- attr(dist.matrix, "Labels")
  if (is.null(labels)) {
    labels <- attr(dist.matrix, "labels")
  }
  if (is.null(labels)) {
    labels <- rownames(as.matrix(dist.matrix))
  }
  labels
}


#' @keywords internal
mStat_validate_dist_object <- function(dist.obj, dist.name = NULL) {
  if (is.null(dist.obj) || !is.list(dist.obj) || length(dist.obj) == 0) {
    stop("dist.obj must be a non-empty list of distance matrices.", call. = FALSE)
  }

  if (is.null(dist.name)) {
    dist.name <- names(dist.obj)
  }

  missing_dist <- setdiff(dist.name, names(dist.obj))
  if (length(missing_dist) != 0) {
    stop(
      "dist.obj is missing the following requested distance matrices: ",
      paste(missing_dist, collapse = ", "),
      call. = FALSE
    )
  }

  reference_labels <- NULL
  reference_key <- NULL

  for (dist_key in dist.name) {
    dist_matrix <- dist.obj[[dist_key]]
    if (is.null(dist_matrix)) {
      stop("dist.obj[['", dist_key, "']] is NULL.", call. = FALSE)
    }

    matrix_view <- as.matrix(dist_matrix)
    labels <- mStat_get_dist_labels(dist_matrix)

    if (!is.matrix(matrix_view) || nrow(matrix_view) != ncol(matrix_view)) {
      stop(
        "dist.obj[['", dist_key, "']] must be coercible to a square matrix.",
        call. = FALSE
      )
    }

    if (is.null(labels) || length(labels) == 0) {
      stop(
        "dist.obj[['", dist_key, "']] must carry sample labels.",
        call. = FALSE
      )
    }

    if (anyNA(labels) || anyDuplicated(labels)) {
      stop(
        "dist.obj[['", dist_key, "']] has missing or duplicated sample labels.",
        call. = FALSE
      )
    }

    if (nrow(matrix_view) != length(labels)) {
      stop(
        "dist.obj[['", dist_key, "']] has inconsistent dimensions and labels.",
        call. = FALSE
      )
    }

    if (is.null(rownames(matrix_view)) || is.null(colnames(matrix_view)) ||
        !identical(rownames(matrix_view), colnames(matrix_view)) ||
        !identical(rownames(matrix_view), labels)) {
      stop(
        "dist.obj[['", dist_key, "']] must have matching row, column, and distance labels.",
        call. = FALSE
      )
    }

    if (is.null(reference_labels)) {
      reference_labels <- labels
      reference_key <- dist_key
    } else if (!identical(labels, reference_labels)) {
      stop(
        "All requested distance matrices must share identical sample labels and order. ",
        "dist.obj[['", dist_key, "']] does not match dist.obj[['", reference_key, "']].",
        call. = FALSE
      )
    }
  }

  invisible(dist.obj)
}


#' @keywords internal
mStat_align_pc_object <- function(pc.obj, dist.obj, dist.name = NULL) {
  if (is.null(pc.obj)) {
    return(pc.obj)
  }

  if (is.null(dist.name)) {
    dist.name <- intersect(names(pc.obj), names(dist.obj))
  }

  for (dist_key in dist.name) {
    dist_labels <- mStat_get_dist_labels(dist.obj[[dist_key]])
    points <- as.matrix(pc.obj[[dist_key]]$points)

    if (!setequal(rownames(points), dist_labels) || nrow(points) != length(dist_labels)) {
      stop(
        "pc.obj[['", dist_key, "']] sample labels must match dist.obj[['", dist_key, "']].",
        call. = FALSE
      )
    }

    pc.obj[[dist_key]]$points <- points[dist_labels, , drop = FALSE]
  }

  pc.obj
}


#' @keywords internal
mStat_validate_pc_object <- function(pc.obj,
                                     dist.name = NULL,
                                     dist.obj = NULL,
                                     required_pc_axes = NULL,
                                     sample_ids = NULL) {
  if (is.null(pc.obj) || !is.list(pc.obj) || length(pc.obj) == 0) {
    stop("pc.obj must be a non-empty list of ordination results.", call. = FALSE)
  }

  if (is.null(dist.name)) {
    dist.name <- names(pc.obj)
  }

  missing_pc <- setdiff(dist.name, names(pc.obj))
  if (length(missing_pc) != 0) {
    stop(
      "pc.obj is missing the following requested ordination results: ",
      paste(missing_pc, collapse = ", "),
      call. = FALSE
    )
  }

  for (dist_key in dist.name) {
    pc_entry <- pc.obj[[dist_key]]
    if (is.null(pc_entry) || is.null(pc_entry$points)) {
      stop(
        "pc.obj[['", dist_key, "']] must contain sample coordinates in `$points`.",
        call. = FALSE
      )
    }

    points <- as.matrix(pc_entry$points)
    point_ids <- rownames(points)

    if (!is.matrix(points) || nrow(points) == 0 || ncol(points) == 0) {
      stop(
        "pc.obj[['", dist_key, "']] must contain a non-empty coordinate matrix.",
        call. = FALSE
      )
    }

    if (is.null(point_ids) || anyNA(point_ids) || anyDuplicated(point_ids)) {
      stop(
        "pc.obj[['", dist_key, "']] must have unique, non-missing sample row names.",
        call. = FALSE
      )
    }

    if (!is.null(required_pc_axes) && ncol(points) < required_pc_axes) {
      stop(
        "pc.obj[['", dist_key, "']] has only ",
        ncol(points),
        " axis/axes, but ",
        required_pc_axes,
        " are required.",
        call. = FALSE
      )
    }

    if (!is.null(dist.obj)) {
      dist_labels <- mStat_get_dist_labels(dist.obj[[dist_key]])
      if (!identical(point_ids, dist_labels)) {
        stop(
          "pc.obj[['", dist_key, "']] sample order must match dist.obj[['", dist_key, "']].",
          call. = FALSE
        )
      }
    }

    if (!is.null(sample_ids)) {
      missing_samples <- setdiff(sample_ids, point_ids)
      if (length(missing_samples) != 0) {
        stop(
          "pc.obj[['", dist_key, "']] is missing the following samples: ",
          paste(missing_samples, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  invisible(pc.obj)
}


#' @keywords internal
mStat_attach_dist_metadata <- function(dist.matrix, meta.dat = NULL) {
  if (is.null(meta.dat)) {
    return(dist.matrix)
  }

  labels <- mStat_get_dist_labels(dist.matrix)
  if (is.null(labels)) {
    return(dist.matrix)
  }

  meta.dat <- as.data.frame(meta.dat, stringsAsFactors = FALSE)
  if (is.null(rownames(meta.dat))) {
    return(dist.matrix)
  }

  missing_rows <- setdiff(labels, rownames(meta.dat))
  if (length(missing_rows) != 0) {
    return(dist.matrix)
  }

  attr(dist.matrix, "metadata") <- meta.dat[labels, , drop = FALSE]
  dist.matrix
}


#' @keywords internal
mStat_extract_dist_metadata <- function(dist.obj,
                                        dist.name,
                                        vars = NULL,
                                        data.obj = NULL) {
  mStat_validate_dist_object(dist.obj, dist.name)
  labels <- mStat_get_dist_labels(dist.obj[[dist.name[1]]])

  if (!is.null(data.obj) && !is.null(data.obj$meta.dat)) {
    meta.tab <- as.data.frame(data.obj$meta.dat, stringsAsFactors = FALSE)
  } else {
    meta.tab <- attr(dist.obj[[dist.name[1]]], "metadata")
  }

  if (is.null(meta.tab)) {
    stop(
      "Metadata is required when using a precomputed dist.obj without data.obj. ",
      "Pass data.obj explicitly or use a dist.obj generated by MicrobiomeStat so metadata is attached.",
      call. = FALSE
    )
  }

  if (is.null(rownames(meta.tab))) {
    stop("Distance metadata must have sample row names.", call. = FALSE)
  }

  missing_rows <- setdiff(labels, rownames(meta.tab))
  if (length(missing_rows) != 0) {
    stop(
      "Metadata is missing the following distance labels: ",
      paste(missing_rows, collapse = ", "),
      call. = FALSE
    )
  }

  meta.tab <- meta.tab[labels, , drop = FALSE]
  if (is.null(vars)) {
    return(meta.tab)
  }
  mStat_select_metadata_columns(meta.tab, vars)
}


#' @keywords internal
mStat_dist_to_tibble <- function(dist.matrix, sample_col = "sample") {
  tibble::rownames_to_column(
    as.data.frame(as.matrix(dist.matrix), stringsAsFactors = FALSE),
    var = sample_col
  )
}


#' @keywords internal
mStat_align_dist_matrix_to_metadata <- function(dist.matrix,
                                                meta.dat) {
  dist.mat <- tryCatch(
    as.matrix(dist.matrix),
    error = function(e) {
      stop(
        "Failed to convert distance object to a matrix: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  meta.tbl <- mStat_meta_with_sample(meta.dat, sample_col = "sample")
  sample_ids <- meta.tbl[["sample"]]
  if (anyNA(sample_ids) || anyDuplicated(sample_ids)) {
    stop("Metadata must have unique, non-missing sample identifiers.", call. = FALSE)
  }

  missing_rows <- setdiff(sample_ids, rownames(dist.mat))
  missing_cols <- setdiff(sample_ids, colnames(dist.mat))
  if (length(missing_rows) > 0 || length(missing_cols) > 0) {
    stop(
      "Distance matrix labels do not match metadata rows.",
      if (length(missing_rows) > 0) {
        paste0(" Missing rows: ", paste(missing_rows, collapse = ", "), ".")
      } else {
        ""
      },
      if (length(missing_cols) > 0) {
        paste0(" Missing columns: ", paste(missing_cols, collapse = ", "), ".")
      } else {
        ""
      },
      call. = FALSE
    )
  }

  dist.mat[sample_ids, sample_ids, drop = FALSE]
}


#' @keywords internal
mStat_prepare_subject_time_distance_long_data <- function(dist.matrix,
                                                          meta.dat,
                                                          subject.var,
                                                          time.var,
                                                          mode = c("baseline", "adjacent"),
                                                          baseline_time = NULL,
                                                          followup_time = NULL,
                                                          sample_col = "sample") {
  mode <- match.arg(mode)

  meta.tbl <- mStat_meta_to_tibble(meta.dat, sample_col = sample_col)
  dist.mat <- mStat_align_dist_matrix_to_metadata(dist.matrix, meta.dat)

  if (!all(c(subject.var, time.var) %in% colnames(meta.tbl))) {
    stop("`subject.var` and `time.var` must both be present in metadata.", call. = FALSE)
  }

  subject_groups <- split(meta.tbl, meta.tbl[[subject.var]])

  long_parts <- lapply(subject_groups, function(subject_meta) {
    if (nrow(subject_meta) == 0) {
      return(NULL)
    }

    ordered_times <- mStat_order_time_labels(subject_meta[[time.var]])
    sample_groups <- split(
      as.character(subject_meta[[sample_col]]),
      as.character(subject_meta[[time.var]])
    )
    subject_id <- subject_meta[[subject.var]][[1]]

    if (mode == "baseline") {
      resolved_baseline <- if (is.null(baseline_time)) {
        ordered_times[[1]]
      } else {
        as.character(baseline_time)
      }

      if (!resolved_baseline %in% ordered_times) {
        return(NULL)
      }

      target_times <- if (is.null(followup_time)) {
        setdiff(ordered_times, resolved_baseline)
      } else {
        intersect(as.character(followup_time), ordered_times)
      }

      if (length(target_times) == 0) {
        return(NULL)
      }

      baseline_samples <- sample_groups[[resolved_baseline]]

      dplyr::bind_rows(lapply(target_times, function(target_time) {
        target_samples <- sample_groups[[target_time]]
        distance_values <- as.vector(dist.mat[target_samples, baseline_samples, drop = FALSE])

        out <- data.frame(
          rep(subject_id, length(distance_values)),
          rep(target_time, length(distance_values)),
          distance = distance_values,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        colnames(out)[1:2] <- c(subject.var, time.var)
        out
      }))
    } else {
      if (length(ordered_times) < 2) {
        return(NULL)
      }

      dplyr::bind_rows(lapply(seq_len(length(ordered_times) - 1L), function(i) {
        time_before <- ordered_times[[i]]
        time_after <- ordered_times[[i + 1L]]
        before_samples <- sample_groups[[time_before]]
        after_samples <- sample_groups[[time_after]]
        distance_values <- as.vector(dist.mat[after_samples, before_samples, drop = FALSE])

        out <- data.frame(
          rep(subject_id, length(distance_values)),
          rep(time_after, length(distance_values)),
          rep(time_before, length(distance_values)),
          distance = distance_values,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        colnames(out)[1:3] <- c(subject.var, time.var, paste0(time.var, ".before"))
        out
      }))
    }
  })

  dplyr::bind_rows(long_parts)
}


#' @keywords internal
mStat_prepare_dist_group_time_long_data <- function(dist.matrix,
                                                    meta.dat,
                                                    group.var,
                                                    time.var,
                                                    strata.var = NULL,
                                                    sample_col = "sample",
                                                    pair_col = "sample2",
                                                    distance_col = "Distance") {
  meta.tbl <- mStat_meta_to_tibble(meta.dat, sample_col = sample_col)
  dist.mat <- mStat_align_dist_matrix_to_metadata(dist.matrix, meta.dat)

  dist_meta_lookup <- meta.tbl %>%
    dplyr::select(
      dplyr::all_of(sample_col),
      Group = dplyr::all_of(group.var),
      Time = dplyr::all_of(time.var),
      dplyr::any_of(strata.var)
    )

  if (!is.null(strata.var)) {
    colnames(dist_meta_lookup)[colnames(dist_meta_lookup) == strata.var] <- "Strata"
  }

  pair_index <- which(upper.tri(dist.mat), arr.ind = TRUE)
  pair_df <- data.frame(
    sample = rownames(dist.mat)[pair_index[, 1]],
    pair = colnames(dist.mat)[pair_index[, 2]],
    distance = dist.mat[pair_index],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(pair_df) <- c(sample_col, pair_col, distance_col)

  pair_df %>%
    dplyr::left_join(dist_meta_lookup, by = sample_col) %>%
    dplyr::left_join(
      dist_meta_lookup,
      by = stats::setNames(sample_col, pair_col),
      suffix = c(".x", ".y")
    )
}


#' @keywords internal
mStat_prepare_beta_change_long_data <- function(dist.matrix,
                                                meta.dat,
                                                subject.var,
                                                time.var,
                                                change.base,
                                                change.after = NULL,
                                                sample_col = "sample") {
  mStat_prepare_subject_time_distance_long_data(
    dist.matrix = dist.matrix,
    meta.dat = meta.dat,
    subject.var = subject.var,
    time.var = time.var,
    mode = "baseline",
    baseline_time = change.base,
    followup_time = change.after,
    sample_col = sample_col
  )
}


#' @keywords internal
mStat_attach_change_metadata <- function(change.df,
                                         meta.dat,
                                         by,
                                         vars = NULL,
                                         sample_col = "sample") {
  if (length(vars) == 0 || is.null(vars)) {
    return(change.df)
  }

  meta.tbl <- mStat_meta_to_tibble(meta.dat, sample_col = sample_col)
  keep_cols <- unique(c(by, vars))

  missing_cols <- setdiff(keep_cols, colnames(meta.tbl))
  if (length(missing_cols) != 0) {
    stop(
      "Change metadata is missing the following columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  meta.tbl <- meta.tbl %>%
    dplyr::select(-dplyr::any_of(sample_col), dplyr::all_of(keep_cols)) %>%
    dplyr::distinct()

  dplyr::left_join(
    change.df,
    meta.tbl,
    by = by,
    relationship = "many-to-one"
  )
}


#' @keywords internal
mStat_prepare_subject_metadata <- function(meta.dat,
                                           subject.var,
                                           vars = NULL,
                                           time.var = NULL,
                                           time.value = NULL,
                                           sample_col = "sample") {
  if (length(vars) == 0 || is.null(vars)) {
    return(NULL)
  }

  meta.tbl <- mStat_meta_to_tibble(meta.dat, sample_col = sample_col)
  keep_cols <- unique(c(subject.var, vars, time.var))

  missing_cols <- setdiff(keep_cols, colnames(meta.tbl))
  if (length(missing_cols) != 0) {
    stop(
      "Subject metadata is missing the following columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(time.var) && !is.null(time.value)) {
    meta.tbl <- meta.tbl[
      mStat_match_metadata_values(meta.tbl[[time.var]], time.value),
      ,
      drop = FALSE
    ]
  }

  selected_cols <- unique(c(subject.var, vars))

  meta.tbl %>%
    dplyr::select(dplyr::all_of(selected_cols)) %>%
    dplyr::distinct(.data[[subject.var]], .keep_all = TRUE)
}


#' @keywords internal
mStat_canonicalize_group_pair_labels <- function(group_x,
                                                 group_y,
                                                 sep = "-") {
  if (length(group_x) != length(group_y)) {
    stop("`group_x` and `group_y` must have the same length.", call. = FALSE)
  }

  vapply(seq_along(group_x), function(i) {
    paste(sort(c(as.character(group_x[[i]]), as.character(group_y[[i]]))), collapse = sep)
  }, character(1))
}


#' @keywords internal
mStat_prepare_pc_pair_change_data <- function(long.df,
                                              subject.var,
                                              time.var,
                                              change.base = NULL,
                                              change.func,
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
      by = c("PC", subject.var),
      suffix = c("_time_1", "_time_2"),
      relationship = "one-to-one"
    ) %>%
    dplyr::mutate(
      value_diff = compute_taxa_change(
        value_after = value_time_2,
        value_before = value_time_1,
        method = change.func,
        verbose = FALSE
      )
    )

  list(
    combined_data = combined_data,
    change.base = pair_slices$change.base,
    change.after = pair_slices$change.after
  )
}


#' @keywords internal
mStat_prepare_beta_trend_long_data <- function(dist.matrix,
                                               meta.dat,
                                               subject.var,
                                               time.var,
                                               vars = NULL,
                                               sample_col = "sample") {
  long.df <- mStat_prepare_subject_time_distance_long_data(
    dist.matrix = dist.matrix,
    meta.dat = meta.dat,
    subject.var = subject.var,
    time.var = time.var,
    mode = "baseline",
    sample_col = sample_col
  )

  mStat_attach_change_metadata(
    change.df = long.df,
    meta.dat = meta.dat,
    by = c(subject.var, time.var),
    vars = vars,
    sample_col = sample_col
  )
}


#' @keywords internal
mStat_prepare_beta_adjacent_long_data <- function(dist.matrix,
                                                  meta.dat,
                                                  subject.var,
                                                  time.var,
                                                  vars = NULL,
                                                  sample_col = "sample") {
  long.df <- mStat_prepare_subject_time_distance_long_data(
    dist.matrix = dist.matrix,
    meta.dat = meta.dat,
    subject.var = subject.var,
    time.var = time.var,
    mode = "adjacent",
    sample_col = sample_col
  )

  mStat_attach_change_metadata(
    change.df = long.df,
    meta.dat = meta.dat,
    by = c(subject.var, time.var),
    vars = vars,
    sample_col = sample_col
  )
}


#' @keywords internal
mStat_meta_to_tibble <- function(meta.dat, sample_col = "sample") {
  meta.df <- as.data.frame(
    meta.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  row_ids <- rownames(meta.df)
  has_explicit_rownames <- !is.null(row_ids) && !identical(row_ids, as.character(seq_len(nrow(meta.df))))

  if (!sample_col %in% colnames(meta.df) || has_explicit_rownames) {
    meta.df[[sample_col]] <- row_ids
  }

  meta.df <- meta.df[, c(sample_col, setdiff(colnames(meta.df), sample_col)), drop = FALSE]
  rownames(meta.df) <- NULL

  tibble::as_tibble(meta.df)
}


#' @keywords internal
mStat_pc_to_tibble <- function(pc.points, sample_col = "sample") {
  pc.points <- as.matrix(pc.points)
  sample_ids <- rownames(pc.points)

  if (is.null(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) {
    stop("PC coordinates must have unique, non-missing sample row names.", call. = FALSE)
  }

  pc.mat <- pc.points
  colnames(pc.mat) <- paste0("PC", seq_len(ncol(pc.mat)))

  tibble::rownames_to_column(
    as.data.frame(pc.mat, stringsAsFactors = FALSE),
    var = sample_col
  )
}


#' @keywords internal
mStat_prepare_pc_long_data <- function(pc.points,
                                       pc.ind,
                                       meta.dat = NULL,
                                       vars = NULL,
                                       sample_col = "sample",
                                       join = c("inner", "left"),
                                       pc_col = "PC",
                                       value_col = "value") {
  join <- match.arg(join)

  pc.tbl <- mStat_pc_to_tibble(pc.points, sample_col = sample_col)
  pc_cols <- paste0("PC", pc.ind)
  missing_cols <- setdiff(pc_cols, colnames(pc.tbl))
  if (length(missing_cols) != 0) {
    stop(
      "PC coordinates are missing the following requested axes: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  pc.tbl <- pc.tbl[, c(sample_col, pc_cols), drop = FALSE]

  if (!is.null(meta.dat)) {
    meta.tbl <- mStat_meta_to_tibble(
      mStat_select_metadata_columns(meta.dat, vars),
      sample_col = sample_col
    )
    join_fun <- if (identical(join, "inner")) dplyr::inner_join else dplyr::left_join
    pc.tbl <- join_fun(pc.tbl, meta.tbl, by = sample_col)
  }

  tidyr::pivot_longer(
    pc.tbl,
    cols = dplyr::all_of(pc_cols),
    names_to = pc_col,
    values_to = value_col
  )
}


#' @keywords internal
mStat_subset_beta_objects <- function(dist.obj, pc.obj = NULL, samIDs = NULL) {
  if (is.null(samIDs)) {
    return(list(dist.obj = dist.obj, pc.obj = pc.obj))
  }

  if (!is.null(dist.obj)) {
    dist.obj <- mStat_subset_dist(dist.obj, samIDs)
  }
  if (!is.null(pc.obj)) {
    pc.obj <- mStat_subset_PC(pc.obj, samIDs)
  }

  list(
    dist.obj = dist.obj,
    pc.obj = pc.obj
  )
}


#' @keywords internal
mStat_prepare_precomputed_beta_context <- function(dist.obj,
                                                   dist.name,
                                                   pc.obj = NULL,
                                                   data.obj = NULL,
                                                   time.var = NULL,
                                                   t.level = NULL,
                                                   t0.level = NULL,
                                                   ts.levels = NULL,
                                                   process_time = FALSE,
                                                   required_pc_axes = NULL) {
  mStat_validate_dist_object(dist.obj, dist.name)

  if (!is.null(pc.obj)) {
    mStat_validate_pc_object(
      pc.obj = pc.obj,
      dist.name = dist.name,
      required_pc_axes = required_pc_axes
    )
    pc.obj <- mStat_align_pc_object(pc.obj, dist.obj, dist.name)
    mStat_validate_pc_object(
      pc.obj = pc.obj,
      dist.name = dist.name,
      dist.obj = dist.obj,
      required_pc_axes = required_pc_axes
    )
  }

  if (!is.null(data.obj) && !is.null(data.obj$meta.dat)) {
    context.obj <- data.obj
  } else {
    context.obj <- list(meta.dat = mStat_extract_dist_metadata(dist.obj, dist.name))
  }

  aligned_meta <- mStat_extract_dist_metadata(
    dist.obj = dist.obj,
    dist.name = dist.name,
    data.obj = context.obj
  )
  if (!is.null(context.obj$feature.tab)) {
    context.obj <- mStat_subset_data(context.obj, samIDs = rownames(aligned_meta))
  } else {
    context.obj$meta.dat <- aligned_meta
  }

  if (process_time) {
    context.obj <- mStat_process_time_variable(context.obj, time.var, t0.level, ts.levels)
  } else if (!is.null(time.var) && !is.null(t.level)) {
    samIDs <- get_sample_ids(context.obj, time.var, t.level)
    if (!is.null(context.obj$feature.tab)) {
      context.obj <- mStat_subset_data(context.obj, samIDs = samIDs)
    } else {
      context.obj$meta.dat <- context.obj$meta.dat[samIDs, , drop = FALSE]
    }
  }

  samIDs <- rownames(context.obj$meta.dat)
  aligned <- mStat_subset_beta_objects(dist.obj, pc.obj, samIDs)

  mStat_validate_dist_object(aligned$dist.obj, dist.name)
  if (!is.null(aligned$pc.obj)) {
    mStat_validate_pc_object(
      pc.obj = aligned$pc.obj,
      dist.name = dist.name,
      dist.obj = aligned$dist.obj,
      required_pc_axes = required_pc_axes,
      sample_ids = samIDs
    )
  }

  list(
    data.obj = context.obj,
    dist.obj = aligned$dist.obj,
    pc.obj = aligned$pc.obj
  )
}
