#' @title Alpha Preparation Utilities
#' @description Internal helpers for preparing alpha diversity inputs.
#'   These utilities centralize rarefaction, time-based subsetting, alpha
#'   calculation, and validation so all alpha workflows share one consistent
#'   preparation path.
#' @name alpha_preparation_utils
#' @keywords internal
NULL


#' @keywords internal
mStat_validate_alpha_object <- function(alpha.obj,
                                        alpha.name,
                                        sample_ids = NULL) {
  if (is.null(alpha.obj) || !is.list(alpha.obj) || length(alpha.obj) == 0) {
    stop("alpha.obj must be a non-empty list of alpha diversity tables.", call. = FALSE)
  }

  missing_alphas <- setdiff(alpha.name, names(alpha.obj))
  if (length(missing_alphas) != 0) {
    stop(
      "The following alpha diversity indices are not available in alpha.obj: ",
      paste(missing_alphas, collapse = ", "),
      call. = FALSE
    )
  }

  for (alpha_key in alpha.name) {
    alpha_entry <- alpha.obj[[alpha_key]]
    if (is.null(alpha_entry) || !(is.data.frame(alpha_entry) || is.matrix(alpha_entry))) {
      stop(
        "alpha.obj[['", alpha_key, "']] must be a data frame or matrix.",
        call. = FALSE
      )
    }

    alpha_entry <- as.data.frame(alpha_entry, stringsAsFactors = FALSE)
    sample_names <- rownames(alpha_entry)

    if (is.null(sample_names) || anyNA(sample_names) || anyDuplicated(sample_names)) {
      stop(
        "alpha.obj[['", alpha_key, "']] must have unique, non-missing sample row names.",
        call. = FALSE
      )
    }

    if (!is.null(sample_ids)) {
      missing_samples <- setdiff(sample_ids, sample_names)
      if (length(missing_samples) != 0) {
        stop(
          "alpha.obj[['", alpha_key, "']] is missing the following samples: ",
          paste(missing_samples, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  invisible(alpha.obj)
}


#' @keywords internal
mStat_alpha_to_tibble <- function(alpha.obj, sample_col = "sample") {
  alpha_names <- names(alpha.obj)
  if (is.null(alpha_names) || anyNA(alpha_names) || anyDuplicated(alpha_names)) {
    stop("alpha.obj must be a named list of alpha diversity tables.", call. = FALSE)
  }

  mStat_validate_alpha_object(alpha.obj = alpha.obj, alpha.name = alpha_names)

  sample_ids <- rownames(as.data.frame(alpha.obj[[alpha_names[[1]]]], stringsAsFactors = FALSE))
  alpha_cols <- lapply(alpha_names, function(alpha_key) {
    alpha_entry <- as.data.frame(alpha.obj[[alpha_key]], stringsAsFactors = FALSE)

    if (!setequal(rownames(alpha_entry), sample_ids)) {
      stop(
        "alpha.obj entries must contain the same sample set. Entry '",
        alpha_key,
        "' does not match the first alpha table.",
        call. = FALSE
      )
    }

    alpha_entry <- alpha_entry[sample_ids, , drop = FALSE]

    if (alpha_key %in% colnames(alpha_entry)) {
      alpha_values <- alpha_entry[[alpha_key]]
    } else if (ncol(alpha_entry) == 1) {
      alpha_values <- alpha_entry[[1]]
    } else {
      stop(
        "alpha.obj[['", alpha_key, "']] must contain a column named '",
        alpha_key,
        "' or exactly one column.",
        call. = FALSE
      )
    }

    stats::setNames(list(unname(alpha_values)), alpha_key)
  })

  alpha_df <- as.data.frame(alpha_cols, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(alpha_df) <- sample_ids

  tibble::rownames_to_column(alpha_df, var = sample_col)
}


#' @keywords internal
mStat_prepare_alpha_data <- function(alpha.obj,
                                     meta.dat = NULL,
                                     vars = NULL,
                                     sample_col = "sample",
                                     join = c("inner", "left")) {
  join <- match.arg(join)

  alpha_df <- mStat_alpha_to_tibble(alpha.obj = alpha.obj, sample_col = sample_col)

  if (is.null(meta.dat)) {
    return(alpha_df)
  }

  meta_subset <- if (is.null(vars)) {
    as.data.frame(meta.dat, stringsAsFactors = FALSE)
  } else {
    mStat_select_metadata_columns(meta.dat, vars)
  }

  meta.tbl <- mStat_meta_to_tibble(
    meta_subset,
    sample_col = sample_col
  )

  join_fun <- if (identical(join, "inner")) dplyr::inner_join else dplyr::left_join
  join_fun(alpha_df, meta.tbl, by = sample_col)
}


#' @keywords internal
mStat_prepare_alpha_inputs <- function(data.obj,
                                       alpha.obj = NULL,
                                       alpha.name,
                                       depth = NULL,
                                       time.var = NULL,
                                       t.level = NULL,
                                       t0.level = NULL,
                                       ts.levels = NULL,
                                       process_time = FALSE) {
  if (process_time) {
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
  }

  if (!is.null(time.var) && !is.null(t.level)) {
    subset.ids <- get_sample_ids(data.obj, time.var, t.level)
    data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)
  }

  if (is.null(alpha.obj)) {
    if (!is.null(depth)) {
      message(
        "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj, depth = depth)
    }

    tree <- NULL
    if ("faith_pd" %in% alpha.name) {
      tree <- data.obj$tree
    }

    alpha.obj <- mStat_calculate_alpha_diversity(
      x = data.obj$feature.tab,
      alpha.name = alpha.name,
      tree = tree
    )
  } else {
    mStat_validate_alpha_object(
      alpha.obj = alpha.obj,
      alpha.name = alpha.name,
      sample_ids = rownames(data.obj$meta.dat)
    )
    alpha.obj <- mStat_subset_alpha(alpha.obj, rownames(data.obj$meta.dat))
  }

  list(
    data.obj = data.obj,
    alpha.obj = alpha.obj
  )
}
