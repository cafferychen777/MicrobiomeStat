#' @title Time Utilities
#' @description Internal helpers for ordering and coercing time variables
#'   consistently across plotting, testing, and interface layers.
#' @name time_utils
#' @keywords internal
NULL


#' @keywords internal
mStat_order_time_values <- function(values) {
  if (is.null(values)) {
    return(values)
  }

  values <- values[!is.na(values)]
  if (length(values) == 0) {
    return(values)
  }

  if (is.factor(values)) {
    values_chr <- as.character(values)
    levels_in_use <- levels(values)
    return(levels_in_use[levels_in_use %in% unique(values_chr)])
  }

  if (is.numeric(values)) {
    return(sort(unique(values)))
  }

  values_chr <- unique(as.character(values))

  if (all(grepl("^[-+]?\\d+(\\.\\d+)?$", values_chr))) {
    ord <- order(as.numeric(values_chr), seq_along(values_chr))
    return(values_chr[ord])
  }

  natural_match <- regexec(
    "^([^0-9]*)([-+]?\\d+(?:\\.\\d+)?)([^0-9]*)$",
    values_chr,
    perl = TRUE
  )
  natural_parts <- regmatches(values_chr, natural_match)
  if (all(lengths(natural_parts) == 4)) {
    prefixes <- vapply(natural_parts, `[`, character(1), 2)
    numbers <- as.numeric(vapply(natural_parts, `[`, character(1), 3))
    suffixes <- vapply(natural_parts, `[`, character(1), 4)
    ord <- order(prefixes, numbers, suffixes, seq_along(values_chr))
    return(values_chr[ord])
  }

  values_chr
}


#' @keywords internal
mStat_coerce_time_to_numeric <- function(values,
                                         time.var,
                                         context = "time-based analysis") {
  if (is.numeric(values)) {
    return(values)
  }

  if (is.factor(values)) {
    values_chr <- as.character(values)
  } else if (is.character(values)) {
    values_chr <- values
  } else {
    stop("`", time.var, "` must be numeric, factor, or character.", call. = FALSE)
  }

  numeric_values <- suppressWarnings(as.numeric(values_chr))
  bad_idx <- !is.na(values_chr) & is.na(numeric_values)
  if (any(bad_idx)) {
    bad_values <- unique(values_chr[bad_idx])
    stop(
      "The `", time.var, "` variable must be numeric or coercible to numeric for ",
      context,
      ". Non-numeric values found: ",
      paste(bad_values, collapse = ", "),
      call. = FALSE
    )
  }

  numeric_values
}


#' @keywords internal
mStat_has_variation <- function(values) {
  non_missing <- stats::na.omit(values)

  if (length(non_missing) == 0) {
    return(FALSE)
  }

  if (is.numeric(non_missing)) {
    return(stats::var(non_missing) > 0)
  }

  length(unique(non_missing)) > 1
}


#' @keywords internal
mStat_validate_time_var_contract <- function(meta.dat,
                                             time.var,
                                             context = "analysis",
                                             require_variation = FALSE) {
  if (is.null(time.var) || !nzchar(time.var)) {
    stop("`time.var` is required for ", context, ".", call. = FALSE)
  }

  if (!time.var %in% names(meta.dat)) {
    stop("`time.var` ('", time.var, "') was not found in metadata for ", context, ".", call. = FALSE)
  }

  observed_time <- stats::na.omit(meta.dat[[time.var]])
  if (length(observed_time) == 0) {
    stop("`time.var` ('", time.var, "') has no observed values for ", context, ".", call. = FALSE)
  }

  if (require_variation && !mStat_has_variation(observed_time)) {
    stop(
      "`time.var` ('", time.var, "') must have at least two observed levels/values for ",
      context,
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' @keywords internal
mStat_validate_group_var_contract <- function(meta.dat,
                                              group.var,
                                              subject.var = NULL,
                                              context = "analysis",
                                              require_variation = TRUE) {
  if (is.null(group.var) || !nzchar(group.var)) {
    stop("`group.var` is required for ", context, ".", call. = FALSE)
  }

  if (!group.var %in% names(meta.dat)) {
    stop("`group.var` ('", group.var, "') was not found in metadata for ", context, ".", call. = FALSE)
  }

  if (require_variation && !mStat_has_variation(meta.dat[[group.var]])) {
    stop(
      "`group.var` ('", group.var, "') must have at least two observed levels/values for ",
      context,
      ".",
      call. = FALSE
    )
  }

  if (!is.null(subject.var)) {
    if (!subject.var %in% names(meta.dat)) {
      stop("`subject.var` ('", subject.var, "') was not found in metadata for ", context, ".", call. = FALSE)
    }

    inconsistent_subjects <- meta.dat %>%
      dplyr::select(all_of(c(subject.var, group.var))) %>%
      dplyr::distinct() %>%
      dplyr::count(!!rlang::sym(subject.var), name = "n_groups") %>%
      dplyr::filter(.data$n_groups > 1) %>%
      dplyr::pull(!!rlang::sym(subject.var))

    if (length(inconsistent_subjects) > 0) {
      stop(
        "`group.var` must be constant within each subject for ",
        context,
        ". Inconsistent subjects: ",
        paste(utils::head(inconsistent_subjects, 10), collapse = ", "),
        if (length(inconsistent_subjects) > 10) " ..." else "",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' @keywords internal
mStat_resolve_pair_timepoints <- function(values,
                                          time.var,
                                          change.base = NULL,
                                          context = "paired change analysis") {
  ordered_levels <- mStat_order_time_values(values)
  ordered_labels <- as.character(ordered_levels)

  if (length(ordered_levels) < 2) {
    stop(
      "The `", time.var, "` variable must contain at least two time levels for ",
      context,
      ".",
      call. = FALSE
    )
  }

  if (is.null(change.base)) {
    change.base_label <- ordered_labels[[1]]
    message(
      "The 'change.base' variable was NULL. It has been set to the first ordered value in the '",
      time.var,
      "' variable: ",
      change.base_label
    )
  } else {
    change.base_label <- as.character(change.base)
  }

  if (!change.base_label %in% ordered_labels) {
    stop(
      "The requested change.base ('",
      change.base_label,
      "') is not present in `",
      time.var,
      "`.",
      call. = FALSE
    )
  }

  change.after_labels <- ordered_labels[ordered_labels != change.base_label]
  if (length(change.after_labels) != 1) {
    stop(
      "The `", time.var, "` variable must contain exactly one follow-up level for ",
      context,
      ". Found baseline '", change.base_label, "' and follow-up level(s): ",
      paste(change.after_labels, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    change.base = change.base_label,
    change.after = change.after_labels[[1]],
    levels = ordered_labels
  )
}
