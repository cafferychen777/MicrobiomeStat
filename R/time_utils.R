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

  if (inherits(values, "Date") || inherits(values, "POSIXt")) {
    return(sort(unique(values)))
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
mStat_order_time_labels <- function(values) {
  as.character(mStat_order_time_values(values))
}


#' @keywords internal
mStat_match_metadata_values <- function(values,
                                        targets) {
  if (is.null(targets)) {
    return(rep(TRUE, length(values)))
  }

  if (length(values) == 0) {
    return(logical())
  }

  target_values <- targets[!is.na(targets)]
  if (length(target_values) == 0) {
    return(rep(FALSE, length(values)))
  }

  if (inherits(values, "Date")) {
    if (!inherits(target_values, "Date")) {
      target_values <- suppressWarnings(
        as.Date(target_values, tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))
      )
    }
    target_values <- target_values[!is.na(target_values)]
    return(!is.na(values) & values %in% target_values)
  }

  if (inherits(values, "POSIXt")) {
    if (!inherits(target_values, "POSIXt")) {
      target_values <- suppressWarnings(
        as.POSIXct(
          target_values,
          tz = "UTC",
          tryFormats = c(
            "%Y-%m-%d %H:%M:%OS",
            "%Y-%m-%d %H:%M",
            "%Y-%m-%dT%H:%M:%OS",
            "%Y-%m-%dT%H:%M",
            "%Y/%m/%d %H:%M:%OS",
            "%Y/%m/%d %H:%M"
          )
        )
      )
    }
    target_values <- target_values[!is.na(target_values)]
    return(!is.na(values) & values %in% target_values)
  }

  if (is.numeric(values)) {
    target_values <- suppressWarnings(as.numeric(target_values))
    target_values <- target_values[!is.na(target_values)]
    return(!is.na(values) & values %in% target_values)
  }

  as.character(values) %in% as.character(target_values)
}


#' @keywords internal
mStat_parse_temporal_labels_to_numeric <- function(values_chr) {
  non_missing <- !is.na(values_chr)

  if (!any(non_missing)) {
    return(rep(NA_real_, length(values_chr)))
  }

  date_values <- tryCatch(
    suppressWarnings(
      as.Date(
        values_chr[non_missing],
        tryFormats = c("%Y-%m-%d", "%Y/%m/%d")
      )
    ),
    error = function(...) rep(as.Date(NA), sum(non_missing))
  )
  if (all(!is.na(date_values))) {
    parsed_values <- rep(NA_real_, length(values_chr))
    parsed_values[non_missing] <- as.numeric(date_values)
    return(parsed_values)
  }

  datetime_values <- tryCatch(
    suppressWarnings(
      as.POSIXct(
        values_chr[non_missing],
        tz = "UTC",
        tryFormats = c(
          "%Y-%m-%d %H:%M:%OS",
          "%Y-%m-%d %H:%M",
          "%Y-%m-%dT%H:%M:%OS",
          "%Y-%m-%dT%H:%M",
          "%Y/%m/%d %H:%M:%OS",
          "%Y/%m/%d %H:%M"
        )
      )
    )
    ,
    error = function(...) rep(as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"), sum(non_missing))
  )
  if (all(!is.na(datetime_values))) {
    parsed_values <- rep(NA_real_, length(values_chr))
    parsed_values[non_missing] <- as.numeric(datetime_values)
    return(parsed_values)
  }

  NULL
}


#' @keywords internal
mStat_format_dropped_time_levels_message <- function(dropped_levels,
                                                     time.var,
                                                     context = "time processing") {
  if (length(dropped_levels) == 0) {
    return(NULL)
  }

  paste0(
    "Subsetting `", time.var, "` for ", context, " excluded these observed levels: ",
    paste(dropped_levels, collapse = ", ")
  )
}


#' @keywords internal
mStat_resolve_time_levels <- function(values,
                                      time.var,
                                      t0.level = NULL,
                                      ts.levels = NULL,
                                      requested_levels = NULL,
                                      context = "time processing") {
  observed_levels <- mStat_order_time_labels(values)

  if (length(observed_levels) == 0) {
    stop(
      "`", time.var, "` must contain at least one observed value for ",
      context,
      ".",
      call. = FALSE
    )
  }

  explicit_levels <- unique(c(
    if (!is.null(requested_levels)) as.character(requested_levels) else character(),
    if (is.null(requested_levels) && !is.null(t0.level)) as.character(t0.level) else character(),
    if (is.null(requested_levels) && !is.null(ts.levels)) as.character(ts.levels) else character()
  ))

  missing_levels <- setdiff(explicit_levels, observed_levels)
  if (length(missing_levels) != 0) {
    stop(
      "Requested time levels not found in `", time.var, "` for ",
      context,
      ": ",
      paste(missing_levels, collapse = ", "),
      ". Available values: ",
      paste(observed_levels, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(requested_levels)) {
    kept_levels <- unique(as.character(requested_levels))
    resolved_t0 <- if (length(kept_levels) > 0) kept_levels[[1]] else NULL
    resolved_ts <- if (length(kept_levels) > 1) kept_levels[-1] else character()
  } else {
    resolved_t0 <- if (!is.null(t0.level)) {
      as.character(t0.level)[[1]]
    } else {
      observed_levels[[1]]
    }

    resolved_ts <- if (!is.null(ts.levels)) {
      unique(as.character(ts.levels))
    } else {
      observed_levels[observed_levels != resolved_t0]
    }

    kept_levels <- unique(c(resolved_t0, resolved_ts))
  }

  dropped_levels <- setdiff(observed_levels, kept_levels)

  list(
    observed_levels = observed_levels,
    kept_levels = kept_levels,
    dropped_levels = dropped_levels,
    dropped_levels_message = mStat_format_dropped_time_levels_message(
      dropped_levels = dropped_levels,
      time.var = time.var,
      context = context
    ),
    t0.level = resolved_t0,
    ts.levels = resolved_ts
  )
}


#' @keywords internal
mStat_inform_numeric_time_requirement <- function(function_name,
                                                  analysis_label,
                                                  conversion_behavior = c("coerce", "require_numeric", "preprocess")) {
  conversion_behavior <- match.arg(conversion_behavior)

  action_text <- switch(
    conversion_behavior,
    coerce = "It will be coerced to numeric inside the function when possible.",
    require_numeric = "Provide a numeric time variable for this analysis.",
    preprocess = "If needed, convert it in metadata before calling this function."
  )

  message(
    "`", function_name, "` uses a numeric time scale for ", analysis_label, ".\n",
    "Non-numeric values will fail if they cannot be coerced consistently.\n",
    action_text
  )
}


#' @keywords internal
mStat_coerce_time_to_numeric <- function(values,
                                         time.var,
                                         context = "time-based analysis") {
  if (inherits(values, "Date") || inherits(values, "POSIXt")) {
    return(as.numeric(values))
  }

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

  if (!any(bad_idx)) {
    return(numeric_values)
  }

  temporal_numeric <- mStat_parse_temporal_labels_to_numeric(values_chr)
  if (!is.null(temporal_numeric)) {
    return(temporal_numeric)
  }

  bad_values <- unique(values_chr[bad_idx])
  stop(
    "The `", time.var, "` variable must be numeric or coercible to numeric for ",
    context,
    ". Non-numeric values found: ",
    paste(bad_values, collapse = ", "),
    call. = FALSE
  )
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
mStat_resolve_variable_terms <- function(data,
                                         terms = NULL) {
  model_terms <- unique(unlist(terms, use.names = FALSE))
  model_terms <- model_terms[!is.na(model_terms)]

  if (length(model_terms) == 0) {
    return(character())
  }

  model_terms[vapply(model_terms, function(term) {
    mStat_has_variation(data[[term]])
  }, logical(1))]
}


#' @keywords internal
mStat_build_formula <- function(response,
                                terms = NULL) {
  model_terms <- unique(unlist(terms, use.names = FALSE))
  model_terms <- model_terms[!is.na(model_terms)]

  if (length(model_terms) == 0) {
    return(stats::as.formula(paste(response, "~ 1")))
  }

  stats::as.formula(paste(response, "~", paste(model_terms, collapse = " + ")))
}


#' @keywords internal
mStat_build_mixed_effects_formula <- function(response.var,
                                              time.var,
                                              group.var = NULL,
                                              subject.var,
                                              adj.vars = NULL,
                                              random_slopes = TRUE) {
  fixed_formula <- mStat_build_formula(
    response = response.var,
    terms = if (!is.null(group.var)) {
      c(adj.vars, paste0(group.var, " * ", time.var))
    } else {
      c(adj.vars, time.var)
    }
  )

  random_effects <- if (random_slopes) {
    paste0("(1 + ", time.var, " | ", subject.var, ")")
  } else {
    paste0("(1 | ", subject.var, ")")
  }

  stats::as.formula(paste(deparse(fixed_formula), random_effects, sep = " + "))
}


#' @keywords internal
create_mixed_effects_formula <- mStat_build_mixed_effects_formula
construct_formula <- function(index,
                              group.var = NULL,
                              time.var,
                              subject.var,
                              adj.vars = NULL,
                              random_slopes = TRUE) {
  fixed_effects <- if (is.null(group.var)) {
    time.var
  } else {
    paste0(group.var, " * ", time.var)
  }

  random_effects <- if (random_slopes) {
    paste0("(1 + ", time.var, " | ", subject.var, ")")
  } else {
    paste0("(1 | ", subject.var, ")")
  }

  terms <- c(fixed_effects, random_effects, adj.vars)
  stats::as.formula(paste(index, "~", paste(terms[!is.na(terms) & nzchar(terms)], collapse = " + ")))
}


#' @keywords internal
mStat_fit_mixed_effects_model <- function(response.var,
                                          time.var,
                                          group.var,
                                          subject.var,
                                          data,
                                          adj.vars = NULL,
                                          ...,
                                          context = "mixed-effects analysis") {
  formula <- mStat_build_mixed_effects_formula(
    response.var = response.var,
    time.var = time.var,
    group.var = group.var,
    subject.var = subject.var,
    adj.vars = adj.vars,
    random_slopes = TRUE
  )

  model_fit <- try(lmer(formula, data = data, ...), silent = TRUE)
  if (!inherits(model_fit, "try-error")) {
    return(model_fit)
  }

  error_message <- conditionMessage(attr(model_fit, "condition"))
  fallback_pattern <- "number of observations.*<=.*number of random effects|singular|number of levels of each grouping factor"
  if (!grepl(fallback_pattern, error_message, ignore.case = TRUE)) {
    stop("Model fitting failed in ", context, ": ", error_message, call. = FALSE)
  }

  message("Simplifying the random-effects structure due to overparameterization.")
  formula <- mStat_build_mixed_effects_formula(
    response.var = response.var,
    time.var = time.var,
    group.var = group.var,
    subject.var = subject.var,
    adj.vars = adj.vars,
    random_slopes = FALSE
  )

  model_fit <- try(lmer(formula, data = data, ...), silent = TRUE)
  if (!inherits(model_fit, "try-error")) {
    return(model_fit)
  }

  simplified_error <- conditionMessage(attr(model_fit, "condition"))
  if (!grepl(fallback_pattern, simplified_error, ignore.case = TRUE)) {
    stop(
      "Model fitting failed even after simplifying the random-effects structure in ",
      context,
      ": ",
      simplified_error,
      call. = FALSE
    )
  }

  message("Falling back to a fixed-effects model due to insufficient repeated observations.")
  fixed_formula <- mStat_build_formula(
    response = response.var,
    terms = if (!is.null(group.var)) {
      c(adj.vars, paste0(group.var, " * ", time.var))
    } else {
      c(adj.vars, time.var)
    }
  )

  stats::lm(fixed_formula, data = data)
}


#' @keywords internal
mStat_extract_group_anova_row <- function(anova_result, group.var) {
  if (is.null(group.var) || !nzchar(group.var)) {
    return(NULL)
  }

  anova_df <- as.data.frame(anova_result)
  terms <- rownames(anova_df)
  if (is.null(terms) || length(terms) == 0) {
    return(NULL)
  }

  term_index <- which(terms == group.var)
  if (length(term_index) == 0) {
    return(NULL)
  }

  selected <- anova_df[term_index[[1]], , drop = FALSE]
  data.frame(
    Term = group.var,
    Estimate = NA_real_,
    Std.Error = NA_real_,
    Statistic = selected$`F value`,
    P.Value = selected$`Pr(>F)`,
    check.names = FALSE
  )
}


#' @keywords internal
mStat_compute_group_anova <- function(model) {
  if (inherits(model, "merMod")) {
    return(anova(model, type = "III"))
  }

  stats::anova(model)
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
  ordered_labels <- mStat_order_time_labels(values)

  if (length(ordered_labels) < 2) {
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


#' @keywords internal
mStat_prepare_pair_time_slices <- function(df,
                                           time.var,
                                           change.base = NULL,
                                           context = "paired change analysis") {
  pair_times <- mStat_resolve_pair_timepoints(
    values = df[[time.var]],
    time.var = time.var,
    change.base = change.base,
    context = context
  )

  base_rows <- mStat_match_metadata_values(df[[time.var]], pair_times$change.base)
  followup_rows <- mStat_match_metadata_values(df[[time.var]], pair_times$change.after)

  list(
    data_time_1 = df[base_rows, , drop = FALSE],
    data_time_2 = df[followup_rows, , drop = FALSE],
    change.base = pair_times$change.base,
    change.after = pair_times$change.after,
    levels = pair_times$levels
  )
}


#' @keywords internal
mStat_resolve_followup_timepoints <- function(values,
                                              time.var,
                                              t0.level = NULL,
                                              ts.levels = NULL,
                                              context = "longitudinal analysis") {
  resolved_time <- mStat_resolve_time_levels(
    values = values,
    time.var = time.var,
    t0.level = t0.level,
    ts.levels = ts.levels,
    context = context
  )

  if (length(resolved_time$ts.levels) == 0) {
    stop(
      "`", time.var, "` must contain at least one follow-up level for ",
      context,
      ".",
      call. = FALSE
    )
  }

  resolved_time
}


#' @keywords internal
mStat_order_named_time_entries <- function(x) {
  if (length(x) == 0) {
    return(x)
  }

  ordered_names <- mStat_order_time_labels(names(x))
  x[ordered_names[ordered_names %in% names(x)]]
}
