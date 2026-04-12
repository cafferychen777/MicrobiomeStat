#' @title Per-Time Utilities
#' @description Internal helpers for time-sliced analyses that need one
#'   consistent contract for reshaping grouped results and skipping invalid
#'   time points.
#' @name per_time_utils
#' @keywords internal
NULL


#' @keywords internal
mStat_restructure_group_term_results <- function(test.list,
                                                 group.var,
                                                 reference_level) {
  if (length(test.list) == 0) {
    return(list())
  }

  term_vectors <- lapply(test.list, function(df) {
    if (is.null(df) || !("Term" %in% colnames(df))) {
      return(character())
    }

    df$Term
  })

  group_prefix_pattern <- paste0("^", mStat_escape_regex(group.var))

  all_terms <- unique(unlist(term_vectors, use.names = FALSE))
  all_terms <- setdiff(all_terms, "(Intercept)")
  all_terms <- all_terms[grepl(group_prefix_pattern, all_terms)]

  if (length(all_terms) == 0) {
    return(list())
  }

  grouped_results <- lapply(all_terms, function(term) {
    term_results <- lapply(names(test.list), function(metric_name) {
      df <- test.list[[metric_name]]
      if (is.null(df) || !("Term" %in% colnames(df))) {
        return(NULL)
      }

      filtered_df <- df[df$Term == term, , drop = FALSE]
      if (nrow(filtered_df) == 0) {
        return(NULL)
      }

      filtered_df %>%
        dplyr::select(-Term) %>%
        dplyr::mutate(Term = metric_name) %>%
        dplyr::select(Term, dplyr::everything()) %>%
        tibble::as_tibble()
    })

    term_results <- term_results[!vapply(term_results, is.null, logical(1))]
    if (length(term_results) == 0) {
      return(NULL)
    }

    dplyr::bind_rows(term_results)
  })

  labels <- vapply(all_terms, function(term) {
    if (identical(term, group.var)) {
      return(group.var)
    }

    modified_term <- stringr::str_replace(term, group_prefix_pattern, "")
    modified_term <- stringr::str_trim(modified_term)
    paste(sprintf("%s vs %s", modified_term, reference_level), "(Reference)")
  }, character(1))

  keep_idx <- !vapply(grouped_results, is.null, logical(1))
  grouped_results <- grouped_results[keep_idx]
  labels <- labels[keep_idx]

  stats::setNames(grouped_results, labels)
}


#' @keywords internal
mStat_has_nonempty_time_result <- function(result) {
  if (is.null(result) || length(result) == 0) {
    return(FALSE)
  }

  if (is.list(result) && !inherits(result, "data.frame")) {
    return(any(vapply(result, function(x) {
      !is.null(x) && NROW(x) > 0
    }, logical(1))))
  }

  NROW(result) > 0
}


#' @keywords internal
mStat_keep_valid_time_results <- function(test.list,
                                          time.levels,
                                          context) {
  has_results <- vapply(test.list, mStat_has_nonempty_time_result, logical(1))

  if (!any(has_results)) {
    stop(
      "No time points could be analyzed for ",
      context,
      ". Check if there are enough observations and group contrasts at each time point.",
      call. = FALSE
    )
  }

  test.list <- test.list[has_results]
  names(test.list) <- time.levels[which(has_results)]
  test.list
}


#' @keywords internal
mStat_collect_time_result_names <- function(test.list,
                                            entry = NULL) {
  result_names <- unlist(lapply(test.list, function(time_result) {
    if (is.null(time_result)) {
      return(character())
    }

    if (is.null(entry)) {
      return(names(time_result))
    }

    if (!entry %in% names(time_result) || is.null(time_result[[entry]])) {
      return(character())
    }

    names(time_result[[entry]])
  }), use.names = FALSE)

  result_names <- unique(result_names)
  result_names[!is.na(result_names) & nzchar(result_names)]
}


#' @keywords internal
mStat_bind_time_results <- function(test.list,
                                    result_names,
                                    time.var,
                                    entry = NULL) {
  bound_results <- lapply(result_names, function(result_name) {
    data_list <- lapply(names(test.list), function(time_name) {
      time_result <- test.list[[time_name]]
      if (is.null(time_result)) {
        return(NULL)
      }

      result_df <- if (is.null(entry)) {
        time_result[[result_name]]
      } else if (entry %in% names(time_result) && !is.null(time_result[[entry]])) {
        time_result[[entry]][[result_name]]
      } else {
        NULL
      }

      if (is.null(result_df) || NROW(result_df) == 0) {
        return(NULL)
      }

      result_df[[time.var]] <- time_name
      result_df
    })

    data_list <- data_list[!vapply(data_list, is.null, logical(1))]
    dplyr::bind_rows(data_list)
  })

  stats::setNames(bound_results, result_names)
}


#' @keywords internal
mStat_run_per_time_analysis <- function(time.levels,
                                        context,
                                        analysis_fn) {
  test.list <- lapply(time.levels, function(t.level) {
    tryCatch(
      analysis_fn(t.level),
      error = function(e) {
        warning(
          "Error analyzing time point ",
          t.level,
          " for ",
          context,
          ": ",
          conditionMessage(e),
          "\nSkipping this time point and continuing with others."
        )
        NULL
      }
    )
  })

  mStat_keep_valid_time_results(
    test.list = test.list,
    time.levels = time.levels,
    context = context
  )
}
