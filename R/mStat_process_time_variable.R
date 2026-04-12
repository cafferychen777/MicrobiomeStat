#' Process Time Variable in Metadata
#'
#' Sets factor levels for time variable and optionally subsets data to specified time points.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @return A MicrobiomeStat data object with processed time variable.
#'
#' @examples
#' data(ecam.obj)
#' mStat_process_time_variable(data.obj = ecam.obj, time.var = "month",
#'   t0.level = "0", ts.levels = c("1","2"))
#'
#' @export
mStat_process_time_variable <-
  function(data.obj,
           time.var,
           t0.level = NULL,
           ts.levels = NULL) {
    time_values <- data.obj$meta.dat[[time.var]]
    resolved_time <- mStat_resolve_time_levels(
      values = time_values,
      time.var = time.var,
      t0.level = t0.level,
      ts.levels = ts.levels,
      context = "time processing"
    )
    target_levels <- resolved_time$kept_levels

    data.obj$meta.dat <- data.obj$meta.dat %>%
      dplyr::mutate(
        !!sym(time.var) := factor(as.character(.data[[time.var]]), levels = target_levels, ordered = TRUE)
      )

    if (!is.null(resolved_time$dropped_levels_message)) {
      message(resolved_time$dropped_levels_message)
      keep_rows <- !is.na(data.obj$meta.dat[[time.var]])
      if (!is.null(data.obj$feature.tab)) {
        data.obj <- mStat_subset_data(
          data.obj,
          samIDs = rownames(data.obj$meta.dat)[keep_rows],
          prune.features = TRUE
        )
      } else {
        data.obj$meta.dat <- data.obj$meta.dat[keep_rows, , drop = FALSE]
      }
    }

    # Return the processed data object.
    # The time variable in this object is now properly formatted for longitudinal analyses.
    return(data.obj)
  }
