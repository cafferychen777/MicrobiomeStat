#' Validate and Adjust a MicrobiomeStat Data Object
#'
#' Checks data object structure and consistency, adjusting components if needed.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @return The validated (and possibly adjusted) data object.
#'
#' @examples
#' \dontrun{
#' # Assume 'data.obj' is your data object
#' data(peerj32.obj)
#' mStat_validate_data(peerj32.obj)
#' mStat_validate_data(mStat_normalize_data(peerj32.obj, "TSS")$data.obj.norm)
#' }
#'
#' @details
#' The function checks if 'meta.dat' is a data frame and converts it to one if it isn't. The function then checks two rules: whether the row names of 'feature.tab' match the row names of 'feature.ann', and whether the column names of 'feature.tab' match the row names of 'meta.dat'. If either of these checks fail, the function adjusts the relevant parts of the data object to meet these rules.
#'
#' @export
mStat_validate_data <- function(data.obj) {
  if (!is.list(data.obj)) {
    stop("Rule 1 failed: data.obj should be a list.")
  }

  if (is.null(data.obj$feature.tab)) {
    stop("Rule 2 failed: data.obj$feature.tab is required.")
  }

  if (is.data.frame(data.obj$feature.tab)) {
    data.obj$feature.tab <- as.matrix(data.obj$feature.tab)
    message("Converted feature.tab to a matrix.")
  }
  if (!is.matrix(data.obj$feature.tab)) {
    stop("Rule 3 failed: feature.tab should be a matrix or data.frame.")
  }

  if (is.null(rownames(data.obj$feature.tab)) || is.null(colnames(data.obj$feature.tab))) {
    stop("Rule 4 failed: feature.tab must have both row names (features) and column names (samples).")
  }
  if (anyDuplicated(rownames(data.obj$feature.tab)) != 0) {
    stop("Rule 5 failed: feature.tab row names must be unique.")
  }
  if (anyDuplicated(colnames(data.obj$feature.tab)) != 0) {
    stop("Rule 6 failed: feature.tab column names must be unique.")
  }

  if (is.null(data.obj$meta.dat)) {
    stop("Rule 7 failed: data.obj$meta.dat is required.")
  }

  if (!is.data.frame(data.obj$meta.dat)) {
    data.obj$meta.dat <- as.data.frame(data.obj$meta.dat)
    message("Converted meta.dat to a data.frame.")
  }

  if (is.null(rownames(data.obj$meta.dat)) && nrow(data.obj$meta.dat) == ncol(data.obj$feature.tab)) {
    rownames(data.obj$meta.dat) <- colnames(data.obj$feature.tab)
    message("meta.dat row names were missing and have been aligned to feature.tab column names.")
  }
  if (is.null(rownames(data.obj$meta.dat))) {
    stop("Rule 8 failed: meta.dat must have row names matching feature.tab column names.")
  }
  if (anyDuplicated(rownames(data.obj$meta.dat)) != 0) {
    stop("Rule 9 failed: meta.dat row names must be unique.")
  }

  missing_meta <- setdiff(colnames(data.obj$feature.tab), rownames(data.obj$meta.dat))
  if (length(missing_meta) != 0) {
    stop(
      "Rule 10 failed: meta.dat is missing sample rows for: ",
      paste(missing_meta, collapse = ", ")
    )
  }
  data.obj$meta.dat <- data.obj$meta.dat[colnames(data.obj$feature.tab), , drop = FALSE]

  if (!is.null(data.obj$feature.ann)) {
    if (is.data.frame(data.obj$feature.ann)) {
      data.obj$feature.ann <- as.matrix(data.obj$feature.ann)
      message("Converted feature.ann to a matrix.")
    }
    if (!is.matrix(data.obj$feature.ann)) {
      stop("Rule 11 failed: feature.ann should be a matrix or data.frame when provided.")
    }
    if (is.null(rownames(data.obj$feature.ann)) && nrow(data.obj$feature.ann) == nrow(data.obj$feature.tab)) {
      rownames(data.obj$feature.ann) <- rownames(data.obj$feature.tab)
      message("feature.ann row names were missing and have been aligned to feature.tab row names.")
    }
    if (is.null(rownames(data.obj$feature.ann))) {
      stop("Rule 12 failed: feature.ann must have row names matching feature.tab row names when provided.")
    }
    if (anyDuplicated(rownames(data.obj$feature.ann)) != 0) {
      stop("Rule 13 failed: feature.ann row names must be unique.")
    }

    missing_ann <- setdiff(rownames(data.obj$feature.tab), rownames(data.obj$feature.ann))
    if (length(missing_ann) != 0) {
      stop(
        "Rule 14 failed: feature.ann is missing feature rows for: ",
        paste(missing_ann, collapse = ", ")
      )
    }
    data.obj$feature.ann <- data.obj$feature.ann[rownames(data.obj$feature.tab), , drop = FALSE]
  }

  message("Validation passed.")
  return(invisible(data.obj))
}
