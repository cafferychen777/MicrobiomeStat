#' Validate and Adjust a MicrobiomeStat Data Object
#'
#' This function is a part of the MicrobiomeStat package. It validates a data object, checks if it meets certain rules, and adjusts it if necessary.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#'
#' @return A list. The validated and possibly adjusted data object.
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
  # Rule 1: Check if data.obj is the correct type
  if (!is.list(data.obj)) {
    stop("Rule 1 failed: data.obj should be a list.")
  } else {
    message("Rule 1 passed: data.obj is a list.")
  }

  # Rule 2: Check if meta.dat is a data.frame
  if (!is.data.frame(data.obj$meta.dat)) {
    # If it is not, convert it to a data.frame
    data.obj$meta.dat <- as.data.frame(data.obj$meta.dat)
    message("Rule 2 passed: meta.dat has been converted to a data.frame.")
  } else {
    message("Rule 2 passed: meta.dat is a data.frame.")
  }

  # Rule 3: Check if the row names of feature.tab match the row names of feature.ann
  if (!identical(rownames(data.obj$feature.tab), rownames(data.obj$feature.ann))) {
    # If they do not match, adjust the order of rows in feature.ann to match feature.tab
    data.obj$feature.ann <- data.obj$feature.ann[match(rownames(data.obj$feature.tab), rownames(data.obj$feature.ann)), ]
    message("Rule 3 passed: The order of rows in feature.ann has been adjusted to match feature.tab.")
  } else {
    message("Rule 3 passed: The row names of feature.tab match the row names of feature.ann.")
  }

  # Rule 4: Check if the column names of feature.tab match the row names of meta.dat
  if (!identical(colnames(data.obj$feature.tab), rownames(data.obj$meta.dat))) {
    # If they do not match, adjust the order of rows in meta.dat to match feature.tab
    data.obj$meta.dat <- data.obj$meta.dat[match(colnames(data.obj$feature.tab), rownames(data.obj$meta.dat)), ]
    message("Rule 4 passed: The order of rows in meta.dat has been adjusted to match feature.tab.")
  } else {
    message("Rule 4 passed: The column names of feature.tab match the row names of meta.dat.")
  }

  # Rule 5: Check if feature.tab is a matrix
  if (!is.matrix(data.obj$feature.tab)) {
    stop("Rule 5 failed: feature.tab should be a matrix.")
  } else {
    message("Rule 5 passed: feature.tab is a matrix.")
  }

  # Rule 6: Check if feature.ann is a matrix
  if (!is.matrix(data.obj$feature.ann)) {
    stop("Rule 6 failed: feature.ann should be a matrix.")
  } else {
    message("Rule 6 passed: feature.ann is a matrix.")
  }

  # Add additional rules here...

  # Add a message to remind the user about using base R data structures
  message("Please note: The data components should follow base R data.frame and matrix structures, not phyloseq's formal class.")

  # If all checks passed, print a validation passed message
  message("Validation passed.")

  # Add a final message to remind the user that validation is not exhaustive
  message("Note: Passing validation does not guarantee the absence of all data issues. Further data exploration may be needed.")

}
