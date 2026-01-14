#' Update Sample Names in Data Object
#'
#' Renames samples across all components of a MicrobiomeStat data object.
#'
#' @inheritParams mStat_data_obj_doc
#' @param new.name Character vector of new sample names (must match number of samples, no duplicates).
#'
#' @return A MicrobiomeStat data object with updated sample names.
#'
#' @examples
#' \dontrun{
#' # Load the required libraries
#' library(MicrobiomeStat)
#'
#' data(peerj32.obj)
#' # Update sample names
#' new.names <- paste0("new-", colnames(peerj32.obj$feature.tab))
#' updated_peerj32.obj <- mStat_update_sample_name(data.obj = peerj32.obj,
#' new.name = new.names)
#' }
#'
#' @export
mStat_update_sample_name <- function (data.obj, new.name) {
  # Check if the number of new sample names matches the number of samples in the metadata.
  # This ensures that we have a new name for each existing sample.
  # If the numbers don't match, the function stops execution and returns an error message.
  if (length(new.name) != nrow(data.obj$meta.dat)) stop('The number of sample do not agree!\n')

  # Check for duplicates in the new sample names.
  # Duplicate names could lead to ambiguity and errors in downstream analyses.
  # If duplicates are found, the function stops execution and returns an error message.
  if (length(new.name) != length(unique(new.name))) stop ('The new names have duplicates!\n')

  # Update the row names of the metadata table with the new sample names.
  # This step is crucial for maintaining consistency between sample identifiers and metadata.
  rownames(data.obj$meta.dat) <- new.name

  # Check if the data object contains a feature aggregation list.
  # This list typically contains aggregated data at different taxonomic levels.
  if ("feature.agg.list" %in% names(data.obj)) {
    # Update the column names of each element in the feature aggregation list.
    # We use lapply to apply the renaming function to each element of the list.
    # This ensures that sample names are consistent across all levels of aggregation.
    data.obj$feature.agg.list <- lapply(data.obj$feature.agg.list, function (x) {
      colnames(x) <- new.name
      return(x)
    })
  }

  # Update the column names of the feature table.
  # The feature table contains the abundance data for each feature (e.g., OTU, ASV) across all samples.
  # Updating these names ensures consistency between the abundance data and sample identifiers.
  colnames(data.obj$feature.tab) <- new.name

  # Return the updated data object.
  # The returned object now has consistent sample names across all its components,
  # which is crucial for correct interpretation and analysis of the microbiome data.
  return(data.obj)
}
