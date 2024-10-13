#' Convert a DESeqDataSet Object to a MicrobiomeStat Data Object
#'
#' This function is a part of the MicrobiomeStat package. It converts a DESeqDataSet object into a list (MicrobiomeStat data object) containing the counts, samples data and feature annotations.
#' The returned list represents a microbiome dataset with features (OTUs or taxa), sample metadata and feature annotations.
#'
#' @name mStat_convert_DESeqDataSet_to_data_obj
#' @param dds.obj A DESeqDataSet object to be converted.
#'
#' @return A MicrobiomeStat data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix containing the counts data with rows as features and columns as samples.
#'   \item meta.dat: A data frame containing the samples data.
#'   \item feature.ann: A matrix containing feature annotations.
#' }
#'
#' @examples
#' \dontrun{
#'   # Load necessary packages
#'   # library(airway)
#'   # library(DESeq2)
#'
#'   # Load dataset
#'   # data("airway")
#'   # dds <- DESeqDataSet(airway, design = ~ cell + dex)
#'
#'   # Convert to MicrobiomeStat data object
#'   # data.obj <- mStat_convert_DESeqDataSet_to_data_obj(dds)
#' }
#'
#' @details
#' The function first checks if each component (counts, samples, feature annotations) of the DESeqDataSet object is not null. If a component is not null, it is converted to the appropriate format and added to the output list. The counts data and feature annotations are converted to a matrix, while the samples data is converted to a data frame. Note that only the rows (features) in the counts data that have a sum > 0 are retained.
#'
#' @author Caffery Yang
#'
#' @export
mStat_convert_DESeqDataSet_to_data_obj <- function (dds.obj) {

  # Initialize an empty list to store the converted data
  # This list will contain three main components: feature table, metadata, and feature annotations
  data.obj <- list()

  # Extract and process the count data (feature table)
  if (!is.null(assay(dds.obj))) {
    # Convert the count matrix to a data frame, then to a matrix
    # This ensures consistent data structure and allows for easier manipulation
    data.obj$feature.tab <- assay(dds.obj) %>%
      as.data.frame() %>%
      as.matrix()

    # Remove features (rows) with zero counts across all samples
    # This step helps to reduce the dimensionality of the data and focus on informative features
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  # Extract and process the sample metadata
  if (!is.null(colData(dds.obj))) {
    # Convert the column data to a data frame
    # This preserves the sample-level information associated with the count data
    data.obj$meta.dat <- colData(dds.obj) %>%
      as.data.frame()
  }

  # Extract and process the feature annotations
  if (!is.null(rowData(dds.obj))) {
    # Convert the row data to a data frame, then to a matrix
    # This ensures consistent structure with the feature table
    data.obj$feature.ann <- rowData(dds.obj) %>%
      as.data.frame() %>%
      as.matrix()

    # Ensure that feature annotations correspond to the features in the count data
    # This step is crucial for maintaining data integrity and consistency
    if (exists("feature.tab", data.obj)) {
      data.obj$feature.ann <- data.obj$feature.ann[rownames(data.obj$feature.ann) %in% rownames(data.obj$feature.tab), ]
    }
  }

  # Return the converted data object
  # This object now contains the count data, sample metadata, and feature annotations in a standardized format
  return(data.obj)
}