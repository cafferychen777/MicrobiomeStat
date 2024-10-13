#' Convert a DGEList Object to a MicrobiomeStat Data Object
#'
#' This function is a part of the MicrobiomeStat package. It converts a DGEList object into a list (MicrobiomeStat data object) containing the counts and samples data.
#' The returned list represents a microbiome dataset with features (OTUs or taxa) and metadata.
#'
#' @name mStat_convert_DGEList_to_data_obj
#' @param dge.obj A DGEList object to be converted.
#'
#' @return A MicrobiomeStat data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix containing the counts data with rows as features and columns as samples.
#'   \item meta.dat: A data frame containing the samples data.
#' }
#'
#' @examples
#' \dontrun{
#'   # Load necessary packages
#'   # library(airway)
#'   # library(DESeq2)
#'   # library(edgeR)
#'
#'   # Load dataset
#'   # data("airway")
#'   # dds <- DESeqDataSet(airway, design = ~ cell + dex)
#'   # dge <- DGEList(counts = counts(dds), group = dds$dex)
#'
#'   # Convert to MicrobiomeStat data object
#'   # data.obj <- mStat_convert_DGEList_to_data_obj(dge)
#' }
#'
#' @details
#' The function first checks if each component (counts, samples) of the DGEList object is not null. If a component is not null, it is converted to the appropriate format and added to the output list. The counts data is converted to a matrix, while the samples data is converted to a data frame. Note that only the rows (features) in the counts data that have a sum > 0 are retained.
#'
#' @author Caffery Yang
#'
#' @export
mStat_convert_DGEList_to_data_obj <- function (dge.obj) {

  # Initialize an empty list to store the converted data
  # This list will contain two main components: feature table and metadata
  data.obj <- list()

  # Extract and process the count data (feature table)
  if (!is.null(dge.obj$counts)) {
    # Convert the count matrix to a data frame, then to a matrix
    # This ensures consistent data structure and allows for easier manipulation
    data.obj$feature.tab <- dge.obj$counts %>%
      as.data.frame() %>%
      as.matrix()

    # Remove features (rows) with zero counts across all samples
    # This step helps to reduce the dimensionality of the data and focus on informative features
    # It's particularly important in RNA-seq analysis to remove genes with no expression
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  # Extract and process the sample metadata
  if (!is.null(dge.obj$samples)) {
    # Convert the sample information to a data frame
    # This preserves the sample-level information associated with the count data
    # Sample metadata typically includes experimental factors and other relevant information
    data.obj$meta.dat <- dge.obj$samples %>%
      as.data.frame()
  }

  # Return the converted data object
  # This object now contains the count data and sample metadata in a standardized format
  # The returned object can be used for further analysis in the MicrobiomeStat pipeline
  return(data.obj)
}
