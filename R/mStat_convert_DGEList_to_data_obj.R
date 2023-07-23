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

  data.obj <- list()

  if (!is.null(dge.obj$counts)) {
    data.obj$feature.tab <- dge.obj$counts %>%
      as.data.frame() %>%
      as.matrix()

    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  if (!is.null(dge.obj$samples)) {
    data.obj$meta.dat <- dge.obj$samples %>%
      as.data.frame()
  }

  return(data.obj)
}
