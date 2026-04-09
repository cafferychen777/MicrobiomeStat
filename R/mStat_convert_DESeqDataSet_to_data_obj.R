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
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop(
      "Package 'SummarizedExperiment' is required to convert DESeqDataSet objects.",
      call. = FALSE
    )
  }

  data.obj <- list()

  assay_data <- SummarizedExperiment::assay(dds.obj)
  if (!is.null(assay_data)) {
    data.obj$feature.tab <- assay_data %>%
      as.data.frame() %>%
      as.matrix()
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, , drop = FALSE]
  }

  col_data <- SummarizedExperiment::colData(dds.obj)
  if (!is.null(col_data)) {
    data.obj$meta.dat <- col_data %>%
      as.data.frame()
  }

  row_data <- SummarizedExperiment::rowData(dds.obj)
  if (!is.null(row_data)) {
    data.obj$feature.ann <- row_data %>%
      as.data.frame() %>%
      as.matrix()

    if (!is.null(data.obj$feature.tab)) {
      data.obj$feature.ann <- data.obj$feature.ann[
        rownames(data.obj$feature.tab),
        ,
        drop = FALSE
      ]
    }
  }

  return(data.obj)
}
