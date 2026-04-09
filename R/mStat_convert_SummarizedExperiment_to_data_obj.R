#' Convert SummarizedExperiment Object to a MicrobiomeStat Data Object
#'
#' This function converts a SummarizedExperiment object into a list (MicrobiomeStat data object) containing the assays data, rowData, and colData.
#' The returned list represents a microbiome dataset with features (OTUs or taxa), metadata, and feature annotations.
#'
#' @name mStat_convert_SummarizedExperiment_to_data_obj
#' @param se.obj A SummarizedExperiment object to be converted.
#'
#' @return A MicrobiomeStat data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix containing the assay data with rows as features and columns as samples.
#'   \item feature.ann: A matrix containing the rowData (feature annotations). Only features present in the assay data are included.
#'   \item meta.dat: A data frame containing the colData (metadata).
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming 'se' is your SummarizedExperiment object
#'   # data_obj <- mStat_convert_SummarizedExperiment_to_data_obj(se)
#'
#'   # If you have airway data available as a SummarizedExperiment object
#'   # you can convert it to a MicrobiomeStat data object using:
#'   # library(airway)
#'   # data(airway)
#'   # airway_obj <- mStat_convert_SummarizedExperiment_to_data_obj(airway)
#' }
#'
#' @details
#' The function first checks if each component (assays, rowData, colData) of the SummarizedExperiment object is not null. If a component is not null, it is converted to the appropriate format and added to the output list. The assays data is converted to a matrix, while the rowData and colData are converted to data frames. Note that only the rows (features) in the assay data that have a sum > 0 are retained.
#'
#' @author Caffery Yang
#'
#' @export
mStat_convert_SummarizedExperiment_to_data_obj <- function (se.obj) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop(
      "Package 'SummarizedExperiment' is required to convert SummarizedExperiment objects.",
      call. = FALSE
    )
  }

  data.obj <- list()

  assay_data <- SummarizedExperiment::assay(se.obj)
  if (!is.null(assay_data)) {
    data.obj$feature.tab <- assay_data %>%
      as.data.frame() %>%
      as.matrix()
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, , drop = FALSE]
  }

  col_data <- SummarizedExperiment::colData(se.obj)
  if (!is.null(col_data)) {
    data.obj$meta.dat <- col_data %>%
      as.data.frame()
  }

  row_data <- SummarizedExperiment::rowData(se.obj)
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
