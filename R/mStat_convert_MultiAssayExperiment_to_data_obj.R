#' Convert MultiAssayExperiment Object to a Data Object
#'
#' This function converts a MultiAssayExperiment object into a list (data object) containing the assays data, and colData.
#' The returned list represents a dataset with features (OTUs or taxa), metadata, and feature annotations.
#'
#' @name mStat_convert_MultiAssayExperiment_to_data_obj
#' @param mae.obj A MultiAssayExperiment object to be converted.
#' @param experiment_name A string to specify which experiment in the MultiAssayExperiment object to be converted. Default is the first experiment.
#'
#' @return A data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix containing the assay data with rows as features and columns as samples.
#'   \item meta.dat: A data frame containing the colData (metadata).
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming 'mae' is your MultiAssayExperiment object
#'   # data_obj <- mStat_convert_MultiAssayExperiment_to_data_obj(mae, "16S")
#'
#'   # If you have an experiment data available as a MultiAssayExperiment object
#'   # you can convert it to a data object using:
#'   # library(MultiAssayExperiment)
#'   # data(mae)
#'   # data_obj <- mStat_convert_MultiAssayExperiment_to_data_obj(mae, "16S")
#' }
#'
#' @details
#' The function first checks if each component (assays, colData) of the MultiAssayExperiment object is not null. If a component is not null, it is converted to the appropriate format and added to the output list. The assays data is converted to a matrix, while the colData is converted to a data frame. Note that only the rows (features) in the assay data that have a sum > 0 are retained.
#'
#' @author Caffery Yang
#'
#' @export
mStat_convert_MultiAssayExperiment_to_data_obj <- function (mae.obj, experiment_name = NULL) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop(
      "Package 'MultiAssayExperiment' is required to convert MultiAssayExperiment objects.",
      call. = FALSE
    )
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop(
      "Package 'SummarizedExperiment' is required to extract assay data from MultiAssayExperiment objects.",
      call. = FALSE
    )
  }

  data.obj <- list()

  experiments_list <- MultiAssayExperiment::experiments(mae.obj)
  if (is.null(experiment_name)) {
    experiment_name <- names(experiments_list)[1]
  }

  if (!experiment_name %in% names(experiments_list)) {
    stop(paste("The experiment", experiment_name, "does not exist in the MultiAssayExperiment object."))
  }

  assay_data <- SummarizedExperiment::assay(experiments_list[[experiment_name]])
  if (!is.null(assay_data)) {
    data.obj$feature.tab <- assay_data %>%
      as.data.frame() %>%
      as.matrix()
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, , drop = FALSE]
  }

  col_data <- MultiAssayExperiment::colData(mae.obj)
  if (!is.null(col_data)) {
    data.obj$meta.dat <- col_data %>%
      as.data.frame()
  }

  if (methods::is(experiments_list[[experiment_name]], "SummarizedExperiment")) {
    row_data <- SummarizedExperiment::rowData(experiments_list[[experiment_name]])
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
  }

  return(data.obj)
}
