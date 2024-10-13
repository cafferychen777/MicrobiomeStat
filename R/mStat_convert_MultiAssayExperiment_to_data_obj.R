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

  # Initialize an empty list to store the converted data
  data.obj <- list()

  # If no experiment name is provided, use the first experiment in the MultiAssayExperiment object
  if (is.null(experiment_name)) {
    experiment_name <- names(experiments(mae.obj))[1]
  }

  # Check if the specified experiment exists in the MultiAssayExperiment object
  # If not, stop execution and provide an informative error message
  if (!experiment_name %in% names(experiments(mae.obj))) {
    stop(paste("The experiment", experiment_name, "does not exist in the MultiAssayExperiment object."))
  }

  # Extract the assay data for the specified experiment
  assay_data <- assays(mae.obj)[[experiment_name]]
  
  # Process the assay data if it exists
  if (!is.null(assay_data)) {
    # Convert the assay data to a matrix format
    # This step ensures compatibility with downstream analyses
    data.obj$feature.tab <- assay_data %>%
      as.data.frame() %>%
      as.matrix()

    # Remove features (rows) with zero counts across all samples
    # This step is crucial for reducing sparsity and improving statistical power in subsequent analyses
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  # Extract and process the column data (metadata) from the MultiAssayExperiment object
  # This metadata typically includes sample-specific information
  if (!is.null(colData(mae.obj))) {
    # Convert the column data to a data frame for easier manipulation
    data.obj$meta.dat <- colData(mae.obj) %>%
      as.data.frame()
  }

  # Return the processed data object
  # This object contains the feature table and metadata, ready for further analysis
  return(data.obj)
}