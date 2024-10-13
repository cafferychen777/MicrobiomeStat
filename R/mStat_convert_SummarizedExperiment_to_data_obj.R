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
  # Initialize an empty list to store the converted data
  # This list will contain various components of the SummarizedExperiment object in a format suitable for MicrobiomeStat
  data.obj <- list()

  # Process the assay data if it exists
  # The assay typically contains the main experimental data, such as gene expression or microbiome abundance
  if (!is.null(assay(se.obj))) {
    # Convert the assay to a matrix format
    # This step ensures compatibility with downstream analyses and improves computational efficiency
    data.obj$feature.tab <- assay(se.obj) %>%
      as.data.frame() %>%
      as.matrix()

    # Remove features (rows) with zero counts across all samples
    # This filtering step is crucial for reducing data sparsity, which can improve statistical power and reduce computational burden in subsequent analyses
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  # Process the column data (sample metadata) if it exists
  # This metadata typically includes important sample-specific information such as experimental conditions or clinical data
  if (!is.null(colData(se.obj))) {
    # Convert the column data to a data frame for easier manipulation in downstream analyses
    data.obj$meta.dat <- colData(se.obj) %>%
      as.data.frame()
  }

  # Process the row data (feature metadata) if it exists
  # This metadata often contains annotations for each feature, such as gene names or taxonomic information for microbiome data
  if (!is.null(rowData(se.obj))) {
    # Convert the row data to a matrix format
    data.obj$feature.ann <- rowData(se.obj) %>%
      as.data.frame() %>%
      as.matrix()

    # Ensure that the feature annotations only include features present in the assay data
    # This step maintains consistency between the feature table and feature annotations, which is crucial for accurate downstream analyses
    if (exists("exp.mat", data.obj)) {
      data.obj$feature.ann <- data.obj$feature.ann[rownames(data.obj$feature.ann) %in% rownames(data.obj$feature.tab), ]
    }
  }

  # Return the processed data object
  # This object contains all the components of the SummarizedExperiment object, reformatted for use with MicrobiomeStat functions
  return(data.obj)
}
