#' Convert a MRExperiment Object to a MicrobiomeStat Data Object
#'
#' This function is a part of the MicrobiomeStat package. It converts a MRExperiment object into a list (MicrobiomeStat data object) containing the expression data, samples data, and feature annotations.
#' The returned list represents a microbiome dataset with features (OTUs or taxa), sample metadata, and feature annotations.
#'
#' @name mStat_convert_MRExperiment_to_data_obj
#' @param mr.obj A MRExperiment object to be converted.
#'
#' @return A MicrobiomeStat data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix containing the expression data with rows as features and columns as samples.
#'   \item meta.dat: A data frame containing the samples data.
#'   \item feature.ann: A matrix containing feature annotations.
#' }
#'
#' @examples
#' \dontrun{
#'   # Load necessary packages
#'   # library(MetagenomeSeq)
#'
#'   # Load dataset
#'   # data("MRexperiment")
#'   # mr <- MRexperiment(counts = counts(MRexperiment), pData = pData(MRexperiment)
#'   # , fData = fData(MRexperiment))
#'   # Convert to MicrobiomeStat data object
#'   # data.obj <- mStat_convert_MRExperiment_to_data_obj(mr)
#' }
#'
#' @details
#' The function first checks if each component (expression data, samples, feature annotations) of the MRExperiment object is not null. If a component is not null, it is converted to the appropriate format and added to the output list. The expression data and feature annotations are converted to a matrix, while the samples data is converted to a data frame. Note that only the rows (features) in the expression data that have a sum > 0 are retained.
#'
#' @author Chen Yang
#' @references Paulson JN, Stine OC, Bravo HC, Pop M. Differential abundance analysis for microbial marker-gene surveys. Nature Methods. 2013;10(12):1200-1202.
#'
#' @export
mStat_convert_MRExperiment_to_data_obj <- function (mr.obj) {

  # Initialize an empty list to store the converted data
  # This list will contain three main components: feature table, metadata, and feature annotations
  # The structure mimics the standard format used in MicrobiomeStat for consistency across analyses
  data.obj <- list()

  # Extract and process the expression data (feature table)
  if (!is.null(exprs(mr.obj))) {
    # Convert the expression matrix to a data frame, then to a matrix
    # This ensures consistent data structure and allows for easier manipulation
    data.obj$feature.tab <- exprs(mr.obj) %>%
      as.data.frame() %>%
      as.matrix()

    # Remove features (rows) with zero counts across all samples
    # This step is crucial in microbiome data analysis to reduce sparsity and focus on informative features
    # It helps to improve statistical power and reduce computational burden in downstream analyses
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  # Extract and process the phenotype data (sample metadata)
  if (!is.null(pData(mr.obj))) {
    # Convert the phenotype data to a data frame
    # This preserves the sample-level information associated with the expression data
    # Phenotype data typically includes experimental factors, clinical variables, and other relevant metadata
    data.obj$meta.dat <- pData(mr.obj) %>%
      as.data.frame()
  }

  # Extract and process the feature data (feature annotations)
  if (!is.null(fData(mr.obj))) {
    # Convert the feature data to a data frame, then to a matrix
    # This ensures consistent structure with the feature table
    data.obj$feature.ann <- fData(mr.obj) %>%
      as.data.frame() %>%
      as.matrix()

    # Ensure that feature annotations correspond to the features in the expression data
    # This step is crucial for maintaining data integrity and consistency
    # It allows for accurate mapping between features and their annotations in downstream analyses
    if (exists("feature.tab", data.obj)) {
      data.obj$feature.ann <- data.obj$feature.ann[rownames(data.obj$feature.ann) %in% rownames(data.obj$feature.tab), ]
    }
  }

  # Return the converted data object
  # This object now contains the expression data, sample metadata, and feature annotations in a standardized format
  # The returned object can be seamlessly integrated into the MicrobiomeStat analysis pipeline
  return(data.obj)
}