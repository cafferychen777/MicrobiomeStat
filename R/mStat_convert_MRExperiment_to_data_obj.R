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
  if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop(
      "Package 'metagenomeSeq' is required to convert MRExperiment objects.",
      call. = FALSE
    )
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop(
      "Package 'Biobase' is required to convert MRExperiment objects.",
      call. = FALSE
    )
  }

  data.obj <- list()

  expr_data <- metagenomeSeq::MRcounts(mr.obj)
  if (!is.null(expr_data)) {
    data.obj$feature.tab <- expr_data %>%
      as.data.frame() %>%
      as.matrix()
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, , drop = FALSE]
  }

  pheno_data <- Biobase::pData(mr.obj)
  if (!is.null(pheno_data)) {
    data.obj$meta.dat <- pheno_data %>%
      as.data.frame()
  }

  feature_data <- Biobase::fData(mr.obj)
  if (!is.null(feature_data)) {
    data.obj$feature.ann <- feature_data %>%
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
