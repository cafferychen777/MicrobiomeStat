#' Filter a Microbiome Data Matrix by Prevalence and Average Abundance
#'
#' This function filters taxa in a microbiome matrix based on specified
#' prevalence and average abundance thresholds.
#'
#' @param x A matrix containing taxa (in rows) and sample (in columns) microbial abundance data.
#' @param prev.filter Numeric, the minimum prevalence threshold for a taxon to be
#' retained. Prevalence is calculated as the proportion of samples where the
#' taxon is present.
#' @param abund.filter Numeric, the minimum average abundance threshold for a taxon
#' to be retained. For data with negative values, set to -Inf to disable abundance filtering.
#'
#' @return A matrix with taxa filtered based on the specified thresholds.
#'
#' @details
#' The function first converts the input data into a long format and then groups by taxa.
#' It computes both the average abundance and prevalence for each taxon. Subsequently, it
#' filters out taxa that do not meet the provided prevalence and average abundance thresholds.
#'
#' @examples
#'
#' # Example with simulated data
#' \dontrun{
#' data_matrix <- matrix(c(0, 3, 4, 0, 2, 7, 8, 9, 10), ncol=3)
#' colnames(data_matrix) <- c("sample1", "sample2", "sample3")
#' rownames(data_matrix) <- c("taxa1", "taxa2", "taxa3")
#'
#' filtered_data_simulated <- mStat_filter(data_matrix, 0.5, 5)
#' print(filtered_data_simulated)
#' }
#'
#' # Example with real dataset: peerj32.obj
#' \dontrun{
#' data(peerj32.obj)
#' data_matrix_real <- peerj32.obj$feature.tab
#'
#' # Assuming the matrix contains counts with taxa in rows and samples in columns
#' filtered_data_real <- mStat_filter(data_matrix_real, 0.5, 5)
#' print(filtered_data_real)
#' }
#' @export
mStat_filter <- function(x, prev.filter, abund.filter){
  x_mat <- as.matrix(x)

  avg_abundance <- rowMeans(x_mat, na.rm = TRUE)
  non_missing <- rowSums(!is.na(x_mat))
  prevalence <- rep(NA_real_, nrow(x_mat))
  has_observed <- non_missing > 0
  prevalence[has_observed] <- rowSums((!is.na(x_mat)) & x_mat != 0)[has_observed] / non_missing[has_observed]

  keep_taxa <- !is.na(prevalence) & prevalence >= prev.filter
  if (!(is.infinite(abund.filter) && abund.filter < 0)) {
    keep_taxa <- keep_taxa & avg_abundance >= abund.filter
  }

  filtered_x <- x_mat[keep_taxa, , drop = FALSE]

  if (is.data.frame(x)) {
    return(as.data.frame(filtered_x, stringsAsFactors = FALSE, check.names = FALSE))
  }

  filtered_x
}
