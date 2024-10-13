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
#' to be retained.
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

  # Convert the input matrix to a long format data frame
  # This transformation facilitates group-wise operations and is crucial for calculating prevalence and abundance
  x_long <- as.data.frame(as.table(as.matrix(x)))

  # Perform filtering based on prevalence and abundance
  # This step is critical for reducing noise and focusing on biologically relevant features
  filtered_taxa <- x_long %>%
    # Group the data by taxa (Var1 represents the row names, i.e., taxa)
    dplyr::group_by(Var1) %>%
    # Calculate summary statistics for each taxon
    dplyr::summarise(
      # Calculate the average abundance across all samples
      # This metric helps identify consistently abundant taxa
      avg_abundance = mean(Freq),
      # Calculate the prevalence (proportion of samples where the taxon is present)
      # Prevalence is a key metric in microbiome studies, indicating how widespread a taxon is
      prevalence = sum(Freq > 0) / dplyr::n()
    ) %>%
    # Apply the filtering criteria
    # This step removes taxa that don't meet the specified thresholds, reducing data dimensionality
    dplyr::filter(prevalence >= prev.filter, avg_abundance >= abund.filter) %>%
    # Extract the names of the taxa that pass the filters
    dplyr::pull("Var1")

  # Subset the original matrix to include only the filtered taxa
  # This creates a new matrix with reduced dimensions, focusing on the most relevant taxa
  filtered_x <- x[filtered_taxa,]

  # Return the filtered matrix
  # The resulting matrix contains only taxa that meet both prevalence and abundance criteria
  return(filtered_x)
}
