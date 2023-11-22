#' Compute Top K Statistics for OTU Data
#'
#' This internal function applies a given statistical function to the OTU data
#' after converting it to a matrix, with rows named by a specified feature level.
#' It is designed to be used internally within the MicrobiomeStat package and is
#' not exported for public use.
#'
#' @param top.k.func A function or a string specifying the statistical function
#'   to apply to the OTU data. If a function, it should take a matrix as input.
#'   If a string, it should be either "mean" or "sd" to compute the row means or
#'   standard deviations, respectively.
#' @param otu_tax_agg A data frame or similar object containing the OTU data,
#'   which will be converted to a matrix for calculation.
#' @param feature.level A character string specifying the column name in
#'   `otu_tax_agg` that contains the feature level to be used as row names in
#'   the matrix.
#'
#' @return A numeric vector containing the computed statistics for each feature.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   # This is an example of how to use the .compute_function internally
#'   # Suppose otu_tax_agg is a data frame and "Species" is the feature level
#'
#'   # Create a simulated otu_tax_agg data frame with random data
#'   set.seed(123) # for reproducibility
#'   otu_tax_agg <- data.frame(
#'     Species = paste("Species", 1:10),
#'     Sample1 = sample(0:100, 10, replace = TRUE),
#'     Sample2 = sample(0:100, 10, replace = TRUE),
#'     Sample3 = sample(0:100, 10, replace = TRUE)
#'   )
#'
#'   # Compute the mean for each species
#'   results_mean <- compute_function("mean", otu_tax_agg, "Species")
#'   print(results_mean)
#'
#'   # Compute the standard deviation for each species
#'   results_sd <- compute_function("sd", otu_tax_agg, "Species")
#'   print(results_sd)
#'
#'   # Compute the prevalence for each species
#'   results_prevalence <- compute_function("prevalence", otu_tax_agg, "Species")
#'   print(results_prevalence)
#' }
#' @importFrom tibble column_to_rownames
#' @importFrom matrixStats rowSds
#'
#' @noRd
compute_function <- function(top.k.func, otu_tax_agg, feature.level) {
  # Convert data to matrix with proper row names
  otu_matrix <- otu_tax_agg %>%
    tibble::column_to_rownames(feature.level) %>%
    as.matrix()

  if (is.function(top.k.func)) {
    # If a function is provided, apply it directly to the matrix
    results <- top.k.func(otu_matrix)
  } else {
    # Use switch to handle string options for top.k.func
    results <- switch(top.k.func,
                      "mean" = {
                        # Calculate row means
                        rowMeans(otu_matrix, na.rm = TRUE)
                      },
                      "sd" = {
                        # Calculate row standard deviations
                        results <- matrixStats::rowSds(otu_matrix, na.rm = TRUE)
                        names(results) <- rownames(otu_matrix)
                        results
                      },
                      "prevalence" = {
                        # Calculate the prevalence (percentage of non-zero entries) for each feature
                        prevalence <- rowSums(otu_matrix > 0) / ncol(otu_matrix)
                        names(prevalence) <- rownames(otu_matrix)
                        prevalence
                      },
                      stop("Invalid function specified")
    )
  }
  return(results)
}
