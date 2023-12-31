#' Subset Alpha Diversity Object
#'
#' This function subsets an alpha diversity object by a specified set of sample IDs.
#'
#' @param alpha.obj Alpha diversity object generated by mStat_calculate_alpha_diversity.
#' @param samIDs A vector of sample IDs to subset by. This can be a logical vector, a numeric vector, or a character vector of sample IDs.
#'
#' @return A list of subsetted alpha diversity indices.
#'
#' @details
#' The function first checks if samIDs is logical or numeric, and if so, converts it to a character vector of sample IDs.
#' Then, it subsets each alpha diversity index in alpha.obj by the sample IDs.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(MicrobiomeStat)
#' library(vegan)
#'
#' # Create example OTU table
#' otu.tab <- matrix(data = rpois(100, 5), nrow = 10, ncol = 10)
#' rownames(otu.tab) <- paste0("Taxon_", 1:10)
#' colnames(otu.tab) <- paste0("Sample_", 1:10)
#'
#' # Calculate alpha diversity indices
#' alpha.obj <- mStat_calculate_alpha_diversity(x = otu.tab, alpha.name = c("shannon", "simpson",
#' "observed_species", "chao1", "ace", "pielou"))
#'
#' # Subset alpha diversity object by sample IDs
#' sample_ids <- sample(colnames(otu.tab), 5)
#' subsetted_alpha.obj <- mStat_subset_alpha(alpha.obj, sample_ids)
#' }
#'
#' @export
mStat_subset_alpha <- function(alpha.obj, samIDs) {
  # If samIDs is logical or numeric, convert it to character form of sample IDs
  if (is.logical(samIDs) | is.numeric(samIDs)) {
    samIDs <- colnames(alpha.obj[[1]])[samIDs]
  }

  # Apply the subsetting to each alpha diversity index in the list
  alpha.obj <- lapply(alpha.obj, function(x) {
    if (!is.null(x)) {
      x <- x[samIDs, , drop = FALSE]
    }
    # Return the processed alpha diversity index
    x
  })

  return(alpha.obj)
}
