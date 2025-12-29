#' Subset Principal Coordinates Analysis (PCoA) Results
#'
#' This function subsets the results of a Principal Coordinates Analysis (PCoA) by specified sample IDs.
#'
#' @param pc.obj A list of PCoA results, usually calculated using the
#' \code{\link[MicrobiomeStat]{mStat_calculate_PC}} function.
#' @param samIDs A vector of sample IDs to subset by. This can be a logical vector, a numeric vector,
#' or a character vector of sample IDs.
#'
#' @return A list of subsetted PCoA results.
#'
#' @details
#' The function first checks if samIDs is logical or numeric, and if so, converts it to a character vector of sample IDs.
#' Then, it subsets each PCoA result in pc.obj by the sample IDs.
#' If a result is `NA` or `NULL`, it remains unchanged.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Calculate PCoA
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC'))
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC'))
#'
#' # Select sample IDs
#' sample_ids <- sample(colnames(peerj32.obj$feature.tab), 10)
#'
#' # Subset PCoA object by sample IDs
#' subsetted_pc.obj <- mStat_subset_PC(pc.obj, sample_ids)
#' }
#'
#' @export
mStat_subset_PC <- function (pc.obj, samIDs) {

  # If samIDs is logical or numeric, convert it to character form of sample IDs
  if (is.logical(samIDs) || is.numeric(samIDs)) {
    samIDs <- rownames(pc.obj[[1]])[samIDs]
  }

  # Apply the subsetting to each PCoA result in the list
  pc.obj <- lapply(pc.obj, function(x) {
    # If x is not NA or NULL, subset the matrix by samIDs
    if(!is.na.null(x)){
      x$points <- x$points[samIDs, ]
      if (!is.null(x$eig)) {
        x$eig <- x$eig[samIDs]
      }
    }
    # Return the processed result
    x
  })

  return(pc.obj)
}
