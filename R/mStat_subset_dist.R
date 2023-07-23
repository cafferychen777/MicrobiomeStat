#' Check if an Object is NA or NULL
#'
#' This function checks if the given object is `NA` or `NULL`.
#'
#' @param x The object to check.
#'
#' @return A logical value indicating whether the input is `NA` or `NULL`.
#'
#' @details
#' The function first checks if the object `x` is `NULL`. If not, it checks if the first element of `x` is `NA`.
#' If either condition is true, it returns `TRUE`; otherwise, it returns `FALSE`.
#'
#' @examples
#' \dontrun{
#' # Check if a variable is NA or NULL
#' is.na.null(NA)  # TRUE
#' is.na.null(NULL)  # TRUE
#' is.na.null(1)  # FALSE
#' }
#'
#' @export
is.na.null <- function (x) {
  if (is.null(x)) {
    return(TRUE)
  } else {
    if (is.na(x)[1]) {
      return(TRUE)
    }  else {
      return(FALSE)
    }
  }
}


#' Subset Distance Matrix
#'
#' This function subsets a list of distance matrices by a specified set of sample IDs.
#'
#' @param dist.obj A list of distance matrices.
#' @param samIDs A vector of sample IDs to subset by. This can be a logical vector, a numeric vector, or a character vector of sample IDs.
#'
#' @return A list of subsetted distance matrices.
#'
#' @details
#' The function first checks if samIDs is logical or numeric, and if so, converts it to a character vector of sample IDs.
#' Then, it subsets each distance matrix in dist.obj by the sample IDs.
#' If a matrix is `NA` or `NULL`, it remains unchanged.
#'
##' @examples
#' \dontrun{
#' # Load required libraries
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Calculate beta diversity
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#'
#' # Select sample IDs
#' sample_ids <- sample(colnames(peerj32.obj$feature.tab), 10)
#'
#' # Subset distance object by sample IDs
#' subsetted_dist.obj <- mStat_subset_dist(dist.obj, sample_ids)
#' }
#'
#' @export
mStat_subset_dist <- function (dist.obj, samIDs) {

  # If samIDs is logical or numeric, convert it to character form of sample IDs
  if (is.logical(samIDs) | is.numeric(samIDs)) {
    samIDs <- rownames(dist.obj[[1]])[samIDs]
  }

  # Apply the subsetting to each distance matrix in the list
  dist.obj <- lapply(dist.obj, function(x) {
    # If x is not NA or NULL, subset the matrix by samIDs
    if(!is.na.null(x)){
      x <- as.matrix(x)
      x <- x[samIDs, samIDs]
      x <- as.dist(x)
    }
    # Return the processed matrix
    x
  })

  return(dist.obj)
}
