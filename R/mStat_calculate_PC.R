#' Calculate Principal Coordinates using MDS, NMDS
#'
#' This function calculates the Principal Coordinates based on
#' different methods such as Metric Multi-Dimensional Scaling (MDS),
#' Non-Metric Multi-Dimensional Scaling (NMDS)
#' @name mStat_calculate_PC
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param method A character vector specifying which methods to use for
#'   calculating PCoA. Supported methods are "mds" (MDS), "nmds" (NMDS)
#' @param k An integer specifying the number of principal coordinates to retain.
#'   Default is 2.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#'   to be used from the dist.obj list. Default is 'BC'.
#'
#' @return A list containing the PCoA results for each specified method and
#'   distance matrix. The results are named with the method's abbreviation.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC'))
#'
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC'))
#' }
#'
#' @export
mStat_calculate_PC <- function(dist.obj, method = c('mds'), k = 2, dist.name = NULL) {

  # Check if dist.name is NULL and return early if so
  # This prevents unnecessary computations if no distance metric is specified
  if (is.null(dist.name)){
    return()
  }

  # Define an inner function to calculate principal coordinates for a single method
  # This function encapsulates the logic for different ordination techniques
  calculate_single_method <- function(m, dist_matrix, k, perplexity = NULL) {
    if (m == 'mds') {
      # Metric Multidimensional Scaling (MDS)
      # MDS aims to preserve the between-object distances in a lower-dimensional space
      # It is useful for visualizing the level of similarity of individual cases in a dataset
      message("Calculating MDS...")
      return(cmdscale(dist_matrix, eig = TRUE, k = k))
    } else if (m == 'nmds') {
      # Non-Metric Multidimensional Scaling (NMDS)
      # NMDS is a rank-based approach that maximizes the correlation between 
      # distances in the original high-dimensional space and distances in the ordination space
      # It's often used when a non-linear relationship between dissimilarities is suspected
      message("Calculating NMDS...")
      return(metaMDS(dist_matrix, distance = "euclidean", k = k))
    } else {
      # If an unsupported method is specified, warn the user and return NULL
      warning(paste("Unsupported method:", m))
      return(NULL)
    }
  }

  # Inform the user that the principal coordinate calculation is starting
  message("Calculating PC...")
  
  # Use lapply to iterate over each specified distance metric
  # This allows for efficient calculation of principal coordinates for multiple distance matrices
  pc.obj <- lapply(dist.name, function(dist_name) {
    # Inform the user about which distance matrix is being processed
    message(paste("Processing", dist_name, "distance..."))
    
    # Extract the distance matrix for the current metric
    dist_matrix <- dist.obj[[dist_name]]
    
    # Calculate principal coordinates using the specified method
    pc_results <- calculate_single_method(method, dist_matrix, k)
    
    # Return the results for this distance metric
    return(pc_results)
  })
  
  # Assign names to the list of results based on the distance metrics used
  names(pc.obj) <- dist.name

  # Inform the user that all calculations are complete
  message("Calculation complete.")
  
  # Return the list of principal coordinate results
  return(pc.obj)
}
