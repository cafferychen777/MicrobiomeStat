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

  if (is.null(dist.name)){
    return()
  }

  calculate_single_method <- function(m, dist_matrix, k, perplexity = NULL) {
    if (m == 'mds') {
      message("Calculating MDS...")
      return(cmdscale(dist_matrix, eig = TRUE, k = k))
    } else if (m == 'nmds') {
      message("Calculating NMDS...")
      return(metaMDS(dist_matrix, distance = "euclidean", k = k))
    } else {
      warning(paste("Unsupported method:", m))
      return(NULL)
    }
  }

  message("Calculating PC...")
  pc.obj <- lapply(dist.name, function(dist_name) {
    message(paste("Processing", dist_name, "distance..."))
    dist_matrix <- dist.obj[[dist_name]]
    pc_results <- calculate_single_method(method, dist_matrix, k)
    return(pc_results)
  })
  names(pc.obj) <- dist.name

  message("Calculation complete.")
  return(pc.obj)
}
