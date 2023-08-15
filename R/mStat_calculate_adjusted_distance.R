#' Calculate Adjusted Distances
#'
#' This function calculates the adjusted distances for a given set of variables,
#' using a model matrix dynamically created from the specified adjustment variables.
#'
#' @param data.obj A MicrobiomeStat data object. The heart of the MicrobiomeStat, consisting of several key components:
#'   \itemize{
#'     \item \strong{Feature.tab (Matrix)}: Meeting point for research objects (OTU/ASV/KEGG/Gene, etc.) and samples.
#'     \item \strong{Meta.dat (Data frame)}: Rows correspond to the samples, and columns serve as annotations, describing the samples.
#'     \item \strong{Feature.ann (Matrix)}: Annotations, carrying classification information like Kingdom, Phylum, etc.
#'     \item \strong{Phylogenetic tree (Optional)}: An evolutionary perspective, illuminating the relationships among various research objects.
#'     \item \strong{Feature.agg.list (Optional)}: Aggregated results based on the feature.tab and feature.ann.
#'   }
#' @param dist.obj A list of distance matrices corresponding to the names in `dist.name`.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param dist.name A character vector of names corresponding to distance matrices in `dist.obj`.
#' @return A list of adjusted distance matrices.
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' dist.obj <- mStat_calculate_beta_diversity(subset_T2D.obj, c("BC","Jaccard"))
#' adj.dist.obj <- mStat_calculate_adjusted_distance(
#' data.obj = subset_T2D.obj,
#' dist.obj = dist.obj,
#' adj.vars = c("subject_gender", "subject_race"),
#' dist.name = c("BC"))
#' data("peerj32.obj")
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, c("BC"))
#' adj.dist.obj <- mStat_calculate_adjusted_distance(
#' data.obj = peerj32.obj,
#' dist.obj = dist.obj,
#' adj.vars = c("sex"),
#' dist.name = c("BC"))
#' }
#' @export
mStat_calculate_adjusted_distance <- function (data.obj,
                                               dist.obj,
                                               adj.vars,
                                               dist.name) {
  # Message to inform the user
  message(
    "Calculating adjusted distances using the provided adjustment variables and distance matrices..."
  )

  # Load metadata and select the required adjustment variables
  meta_tab <-
    load_data_obj_metadata(data.obj) %>% select(all_of(c(adj.vars)))

  # Create a formula string from the adjustment variables and convert it into a formula object
  formula_str <- paste("xx ~", paste(sprintf("meta_tab$%s", adj.vars), collapse = "+"))
  dynamic_formula <- as.formula(formula_str)

  # Iterate over the distance names to calculate adjusted distances
  adj.dist.obj <- lapply(dist.name, function(sub_dist.name) {

    D <- as.matrix(dist.obj[[sub_dist.name]])

    obj <- suppressWarnings(cmdscale(D, k = nrow(D) - 1, eig = TRUE))

    eig <- obj$eig

    s <- sign(eig)

    eig <- abs(eig)

    xx <- obj$points

    res <- residuals(lm(dynamic_formula))

    res_positive <- suppressWarnings(sweep(res, 2, sqrt(eig)[s >= 0], "*"))
    res_negative <- suppressWarnings(sweep(res, 2, sqrt(eig)[s < 0], "*"))

    dist_positive <- as.matrix(dist(res_positive))
    dist_negative <- as.matrix(dist(res_negative))

    D.adj <- dist_positive - dist_negative

    return(as.dist(D.adj))
  })

  names(adj.dist.obj) <- dist.name

  return(adj.dist.obj)
}
