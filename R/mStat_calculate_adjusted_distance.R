#' Calculate Adjusted Distances Based on Specified Adjustment Variables
#'
#' This function provides a way to calculate distances adjusted for specified variables in a dataset.
#' It incorporates the results of multidimensional scaling (MDS) and adjusts distances based on a linear model.
#'
#' @title Calculate Adjusted Distances for Specified Variables
#'
#' @description Computes adjusted distance matrices for a set of variables using a dynamic model matrix.
#'
#' @param data.obj A MicrobiomeStat data object. This central component contains:
#'   \itemize{
#'     \item \strong{feature.tab}: Matrix representation linking research objects and samples.
#'     \item \strong{meta.dat}: Data frame with rows representing samples and columns as annotations.
#'     \item \strong{feature.ann}: Matrix annotations with classification data.
#'     \item \strong{phylogenetic tree}: (Optional) Tree depicting evolutionary relationships.
#'     \item \strong{feature.agg.list}: (Optional) Aggregated results from feature.tab and feature.ann.
#'   }
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param adj.vars Character vector listing variable names intended for adjustment.
#' @param dist.name Character vector specifying names of distance matrices present in `dist.obj`.
#'
#' @return A list of adjusted distance matrices. Each matrix is named according to the `dist.name` parameter.
#'
#' @details The function uses cmdscale for multidimensional scaling and then adjusts the resultant distances based on a linear model using the variables provided in `adj.vars`.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' subset_T2D.dist.obj <- mStat_calculate_beta_diversity(subset_T2D.obj, c("BC","Jaccard"))
#' adj.dist.obj <- mStat_calculate_adjusted_distance(
#' data.obj = subset_T2D.obj,
#' dist.obj = subset_T2D.dist.obj,
#' adj.vars = c("subject_gender", "subject_race"),
#' dist.name = c("BC"))
#' }
#' data("peerj32.obj")
#' peerj32.dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, c("BC"))
#' adj.dist.obj <- mStat_calculate_adjusted_distance(
#' data.obj = peerj32.obj,
#' dist.obj = peerj32.dist.obj,
#' adj.vars = c("sex"),
#' dist.name = c("BC"))
#'
#' @export
mStat_calculate_adjusted_distance <- function (data.obj,
                                               dist.obj,
                                               adj.vars,
                                               dist.name) {

  if (is.null(dist.name)){
    return()
  }

  # Message to inform the user
  message(
    "Calculating adjusted distances using the provided adjustment variables and distance matrices..."
  )

  # Load metadata and select the required adjustment variables
  meta_tab <-
    data.obj$meta.dat %>% select(all_of(c(adj.vars)))

  # Create a formula string from the adjustment variables and convert it into a formula object
  formula_str <- paste("xx ~", paste(adj.vars, collapse = "+"))
  dynamic_formula <- as.formula(formula_str)

  # Iterate over the distance names to calculate adjusted distances
  adj.dist.obj <- lapply(dist.name, function(sub_dist.name) {

    D <- as.matrix(dist.obj[[sub_dist.name]])

    obj <- suppressWarnings(cmdscale(D, k = nrow(D) - 1, eig = TRUE))

    if (!all(obj$eig > 0)){
      message("Additive constant c* is being added to the non-diagonal dissimilarities to ensure they are Euclidean.")
      obj <- suppressWarnings(cmdscale(D, k = nrow(D) - 1, eig = TRUE, add = TRUE))
    }

    eig <- obj$eig

    s <- sign(eig)

    eig <- abs(eig)

    xx <- obj$points

    # Create an empty matrix to store residuals
    res_matrix <- matrix(0, nrow = nrow(xx), ncol = ncol(xx))

    rownames(res_matrix) <- rownames(xx)

    # Loop through each column of xx to calculate residuals
    for (i in 1:ncol(xx)) {
      column_name <- paste0("MDS", i)
      formula_str <- paste(column_name, "~", paste(adj.vars, collapse = "+"))
      dynamic_formula <- as.formula(formula_str)

      # Bind the current MDS column to meta_tab
      current_data <- cbind(meta_tab, xx[, i, drop = FALSE])
      colnames(current_data)[ncol(current_data)] <- column_name

      model_res <- residuals(lm(dynamic_formula, data = current_data))
      res_matrix[, i] <- model_res
    }

    D.adj <- sqrt(stats::dist(t(t(res_matrix) * sqrt(eig)))^2)

    return(as.dist(D.adj))
  })

  names(adj.dist.obj) <- dist.name

  return(adj.dist.obj)
}
