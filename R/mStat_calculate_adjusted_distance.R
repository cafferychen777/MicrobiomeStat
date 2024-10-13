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

  # Check if dist.name is provided, if not, return early
  if (is.null(dist.name)){
    return()
  }

  # Inform the user about the calculation process
  message(
    "Calculating adjusted distances using the provided adjustment variables and distance matrices..."
  )

  # Extract metadata and select only the specified adjustment variables
  # This step ensures we only work with the relevant variables for adjustment
  meta_tab <-
    data.obj$meta.dat %>% select(all_of(c(adj.vars)))

  # Create a formula string for the linear model
  # This formula will be used to adjust the distances based on the specified variables
  formula_str <- paste("xx ~", paste(adj.vars, collapse = "+"))
  dynamic_formula <- as.formula(formula_str)

  # Iterate over each specified distance matrix to calculate adjusted distances
  adj.dist.obj <- lapply(dist.name, function(sub_dist.name) {

    # Convert the distance object to a matrix for easier manipulation
    distance_matrix <- as.matrix(dist.obj[[sub_dist.name]])

    # Perform classical multidimensional scaling (MDS) on the distance matrix
    # This step transforms the distances into a set of coordinates in Euclidean space
    mds_result <- suppressWarnings(cmdscale(distance_matrix, k = nrow(distance_matrix) - 1, eig = TRUE))

    # Check if the eigenvalues are all positive
    # If not, add a constant to make the distances Euclidean
    if (!all(mds_result$eig > 0)){
      message("Additive constant c* is being added to the non-diagonal dissimilarities to ensure they are Euclidean.")
      mds_result <- suppressWarnings(cmdscale(distance_matrix, k = nrow(distance_matrix) - 1, eig = TRUE, add = TRUE))
    }

    # Extract eigenvalues from the MDS result
    eigenvalues <- mds_result$eig

    # Store the sign of eigenvalues (not used in current implementation)
    eigenvalue_signs <- sign(eigenvalues)

    # Take the absolute value of eigenvalues
    eigenvalues <- abs(eigenvalues)

    # Extract the MDS coordinates
    mds_coordinates <- mds_result$points

    # Initialize a matrix to store residuals
    # These residuals will represent the adjusted distances
    res_matrix <- matrix(0, nrow = nrow(mds_coordinates), ncol = ncol(mds_coordinates))
    rownames(res_matrix) <- rownames(mds_coordinates)

    # Iterate through each MDS dimension
    for (i in 1:ncol(mds_coordinates)) {
      # Create a column name for the current MDS dimension
      column_name <- paste0("MDS", i)
      
      # Create a formula for the linear model
      # This model will adjust the MDS coordinates based on the specified variables
      formula_str <- paste(column_name, "~", paste(adj.vars, collapse = "+"))
      dynamic_formula <- as.formula(formula_str)

      # Combine the metadata with the current MDS dimension
      current_data <- cbind(meta_tab, mds_coordinates[, i, drop = FALSE])
      colnames(current_data)[ncol(current_data)] <- column_name

      # Fit a linear model and extract residuals
      # These residuals represent the MDS coordinates after adjusting for the specified variables
      model_res <- residuals(lm(dynamic_formula, data = current_data))
      res_matrix[, i] <- model_res
    }

    # Calculate adjusted distances using the residuals
    # This step transforms the adjusted coordinates back into distances
    # The square root of eigenvalues is used to weight the dimensions appropriately
    D.adj <- sqrt(stats::dist(t(t(res_matrix) * sqrt(eigenvalues)))^2)

    # Return the adjusted distance matrix
    return(as.dist(D.adj))
  })

  # Name the list of adjusted distance matrices
  names(adj.dist.obj) <- dist.name

  # Return the list of adjusted distance matrices
  return(adj.dist.obj)
}
