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

  message(
    "Calculating adjusted distances using the provided adjustment variables and distance matrices..."
  )

  meta_tab <- mStat_extract_dist_metadata(
    dist.obj = dist.obj,
    dist.name = dist.name,
    vars = adj.vars,
    data.obj = data.obj
  )

  if (anyNA(meta_tab)) {
    stop(
      "Adjustment variables contain missing values. Please remove or impute missing covariates before adjusting distances.",
      call. = FALSE
    )
  }

  adj.dist.obj <- lapply(dist.name, function(sub_dist.name) {
    distance_matrix <- as.matrix(dist.obj[[sub_dist.name]])

    mds_result <- withCallingHandlers(
      stats::cmdscale(distance_matrix, k = nrow(distance_matrix) - 1, eig = TRUE),
      warning = function(cond) {
        if (grepl("eigenvalue", conditionMessage(cond), ignore.case = TRUE))
          invokeRestart("muffleWarning")
      }
    )

    if (!all(mds_result$eig > 0)){
      message("Additive constant c* is being added to the non-diagonal dissimilarities to ensure they are Euclidean.")
      mds_result <- withCallingHandlers(
        stats::cmdscale(distance_matrix, k = nrow(distance_matrix) - 1, eig = TRUE, add = TRUE),
        warning = function(cond) {
          if (grepl("eigenvalue", conditionMessage(cond), ignore.case = TRUE))
            invokeRestart("muffleWarning")
        }
      )
    }

    mds_coordinates <- mds_result$points
    meta_current <- meta_tab[rownames(mds_coordinates), , drop = FALSE]

    positive_axes <- mds_result$eig[seq_len(ncol(mds_coordinates))] > 0
    if (!any(positive_axes)) {
      stop(
        "Unable to compute adjusted distances because no positive MDS axes were retained for ",
        sub_dist.name,
        ".",
        call. = FALSE
      )
    }

    mds_coordinates <- mds_coordinates[, positive_axes, drop = FALSE]
    axis_weights <- sqrt(mds_result$eig[seq_len(ncol(mds_result$points))][positive_axes])

    res_matrix <- matrix(NA_real_, nrow = nrow(mds_coordinates), ncol = ncol(mds_coordinates))
    rownames(res_matrix) <- rownames(mds_coordinates)

    for (i in seq_len(ncol(mds_coordinates))) {
      column_name <- paste0("MDS", i)
      dynamic_formula <- stats::as.formula(paste(column_name, "~", paste(adj.vars, collapse = "+")))

      current_data <- cbind(meta_current, mds_coordinates[, i, drop = FALSE])
      colnames(current_data)[ncol(current_data)] <- column_name
      model_frame <- stats::model.frame(dynamic_formula, data = current_data, na.action = stats::na.exclude)
      model_fit <- stats::lm(dynamic_formula, data = model_frame, na.action = stats::na.exclude)
      res_matrix[, i] <- stats::residuals(model_fit)
    }

    weighted_coordinates <- sweep(res_matrix, 2, axis_weights, `*`)
    D.adj <- stats::dist(weighted_coordinates)
    mStat_attach_dist_metadata(D.adj, data.obj$meta.dat)
  })

  names(adj.dist.obj) <- dist.name
  adj.dist.obj
}
