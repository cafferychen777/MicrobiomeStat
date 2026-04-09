#' Adjust Alpha Diversity Indices Based on Covariates
#'
#' This function adjusts alpha diversity indices for specified covariates using linear regression. The adjusted alpha diversity values are computed for each index by regressing the original values on the covariates and then using the residuals to adjust the original values.
#'
#' @param alpha.obj A list where each element is a data frame containing alpha diversity indices.
#'        Each data frame should have samples as rows and a single column named after the alpha diversity index.
#' @param meta.dat A data frame containing metadata for the samples.
#'        This data frame should include all the covariates specified in `adj.vars` as well as any other variables needed for the analysis.
#' @param adj.vars A character vector specifying the names of the covariates in `meta.dat`
#'        that should be used to adjust the alpha diversity indices. These variables are used in a linear model to adjust the diversity values.
#' @return A list of data frames, structured identically to `alpha.obj`, where each alpha diversity
#'         index has been adjusted based on the specified covariates. The structure and naming of the
#'         list and its elements are preserved.
#' @details The adjustment is done by fitting a linear model with the alpha diversity index as the
#'          response and the covariates as predictors. The adjusted values are the original values
#'          minus the fitted values plus the intercept. This ensures that the adjusted values are on
#'          the same scale as the original values but have the effect of the covariates removed.
#'
#'          If the covariates include categorical variables, these are converted to factors before
#'          creating the model matrix. The function will stop with an error if any of the specified
#'          covariates are not found in the metadata.
#'
#' @examples
#' \dontrun{
#'   # Assuming you have already calculated alpha diversity and have the metadata ready
#'   alpha.obj <- mStat_calculate_alpha_diversity(peerj32.obj$feature.tab, c("shannon", "simpson"))
#'   # Adjust alpha diversity based on treatment and day
#'   adjusted.alpha.obj <- mStat_calculate_adjusted_alpha_diversity(
#'     alpha.obj = alpha.obj,
#'     meta.dat = peerj32.obj$meta.dat,
#'     adj.vars = c("sex")
#'   )
#' }
#'
#' @export
mStat_calculate_adjusted_alpha_diversity <- function(alpha.obj, meta.dat, adj.vars) {
  if (!is.list(alpha.obj)) {
    stop("`alpha.obj` should be a list of data frames.")
  }

  if (!is.data.frame(meta.dat)) {
    stop("`meta.dat` should be a data frame.")
  }

  if (is.null(adj.vars) || !is.character(adj.vars) || length(adj.vars) == 0) {
    stop("`adj.vars` should be a non-empty character vector of covariate names.")
  }

  alpha_names <- names(alpha.obj)
  mStat_validate_alpha_object(
    alpha.obj = alpha.obj,
    alpha.name = alpha_names
  )

  prepared_alpha <- mStat_prepare_alpha_data(
    alpha.obj = alpha.obj,
    meta.dat = meta.dat,
    vars = adj.vars,
    sample_col = "sample",
    join = "inner"
  )

  if (anyNA(prepared_alpha[, adj.vars, drop = FALSE])) {
    stop(
      "Adjustment variables contain missing values. Please remove or impute missing covariates before adjusting alpha diversity.",
      call. = FALSE
    )
  }

  adjusted.alpha.obj <- lapply(alpha_names, function(index_name) {
    covariate_data <- prepared_alpha[, adj.vars, drop = FALSE]
    covariate_data[] <- lapply(covariate_data, function(column) {
      if (is.character(column) && !is.factor(column)) {
        return(as.factor(column))
      }
      column
    })

    M <- stats::model.matrix(~ 0 + ., data = covariate_data)
    M_centered <- scale(M, scale = FALSE)
    fit <- stats::lm(prepared_alpha[[index_name]] ~ M_centered)
    adjusted_values <- unname(stats::coef(fit)[1] + stats::residuals(fit))

    adjusted_alpha <- data.frame(
      adjusted_values,
      row.names = prepared_alpha$sample,
      check.names = FALSE
    )
    colnames(adjusted_alpha) <- index_name
    adjusted_alpha
  })

  names(adjusted.alpha.obj) <- alpha_names
  adjusted.alpha.obj
}
