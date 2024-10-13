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
  # Perform input validation checks
  # Ensure alpha.obj is a list of data frames containing alpha diversity indices
  if (!is.list(alpha.obj)) {
    stop("`alpha.obj` should be a list of data frames.")
  }

  # Verify that meta.dat is a data frame containing metadata for samples
  if (!is.data.frame(meta.dat)) {
    stop("`meta.dat` should be a data frame.")
  }

  # Check that adj.vars is a non-null character vector specifying covariate names
  if (is.null(adj.vars) || !is.character(adj.vars)) {
    stop("`adj.vars` should be a character vector of covariate names.")
  }

  # Initialize an empty list to store the adjusted alpha diversity data frames
  adjusted.alpha.obj <- list()

  # Iterate over each alpha diversity index and adjust it
  adjusted.alpha.obj <- lapply(names(alpha.obj), function(index_name) {
    # Extract the current alpha diversity data frame
    alpha_df <- alpha.obj[[index_name]]

    # Combine alpha diversity data with metadata
    # This step ensures that we have all necessary information for adjustment
    alpha_df <- dplyr::bind_cols(alpha_df, meta.dat)

    # Verify that all specified covariates are present in the combined data
    if (!all(adj.vars %in% names(alpha_df))) {
      stop("Not all adjustment variables found in the metadata.")
    }

    # Prepare the model matrix for covariates
    # This step handles both continuous and categorical variables
    data_subset <- alpha_df %>%
      dplyr::select(all_of(adj.vars)) %>%
      # Convert character variables to factors for proper handling in the model
      dplyr::mutate(dplyr::across(dplyr::where(is.character) & !is.factor, as.factor))

    # Create a model matrix without intercept
    # This allows for proper handling of categorical variables with multiple levels
    M <- model.matrix(~ 0 + ., data = data_subset)

    # Center the covariates by subtracting their means
    # This step is crucial for interpretability of the adjusted values
    M_centered <- scale(M, scale = FALSE)

    # Fit a linear model to adjust the alpha diversity based on the centered covariates
    # The alpha diversity index is the response variable, and centered covariates are predictors
    fit <- lm(alpha_df[[index_name]] ~ M_centered)

    # Compute adjusted values
    # We add the intercept back to ensure the adjusted values are on the same scale as the original
    # The residuals represent the variation in alpha diversity not explained by the covariates
    adjusted_values <- fit$coefficients[1] + residuals(fit)

    # Update the alpha diversity values in the data frame with the adjusted values
    alpha_df[[index_name]] <- adjusted_values

    # Select only the adjusted alpha diversity column
    # This ensures the output structure matches the input
    alpha_df <- alpha_df %>% select(!!sym(index_name))

    # Return the adjusted data frame for this alpha diversity index
    return(alpha_df)
  })

  # Preserve the original names of the alpha diversity indices in the adjusted list
  names(adjusted.alpha.obj) <- names(alpha.obj)

  # Return the list of adjusted alpha diversity data frames
  return(adjusted.alpha.obj)
}