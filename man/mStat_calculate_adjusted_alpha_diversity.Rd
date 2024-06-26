% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mStat_calculate_adjusted_alpha_diversity.R
\name{mStat_calculate_adjusted_alpha_diversity}
\alias{mStat_calculate_adjusted_alpha_diversity}
\title{Adjust Alpha Diversity Indices Based on Covariates}
\usage{
mStat_calculate_adjusted_alpha_diversity(alpha.obj, meta.dat, adj.vars)
}
\arguments{
\item{alpha.obj}{A list where each element is a data frame containing alpha diversity indices.
Each data frame should have samples as rows and a single column named after the alpha diversity index.}

\item{meta.dat}{A data frame containing metadata for the samples.
This data frame should include all the covariates specified in `adj.vars` as well as any other variables needed for the analysis.}

\item{adj.vars}{A character vector specifying the names of the covariates in `meta.dat`
that should be used to adjust the alpha diversity indices. These variables are used in a linear model to adjust the diversity values.}
}
\value{
A list of data frames, structured identically to `alpha.obj`, where each alpha diversity
        index has been adjusted based on the specified covariates. The structure and naming of the
        list and its elements are preserved.
}
\description{
This function adjusts alpha diversity indices for specified covariates using linear regression. The adjusted alpha diversity values are computed for each index by regressing the original values on the covariates and then using the residuals to adjust the original values.
}
\details{
The adjustment is done by fitting a linear model with the alpha diversity index as the
         response and the covariates as predictors. The adjusted values are the original values
         minus the fitted values plus the intercept. This ensures that the adjusted values are on
         the same scale as the original values but have the effect of the covariates removed.

         If the covariates include categorical variables, these are converted to factors before
         creating the model matrix. The function will stop with an error if any of the specified
         covariates are not found in the metadata.
}
\examples{
\dontrun{
  # Assuming you have already calculated alpha diversity and have the metadata ready
  alpha.obj <- mStat_calculate_alpha_diversity(peerj32.obj$feature.tab, c("shannon", "simpson"))
  # Adjust alpha diversity based on treatment and day
  adjusted.alpha.obj <- mStat_calculate_adjusted_alpha_diversity(
    alpha.obj = alpha.obj,
    meta.dat = peerj32.obj$meta.dat,
    adj.vars = c("sex")
  )
}

}
