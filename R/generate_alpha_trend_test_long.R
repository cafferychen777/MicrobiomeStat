#' Construct a Linear Mixed Effects Model Formula
#'
#' This function constructs a formula for a linear mixed effects model based on the specified
#' index, group variable, time variable, subject variable, and additional covariates.
#' When the `group.var` is NULL, the function will test the slope=0 in the mixed effects model.
#'
#' @param index A character string representing the dependent variable in the model.
#' @param group.var A character string representing the group variable in the model, or NULL if no group variable is included.
#' @param time.var A character string representing the time variable in the model.
#' @param subject.var A character string representing the subject variable in the model.
#' @param adj.vars A character string or vector representing additional covariates to be included in the model, or NULL if no additional covariates are included.
#'
#' @return A formula object suitable for use in functions requiring a linear mixed effects model formula.
#'
#' @noRd
construct_formula <- function(index, group.var, time.var, subject.var, adj.vars) {
  if (!is.null(group.var)) {
    formula_part <- paste(index, "~", group.var, "*", time.var, " + (1 +", time.var, "|", subject.var, ")")
  } else {
    formula_part <- paste(index, "~", time.var, " + (1|", subject.var, ")")
  }

  if (!is.null(adj.vars)) {
    formula_str <- paste(formula_part, "+", adj.vars)
  } else {
    formula_str <- formula_part
  }
  return(as.formula(formula_str))
}

#' Longitudinal Alpha Diversity Volatility Test in MicrobiomeStat
#'
#' This function, part of the MicrobiomeStat package, calculates the volatility
#' of alpha diversity measures in longitudinal data and tests the association
#' between the volatility and a group variable. Volatility is calculated as the mean
#' of absolute differences between consecutive alpha diversity measures, normalized
#' by the time difference (mean(abs(alpha_j - alpha_{j-1}) / (t_j - t_{j-1}))).
#'
#' The function first obtains the residuals by fitting a linear model with alpha
#' diversity as the response and adjustment variables as predictors. It then calculates
#' volatility for each subject and tests its association with a group variable using a
#' linear model. If the group variable specified in `group.var` has more than two levels,
#' an ANOVA is performed to test the association between alpha diversity volatility
#' and the group variable.
#'
#' @param data.obj A list object that includes feature.tab (an OTU table with
#' taxa as rows and samples as columns) and meta.dat (a metadata table with
#' samples as rows and variables as columns).
#' @param alpha.obj An object containing precomputed alpha diversity indices,
#' typically calculated by the function `mStat_calculate_alpha_diversity`.
#' If NULL (the default), alpha diversity will be calculated from the data.obj using
#' `mStat_calculate_alpha_diversity`.
#' @param alpha.name A string with the name of the alpha diversity index to compute.
#' Options could include: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param time.var A string representing the time variable's name in the
#' metadata. The default is NULL.
#' @param subject.var A string indicating the variable for subject identifiers.
#' @param group.var A string representing the group variable's name in the
#' metadata.
#' @param adj.vars A character vector with the names of adjustment variables in
#' the metadata, used in fitting the linear model for residuals.
#' @return A summary object containing the results of the linear model testing
#' the association between alpha diversity volatility and the group variable.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' subset_T2D.obj2 <- subset_T2D.obj
#' subset_T2D.obj2$meta.dat$visit_number <- as.numeric(subset_T2D.obj2$meta.dat$visit_number)
#' generate_alpha_trend_test_long(
#' data.obj = subset_T2D.obj2,
#' alpha.name = c("shannon","simpson"),
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_gender",
#' adj.vars = NULL
#' )
#' }
#' @export
generate_alpha_trend_test_long <- function(data.obj,
                                                alpha.obj = NULL,
                                                alpha.name,
                                                time.var,
                                                subject.var,
                                                group.var,
                                                adj.vars = NULL) {

  if (is.null(alpha.obj)) {
    if (!is_rarefied(data.obj)) {
      message(
        "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj)
    }
    otu_tab <- load_data_obj_count(data.obj)
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
  } else {
    # Verify that all alpha.name are present in alpha.obj
    if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x)
      colnames(x))))) {
      missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
      stop(
        "The following alpha diversity indices are not available in alpha.obj: ",
        paste(missing_alphas, collapse = ", "),
        call. = FALSE
      )
    }
  }

  message(
    "The trend test in 'generate_alpha_trend_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
    "The time variable will be converted to numeric within the function if needed."
  )

  data.obj$meta.dat <- data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
      subject.var, group.var, time.var, adj.vars
    )))

  alpha_df <-
    dplyr::bind_cols(alpha.obj) %>% tibble::rownames_to_column("sample") %>%
    dplyr::inner_join(meta_tab %>% rownames_to_column("sample"),
                      by = c("sample"))

  test.list <- lapply(alpha.name, function(index) {

    formula <- construct_formula(index, group.var, time.var, subject.var, adj.vars)

    model <- lmer(formula, data = alpha_df)

    coef.tab <- extract_coef(model)

    return(as_tibble(coef.tab))
  })

  # Assign names to the elements of test.list
  names(test.list) <- alpha.name

  return(test.list)
}
