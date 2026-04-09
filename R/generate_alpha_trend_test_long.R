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
    adj_str <- paste(adj.vars, collapse = " + ")
    formula_str <- paste(formula_part, "+", adj_str)
  } else {
    formula_str <- formula_part
  }
  return(as.formula(formula_str))
}

#' @title Alpha Diversity Trend Test (Longitudinal)
#'
#' @description Tests for temporal trends in alpha diversity using linear
#'   mixed-effects models with a numeric time variable.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @return A list object containing the results of the trend test for each alpha
#'   diversity index specified in `alpha.name`. Each element of the list contains a
#'   summary table of the trend test results for a specific alpha diversity index.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("subset_T2D.obj")
#'
#' # Example 1: Test trend for multiple alpha indices with group and no adjustment
#' result1 <- generate_alpha_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "observed_species"),
#'   time.var = "visit_number_num",
#'   subject.var = "subject_id",
#'   group.var = "subject_race"
#' )
#'
#' # Example 2: Test trend for multiple alpha indices with group and adjustment for one covariate
#' result2 <- generate_alpha_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "observed_species"),
#'   time.var = "visit_number_num",
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = "subject_gender"
#' )
#'
#' # Example 3: Test trend for multiple alpha indices with group and adjustment for multiple covariates
#' result3 <- generate_alpha_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "observed_species"),
#'   time.var = "visit_number_num",
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = c("subject_gender", "sample_body_site")
#' )
#'
#' # Example 4: Test trend for multiple alpha indices without group
#' result4 <- generate_alpha_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "observed_species"),
#'   time.var = "visit_number_num",
#'   subject.var = "subject_id"
#' )
#'
#' # Example 5: Test trend for a single alpha index with group
#' result5 <- generate_alpha_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = "chao1",
#'   time.var = "visit_number_num",
#'   subject.var = "subject_id",
#'   group.var = "subject_race"
#' )
#'
#' # Example 6: Test trend using pre-calculated alpha diversity
#' alpha.obj <- mStat_calculate_alpha_diversity(subset_T2D.obj$feature.tab,
#' c("shannon", "observed_species"))
#' result6 <- generate_alpha_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = alpha.obj,
#'   time.var = "visit_number_num",
#'   subject.var = "subject_id",
#'   group.var = "subject_race"
#' )
#' }
#' @export
generate_alpha_trend_test_long <- function(data.obj,
                                           alpha.obj = NULL,
                                           alpha.name = c("shannon","observed_species"),
                                           depth = NULL,
                                           time.var,
                                           subject.var,
                                           group.var = NULL,
                                           adj.vars = NULL) {

  # Exit the function if no alpha diversity indices are specified
  if (is.null(alpha.name)){
    return()
  }

  # Calculate alpha diversity if not provided
  prepared <- mStat_prepare_alpha_inputs(
    data.obj = data.obj,
    alpha.obj = alpha.obj,
    alpha.name = alpha.name,
    depth = depth
  )
  data.obj <- prepared$data.obj
  alpha.obj <- prepared$alpha.obj

  # Inform the user about the importance of numeric time variable for trend test
  # This message ensures that the user understands the requirements for proper analysis
  message(
    "The trend test in 'generate_alpha_trend_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
    "The time variable will be converted to numeric within the function if needed."
  )

  # Convert the time variable to numeric
  # This step is crucial for performing trend analysis
  data.obj$meta.dat[[time.var]] <- mStat_coerce_time_to_numeric(
    data.obj$meta.dat[[time.var]],
    time.var = time.var,
    context = "alpha trend analysis"
  )

  # Extract relevant metadata for the analysis
  meta_tab <-
    data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
      subject.var, group.var, time.var, adj.vars
    )))

  # Combine alpha diversity data with metadata
  # This creates a comprehensive dataset for our analysis
  alpha_df <- mStat_prepare_alpha_data(
    alpha.obj = alpha.obj,
    meta.dat = meta_tab,
    sample_col = "sample",
    join = "inner"
  )

  # Perform statistical tests for each alpha diversity index
  test.list <- lapply(alpha.name, function(index) {

    # Construct the formula for the linear mixed-effects model
    # This model accounts for repeated measures and potential group differences
    formula <- construct_formula(index, group.var, time.var, subject.var, adj.vars)

    # Fit the linear mixed-effects model
    # This model allows for the analysis of longitudinal data with potential confounders
    model <- lmer(formula, data = alpha_df)

    if (!is.null(group.var)){
      # Check if the grouping variable has more than two categories
      if (length(unique(alpha_df[[group.var]])) > 2) {
        # Perform Type III ANOVA for multi-category grouping variables
        # This tests for overall differences among groups, accounting for other variables
        anova_result <- anova(model, type = "III")

        # Extract coefficients from the model
        coef.tab <- extract_coef(model)
        
        # Append the ANOVA result for the grouping variable to the coefficient table
        # This provides both individual coefficient estimates and overall group effects
        last_row <- utils::tail(anova_result, 1)
        var_name <- rownames(last_row)[1]

        adjusted_last_row <- data.frame(
          Term = var_name,
          Estimate = NA,
          Std.Error = NA,
          Statistic = last_row$`F value`,
          P.Value = last_row$`Pr(>F)`
        )

        coef.tab <- rbind(coef.tab, adjusted_last_row)
      } else {
        # For binary grouping variables, extract coefficients directly
        coef.tab <- extract_coef(model)
      }
    } else {
      # If no grouping variable is specified, extract coefficients directly
      coef.tab <- extract_coef(model)
    }

    return(tibble::as_tibble(coef.tab))
  })

  # Assign names to the elements of test.list based on the alpha diversity indices
  names(test.list) <- alpha.name

  return(test.list)
}
