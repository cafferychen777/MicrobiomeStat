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

  mStat_inform_numeric_time_requirement(
    function_name = "generate_alpha_trend_test_long",
    analysis_label = "trend analysis",
    conversion_behavior = "coerce"
  )

  # Convert the time variable to numeric
  # This step is crucial for performing trend analysis
  data.obj$meta.dat[[time.var]] <- mStat_coerce_time_to_numeric(
    data.obj$meta.dat[[time.var]],
    time.var = time.var,
    context = "alpha trend analysis"
  )

  # Extract relevant metadata for the analysis
  meta_tab <- mStat_prepare_alpha_meta_tab(
    data.obj = data.obj,
    vars = c(subject.var, group.var, time.var, adj.vars)
  )

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

    model <- mStat_fit_mixed_effects_model(
      response.var = index,
      time.var = time.var,
      group.var = group.var,
      subject.var = subject.var,
      data = alpha_df,
      adj.vars = adj.vars,
      context = "alpha trend analysis"
    )

    coef.tab <- extract_coef(model)
    if (!is.null(group.var) && length(unique(stats::na.omit(alpha_df[[group.var]]))) > 2) {
      anova_result <- mStat_compute_group_anova(model)
      group_row <- mStat_extract_group_anova_row(anova_result, group.var)
      if (!is.null(group_row)) {
        coef.tab <- rbind(coef.tab, group_row)
      }
    }

    return(tibble::as_tibble(coef.tab))
  })

  # Assign names to the elements of test.list based on the alpha diversity indices
  names(test.list) <- alpha.name

  return(test.list)
}
