#' @title Alpha Diversity Test Per Time Point (Longitudinal)
#'
#' @description Performs alpha diversity testing at each time point in longitudinal
#'   data using linear models to compare groups.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @return A list where each element corresponds to a different time point and contains the results of alpha diversity tests for that time point.
#'
#' @param group.var Required. Character string specifying the grouping variable in metadata.
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset (without providing alpha.obj)
#' data("ecam.obj")
#'
#' # Perform the longitudinal alpha diversity test for the ECAM dataset
#' alpha_test_results_ecam <- generate_alpha_per_time_test_long(
#'   data.obj = ecam.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "pielou"),
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "delivery",
#'   adj.vars = c("diet")
#' )
#'
#' # Generate dot plots for the ECAM dataset results
#' dot_plots_ecam <- generate_alpha_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = alpha_test_results_ecam,
#'   group.var = "delivery",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#'
#' # Example 2: Analyzing the Type 2 Diabetes (T2D) dataset (with providing alpha.obj)
#' data("subset_T2D.obj")
#'
#' # Calculate alpha diversity indices
#' alpha.obj <- mStat_calculate_alpha_diversity(
#'   x = subset_T2D.obj$feature.tab,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou")
#' )
#'
#' # Perform the longitudinal alpha diversity test for the T2D dataset
#' alpha_test_results_T2D <- generate_alpha_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = alpha.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   group.var = "subject_race",
#'   adj.vars = c("sample_body_site")
#' )
#'
#' # Generate dot plots for the T2D dataset results
#' dot_plots_T2D <- generate_alpha_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = alpha_test_results_T2D,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#'
#' # Example 3: Analyzing the ECAM dataset without adjustment variables
#' alpha_test_results_ecam_no_adj <- generate_alpha_per_time_test_long(
#'   data.obj = ecam.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "pielou"),
#'   depth = 1000,
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "delivery",
#'   adj.vars = NULL
#' )
#'
#' # Generate dot plots for the ECAM dataset results without adjustment variables
#' dot_plots_ecam_no_adj <- generate_alpha_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = alpha_test_results_ecam_no_adj,
#'   group.var = "delivery",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 16,
#'   theme.choice = "minimal"
#' )
#' }
#' @export
generate_alpha_per_time_test_long <- function(data.obj,
                                     alpha.obj = NULL,
                                     alpha.name = NULL,
                                     depth = NULL,
                                     time.var,
                                     t0.level,
                                     ts.levels,
                                     group.var,
                                     adj.vars = NULL) {
  # Validate the input data object to ensure it meets the required format.
  data.obj <- mStat_validate_data(data.obj)

  # If no alpha diversity indices are specified, exit the function.
  if (is.null(alpha.name)){
    return()
  }

  mStat_validate_time_var_contract(
    meta.dat = data.obj$meta.dat,
    time.var = time.var,
    context = "alpha per-time testing"
  )

  prepared <- mStat_prepare_alpha_inputs(
    data.obj = data.obj,
    alpha.obj = alpha.obj,
    alpha.name = alpha.name,
    depth = depth,
    time.var = time.var,
    t0.level = t0.level,
    ts.levels = ts.levels,
    process_time = TRUE
  )
  data.obj <- prepared$data.obj
  alpha.obj <- prepared$alpha.obj

  # Extract relevant metadata for the analysis.
  meta_tab <- mStat_prepare_alpha_meta_tab(
    data.obj = data.obj,
    vars = c(group.var, time.var, adj.vars)
  )

  mStat_validate_group_var_contract(
    meta.dat = meta_tab,
    group.var = group.var,
    context = "alpha per-time testing",
    require_variation = FALSE
  )

  # Determine the reference level for the group variable.
  # This will be used as the baseline for comparisons.
  reference_level <- .mStat_get_group_reference_level(
    meta.dat = meta_tab,
    group.var = group.var
  )

  resolved_time <- mStat_resolve_followup_timepoints(
    values = meta_tab[[time.var]],
    time.var = time.var,
    t0.level = t0.level,
    ts.levels = ts.levels,
    context = "alpha per-time testing"
  )
  t0.level <- resolved_time$t0.level
  ts.levels <- resolved_time$ts.levels

  # Combine all time levels for analysis.
  time.levels <- resolved_time$kept_levels

  # Perform alpha diversity tests for each time point.
  mStat_run_per_time_analysis(
    time.levels = time.levels,
    context = "alpha per-time testing",
    analysis_fn = function(t.level) {
      # Subset the data for the specific time level.
      # This allows for time-specific analysis of alpha diversity.
      subset_inputs <- mStat_subset_analysis_inputs_by_meta_values(
        data.obj = data.obj,
        var = time.var,
        values = t.level,
        alpha.obj = alpha.obj,
        prune.features = FALSE
      )
      subset_data.obj <- subset_inputs$data.obj
      subset_alpha.obj <- subset_inputs$alpha.obj

      subset_meta_tab <- mStat_prepare_alpha_meta_tab(
        data.obj = subset_data.obj,
        vars = c(group.var, time.var, adj.vars)
      )

      mStat_validate_group_var_contract(
        meta.dat = subset_meta_tab,
        group.var = group.var,
        context = paste0("alpha per-time testing at ", t.level)
      )

      # Perform alpha diversity test for the subset data.
      subset.test.list <- generate_alpha_test_single(
        data.obj = subset_data.obj,
        alpha.obj = subset_alpha.obj,
        alpha.name = alpha.name,
        time.var = NULL,
        t.level = NULL,
        group.var = group.var,
        adj.vars = adj.vars
      )

      mStat_restructure_group_term_results(
        test.list = subset.test.list,
        group.var = group.var,
        reference_level = reference_level
      )
    }
  )
}
