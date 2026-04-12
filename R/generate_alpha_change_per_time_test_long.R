#' @title Alpha Diversity Change Test Per Time Point (Longitudinal)
#'
#' @description Performs paired tests comparing alpha diversity changes between
#'   baseline and each follow-up time point in longitudinal data.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. Options include 'log fold change', 'absolute change',
#'   or a custom function.
#'
#' @return A list of tables, one for each alpha diversity metric, summarizing the results of the statistical tests.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data("subset_pairs.obj")
#'
#' alpha_test_results <- generate_alpha_change_per_time_test_long(
#'   data.obj = subset_pairs.obj,
#'   alpha.name = c("shannon", "simpson"),
#'   time.var = "Antibiotic",
#'   t0.level = unique(subset_pairs.obj$meta.dat$Antibiotic)[1],
#'   ts.levels = unique(subset_pairs.obj$meta.dat$Antibiotic)[-1],
#'   subject.var = "MouseID",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   alpha.change.func = "absolute change"
#' )
#'
#' generate_alpha_per_time_dotplot_long(
#'   data.obj = subset_pairs.obj,
#'   test.list = alpha_test_results,
#'   group.var = "Sex",
#'   time.var = "Antibiotic",
#'   t0.level = unique(subset_pairs.obj$meta.dat$Antibiotic)[1],
#'   ts.levels = unique(subset_pairs.obj$meta.dat$Antibiotic)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#' }
#' @export
generate_alpha_change_per_time_test_long <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = NULL,
           depth = NULL,
           time.var,
           t0.level,
           ts.levels,
           subject.var,
           group.var,
           adj.vars = NULL,
           alpha.change.func = "log fold change") {

    data.obj <- mStat_validate_data(data.obj)

    if (is.null(alpha.name)){
      return()
    }

    mStat_validate_time_var_contract(
      meta.dat = data.obj$meta.dat,
      time.var = time.var,
      context = "alpha change per-time testing"
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

    meta_tab <- mStat_prepare_alpha_meta_tab(
      data.obj = data.obj,
      vars = c(subject.var, group.var, time.var, adj.vars)
    )

    mStat_validate_group_var_contract(
      meta.dat = meta_tab,
      group.var = group.var,
      subject.var = subject.var,
      context = "alpha change per-time testing",
      require_variation = FALSE
    )

    reference_level <- .mStat_get_group_reference_level(
      meta.dat = meta_tab,
      group.var = group.var
    )

    resolved_time <- mStat_resolve_followup_timepoints(
      values = meta_tab[[time.var]],
      time.var = time.var,
      t0.level = t0.level,
      ts.levels = ts.levels,
      context = "alpha change per-time testing"
    )
    t0.level <- resolved_time$t0.level
    ts.levels <- resolved_time$ts.levels

    # Generate tests
    mStat_run_per_time_analysis(
      time.levels = ts.levels,
      context = "alpha change per-time testing",
      analysis_fn = function(ts.level) {
        subset_inputs <- mStat_subset_analysis_inputs_by_meta_values(
          data.obj = data.obj,
          var = time.var,
          values = c(t0.level, ts.level),
          alpha.obj = alpha.obj,
          prune.features = TRUE
        )
        subset_data.obj <- subset_inputs$data.obj
        subset_alpha.obj <- subset_inputs$alpha.obj

        subset.test.list <- generate_alpha_change_test_pair(
          data.obj = subset_data.obj,
          alpha.obj = subset_alpha.obj,
          alpha.name = alpha.name,
          depth = depth,
          time.var = time.var,
          subject.var = subject.var,
          group.var = group.var,
          adj.vars = adj.vars,
          change.base = t0.level,
          alpha.change.func = alpha.change.func
        )

        mStat_restructure_group_term_results(
          test.list = subset.test.list,
          group.var = group.var,
          reference_level = reference_level
        )
      }
    )
  }
