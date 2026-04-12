#' Beta Diversity Change Tests Per Time Point
#'
#' Performs beta diversity change tests at each time point relative to baseline.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param ... Additional arguments passed to underlying functions.
#' @param group.var Required. Character string specifying the grouping variable in metadata.
#'
#' @return A list of test results, structured to facilitate further analysis and visualization.
#'
#' @seealso Related functions in MicrobiomeStat for data preparation and beta diversity calculation,
#'          such as \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}},
#'          \code{\link[MicrobiomeStat]{mStat_calculate_PC}}, and data conversion functions like
#'          \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}.
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' result1 <- generate_beta_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_beta_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result1,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#' )
#'
#' generate_beta_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = c("sample_body_site", "subject_gender"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, c('BC', 'Jaccard'))
#' result2 <- generate_beta_change_per_time_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   subject.var = "subject.id",
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_beta_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = result2,
#'   group.var = "delivery",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 15
#' )
#' }
#' @export
generate_beta_change_per_time_test_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           subject.var,
           group.var,
           adj.vars = NULL,
           dist.name = c('BC', 'Jaccard'),
           ...) {

    # Check if distance names are provided
    if (is.null(dist.name)){
      return()
    }

    if (is.null(data.obj) && is.null(dist.obj)) {
      stop("Either `data.obj` or `dist.obj` must be provided.", call. = FALSE)
    }

    if (!is.null(data.obj)) {
      data.obj <- mStat_validate_data(data.obj)
    }

    # Calculate beta diversity if not provided
    if (is.null(dist.obj)) {
      mStat_validate_time_var_contract(
        meta.dat = data.obj$meta.dat,
        time.var = time.var,
        context = "beta change per-time testing"
      )

      # Process time variable to ensure proper ordering
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      # Extract relevant metadata
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var)))
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Adjust distance matrix if adjustment variables are provided
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      prepared_context <- mStat_prepare_precomputed_beta_context(
        dist.obj = dist.obj,
        dist.name = dist.name,
        data.obj = data.obj,
        time.var = time.var,
        t0.level = t0.level,
        ts.levels = ts.levels,
        process_time = TRUE
      )
      data.obj <- prepared_context$data.obj
      dist.obj <- prepared_context$dist.obj
      meta_tab <- mStat_extract_dist_metadata(
        dist.obj = dist.obj,
        dist.name = dist.name,
        vars = c(subject.var, time.var, group.var),
        data.obj = data.obj
      )
    }

    mStat_validate_time_var_contract(
      meta.dat = meta_tab,
      time.var = time.var,
      context = "beta change per-time testing"
    )

    mStat_validate_group_var_contract(
      meta.dat = meta_tab,
      group.var = group.var,
      subject.var = subject.var,
      context = "beta change per-time testing",
      require_variation = FALSE
    )

    # Determine the reference level for the grouping variable
    reference_level <- .mStat_get_group_reference_level(
      meta.dat = meta_tab,
      group.var = group.var
    )

    resolved_time <- mStat_resolve_followup_timepoints(
      values = meta_tab[[time.var]],
      time.var = time.var,
      t0.level = t0.level,
      ts.levels = ts.levels,
      context = "beta change per-time testing"
    )
    t0.level <- resolved_time$t0.level
    ts.levels <- resolved_time$ts.levels

    # Perform beta diversity change tests for each follow-up time point
    mStat_run_per_time_analysis(
      time.levels = ts.levels,
      context = "beta change per-time testing",
      analysis_fn = function(ts.level) {
        # Create subsets of data and distance objects
        subset_inputs <- mStat_subset_analysis_inputs_by_meta_values(
          data.obj = data.obj,
          var = time.var,
          values = c(t0.level, ts.level),
          dist.obj = dist.obj,
          prune.features = TRUE
        )
        subset_data.obj <- subset_inputs$data.obj
        subset_dist.obj <- subset_inputs$dist.obj

        # Perform beta diversity change test for the current time point pair
        subset.test.list <- generate_beta_change_test_pair(
          data.obj = subset_data.obj,
          dist.obj = subset_dist.obj,
          time.var = time.var,
          subject.var = subject.var,
          group.var = group.var,
          adj.vars = adj.vars,
          change.base = t0.level,
          dist.name = dist.name
        )

        mStat_restructure_group_term_results(
          test.list = subset.test.list,
          group.var = group.var,
          reference_level = reference_level
        )
      }
    )
  }
