#' Test Beta Diversity Change Between Time Points
#'
#' Tests within-subject beta diversity changes between two time points using
#' linear models. Supports group comparisons and covariate adjustment.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#'   If NULL, uses the first time point.
#'
#' @examples
#' \dontrun{
#'
#' # Load packages
#' library(vegan)
#'
#' # Load data
#' data(peerj32.obj)
#' generate_beta_change_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = NULL,
#'   change.base = "1",
#'   dist.name = c('BC', 'Jaccard')
#' )
#' generate_beta_change_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' data("subset_pairs.obj")
#' generate_beta_change_test_pair(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   dist.name = c('BC', 'Jaccard')
#' )
#' }
#' @return A named list containing linear modeling results for each beta diversity metric.

#' Each list element corresponds to one of the distance metrics specified in \code{dist.name}.

#' It contains a coefficient table from fitting a linear model with the beta diversity change as
#' response and the \code{group_var} and \code{adj_vars} as predictors.

#' If \code{group_var} has multiple levels, ANOVA results are also included after the coefficients.

#' Column names include:
#' Term, Estimate, Std.Error, Statistic, P.Value
#' @export
#' @name generate_beta_change_test_pair
generate_beta_change_test_pair <-
  function(data.obj,
           dist.obj = NULL,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars = NULL,
           change.base = NULL,
           dist.name = c('BC', 'Jaccard')) {

    # Check if dist.name is provided, if not, return early
    if (is.null(dist.name)){
      return()
    }

    meta_vars <- c(subject.var, group.var, time.var)
    if (!is.null(adj.vars)) {
      meta_vars <- c(meta_vars, adj.vars)
    }

    if (is.null(dist.obj) && !is.null(data.obj)) {
      dist.obj <- mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
    } else {
      prepared_context <- mStat_prepare_precomputed_beta_context(
        dist.obj = dist.obj,
        dist.name = dist.name,
        data.obj = data.obj
      )
      data.obj <- prepared_context$data.obj
      dist.obj <- prepared_context$dist.obj
    }

    meta_tab <- mStat_extract_dist_metadata(
      dist.obj = dist.obj,
      dist.name = dist.name,
      vars = meta_vars,
      data.obj = data.obj
    ) %>% mStat_meta_to_tibble(sample_col = "sample")

    mStat_validate_group_var_contract(
      meta.dat = tibble::column_to_rownames(meta_tab, "sample"),
      group.var = group.var,
      subject.var = subject.var,
      context = "beta change testing"
    )

    # Initialize variable to store time-varying information
    time_varying_info <- list(
      time_varying_vars = character(),
      non_time_varying_vars = character()
    )

    # Handle time-varying covariates
    if (!is.null(adj.vars)){
      # Identify time-varying and non-time-varying variables
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # Adjust distances for non-time-varying variables
      if (length(time_varying_info$non_time_varying_vars) > 0){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj, dist.obj, time_varying_info$non_time_varying_vars, dist.name)
      }
    }

    pair_times <- mStat_resolve_pair_timepoints(
      values = meta_tab[[time.var]],
      time.var = time.var,
      change.base = change.base,
      context = "beta change testing"
    )
    change.base <- pair_times$change.base
    change.after <- pair_times$change.after

    # Perform analysis for each distance metric
    test.list <- lapply(dist.name, function(dist.name){

      long.df <- mStat_prepare_beta_change_long_data(
        dist.matrix = dist.obj[[dist.name]],
        meta.dat = meta_tab,
        subject.var = subject.var,
        time.var = time.var,
        change.base = change.base,
        change.after = change.after
      )

      long.df <- mStat_attach_change_metadata(
        change.df = long.df,
        meta.dat = meta_tab,
        by = c(subject.var, time.var),
        vars = c(group.var, time_varying_info$time_varying_vars)
      )

      predictors <- c(time_varying_info$time_varying_vars, group.var)

      if (!is.null(group.var) && group.var %in% names(long.df)) {
        remaining_group_levels <- unique(stats::na.omit(long.df[[group.var]]))
        if (length(remaining_group_levels) < 2) {
          warning(
            "After data filtering, group variable '", group.var, "' has only ",
            length(remaining_group_levels), " level(s): ",
            paste(remaining_group_levels, collapse = ", "),
            ". Removing it from the model."
          )
          predictors <- setdiff(predictors, group.var)
        }
      }

      valid_terms <- mStat_resolve_variable_terms(
        data = long.df,
        terms = predictors
      )
      model_formula <- mStat_build_formula(
        response = "distance",
        terms = valid_terms
      )

      lm.model <- stats::lm(model_formula, data = long.df)
      coef.tab <- extract_coef(lm.model)

      if (group.var %in% valid_terms && length(unique(stats::na.omit(long.df[[group.var]]))) > 2) {
        group_row <- mStat_extract_group_anova_row(stats::anova(lm.model), group.var)
        if (!is.null(group_row)) {
          coef.tab <- rbind(coef.tab, group_row)
        }
      }

      return(coef.tab)
    })

    names(test.list) <- dist.name

    return(test.list)

  }
