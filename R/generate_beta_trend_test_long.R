#' Beta Diversity Trend Test for Longitudinal Data
#'
#' Performs linear mixed effects models to test longitudinal trends in beta diversity.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing the result of the trend test for each specified beta diversity index. The result includes a tibble with the coefficients extracted from the mixed-effects model fitted for each distance.
#'
#' @details The function starts by validating the input data, followed by processing the time variable and calculating the beta diversity if necessary. Adjustments are made based on the provided adjusting variables, and the mixed-effects model is fitted to the long-format data. The coefficients of the model are extracted and returned for each beta diversity index specified.
#'
#' @note A warning message will be displayed to ensure that the time variable is coded as numeric. Non-numeric coding may lead to issues in the trend test computation.
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_beta_trend_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   group.var = "diet",
#'   adj.vars = c("antiexposedall","delivery"),
#'   dist.name = c("BC", "Jaccard")
#'   )
#'
#' generate_beta_trend_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   group.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = c("BC", "Jaccard")
#'   )
#'
#' data(subset_T2D.obj)
#' generate_beta_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   adj.vars = c("subject_gender"),
#'   dist.name = c("BC", "Jaccard")
#' )
#'
#' generate_beta_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = c("BC", "Jaccard")
#' )
#' }
#' @export
generate_beta_trend_test_long <-
  function(data.obj,
           dist.obj = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           adj.vars = NULL,
           dist.name = c("BC"),
           ...) {

    # Check if distance metrics are provided
    if (is.null(dist.name)){
      return()
    }

    if (is.null(data.obj) && is.null(dist.obj)) {
      stop("Either `data.obj` or `dist.obj` must be provided.", call. = FALSE)
    }

    # Validate the input data object only when it is actually used.
    if (!is.null(data.obj)) {
      data.obj <- mStat_validate_data(data.obj)
    }

    mStat_inform_numeric_time_requirement(
      function_name = "generate_beta_trend_test_long",
      analysis_label = "trend analysis",
      conversion_behavior = "coerce"
    )

    # Convert time variable to numeric when raw metadata is available locally.
    if (!is.null(data.obj)) {
      data.obj$meta.dat[[time.var]] <- mStat_coerce_time_to_numeric(
        data.obj$meta.dat[[time.var]],
        time.var = time.var,
        context = "beta trend analysis"
      )
    }

    # If distance object is not provided, calculate it from the data object
    if (is.null(dist.obj)) {
      # Extract relevant metadata
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # If adjustment variables are provided, calculate adjusted distances
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      prepared_context <- mStat_prepare_precomputed_beta_context(
        dist.obj = dist.obj,
        dist.name = dist.name,
        data.obj = data.obj
      )
      data.obj <- prepared_context$data.obj
      dist.obj <- prepared_context$dist.obj
      meta_tab <- mStat_extract_dist_metadata(
        dist.obj = dist.obj,
        dist.name = dist.name,
        vars = c(subject.var, time.var, group.var, adj.vars),
        data.obj = data.obj
      )
    }

    # Ensure the distance object and metadata have matching dimensions
    if (nrow(as.matrix(dist.obj[[dist.name[1]]])) > nrow(meta_tab)){
      samIDs <- rownames(meta_tab)
      dist.obj <- mStat_subset_dist(dist.obj = dist.obj, samIDs = samIDs)
    }

    # Perform trend test for each distance metric
    test.list <- lapply(dist.name, function(dist.name){

      # Convert baseline-to-follow-up distances to long format and attach
      # follow-up metadata for model covariates.
      long.df <- mStat_prepare_beta_trend_long_data(
        dist.matrix = dist.obj[[dist.name]],
        meta.dat = meta_tab,
        subject.var = subject.var,
        time.var = time.var,
        vars = c(group.var, adj.vars)
      )

      # Ensure time variable is numeric
      long.df[[time.var]] <- mStat_coerce_time_to_numeric(
        long.df[[time.var]],
        time.var = time.var,
        context = "beta trend analysis"
      )

      model <- mStat_fit_mixed_effects_model(
        response.var = "distance",
        time.var = time.var,
        group.var = group.var,
        subject.var = subject.var,
        data = long.df,
        adj.vars = adj.vars,
        context = "beta trend analysis"
      )

      # Extract and format model results
      coef.tab <- extract_coef(model)
      if (!is.null(group.var) && length(unique(stats::na.omit(long.df[[group.var]]))) > 2) {
        anova_result <- mStat_compute_group_anova(model)
        group_row <- mStat_extract_group_anova_row(anova_result, group.var)
        if (!is.null(group_row)) {
          coef.tab <- rbind(coef.tab, group_row)
        }
      }

      return(tibble::as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list
    names(test.list) <- dist.name

    return(test.list)
  }
