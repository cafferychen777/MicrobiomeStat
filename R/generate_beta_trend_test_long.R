#' Create Mixed Effects Formula
#'
#' Constructs a mixed-effects model formula based on the provided response, time, group, and subject variables.
#'
#' @param response.var The response variable for the formula.
#' @param time.var The time variable for the formula.
#' @param group.var (Optional) The group variable for the formula.
#' @param subject.var The subject variable for the formula.
#'
#' @return A formula object representing the mixed-effects model, which includes both fixed and random effects parts.
#'
#' @details This function creates a formula object for a mixed-effects model by considering fixed and random effects. The fixed effects part includes the response, time, and optional group variables, while the random effects part includes the time and subject variables.
#' The resulting formula object can be used as an input to mixed-effects modeling functions such as 'lmer' from the lme4 package.
#'
#' @noRd
create_mixed_effects_formula <- function(response.var,
                                        time.var,
                                        group.var = NULL,
                                        subject.var,
                                        adj.vars = NULL,
                                        random_slopes = TRUE) {
  fixed_effects <- response.var
  if (!is.null(group.var)) {
    fixed_effects <- paste(fixed_effects, group.var, sep = " ~ ")
    fixed_effects <- paste(fixed_effects, " * ", time.var, sep = "")
  } else {
    fixed_effects <- paste(fixed_effects, "~", time.var, sep = "")
  }

  if (!is.null(adj.vars) && length(adj.vars) > 0) {
    fixed_effects <- paste(fixed_effects, paste(adj.vars, collapse = " + "), sep = " + ")
  }

  if (random_slopes) {
    random_effects <- paste("(1 +", time.var, "|", subject.var, ")", sep = "")
  } else {
    random_effects <- paste("(1|", subject.var, ")", sep = "")
  }

  as.formula(paste(fixed_effects, random_effects, sep = " + "))
}


#' @noRd
create_fixed_effects_formula <- function(response.var,
                                         time.var,
                                         group.var = NULL,
                                         adj.vars = NULL) {
  fixed_effects <- response.var
  if (!is.null(group.var)) {
    fixed_effects <- paste(fixed_effects, group.var, sep = " ~ ")
    fixed_effects <- paste(fixed_effects, " * ", time.var, sep = "")
  } else {
    fixed_effects <- paste(fixed_effects, "~", time.var, sep = "")
  }

  if (!is.null(adj.vars) && length(adj.vars) > 0) {
    fixed_effects <- paste(fixed_effects, paste(adj.vars, collapse = " + "), sep = " + ")
  }

  stats::as.formula(fixed_effects)
}


#' @noRd
mStat_fit_mixed_effects_model <- function(response.var,
                                          time.var,
                                          group.var,
                                          subject.var,
                                          data,
                                          adj.vars = NULL,
                                          ...,
                                          context = "mixed-effects analysis") {
  formula <- create_mixed_effects_formula(
    response.var = response.var,
    time.var = time.var,
    group.var = group.var,
    subject.var = subject.var,
    adj.vars = adj.vars,
    random_slopes = TRUE
  )

  model_fit <- try(lmer(formula, data = data, ...), silent = TRUE)
  if (!inherits(model_fit, "try-error")) {
    return(model_fit)
  }

  error_message <- conditionMessage(attr(model_fit, "condition"))
  fallback_pattern <- "number of observations.*<=.*number of random effects|singular|number of levels of each grouping factor"
  if (!grepl(fallback_pattern, error_message, ignore.case = TRUE)) {
    stop("Model fitting failed in ", context, ": ", error_message, call. = FALSE)
  }

  message("Simplifying the random-effects structure due to overparameterization.")
  formula <- create_mixed_effects_formula(
    response.var = response.var,
    time.var = time.var,
    group.var = group.var,
    subject.var = subject.var,
    adj.vars = adj.vars,
    random_slopes = FALSE
  )

  model_fit <- try(lmer(formula, data = data, ...), silent = TRUE)
  if (!inherits(model_fit, "try-error")) {
    return(model_fit)
  }

  simplified_error <- conditionMessage(attr(model_fit, "condition"))
  if (!grepl(fallback_pattern, simplified_error, ignore.case = TRUE)) {
    stop(
      "Model fitting failed even after simplifying the random-effects structure in ",
      context,
      ": ",
      simplified_error,
      call. = FALSE
    )
  }

  message("Falling back to a fixed-effects model due to insufficient repeated observations.")
  fixed_formula <- create_fixed_effects_formula(
    response.var = response.var,
    time.var = time.var,
    group.var = group.var,
    adj.vars = adj.vars
  )

  stats::lm(fixed_formula, data = data)
}


#' @noRd
mStat_extract_group_anova_row <- function(anova_result, group.var) {
  if (is.null(group.var) || !nzchar(group.var)) {
    return(NULL)
  }

  anova_df <- as.data.frame(anova_result)
  terms <- rownames(anova_df)
  if (is.null(terms) || length(terms) == 0) {
    return(NULL)
  }

  term_index <- which(terms == group.var)
  if (length(term_index) == 0) {
    return(NULL)
  }

  selected <- anova_df[term_index[[1]], , drop = FALSE]
  data.frame(
    Term = group.var,
    Estimate = NA_real_,
    Std.Error = NA_real_,
    Statistic = selected$`F value`,
    P.Value = selected$`Pr(>F)`,
    check.names = FALSE
  )
}


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
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   adj.vars = c("subject_gender","sample_body_site"),
#'   dist.name = c("BC", "Jaccard")
#' )
#'
#' generate_beta_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
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

    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Inform the user about the importance of numeric time variable
    message(
      "The trend test in 'generate_beta_trend_test_long' relies on a numeric time variable.\n",
      "Please ensure that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
      "The time variable will be converted to numeric within the function if needed."
    )

    # Convert time variable to numeric
    data.obj$meta.dat[[time.var]] <- mStat_coerce_time_to_numeric(
      data.obj$meta.dat[[time.var]],
      time.var = time.var,
      context = "beta trend analysis"
    )

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

      # Convert distance matrix to long format
      dist.df <- as.matrix(dist.obj[[dist.name]])
      dist.df <- mStat_dist_to_tibble(dist.df, sample_col = "sample")
      meta_tab <- mStat_meta_to_tibble(meta_tab, sample_col = "sample")

      # Prepare data for longitudinal analysis
      # This step calculates the distance from each time point to the baseline for each subject
      selected_cols <- c(
        paste0(subject.var, ".subject"),
        paste0(time.var, ".subject"),
        "distance",
        paste0(adj.vars, ".subject")
      )

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var, ".sample")) == min(!!sym(paste0(time.var, ".sample")))) %>%
        filter(!!sym(paste0(time.var, ".subject")) != !!sym(paste0(time.var, ".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(all_of(selected_cols))

      colnames(long.df) <- c(subject.var, time.var, "distance", adj.vars)

      # Add group information to the long format data
      if (!is.null(group.var)) {
        long.df <- long.df %>%
          dplyr::left_join(
            meta_tab %>% dplyr::select(all_of(c(subject.var, group.var))) %>% dplyr::distinct(),
            by = subject.var,
            relationship = "many-to-many"
          )
      }

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
      if (!is.null(group.var) && length(unique(long.df[[group.var]])) > 2) {
        anova_result <- anova(model, type = "III")
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
