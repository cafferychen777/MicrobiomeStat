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
create_mixed_effects_formula <- function(response.var, time.var, group.var = NULL, subject.var, random_slopes = TRUE) {
  # Formula for the fixed effects part
  fixed_effects <- response.var
  if (!is.null(group.var)) {
    fixed_effects <- paste(fixed_effects, group.var, sep = " ~ ")
    fixed_effects <- paste(fixed_effects, " * ", time.var, sep = "")
  } else {
    fixed_effects <- paste(fixed_effects, "~", time.var, sep = "")
  }

  # Formula for constructing the random effects part
  if (random_slopes) {
    random_effects <- paste("(1 +", time.var, "|", subject.var, ")", sep = "")
  } else {
    random_effects <- paste("(1|", subject.var, ")", sep = "")
  }

  # Merge fixed effects and random effects parts
  formula_str <- paste(fixed_effects, random_effects, sep = " + ")

  # Convert to formula object using as.formula
  formula_obj <- as.formula(formula_str)

  # Return formula object
  return(formula_obj)
}


#' @title Generate Beta Diversity Trend Test for Longitudinal Data
#'
#' @description The function `generate_beta_trend_test_long` performs a linear mixed effects model to test the longitudinal trend in beta diversity. It models the distance matrix as the response, with time as a fixed effect and subject as a random effect. An optional grouping variable can be included as an interaction with time. Covariates for adjustment can also be added to the model. The time variable should be numeric, and the function coerces it as needed. The linear mixed model is fitted using lmer, and the estimated coefficients are extracted and returned, representing the trends over time and differences between groups.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param subject.var The variable name in metadata that represents the subject identifiers. Should match the subject identifiers used when calculating distance matrix.
#' @param time.var The variable name in metadata that represents the time points. Should be a numeric vector coded as continuous values. The function will try to coerce to numeric if needed.
#' @param group.var (Optional) The variable name in metadata that represents the grouping factor of interest. Default is NULL.
#' @param adj.vars (Optional) A character vector of variable names in metadata to be used for adjustment in the model. Default is NULL.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Default is "BC" (Bray-Curtis). Supported indices include "BC", "Jaccard", "UniFrac", "GUniFrac", "WUniFrac", and "JS".
#' @param ... (Optional) Additional arguments to pass to internal functions.
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
    mStat_validate_data(data.obj)

    # Inform the user about the importance of numeric time variable
    message(
      "The trend test in 'generate_alpha_trend_test_long' relies on a numeric time variable.\n",
      "Please ensure that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
      "The time variable will be converted to numeric within the function if needed."
    )

    # Convert time variable to numeric
    data.obj$meta.dat <- data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

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
      # If distance object is provided, extract metadata from the appropriate source
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      }
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
      dist.df <- dist.df %>%
        as.data.frame() %>%
        rownames_to_column("sample")
      meta_tab <- meta_tab %>% rownames_to_column("sample")

      # Prepare data for longitudinal analysis
      # This step calculates the distance from each time point to the baseline for each subject
      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == min(!!sym(paste0(time.var,".sample")))) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      # Add group information to the long format data
      long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(all_of(c(subject.var, group.var))) %>% dplyr::distinct(), by = subject.var, relationship = "many-to-many")

      # Ensure time variable is numeric
      long.df <- long.df %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

      # Attempt to fit the mixed effects model with random slopes
      # This model accounts for individual variations in the trend over time
      formula <- create_mixed_effects_formula(
        response.var = "distance",
        time.var = time.var,
        group.var = group.var,
        subject.var = subject.var,
        random_slopes = TRUE
      )

      model_fit <- try(lmer(formula, data = long.df), silent = TRUE)

      # If the model fitting fails, attempt to simplify the random effects structure
      if (inherits(model_fit, "try-error")) {
        error_message <- attr(model_fit, "condition")$message
        if (grepl("number of observations.*<=.*number of random effects", error_message, ignore.case = TRUE)) {
          message("Simplifying the random-effects structure due to overparameterization.")
          # Simplify the random effects structure by removing random slopes
          formula <- create_mixed_effects_formula(
            response.var = "distance",
            time.var = time.var,
            group.var = group.var,
            subject.var = subject.var,
            random_slopes = FALSE
          )
          model_fit <- try(lmer(formula, data = long.df), silent = TRUE)
          if (inherits(model_fit, "try-error")) {
            stop("Model fitting failed even after simplifying the random-effects structure: ", model_fit)
          }
        } else {
          stop("Model fitting failed due to an error: ", model_fit)
        }
      }

      model <- model_fit

      # Extract and format model results
      if (!is.null(group.var)){
        # For multi-category group variables, perform Type III ANOVA
        if (length(unique(long.df[[group.var]])) > 2) {
          anova_result <- anova(model, type = "III")

          # Extract coefficients and ANOVA results
          coef.tab <- extract_coef(model)
          last_row <- utils::tail(anova_result, 1)
          var_name <- rownames(last_row)[1]

          # Adjust ANOVA results to match coefficient table format
          adjusted_last_row <- data.frame(
            Term = var_name,
            Estimate = NA,
            Std.Error = NA,
            Statistic = last_row$`F value`,
            P.Value = last_row$`Pr(>F)`
          )

          # Combine coefficient table with ANOVA results
          coef.tab <- rbind(coef.tab, adjusted_last_row)

        } else {
          # For binary group variables, extract coefficients directly
          coef.tab <- extract_coef(model)
        }
      } else {
        # If no group variable, extract coefficients
        coef.tab <- extract_coef(model)
      }

      return(as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list
    names(test.list) <- dist.name

    return(test.list)
  }