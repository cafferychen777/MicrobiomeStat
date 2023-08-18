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
#' generate_alpha_volatility_test_long(
#' data.obj = subset_T2D.obj2,
#' alpha.obj = NULL,
#' alpha.name = c("shannon","simpson"),
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "sample_body_site"
#' )
#' }
#' @export
generate_alpha_volatility_test_long <- function(data.obj,
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
    "The volatility calculation in generate_alpha_volatility_test_long relies on a numeric time variable.\n",
    "Please check that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
    "You can ensure the time variable is numeric by mutating it in the metadata."
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

  if (!is.null(adj.vars)){
    # Obtain the residuals by fitting lm(alpha ~ adj.vars)
    residuals_lm <- stats::lm(as.formula(paste0(index, "~", paste(adj.vars, collapse = "+"))), data = alpha_df)$residuals
  } else {
    residuals_lm <- alpha_df %>% dplyr::select(all_of(c(index))) %>% dplyr::pull()
  }

  # Calculate the volatility for each subject
  # Add residuals to alpha_df
  alpha_df$residuals <- residuals_lm

  # Group data by subject and calculate volatility
  volatility_df <- alpha_df %>%
    dplyr::group_by(!!sym(subject.var)) %>%
    dplyr::arrange(!!sym(time.var)) %>%
    dplyr::mutate(diff_residuals = abs(residuals - dplyr::lag(residuals)),
           diff_time = !!sym(time.var) - dplyr::lag(!!sym(time.var))) %>%
    dplyr::filter(!is.na(diff_residuals), !is.na(diff_time)) %>%
    dplyr::filter(diff_time != 0) %>%
    dplyr::summarize(volatility = mean(diff_residuals / diff_time), .groups = 'drop')

  test_df <- volatility_df %>%
    dplyr::left_join(meta_tab %>%
                select(all_of(c(subject.var, group.var))) %>%
                  dplyr::distinct(),
              by = subject.var)

  # Test the association between the volatility and the grp.var
  formula <- as.formula(paste("volatility ~", group.var))

  test_result <- stats::lm(formula, data = test_df)

  coef.tab <- extract_coef(test_result)

  # Run ANOVA on the model if group.var is multi-categorical
  if (length(unique(alpha_df[[group.var]])) > 2) {
    anova.tab <- broom::tidy(anova(test_result))

    # Rearrange the table and add missing columns
    anova.tab <- anova.tab %>%
      select(
        term = term,
        Statistic = statistic,
        df = df,
        P.Value = p.value
      ) %>%
      dplyr::mutate(Estimate = NA, Std.Error = NA)

    # Reorder the columns to match coef.tab
    anova.tab <- anova.tab %>%
      select(
        Term = term,
        Estimate = Estimate,
        Std.Error = Std.Error,
        Statistic = Statistic,
        P.Value = P.Value
      )

    coef.tab <-
      rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
  }

  return(as_tibble(coef.tab))
  })

  # Assign names to the elements of test.list
  names(test.list) <- alpha.name

  return(test.list)
}
