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
    formula_str <- paste(formula_part, "+", adj.vars)
  } else {
    formula_str <- formula_part
  }
  return(as.formula(formula_str))
}

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
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name A string with the name of the alpha diversity index to compute.
#' Options could include: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param time.var Character string specifying the column name in metadata containing the
#'                numeric time variable. Should contain ordered time points for each subject.
#'                Required to calculate volatility over time.
#' @param subject.var Character string specifying the column name in metadata containing
#'                    unique subject IDs. Required to calculate volatility within subjects.
#' @param group.var Character string specifying the column name in metadata containing grouping
#'                 categories. Volatility will be compared between groups using linear models.
#' @param adj.vars Character vector specifying column names in metadata containing covariates
#'                to adjust for when calculating residuals.
#' @return A summary object containing the results of the linear model testing
#' the association between alpha diversity volatility and the group variable.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' subset_T2D.obj2 <- subset_T2D.obj
#' subset_T2D.obj2$meta.dat$visit_number <- as.numeric(subset_T2D.obj2$meta.dat$visit_number)
#' generate_alpha_trend_test_long(
#' data.obj = subset_T2D.obj2,
#' alpha.name = c("shannon","simpson"),
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = NULL
#' )
#' }
#' @export
generate_alpha_trend_test_long <- function(data.obj,
                                           alpha.obj = NULL,
                                           alpha.name = NULL,
                                           depth = NULL,
                                           time.var,
                                           subject.var,
                                           group.var,
                                           adj.vars = NULL) {

  if (is.null(alpha.obj)) {
    if (!is_rarefied(data.obj)) {
      message(
        "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj, depth = depth)
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
    "The trend test in 'generate_alpha_trend_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
    "The time variable will be converted to numeric within the function if needed."
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

    formula <- construct_formula(index, group.var, time.var, subject.var, adj.vars)

    model <- lmer(formula, data = alpha_df)

    # Check if group.var is multi-category
    if (length(unique(alpha_df[[group.var]])) > 2) {
      anova_result <- anova(model, type = "III")

      # Here, I assume you want to append this p-value to the result.
      # Adjust the way of appending the p-value based on your desired output.
      coef.tab <- extract_coef(model)
      # Append the last row of the anova_result to the coef.tab
      last_row <- utils::tail(anova_result, 1)
      # 获取last_row的列名
      var_name <- rownames(last_row)[1]

      # 调整last_row以匹配coef.tab的格式
      adjusted_last_row <- data.frame(
        Term = var_name,
        Estimate = NA,  # 你可以根据需求进行更改
        Std.Error = NA,  # 你可以根据需求进行更改
        Statistic = last_row$`F value`,
        P.Value = last_row$`Pr(>F)`
      )

      # 合并coef.tab和adjusted_last_row
      coef.tab <- rbind(coef.tab, adjusted_last_row)

    } else {
      coef.tab <- extract_coef(model)
    }

    return(as_tibble(coef.tab))
  })

  # Assign names to the elements of test.list
  names(test.list) <- alpha.name

  return(test.list)
}
