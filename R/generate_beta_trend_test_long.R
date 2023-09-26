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
create_mixed_effects_formula <- function(response.var, time.var, group.var = NULL, subject.var) {
  # 构建固定效应部分的公式
  fixed_effects <- response.var
  if (!is.null(group.var)) {
    fixed_effects <- paste(fixed_effects, group.var, sep = " ~ ")
    fixed_effects <- paste(fixed_effects, " * ", time.var, sep = "")
  } else {
    fixed_effects <- paste(fixed_effects, "~", time.var, sep = "")
  }

  # 构建随机效应部分的公式
  random_effects <- paste("(1 +", time.var, "|", subject.var, ")", sep = "")

  # 合并固定效应和随机效应部分
  formula_str <- paste(fixed_effects, random_effects, sep = " + ")

  # 用as.formula转换为公式对象
  formula_obj <- as.formula(formula_str)

  # 返回公式对象
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
#' data(esubset_T2D.obj)
#' generate_beta_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   adj.vars = c("subject_gender","sample_body_site"),
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

    mStat_validate_data(data.obj)

    message(
      "The trend test in 'generate_beta_trend_test_long' relies on a numeric time variable.\n",
      "Please ensure that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
      "The time variable will be processed within the function if needed."
    )

    if (is.null(dist.obj)) {
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      }
    }

    if (nrow(as.matrix(dist.obj[[dist.name[1]]])) > nrow(meta_tab)){
      samIDs <- rownames(meta_tab)
      dist.obj <- mStat_subset_dist(dist.obj = dist.obj, samIDs = samIDs)
    }

    test.list <- lapply(dist.name,function(dist.name){

      dist.df <- as.matrix(dist.obj[[dist.name]])

      dist.df <- dist.df %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      meta_tab <- meta_tab %>% rownames_to_column("sample")

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == min(meta_tab[,time.var])) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(all_of(c(subject.var, group.var))) %>% dplyr::distinct(), by = subject.var, relationship = "many-to-many")

      formula <- create_mixed_effects_formula(response.var = "distance",
                                              time.var = time.var,
                                              group.var = group.var,
                                              subject.var = subject.var)

      long.df <- long.df %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

      model <- lmer(formula, data = long.df)

      # Check if group.var is multi-category
      if (length(unique(long.df[[group.var]])) > 2) {
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
          Estimate = NA,
          Std.Error = NA,
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
    names(test.list) <- dist.name

    return(test.list)
  }
