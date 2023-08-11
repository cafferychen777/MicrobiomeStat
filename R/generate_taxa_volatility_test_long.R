#' Longitudinal Taxa Abundance Volatility Test
#'
#' This function calculates the volatility of taxa abundances in longitudinal microbiome data.
#' It tests for association between abundance volatility and a grouping variable.
#'
#' Volatility is calculated as the mean absolute difference in abundance between
#' consecutive time points, normalized by time difference:
#'   mean(|abundance(t+1) - abundance(t)| / (time(t+1) - time(t)))
#'
#' The function transforms the abundance data first before volatility calculation.
#' Default transform is 'CLR' for count and proportion data. No transform for other types.
#'
#' For count data, a pseudocount of 0.5 is added before CLR transform.
#' For proportion data, zeros are replaced with 1/2 of the minimum positive value before CLR.
#'
#' It then calculates volatility within each subject, and tests for association with
#' the grouping variable using linear models. If the grouping variable has multiple levels,
#' an ANOVA is performed.
#'
#' @param data.obj List containing OTU table and metadata
#' @param feature.level Taxonomic level for aggregation
#' @param time.var Time variable name in metadata
#' @param subject.var Subject variable name in metadata
#' @param group.var Grouping variable to test
#' @param adj.vars Covariates for adjustment
#' @param feature.dat.type Either 'count' or 'proportion'
#' @param transform Transform to apply before volatility calculation
#'
#' @return Tidy dataframe of test results
#'
#' @examples
#'
#' generate_alpha_volatility_test_long(
#' data.obj = subset_T2D.obj2,
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "subject_gender",
#' feature.level = "Phylum",
#' feature.dat.type = "count"
#' )
#'
#' @export
generate_taxa_volatility_test_long <- function(data.obj,
                                               time.var,
                                               subject.var,
                                               group.var,
                                               adj.vars = NULL,
                                               feature.level,
                                               feature.dat.type = c("count", "proportion","other")
                                               ) {

  # Validate and extract data
  mStat_validate_data(data.obj)

  data.obj <- mStat_process_time_variable(data.obj, time.var)

  meta_tab <- load_data_obj_metadata(data.obj) %>%
    dplyr::select(all_of(c(time.var, subject.var, group.var, adj.vars)))

  if (feature.dat.type %in% c("count","proportion")){
    otu_tab <-
      load_data_obj_count(mStat_normalize_data(data.obj, method = "CLR")$data.obj.norm)
  } else {
    otu_tab <- load_data_obj_count(data.obj)
  }

  tax_tab <- load_data_obj_taxonomy(data.obj) %>%
    as.data.frame() %>%
    {
      if ("original" %in% feature.level)
        dplyr::mutate(., original = rownames(.))
      else
        .
    } %>%
    select(all_of(feature.level))

  otu_tax <-
    cbind(otu_tab, tax_tab)

  test.list <- lapply(feature.level, function(feature.level) {

    # 聚合 OTU 表
    otu_tax_agg <- otu_tax %>%
      tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(sample, !!sym(feature.level))) %>%
      dplyr::summarise(count = sum(count)) %>%
      tidyr::spread(key = "sample", value = "count")

    # 转换计数为数值类型
    otu_tax_agg_numeric <-
      otu_tax_agg %>%
      dplyr::mutate_at(vars(-!!sym(feature.level)), as.numeric) %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "unclassified")) %>% column_to_rownames(feature.level)

    lapply(rownames(otu_tax_agg_numeric), function(taxon){

    })
  # Calculate volatility
  volatility_df <- otu_tab %>%
    as.data.frame() %>%
    rownames_to_column(feature.level) %>%
    tidyr::gather(key = "sample", value = "abundance", -!!sym(feature.level)) %>%
    dplyr::inner_join(meta_tab, by = "sample") %>%
    arrange(!!sym(subject.var), !!sym(time.var)) %>%
    group_by(!!sym(subject.var), !!sym(feature.level)) %>%
    mutate(diff_abundance = abs(abundance - lag(abundance)),
           diff_time = !!sym(time.var) - lag(!!sym(time.var))) %>%
    filter(!is.na(diff_abundance), !is.na(diff_time)) %>%
    filter(diff_time != 0) %>%
    summarize(volatility = mean(diff_abundance/diff_time), .groups = 'drop')

  # Test association with group variable
  test_df <- volatility_df %>%
    left_join(meta_tab %>%
                dplyr::select(!!sym(subject.var), !!sym(group.var)) %>%
                distinct(),
              by = subject.var)

  formula <- as.formula(paste("volatility ~", group.var))

  test_result <- lm(formula, data = test_df)

  return(broom::tidy(test_result))
}
