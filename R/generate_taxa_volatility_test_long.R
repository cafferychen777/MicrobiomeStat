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
#' @param data.obj A MicrobiomeStat data object. The heart of the MicrobiomeStat, consisting of several key components:
#'   \itemize{
#'     \item \strong{Feature.tab (Matrix)}: Meeting point for research objects (OTU/ASV/KEGG/Gene, etc.) and samples.
#'     \item \strong{Meta.dat (Data frame)}: Rows correspond to the samples, and columns serve as annotations, describing the samples.
#'     \item \strong{Feature.ann (Matrix)}: Annotations, carrying classification information like Kingdom, Phylum, etc.
#'     \item \strong{Phylogenetic tree (Optional)}: An evolutionary perspective, illuminating the relationships among various research objects.
#'     \item \strong{Feature.agg.list (Optional)}: Aggregated results based on the feature.tab and feature.ann.
#'   }
#' @param time.var A string. The name of the time variable in the metadata.
#' @param subject.var A string. The name of the subject variable in the metadata.
#' @param group.var A string. The grouping variable to test, found in metadata.
#' @param adj.vars A vector of strings. Covariates for adjustment in the linear models. Defaults to NULL.
#' @param prev.filter A numeric value indicating the minimum prevalence for a feature to be tested.
#' @param abund.filter A numeric value indicating the minimum abundance for a feature to be tested.
#' @param feature.level A vector of strings. Taxonomic level(s) for aggregation, e.g. "Phylum", "Class".
#' @param feature.dat.type A string. Either 'count', 'proportion', or 'other'. Specifies the data type of feature for appropriate transformation.
#' @return A list of test results. The results are returned in a tidy dataframe format, including coefficients, standard errors, statistics, and p-values from linear models and ANOVA tests.
#'
#' @examples
#' data("subset_T2D.obj")
#' subset_T2D.obj2 <- subset_T2D.obj
#' subset_T2D.obj2$meta.dat$visit_number <- as.numeric(subset_T2D.obj2$meta.dat$visit_number)
#' generate_taxa_volatility_test_long(
#' data.obj = subset_T2D.obj,
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_gender",
#' adj.vars = NULL,
#' prev.filter = 1e-7,
#' abund.filter = 0,
#' feature.level = c("Phylum"),
#' feature.dat.type = "count"
#' )
#' @export
generate_taxa_volatility_test_long <- function(data.obj,
                                               time.var,
                                               subject.var,
                                               group.var,
                                               adj.vars = NULL,
                                               prev.filter = 0,
                                               abund.filter = 0,
                                               feature.level,
                                               feature.dat.type = c("count", "proportion", "other"),
                                               Transform = "CLR",
                                               ...) {
  # Validate and extract data
  mStat_validate_data(data.obj)

  message(
    "The volatility calculation in generate_taxa_volatility_test_long relies on a numeric time variable.\n",
    "Please check that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
    "You can ensure the time variable is numeric by mutating it in the metadata."
  )

  data.obj$meta.dat <-
    data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

  meta_tab <- load_data_obj_metadata(data.obj) %>%
    dplyr::select(all_of(c(
      time.var, subject.var, group.var, adj.vars
    ))) %>% rownames_to_column("sample")

  if (feature.dat.type %in% c("count", "proportion")) {
    otu_tab <-
      load_data_obj_count(mStat_normalize_data(data.obj, method = Transform)$data.obj.norm)
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

  if (Transform == "CLR"){
    abund.filter <- 0
  }

  test.list <- lapply(feature.level, function(feature.level) {

    otu_tax <-
      cbind(otu_tab,
            tax_tab %>% select(all_of(feature.level)))

    # 聚合 OTU 表
    otu_tax_filtered <- otu_tax %>%
      tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
      dplyr::group_by_at(vars(!!sym(feature.level))) %>%
      dplyr::summarise(total_count = mean(value),
                       prevalence = sum(value > 0) / dplyr::n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count,-prevalence) %>%
      dplyr::left_join(otu_tax, by = feature.level)

    otu_tax_agg <- otu_tax_filtered %>%
      tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
      dplyr::group_by_at(vars(sample,!!sym(feature.level))) %>%
      dplyr::summarise(value = sum(value)) %>%
      tidyr::spread(key = "sample", value = "value")

    # 转换计数为数值类型
    otu_tax_agg_numeric <-
      otu_tax_agg %>%
      dplyr::mutate_at(vars(-!!sym(feature.level)), as.numeric) %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "unclassified")) %>% column_to_rownames(feature.level)

    sub_test.list <-
      lapply(rownames(otu_tax_agg_numeric), function(taxon) {
        taxa_df <- otu_tax_agg_numeric %>%
          rownames_to_column("taxa") %>%
          tidyr::gather(key = "sample", value = "value",-taxa) %>%
          dplyr::filter(taxa == taxon) %>%
          dplyr::left_join(meta_tab, by = "sample")

        # Group data by subject and calculate volatility
        volatility_df <- taxa_df %>%
          dplyr::arrange(!!sym(subject.var),!!sym(time.var)) %>%
          dplyr::group_by(!!sym(subject.var)) %>%
          dplyr::mutate(
            diff_value = abs(value - dplyr::lag(value)),
            diff_time = !!sym(time.var) - dplyr::lag(!!sym(time.var))
          ) %>%
          dplyr::filter(!is.na(diff_value),!is.na(diff_time)) %>%
          dplyr::filter(diff_time != 0) %>%
          dplyr::summarize(volatility = mean(diff_value / diff_time),
                    .groups = 'drop')

        test_df <- volatility_df %>%
          dplyr::left_join(meta_tab %>%
                             select(all_of(c(subject.var, group.var, adj.vars))) %>%
                             dplyr::distinct(),
                           by = subject.var)

        # Create a formula including the group variable and adjustment variables (if any)
        formula_str <- paste("volatility ~", group.var)
        if (!is.null(adj.vars)) {
          formula_str <- paste(formula_str, "+", paste(adj.vars, collapse = " + "))
        }
        formula <- as.formula(formula_str)

        # Run the linear model
        test_result <- lm(formula, data = test_df, ...)

        coef.tab <- extract_coef(test_result)

        # Run ANOVA on the model if group.var is multi-categorical
        if (length(unique(taxa_df[[group.var]])) > 2) {
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

          coef.tab <- rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
        }
        return(as_tibble(coef.tab))
      })

    # Assign names to the elements of test.list
    names(sub_test.list) <- rownames(otu_tax_agg_numeric)
    return(sub_test.list)
  })

  names(test.list) <- feature.level

  return(test.list)
}
