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
#' @param feature.level A vector of strings. Taxonomic level(s) for aggregation, e.g. "Phylum", "Class".
#' @param feature.dat.type A string. Either 'count', 'proportion', or 'other'. Specifies the data type of feature for appropriate transformation.
#' @return A list of test results. The results are returned in a tidy dataframe format, including coefficients, standard errors, statistics, and p-values from linear models and ANOVA tests.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' subset_T2D.obj2 <- subset_T2D.obj
#' subset_T2D.obj2$meta.dat$visit_number <- as.numeric(subset_T2D.obj2$meta.dat$visit_number)
#' generate_taxa_volatility_test_long(
#' data.obj = subset_T2D.obj2,
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "subject_gender",
#' feature.level = c("Phylum","Class"),
#' feature.dat.type = "count"
#' )
#' }
#' @export
generate_taxa_volatility_test_long <- function(data.obj,
                                               time.var,
                                               subject.var,
                                               group.var,
                                               adj.vars = NULL,
                                               feature.level,
                                               feature.dat.type = c("count", "proportion", "other")) {
  # Validate and extract data
  mStat_validate_data(data.obj)

  data.obj <- mStat_process_time_variable(data.obj, time.var)

  message(
    "The volatility calculation in generate_alpha_volatility_test_long relies on a numeric time variable.
         Please check that your time variable is coded as numeric.
         If the time variable is not numeric, it may cause issues in computing the results of the volatility test.
         You can ensure the time variable is numeric by mutating it in the metadata."
  )

  data.obj$meta.dat <-
    data.obj$meta.dat %>% mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

  meta_tab <- load_data_obj_metadata(data.obj) %>%
    dplyr::select(all_of(c(
      time.var, subject.var, group.var, adj.vars
    ))) %>% rownames_to_column("sample")

  if (feature.dat.type %in% c("count", "proportion")) {
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

    sub_test.list <-
      lapply(rownames(otu_tax_agg_numeric), function(taxon) {
        taxa_df <- otu_tax_agg_numeric %>%
          rownames_to_column("taxa") %>%
          gather(key = "sample", value = "value",-taxa) %>%
          filter(taxa == taxon) %>%
          left_join(meta_tab, by = "sample")

        if (!is.null(adj.vars)) {
          # Obtain the residuals by fitting lm(alpha ~ adj.vars)
          residuals_lm <-
            lm(as.formula(paste0(
              "value", "~", paste(adj.vars, collapse = "+")
            )), data = taxa_df)$residuals
        }

        taxa_df$residuals <- residuals_lm

        # Group data by subject and calculate volatility
        volatility_df <- taxa_df %>%
          arrange(!!sym(subject.var),!!sym(time.var)) %>%
          group_by(!!sym(subject.var)) %>%
          mutate(
            diff_residuals = abs(residuals - lag(residuals)),
            diff_time = !!sym(time.var) - lag(!!sym(time.var))
          ) %>%
          filter(!is.na(diff_residuals),!is.na(diff_time)) %>%
          filter(diff_time != 0) %>%
          summarize(volatility = mean(diff_residuals / diff_time),
                    .groups = 'drop')

        test_df <- volatility_df %>%
          left_join(meta_tab %>%
                      select(all_of(c(
                        subject.var, group.var
                      ))) %>%
                      distinct(),
                    by = subject.var)

        # Test the association between the volatility and the grp.var
        formula <- as.formula(paste("volatility ~", group.var))
        test_result <- lm(formula, data = test_df)

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

          coef.tab <-
            rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
        }
        return(as_tibble(coef.tab))
      })

    # Assign names to the elements of test.list
    names(sub_test.list) <- rownames(otu_tax_agg_numeric)
    return(sub_test.list)
  })

  # Assign names to the elements of test.list
  names(test.list) <- feature.level

  return(test.list)
}
