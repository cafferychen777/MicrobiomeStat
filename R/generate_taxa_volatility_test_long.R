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
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param time.var Character string specifying the column name in metadata containing
#'                the numeric time variable. Should contain ordered time points for each
#'                subject. Required to calculate volatility over time.
#' @param subject.var Character string specifying the column name in metadata containing
#'                    unique subject IDs. Required to calculate volatility within subjects
#'                    over time.
#' @param group.var Character string specifying the column name in metadata containing
#'                 grouping categories. Volatility will be compared between groups using
#'                 linear models. Required.
#' @param adj.vars Character vector specifying column names in metadata containing covariates
#'                to adjust for in linear models. Optional, can be NULL.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' @param feature.level Character vector specifying taxonomic level(s) to aggregate abundance data to
#'                     before volatility calculation, e.g. c("Phylum", "Genus"). The special value
#'                     "original" can also be provided, which will use the original taxon identifiers.
#' @param feature.dat.type Character string specifying the data type of the abundance data. Should be
#'                        one of "count", "proportion", or "other". Determines transform. This should
#'                        match the units of data used in feature.level.
#' @param transform Character string specifying transformation method. If "CLR", count and
#'                 proportion data will be CLR transformed before volatility calculation.
#'                 Default "CLR".
#' @param ... Additional arguments passed to other methods.
#' @return A list of test results. The results are returned in a tidy dataframe format, including coefficients, standard errors, statistics, and p-values from linear models and ANOVA tests.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' test.list <- generate_taxa_volatility_test_long(
#' data.obj = subset_T2D.obj,
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "sample_body_site",
#' prev.filter = 0.1,
#' abund.filter = 0.0001,
#' feature.level = c("Genus"),
#' feature.dat.type = "count",
#' transform = "CLR"
#' )
#' plot.list <- generate_taxa_volatility_volcano_long(data.obj = subset_T2D.obj,
#'                                                    group.var = "subject_race",
#'                                                    test.list = test.list,
#'                                                    feature.sig.level = 0.1,
#'                                                    feature.mt.method = "none")
#'
#' data("ecam.obj")
#' test.list <- generate_taxa_volatility_test_long(
#'   data.obj = ecam.obj,
#'   time.var = "month_num",
#'   subject.var = "subject.id",
#'   group.var = "antiexposedall",
#'   adj.vars = "delivery",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.level = c("Order", "Family", "Genus"),
#'   feature.dat.type = "proportion",
#'   transform = "CLR"
#' )
#' plot.list <- generate_taxa_volatility_volcano_long(
#'   data.obj = ecam.obj,
#'   group.var = "antiexposedall",
#'   test.list = test.list,
#'   feature.sig.level = 0.2,
#'   feature.mt.method = "none"
#' )
#' }
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
                                               transform = "CLR",
                                               ...) {
  # Validate and extract data
  mStat_validate_data(data.obj)

  feature.dat.type <- match.arg(feature.dat.type)

  message(
    "The volatility calculation relies on a numeric time variable.\n",
    "Please check that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
    "You can ensure the time variable is numeric by mutating it in the metadata."
  )

  data.obj$meta.dat <-
    data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

  meta_tab <- data.obj$meta.dat %>%
    dplyr::select(all_of(c(
      time.var, subject.var, group.var, adj.vars
    ))) %>% rownames_to_column("sample")

  group_level <- meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

  reference_level <- group_level[1]

  if (transform == "CLR"){
    abund.filter <- 0
  }

  # Create a formula including the group variable and adjustment variables (if any)
  formula_str <- paste("volatility ~", group.var)
  if (!is.null(adj.vars)) {
    formula_str <- paste(formula_str, "+", paste(adj.vars, collapse = " + "))
  }
  formula <- as.formula(formula_str)

  test.list <- lapply(feature.level, function(feature.level) {

    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
      data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    if (feature.level != "original"){
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- data.obj$feature.tab
    }

    # Calculate the average abundance of each group
    prop_prev_data <-
      otu_tax_agg %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::group_by(Var1) %>%
      dplyr::summarise(
        avg_abundance = mean(Freq),
        prevalence = sum(Freq > 0) / dplyr::n()
      ) %>% column_to_rownames("Var1") %>%
      rownames_to_column(feature.level)

    # Function for performing CLR transformation on data
    clr_transform <- function(x) {
      gm = exp(mean(log(x)))
      return(log(x / gm))
    }

    impute_zeros_rowwise <- function(row) {
      min_nonzero <- min(row[row > 0])
      row[row == 0] <- min_nonzero / 2
      return(row)
    }

    otu_tax_agg_imputed_temp <- apply(otu_tax_agg, 1, impute_zeros_rowwise)

    # Transpose the matrix to return the rows and columns to their original positions
    otu_tax_agg_imputed <- t(otu_tax_agg_imputed_temp)

    # Application CLR conversion
    otu_tax_agg_clr <- apply(otu_tax_agg_imputed, 1, clr_transform)

    # Transpose back
    otu_tax_agg_clr <- t(otu_tax_agg_clr)

    otu_tax_agg_clr_long <- otu_tax_agg_clr %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter) %>%
      rownames_to_column(feature.level) %>%
      tidyr::gather(key = "sample", value = "value",-feature.level)

    taxa.levels <- otu_tax_agg_clr_long %>% select(all_of(feature.level)) %>% pull() %>% unique()

    sub_test.list <-
      lapply(taxa.levels, function(taxon) {

        taxa_df <- otu_tax_agg_clr_long %>%
          dplyr::filter(!!sym(feature.level) == taxon) %>%
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

        # Run the linear model
        test_result <- lm(formula, data = test_df)

        coef.tab <- extract_coef(test_result)

        # Run ANOVA on the model if group.var is multi-categorical
        if (length(unique(taxa_df[[group.var]])) > 2) {
          anova <- anova(test_result)
          anova.tab <- anova %>% as.data.frame() %>%
            rownames_to_column("Term") %>%
            dplyr::select(
                   Term,
                   Statistic = `F value`,
                   P.Value = `Pr(>F)`) %>%
            as_tibble() %>%
            dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
            select(
              Term,
              Estimate,
              Std.Error,
              Statistic,
              P.Value
            )

          coef.tab <- rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
        }
        return(as_tibble(coef.tab))
      })

    # Assign names to the elements of test.list
    names(sub_test.list) <- otu_tax_agg_clr_long %>% select(all_of(feature.level)) %>% pull() %>% unique()

    # Find all unique terms
    unique_terms <- grep(paste0("^", group.var, "$|^", group.var, ".*"), unique(unlist(lapply(sub_test.list, function(df) unique(df$Term)))), value = TRUE)

    # Extract data for each Term and store it in a new list
    result_list <- lapply(unique_terms, function(term) {
      do.call(rbind, lapply(sub_test.list, function(df) {
        df %>% dplyr::filter(Term == term)
      })) %>%
        dplyr::mutate(!!sym(feature.level) := names(sub_test.list)) %>%
        dplyr::left_join(prop_prev_data, by = feature.level) %>%
        dplyr::mutate(Adjusted.P.Value = p.adjust(P.Value, method = "fdr")) %>%
        dplyr::select(all_of(c(feature.level, "Estimate", "Std.Error", "P.Value", "Adjusted.P.Value", "avg_abundance", "prevalence"))) %>%
        dplyr::rename(Coefficient = Estimate,
                      SE = Std.Error,
                      Variable = feature.level,
                      Mean.Abundance = avg_abundance,
                      Prevalence = prevalence)
    })

    # Name the new list
    names(result_list) <- unique_terms

    new_names <- sapply(names(result_list), function(name) {
      # Check whether the name matches the specified pattern and is not the result of ANOVA
      if (grepl(paste0("^", group.var), name) && !grepl(paste0("^", group.var, "$"), name)) {
        sub_name <- sub(paste0(group.var), "", name)
        return(paste(sub_name, "vs", reference_level, "(Reference)"))
      }
      return(name)
    })

    names(result_list) <- new_names

    return(result_list)
  })

  names(test.list) <- feature.level

  return(test.list)
}
