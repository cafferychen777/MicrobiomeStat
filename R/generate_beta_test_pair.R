#' Perform pairwise beta diversity tests using PERMANOVA
#'
#' This function calculates various beta diversity indices for a microbiome dataset and performs pairwise PERMANOVA tests on them to assess the effect of a grouping variable and any additional covariates on the beta diversity. It returns the PERMANOVA results for each of the specified beta diversity indices.
#'
#' @name generate_beta_test_pair
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param time.var Character string specifying the column name in metadata containing a time variable.
#'                Used for subsetting data before beta diversity calculation if needed. Optional.
#' @param subject.var Character string specifying the column name in metadata containing unique subject
#'                    IDs. Used as a stratifying variable in PERMANOVA tests. Required.
#' @param group.var Character string specifying the column name in metadata containing grouping
#'                 categories. Used as the main predictor in PERMANOVA tests. Required.
#' @param adj.vars Character vector specifying column names in metadata containing covariates
#'                to include in PERMANOVA model. Optional, can be NULL if no adjustment needed.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#'
#' # Perform pairwise beta diversity tests using PERMANOVA
#' generate_beta_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   subject.var = "subject",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#' }
#'
#' @return A list containing the PERMANOVA results for each beta diversity index. The list includes two elements: "p.tab" - a table of p-values for the PERMANOVA tests dplyr::across all indices, and "aov.tab" - a table containing detailed PERMANOVA results for each index. The p.tab and aov.tab tables include columns for the terms in the PERMANOVA model, the degrees of freedom, sums of squares, mean squares, F statistics, R-squared values, and p-values.
#' @export
#' @seealso
#' \code{\link[GUniFrac]{PermanovaG2}}
#' \code{\link[vegan]{adonis2}}
generate_beta_test_pair <- function(data.obj,
                                    dist.obj = NULL,
                                    time.var,
                                    subject.var,
                                    group.var,
                                    adj.vars = NULL,
                                    dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

  if (is.null(dist.name)){
    return()
  }

  # Calculate beta diversity indices if not provided
  if (is.null(dist.obj)) {
    dist.obj <- mStat_calculate_beta_diversity(data.obj, dist.name)
  } else {
    message("Using provided dist.obj...")
  }

  # Run PermanovaG2 for the entire dist.obj
  message("Running PermanovaG2 for all distances...")

  # Create the formula for PermanovaG2
  if (is.null(adj.vars)) {
    formula_str <- paste0("~ ", group.var)
  } else {
    adj_vars_terms <- paste0(adj.vars, collapse = " + ")
    formula_str <- paste0("~ ", adj_vars_terms, " + ", group.var)
  }
  # Add time.var to the formula
  formula_str <- paste0("dist.obj", formula_str, " + ", time.var)
  formula <- as.formula(formula_str)

  # Run PermanovaG2 for the entire dist.obj
  message("Running PermanovaG2 for all distances...")
  result <- GUniFrac::PermanovaG2(formula, data = data.obj$meta.dat, strata = data.obj$meta.dat[[subject.var]])

  # Format p.tab
  p.tab <- result$p.tab %>%
    as_tibble(rownames = "Term") %>%
    dplyr::mutate(
      Term = gsub("data.obj\\$meta.dat\\[\\[\"", "", Term),
      Term = gsub("\"\\]\\]", "", Term),
      Term = ifelse(Term == group.var, group.var, Term),
      dplyr::across(where(is.numeric), ~ round(., 3))  # round all numeric values to 3 decimal places
    )

  # Format aov.tab
  aov.tab <- dplyr::bind_rows(result$aov.tab.list) %>%
    as_tibble(rownames = "Variable") %>%
    dplyr::mutate(
      Variable = gsub("data.obj\\$meta.dat\\[\\[\"", "", Variable),
      Variable = gsub("\"\\]\\]", "", Variable),
      Variable = gsub("\\...\\d+$", "", Variable), # Remove "...1", "...2", etc.
      Variable = ifelse(Variable == group.var, group.var, Variable),
      Distance = rep(dist.name, each = dplyr::n() / length(dist.name))
    ) %>%
    dplyr::rename(
      `DF` = Df,
      `Sum.Sq` = SumsOfSqs,
      `Mean.Sq` = MeanSqs,
      `F.Statistic` = F.Model,
      `R.Squared` = R2,
      `P.Value` = `Pr(>F)`
    ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.), "NA", round(., 3))))

  aov.tab <- aov.tab %>% select(Distance, everything())

  return(list("p.tab" = p.tab, "aov.tab" = aov.tab))
}
