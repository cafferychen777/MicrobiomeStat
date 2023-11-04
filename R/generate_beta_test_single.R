#' Calculate Beta Diversity Indices and Perform PERMANOVA Tests
#'
#' This function computes a range of beta diversity indices for a microbiome dataset and evaluates them using
#' PERMANOVA tests. The goal is to assess the effect of a specific grouping variable along with potential covariates
#' on beta diversity. The function outputs the PERMANOVA results for the selected beta diversity indices.
#'
#' @name generate_beta_test_single
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param time.var The column name in metadata containing the time variable for subsetting data before analysis. Default is NULL.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var The column name in metadata containing the grouping variable of interest for PERMANOVA test.
#' @param adj.vars Character vector of column names in metadata to be included as covariates in PERMANOVA test. Default is NULL.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Perform beta diversity tests using PERMANOVA
#' generate_beta_test_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   t.level = NULL,
#'   group.var = "group",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Perform beta diversity tests using PERMANOVA
#' generate_beta_test_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   t.level = NULL,
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Perform beta diversity tests using PERMANOVA
#' generate_beta_test_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Perform beta diversity tests using PERMANOVA
#' generate_beta_test_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   t.level = "2",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' data(ecam.obj)
#'
#' # Perform beta diversity tests using PERMANOVA
#' generate_beta_test_single(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   time.var = "month",
#'   t.level = "0",
#'   group.var = "delivery",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Perform beta diversity tests using PERMANOVA
#' generate_beta_test_single(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   time.var = "month",
#'   t.level = "0",
#'   group.var = "delivery",
#'   adj.vars = "diet",
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, dist.name = c('BC', 'Jaccard'))
#' generate_beta_test_single(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = "delivery",
#'   adj.vars = "diet",
#'   dist.name = c('BC', 'Jaccard')
#' )
#' }
#'
#' @return A list containing the PERMANOVA results for each beta diversity index and an omnibus p-value. The list includes two elements: "p.tab" - a table of p-values for the PERMANOVA tests dplyr::across all indices, and "aov.tab" - a table containing detailed PERMANOVA results for each index. The p.tab and aov.tab tables include columns for the terms in the PERMANOVA model, the degrees of freedom, sums of squares, mean squares, F statistics, R-squared values, and p-values.
#' @export
generate_beta_test_single <- function(data.obj,
                                      dist.obj = NULL,
                                      time.var = NULL,
                                      t.level = NULL,
                                      group.var,
                                      adj.vars = NULL,
                                      dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

  if (is.null(dist.name)){
    return()
  }

  # Calculate beta diversity indices if not provided
  if (is.null(dist.obj)) {
    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
    dist.obj <- mStat_calculate_beta_diversity(data.obj, dist.name)
  } else {
    message("Using provided dist.obj...")
    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
    dist.obj <- mStat_subset_dist(dist.obj, rownames(data.obj$meta.dat))
  }

  # Run PermanovaG2 for the entire dist.obj
  message("Running PermanovaG2 for all distances...")

  # Create the formula for PermanovaG2
  if (is.null(adj.vars)) {
    formula_str <- paste0("dist.obj ~ ", group.var)
  } else {
    adj_vars_terms <- paste0(adj.vars, collapse = " + ")
    formula_str <- paste0("dist.obj ~ ", adj_vars_terms, " + ", group.var)
  }

  formula <- as.formula(formula_str)

  # Run PermanovaG2 for the entire dist.obj
  message("Running PermanovaG2 for all distances...")
  result <- GUniFrac::PermanovaG2(formula, data = data.obj$meta.dat)

  # Initialize a list to store the PERMANOVA results in tibble format
  permanova.results <- list()

  # Convert the p.tab results to tibble
  permanova.results$p.tab <- result$p.tab %>%
    as_tibble(rownames = "Term") %>%
    dplyr::mutate(Term = ifelse(Term == "data.obj$meta.dat[[group.var]]", group.var, Term))

  # Convert each aov.tab result to tibble and store them in the list
  for (i in seq_along(result$aov.tab.list)) {
    permanova.results$aov.tab[[dist.name[i]]] <- result$aov.tab.list[[i]] %>%
      as_tibble(rownames = "Variable") %>%
      dplyr::mutate(
        Variable = gsub("data.obj\\$meta.dat\\[\\[\"", "", Variable),
        Variable = gsub("\"\\]\\]", "", Variable),
        Variable = ifelse(Variable == paste0("data.obj$meta.dat[[group.var]]"), group.var, Variable),
        Distance = dist.name[i]
      ) %>%
      dplyr::rename(
        `DF` = Df,
        `Sum.Sq` = SumsOfSqs,
        `Mean.Sq` = MeanSqs,
        `F.Statistic` = F.Model,
        `R.Squared` = R2,
        `P.Value` = `Pr(>F)`
      )
  }

  # Format p.tab
  p.tab <- permanova.results$p.tab %>%
    dplyr::mutate(
      Term = gsub("data.obj\\$meta.dat\\[\\[\"", "", Term),
      Term = gsub("\"\\]\\]", "", Term),
      dplyr::across(where(is.numeric), ~ round(., 3))  # round all numeric values to 3 decimal places
    )

  # Format aov.tab
  aov.tab <- dplyr::bind_rows(permanova.results$aov.tab) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.), "NA", round(., 3))))  # round all numeric values to 3 decimal places, replace NA with "NA"

  aov.tab <- aov.tab %>% dplyr::select(Distance, everything())

  return(list("p.tab" = p.tab, "aov.tab" = aov.tab))
}
