#' Perform pairwise beta diversity tests using PERMANOVA
#'
#' This function calculates various beta diversity indices for a microbiome dataset and performs pairwise PERMANOVA tests on them to assess the effect of a grouping variable and any additional covariates on the beta diversity. It returns the PERMANOVA results for each of the specified beta diversity indices.
#'
#' @name generate_beta_test_pair
#' @param data.obj A list object in a format specific to MicrobiomeStat, containing microbiome data and metadata.
#' @param dist.obj A distance object, if NULL it will be calculated based on data.obj and dist.name. This distance object should be a dissimilarity matrix or an object that can be coerced to a dissimilarity matrix.
#' @param time.var The name of the column in data.obj$meta.dat that contains the time variable for subsetting data, if applicable. This should be a character string representing the column name.
#' @param subject.var The name of the column in data.obj$meta.dat that contains the subject variable for the PERMANOVA tests. This should be a character string representing the column name.
#' @param group.var The name of the column in data.obj$meta.dat that contains the grouping variable for the PERMANOVA tests. This should be a character string representing the column name.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param dist.name A character vector specifying which beta diversity indices to calculate and test (e.g., 'BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'). Default is c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS').
#' @examples
#'
#' library(vegan)
#' library(GUniFrac)
#' library(ape)
#' library(philentropy)
#' library(MicrobiomeStat)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#'
#' # Perform pairwise beta diversity tests using PERMANOVA
#' beta_test_pair_results <- generate_beta_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   subject.var = "subject",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#'
#' @return A list containing the PERMANOVA results for each beta diversity index. The list includes two elements: "p.tab" - a table of p-values for the PERMANOVA tests across all indices, and "aov.tab" - a table containing detailed PERMANOVA results for each index. The p.tab and aov.tab tables include columns for the terms in the PERMANOVA model, the degrees of freedom, sums of squares, mean squares, F statistics, R-squared values, and p-values.
#' @export
#' @seealso
#' \code{\link[GUniFrac]{PermanovaG2}}
#' \code{\link[microbiome]{peerj32}}
#' \code{\link[vegan]{adonis2}}
generate_beta_test_pair <- function(data.obj,
                                    dist.obj = NULL,
                                    time.var,
                                    subject.var,
                                    group.var,
                                    adj.vars = NULL,
                                    dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

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
      across(where(is.numeric), ~ round(., 3))  # round all numeric values to 3 decimal places
    )

  # Format aov.tab
  aov.tab <- bind_rows(result$aov.tab.list) %>%
    as_tibble(rownames = "Variable") %>%
    dplyr::mutate(
      Variable = gsub("data.obj\\$meta.dat\\[\\[\"", "", Variable),
      Variable = gsub("\"\\]\\]", "", Variable),
      Variable = gsub("\\...\\d+$", "", Variable), # Remove "...1", "...2", etc.
      Variable = ifelse(Variable == group.var, group.var, Variable),
      Distance = rep(dist.name, each = n() / length(dist.name))
    ) %>%
    dplyr::rename(
      `DF` = Df,
      `Sum_Sq` = SumsOfSqs,
      `Mean_Sq` = MeanSqs,
      `F_Statistic` = F.Model,
      `R_Squared` = R2,
      `P_Value` = `Pr(>F)`
    ) %>%
    dplyr::mutate(across(where(is.numeric), ~ ifelse(is.na(.), "NA", round(., 3))))

  return(list("p.tab" = p.tab, "aov.tab" = aov.tab))
}
