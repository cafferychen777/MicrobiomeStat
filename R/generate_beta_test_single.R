#' Perform beta diversity tests using PERMANOVA
#'
#' This function calculates various beta diversity indices for a microbiome dataset and performs PERMANOVA tests on them to assess the effect of a grouping variable and any additional covariates on the beta diversity. It returns the PERMANOVA results for each of the specified beta diversity indices.
#'
#' @name generate_beta_test_single
#' @param data.obj A list object in a format specific to MicrobiomeStat, containing a 'meta.dat' component with metadata and a 'feature.tab' component with taxa in rows and samples in columns.
#' @param dist.obj A distance object, usually a dissimilarity matrix computed based on the microbiome data in 'data.obj'. If NULL, it will be calculated using 'dist.name'.
#' @param time.var The name of the column in the 'meta.dat' component of 'data.obj' that contains the time variable for subsetting data. Should be NULL if no time variable is needed.
#' @param t.level The specific level of the time variable to use when subsetting data. Should be NULL if no time variable is used.
#' @param group.var The name of the column in the 'meta.dat' component of 'data.obj' that contains the grouping variable for the PERMANOVA tests.
#' @param adj.vars A character vector with the names of columns in the 'meta.dat' component of 'data.obj' that are used as covariates in the PERMANOVA tests. Should be NULL if no covariates are used.
#' @param dist.name A character vector indicating the types of beta diversity indices to calculate and test. Possible values are: 'BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'.
#' @examples
#'
#' library(vegan)
#' library(GUniFrac)
#' library(ape)
#' library(philentropy)
#' library(MicrobiomeStat)
#' data(peerj32.obj)
#'
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#'
#' # Perform beta diversity tests using PERMANOVA
#' beta_diversity_test <- generate_beta_test_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
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

  # Calculate beta diversity indices if not provided
  if (is.null(dist.obj)) {
    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
    dist.obj <- mStat_calculate_beta_diversity(data.obj, dist.name)
  } else {
    message("Using provided dist.obj...")
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
        `Sum_Sq` = SumsOfSqs,
        `Mean_Sq` = MeanSqs,
        `F_Statistic` = F.Model,
        `R_Squared` = R2,
        `P_Value` = `Pr(>F)`
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
  aov.tab <- bind_rows(permanova.results$aov.tab) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.), "NA", round(., 3))))  # round all numeric values to 3 decimal places, replace NA with "NA"

  return(list("p.tab" = p.tab, "aov.tab" = aov.tab))
}
