#' Perform longitudinal beta diversity tests using PERMANOVA
#'
#' This function, part of the MicrobiomeStat package, calculates various beta diversity indices for a microbiome dataset. It then performs longitudinal PERMANOVA tests on these indices to assess the effect of a grouping variable and any additional covariates on the beta diversity. The function returns the PERMANOVA results for each of the specified beta diversity indices.
#'
#' @name generate_beta_test_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, containing microbiome data and metadata.
#' @param dist.obj A distance object, if NULL it will be calculated based on data.obj and dist.name. This distance object should be a dissimilarity matrix or an object that can be coerced to a dissimilarity matrix.
#' @param time.var The name of the column in data.obj$meta.dat that contains the time variable for subsetting data, if applicable. This should be a character string representing the column name.
#' @param t0.level A character or numeric value indicating the baseline level of the time variable for subsetting data.
#' @param ts.levels A character or numeric vector indicating the subsequent time points of interest for subsetting data.
#' @param subject.var The name of the column in data.obj$meta.dat that contains the subject variable for the PERMANOVA tests. This should be a character string representing the column name.
#' @param group.var The name of the column in data.obj$meta.dat that contains the grouping variable for the PERMANOVA tests. This should be a character string representing the column name.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param dist.name A character vector specifying which beta diversity indices to calculate and test (e.g., 'BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'). Default is c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS').
#' @examples
#' data(subset_T2D.obj)
#'
#' # Perform pairwise beta diversity tests using PERMANOVA
#' beta_test_long_results <- generate_beta_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   time.var = "visit_number",
#'   t0.level = sort(unique(subset_T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = sort(unique(subset_T2D.obj$meta.dat$visit_number))[2:6],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = "subject_gender",
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
generate_beta_test_long <- function(data.obj,
                                    dist.obj = NULL,
                                    time.var,
                                    t0.level,
                                    ts.levels,
                                    subject.var,
                                    group.var,
                                    adj.vars = NULL,
                                    dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

  if (is.null(dist.obj)) {
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    dist.obj <-
      mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
  } else {
    if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    } else {
      metadata <- attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      data.obj <- list(meta.dat = metadata)
      data.obj <- mStat_process_time_variable(metadata, time.var, t0.level, ts.levels)
    }
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
  formula_str <- paste0(formula_str, " + ", time.var)
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
