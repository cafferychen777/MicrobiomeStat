#' Longitudinal Alpha Diversity Test in Microbiome Data
#'
#' This function performs alpha diversity testing at multiple time points in longitudinal microbiome data.
#' It applies linear model analyses to assess the diversity of microbial communities over time within different groups or under various conditions.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name character vector containing the names of alpha diversity indices to calculate.
#'                Possible values are: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param time.var A string representing the time variable in the meta.dat.
#' @param t0.level A baseline time point for longitudinal analysis, specified as a character
#'                 or numeric value, e.g., "week_0" or 0.
#' @param ts.levels A vector of character strings indicating follow-up time points, such as
#'                  c("week_4", "week_8").
#' @param group.var Optional; a string specifying the group variable in meta.dat for between-group comparisons.
#' @param adj.vars Optional; a vector of strings representing covariates in meta.dat for adjustment in the analysis.
#' @return A list where each element corresponds to a different time point and contains the results of alpha diversity tests for that time point.
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset (without providing alpha.obj)
#' data("ecam.obj")
#'
#' # Perform the longitudinal alpha diversity test for the ECAM dataset
#' alpha_test_results_ecam <- generate_alpha_test_long(
#'   data.obj = ecam.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "pielou"),
#'   depth = 1000,
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "delivery",
#'   adj.vars = c("diet")
#' )
#'
#' # Generate dot plots for the ECAM dataset results
#' dot_plots_ecam <- generate_alpha_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = alpha_test_results_ecam,
#'   group.var = "delivery",
#'   time.var = "month_num",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#'
#' # Example 2: Analyzing the Type 2 Diabetes (T2D) dataset (with providing alpha.obj)
#' data("subset_T2D.obj")
#'
#' # Calculate alpha diversity indices
#' alpha.obj <- mStat_calculate_alpha_diversity(
#'   x = subset_T2D.obj$feature.tab,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou")
#' )
#'
#' # Perform the longitudinal alpha diversity test for the T2D dataset
#' alpha_test_results_T2D <- generate_alpha_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = alpha.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   adj.vars = c("sample_body_site")
#' )
#'
#' # Generate dot plots for the T2D dataset results
#' dot_plots_T2D <- generate_alpha_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = alpha_test_results_T2D,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#' }
#' @export
generate_alpha_test_long <- function(data.obj,
                                     alpha.obj = NULL,
                                     alpha.name = NULL,
                                     depth = NULL,
                                     time.var,
                                     t0.level,
                                     ts.levels,
                                     group.var = NULL,
                                     adj.vars = NULL) {
  # Validate input data
  mStat_validate_data(data.obj)

  if (is.null(alpha.name)){
    return()
  }

  if (is.null(alpha.obj)) {
    if (!is_rarefied(data.obj)) {
      message(
        "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj, depth = depth)
    }
    otu_tab <- data.obj$feature.tab
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
  }

  meta_tab <-
    data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
      group.var, time.var, adj.vars
    )))

  reference_level <- levels(as.factor(meta_tab[,group.var]))[1]

  if (is.null(t0.level)) {
    if (is.numeric(meta_tab[, time.var])) {
      t0.level <- sort(unique(meta_tab[, time.var]))[1]
    } else {
      t0.level <- levels(meta_tab[, time.var])[1]
    }
  }

  if (is.null(ts.levels)) {
    if (is.numeric(meta_tab[, time.var])) {
      ts.levels <- sort(unique(meta_tab[, time.var]))[-1]
    } else {
      ts.levels <- levels(meta_tab[, time.var])[-1]
    }
  }

  # Get unique time levels
  time.levels <- c(t0.level, ts.levels)

  test.list <- lapply(time.levels, function(t.level){
    # Subset the data for the specific time level
    condition <- paste(time.var, "== '", t.level, "'", sep = "")
    subset_data.obj <- mStat_subset_data(data.obj, condition = condition)

    # Subset the alpha.obj to match the subsetted data
    subset_alpha.obj <- mStat_subset_alpha(alpha.obj, rownames(subset_data.obj$meta.dat))

    # Perform alpha diversity test for the subset data
    subset.test.list <- generate_alpha_test_single(
      data.obj = subset_data.obj,
      alpha.obj = subset_alpha.obj,
      alpha.name = alpha.name,
      time.var = time.var,
      t.level = t.level,
      group.var = group.var,
      adj.vars = adj.vars
    )

    all_terms <- unique(unlist(lapply(subset.test.list, \(df) df$Term)))
    all_terms <- setdiff(all_terms, "(Intercept)")
    all_terms <- all_terms[grepl(paste0("^", group.var), all_terms)]

    new_list <- lapply(all_terms, \(term) {
      alpha_dfs <- lapply(names(subset.test.list), \(alpha_name) {
        df <- subset.test.list[[alpha_name]]
        filtered_df <- filter(df, Term == term) %>%
          dplyr::select(-Term) %>%
          dplyr::mutate(Term = alpha_name) %>%
          dplyr::select(Term, everything())

        if(nrow(filtered_df) == 0) {
          return(NULL)
        }

        return(filtered_df)
      })
      do.call(rbind, alpha_dfs)
    })

    all_terms <- lapply(all_terms, function(term) {
      if (term != group.var) {
        modified_term <- stringr::str_replace(term, paste0("^", group.var), "")
        modified_term <- stringr::str_trim(modified_term)
        paste(sprintf("%s vs %s", modified_term, reference_level), "(Reference)")
      } else {
        term
      }
    })

    new_subset_test_list <- setNames(new_list, all_terms)

    return(new_subset_test_list)
  })

  names(test.list) <- time.levels

  return(test.list)
}
