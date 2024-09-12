#' Generate Beta Change Tests for Longitudinal Data
#'
#' This function is designed to perform beta change tests on longitudinal microbiome data.
#' It integrates multiple components of the MicrobiomeStat package to process and analyze
#' the data, focusing on beta diversity changes over time within subjects.
#'
#' @param data.obj A list object formatted specifically for MicrobiomeStat, potentially
#'                 including components like feature.tab (matrix), feature.ann (matrix),
#'                 meta.dat (data.frame), tree, and feature.agg.list (list). This object can
#'                 be derived from various formats using MicrobiomeStat functions such as
#'                 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj',
#'                 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj',
#'                 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj',
#'                 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'.
#' @param dist.obj A distance matrix between samples, typically generated using
#'                 \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}}. If NULL,
#'                 the function will compute beta diversity from \code{data.obj}.
#' @param time.var A character string identifying the column in metadata that contains
#'                 time information. This column should list time points for each sample,
#'                 necessary for ordering and linking samples over time for the same subject.
#' @param t0.level A baseline time point for longitudinal analysis, specified as a character
#'                 or numeric value, e.g., "week_0" or 0.
#' @param ts.levels A vector of character strings indicating follow-up time points, such as
#'                  c("week_4", "week_8").
#' @param subject.var The name of the column in metadata containing the subject IDs.
#'                    This should uniquely identify each subject in the study. Required
#'                    to identify samples that belong to the same subject.
#' @param group.var A character string naming the column in metadata that contains the
#'                  grouping variable, which can be used for color mapping in plots.
#'                  Optional; can be NULL.
#' @param adj.vars A vector of character strings naming columns in metadata to include as
#'                 covariates for adjusting the distance matrix before ordination.
#'                 Can be empty or NULL if no adjustment is needed.
#' @param dist.name A vector of character strings specifying which beta diversity indices to
#'                  calculate. Supported indices include "BC" (Bray-Curtis), "Jaccard",
#'                  "UniFrac" (unweighted), "GUniFrac" (generalized UniFrac), "WUniFrac"
#'                  (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a specified
#'                  index is not present in dist.obj, it will be computed internally.
#'                  An error is returned if an unsupported index is specified.
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return A list of test results, structured to facilitate further analysis and visualization.
#'
#' @seealso Related functions in MicrobiomeStat for data preparation and beta diversity calculation,
#'          such as \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}},
#'          \code{\link[MicrobiomeStat]{mStat_calculate_PC}}, and data conversion functions like
#'          \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}.
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' result1 <- generate_beta_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_beta_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result1,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#' )
#'
#' generate_beta_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = c("sample_body_site", "subject_gender"),
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, c('BC', 'Jaccard'))
#' result2 <- generate_beta_change_per_time_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   subject.var = "subject.id",
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_beta_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = result2,
#'   group.var = "delivery",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 15
#' )
#' }
#' @export
generate_beta_change_per_time_test_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           subject.var,
           group.var = NULL,
           adj.vars = NULL,
           dist.name = c('BC', 'Jaccard'),
           ...) {

    if (is.null(dist.name)){
      return()
    }

    if (is.null(dist.obj)) {
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(time.var, group.var)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(time.var, group.var)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(time.var, group.var)))
        data.obj <- list(meta.dat = meta_tab)
        data.obj <- mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat
      }
    }

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

    time.levels <- c(t0.level, ts.levels)

    test.list <- lapply(ts.levels, function(ts.level){
      # Subset the data for the specific time level
      subset.ids <- rownames(data.obj$meta.dat %>%
                               filter(!!sym(time.var) %in% c(t0.level, ts.level)))

      subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)

      subset_dist.obj <- mStat_subset_dist(dist.obj, samIDs = subset.ids)

      subset.test.list <- generate_beta_change_test_pair(
        data.obj = subset_data.obj,
        dist.obj = subset_dist.obj,
        time.var = time.var,
        subject.var = subject.var,
        group.var = group.var,
        adj.vars = adj.vars,
        change.base = t0.level,
        dist.name = dist.name
      )

      all_terms <- unique(unlist(lapply(subset.test.list, \(df) df$Term)))
      all_terms <- setdiff(all_terms, "(Intercept)")
      all_terms <- all_terms[grepl(paste0("^", group.var), all_terms)]

      new_list <- lapply(all_terms, \(term) {
        alpha_dfs <- lapply(names(subset.test.list), \(alpha_name) {
          df <- subset.test.list[[alpha_name]]
          filtered_df <- filter(df, Term == term) %>%
            select(-Term) %>%
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

    names(test.list) <- ts.levels

    return(test.list)
  }
