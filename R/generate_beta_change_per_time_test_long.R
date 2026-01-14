#' Beta Diversity Change Tests Per Time Point
#'
#' Performs beta diversity change tests at each time point relative to baseline.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
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

    # Check if distance names are provided
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    if (is.null(dist.obj)) {
      # Process time variable to ensure proper ordering
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      # Extract relevant metadata
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(time.var, group.var)))
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Adjust distance matrix if adjustment variables are provided
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # If distance object is provided, extract metadata
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

    # Determine the reference level for the grouping variable
    reference_level <- levels(as.factor(meta_tab[,group.var]))[1]

    # Set default baseline time point if not provided
    if (is.null(t0.level)) {
      if (is.numeric(meta_tab[, time.var])) {
        t0.level <- sort(unique(meta_tab[, time.var]))[1]
      } else {
        t0.level <- levels(meta_tab[, time.var])[1]
      }
    }

    # Set default follow-up time points if not provided
    if (is.null(ts.levels)) {
      if (is.numeric(meta_tab[, time.var])) {
        ts.levels <- sort(unique(meta_tab[, time.var]))[-1]
      } else {
        ts.levels <- levels(meta_tab[, time.var])[-1]
      }
    }

    # Combine baseline and follow-up time points
    time.levels <- c(t0.level, ts.levels)

    # Perform beta diversity change tests for each follow-up time point
    test.list <- lapply(ts.levels, function(ts.level){
      # Subset the data for the specific time level
      subset.ids <- get_sample_ids(data.obj, time.var, c(t0.level, ts.level))

      # Create subsets of data and distance objects
      subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)
      subset_dist.obj <- mStat_subset_dist(dist.obj, samIDs = subset.ids)

      # Perform beta diversity change test for the current time point pair
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

      # Extract all terms from the test results
      all_terms <- unique(unlist(lapply(subset.test.list, \(df) df$Term)))
      all_terms <- setdiff(all_terms, "(Intercept)")
      all_terms <- all_terms[grepl(paste0("^", group.var), all_terms)]

      # Reorganize test results by term
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

      # Modify term names for better readability
      all_terms <- lapply(all_terms, function(term) {
        if (term != group.var) {
          modified_term <- stringr::str_replace(term, paste0("^", group.var), "")
          modified_term <- stringr::str_trim(modified_term)
          paste(sprintf("%s vs %s", modified_term, reference_level), "(Reference)")
        } else {
          term
        }
      })

      # Set names for the reorganized test results
      new_subset_test_list <- setNames(new_list, all_terms)

      return(new_subset_test_list)
    })

    # Set names for the test results list using follow-up time points
    names(test.list) <- ts.levels

    return(test.list)
  }