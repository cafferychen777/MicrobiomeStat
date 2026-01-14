#' @title Alpha Diversity Change Test Per Time Point (Longitudinal)
#'
#' @description Performs paired tests comparing alpha diversity changes between
#'   baseline and each follow-up time point in longitudinal data.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. Options include 'log fold change', 'absolute change',
#'   or a custom function.
#'
#' @return A list of tables, one for each alpha diversity metric, summarizing the results of the statistical tests.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#' # Perform the longitudinal alpha diversity test for the peerj32 dataset
#' alpha_test_results_peerj32 <- generate_alpha_change_per_time_test_long(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "time",
#'   t0.level = unique(peerj32.obj$meta.dat$time)[1],
#'   ts.levels = unique(peerj32.obj$meta.dat$time)[-1],
#'   subject.var = "subject",
#'   group.var = "sex",
#'   adj.vars = "group",
#'   alpha.change.func = "log fold change"
#' )
#' # Generate dot plots for the peerj32 dataset results
#' generate_alpha_per_time_dotplot_long(
#'   data.obj = peerj32.obj,
#'   test.list = alpha_test_results_peerj32,
#'   group.var = "sex",
#'   time.var = "time",
#'   t0.level = unique(peerj32.obj$meta.dat$time)[1],
#'   ts.levels = unique(peerj32.obj$meta.dat$time)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#'
#' # Perform the longitudinal alpha diversity test for the peerj32 dataset without adjustment variables
#' alpha_test_results_peerj32_no_adj <- generate_alpha_change_per_time_test_long(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "time",
#'   t0.level = unique(peerj32.obj$meta.dat$time)[1],
#'   ts.levels = unique(peerj32.obj$meta.dat$time)[-1],
#'   subject.var = "subject",
#'   group.var = "sex",
#'   adj.vars = NULL,
#'   alpha.change.func = "log fold change"
#' )
#' # Generate dot plots for the peerj32 dataset results without adjustment variables
#' generate_alpha_per_time_dotplot_long(
#'   data.obj = peerj32.obj,
#'   test.list = alpha_test_results_peerj32_no_adj,
#'   group.var = "sex",
#'   time.var = "time",
#'   t0.level = unique(peerj32.obj$meta.dat$time)[1],
#'   ts.levels = unique(peerj32.obj$meta.dat$time)[-1],
#'   base.size = 16,
#'   theme.choice = "minimal"
#' )
#'
#' data("subset_T2D.obj")
#' # Perform the longitudinal alpha diversity test for the T2D dataset
#' alpha_test_results_T2D <- generate_alpha_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = c("sample_body_site"),
#'   alpha.change.func = "log fold change"
#' )
#' # Generate dot plots for the T2D dataset results
#' dot_plots_T2D <- generate_alpha_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = alpha_test_results_T2D,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#'
#' # Perform the longitudinal alpha diversity test for the T2D dataset without adjustment variables
#' alpha_test_results_T2D_no_adj <- generate_alpha_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   alpha.change.func = "log fold change"
#' )
#' # Generate dot plots for the T2D dataset results without adjustment variables
#' dot_plots_T2D_no_adj <- generate_alpha_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = alpha_test_results_T2D_no_adj,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   base.size = 16,
#'   theme.choice = "bw"
#' )
#' }
#' @export
generate_alpha_change_per_time_test_long <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = NULL,
           depth = NULL,
           time.var,
           t0.level,
           ts.levels,
           subject.var,
           group.var,
           adj.vars = NULL,
           alpha.change.func = "log fold change") {

    if (is.null(alpha.name)){
      return()
    }

    if (is.null(alpha.obj)) {
      if (!is.null(depth)) {
        message(
          "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj, depth = depth)
      }
      otu_tab <- data.obj$feature.tab
      
      # Extract tree if faith_pd is requested
      tree <- NULL
      if ("faith_pd" %in% alpha.name) {
        tree <- data.obj$tree
      }
      
      alpha.obj <- mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name, tree = tree)
    } else {
      # Verify that all alpha.name are present in alpha.obj
      if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x)
        colnames(x))))) {
        missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
        stop(
          "The following alpha diversity indices are not available in alpha.obj: ",
          paste(missing_alphas, collapse = ", "),
          call. = FALSE
        )
      }
    }

    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
        subject.var, group.var, time.var, adj.vars
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

    # Generate tests
    test.list <- lapply(ts.levels, function(ts.level) {

      subset.ids <- get_sample_ids(data.obj, time.var, c(t0.level, ts.level))

      subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)

      subset_alpha.obj <- mStat_subset_alpha(alpha.obj, samIDs = subset.ids)

      subset.test.list <- generate_alpha_change_test_pair(
        data.obj = subset_data.obj,
        alpha.obj = subset_alpha.obj,
        alpha.name = alpha.name,
        depth = depth,
        time.var = time.var,
        subject.var = subject.var,
        group.var = group.var,
        adj.vars = adj.vars,
        change.base = t0.level,
        alpha.change.func = alpha.change.func
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

    # Assign names to the elements of test.list
    names(test.list) <- ts.levels

    return(test.list)
  }
