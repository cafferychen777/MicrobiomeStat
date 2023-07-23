#' Longitudinal Alpha Diversity Association Test in MicrobiomeStat
#'
#' This function, part of the MicrobiomeStat package, implements a method for
#' performing an association test for multiple alpha diversity measures in
#' longitudinal data. The test is based on a mixed-effects model, which is more
#' robust and efficient than previous implementations. The function is especially
#' suitable for microbiome studies involving repeated measurements over time.
#'
#' The mixed-effects model includes a temporal component, accounting for the
#' correlation among repeated measurements on the same subject over time.
#' The function also accepts additional adjustment variables for greater
#' flexibility in modeling.
#'
#' The output is a list of coefficient tables, one for each alpha diversity index.
#' Each table includes the term, estimate, standard error, t value, and p-value
#' for each fixed effect in the model.
#'
#' @param data.obj A list object that includes feature.tab (an OTU table with
#' taxa as rows and samples as columns) and meta.dat (a metadata table with
#' samples as rows and variables as columns).
#' @param alpha.obj A pre-calculated alpha diversity object. If NULL, the
#' function will calculate alpha diversity based on data.obj.
#' @param time.var A string representing the time variable's name in the
#' metadata. The default is NULL.
#' @param t0.level The level in time.var considered as the reference level in
#' the longitudinal study. Typically the first time point.
#' @param ts.levels The levels in time.var that are to be considered in
#' the analysis, excluding the reference level.
#' @param alpha.name A character vector with the names of alpha diversity
#' indices to compute. Options include: "shannon", "simpson",
#' "observed_species", "chao1", "ace", and "pielou".
#' @param group.var A string representing the group variable's name in the
#' metadata.
#' @param adj.vars A character vector with the names of adjustment variables in
#' the metadata.
#' @return A list containing the association tests for each alpha diversity
#' index.
#'
#' @examples
#' library("HMP2Data")
#' T2D <- T2D16S()
#' T2D.obj <- mStat_convert_phyloseq_to_data_obj(T2D.phy)
#'
#' subset_T2D.obj <- mStat_subset_data(T2D.obj,colnames(T2D.obj$feature.tab
#' [,colSums(T2D.obj$feature.tab) >= 2000]))
#' alpha_test_results <- generate_alpha_test_long(
#' data.obj = subset_T2D.obj,
#' alpha.obj = NULL,
#' time.var = "visit_number",
#' t0.level = sort(unique(T2D.obj$meta.dat$visit_number))[1],
#' ts.levels = sort(unique(T2D.obj$meta.dat$visit_number))[2:6],
#' alpha.name = c("shannon", "simpson"),
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "subject_gender"
#' )
#'
#' @export
generate_alpha_test_long <-
  function(data.obj,
           alpha.obj = NULL,
           time.var,
           t0.level,
           ts.levels,
           alpha.name,
           subject.var,
           group.var,
           adj.vars) {
    if (is.null(alpha.obj)) {
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      if (!is_rarefied(data.obj)) {
        message(
          "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj)
      }
      otu_tab <- load_data_obj_count(data.obj)
      alpha.obj <-
        mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
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
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    }

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
        subject.var, group.var, time.var, adj.vars
      )))

    # Convert the alpha.obj list to a data frame
    alpha_df <-
      bind_cols(alpha.obj) %>% bind_cols(tibble("sample" = colnames(otu_tab))) %>%
      inner_join(meta_tab %>% rownames_to_column("sample"),
                 by = c("sample"))

    test.list <- lapply(alpha.name, function(index) {
      formula_str <- paste0(index, "~", time.var)
      if (!is.null(adj.vars)) {
        formula_str <-
          paste0(formula_str, "+", paste(adj.vars, collapse = "+"))
      }
      formula_str <-
        paste0(formula_str, "+", paste(group.var, collapse = "+"))
      formula_str <- paste0(formula_str, " + (1|", subject.var, ")")
      formula <- as.formula(formula_str)

      lme.model <- lmerTest::lmer(formula, data = alpha_df)
      coef.tab <- extract_coef(lme.model)

      return(as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list
    names(test.list) <- alpha.name

    return(test.list)
  }
