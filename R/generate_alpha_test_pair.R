#' Check if the data has been rarefied
#'
#' This function checks if the data has been rarefied by inspecting if the sum
#' of each column (which represents each sample in the OTU table) is equal.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#'
#' @return A boolean value indicating whether the data is rarefied. It returns
#'   TRUE if the data is rarefied, FALSE otherwise.
#'
#' @examples
#' # Assuming peerj32.obj is a data object with OTU and metadata tables
#' data(peerj32.obj)
#' is_rarefied(peerj32.obj)
#'
#' @export
is_rarefied <- function(data.obj) {
  unique_colsums <-
    unique(round(colSums(data.obj$feature.tab), 5))
  return(length(unique_colsums) == 1)
}


#' Extract coefficient table from a model
#'
#' This is a helper function for the `generate_alpha_test_pair` function. It
#' extracts a coefficient table from a mixed-effects model fitted by the `lmer`
#' function from the `lmerTest` package. The coefficient table includes the
#' term, estimate, standard error, t value (statistic), and p-value.
#'
#' @param model A mixed-effects model object.
#' @return A data frame containing the coefficient table.
#'
#'
#' @noRd
extract_coef <- function(model) {
  s <- summary(model)
  fixed_eff <- s$coefficients
  coef_tab <- data.frame(
    Term = rownames(fixed_eff),
    Estimate = fixed_eff[, "Estimate"],
    Std.Error = fixed_eff[, "Std. Error"],
    Statistic = fixed_eff[, "t value"],
    P.Value = fixed_eff[, "Pr(>|t|)"]
  )
  return(coef_tab)
}

#' Alpha Diversity Association Test
#'
#' This function implements an association test for multiple alpha diversity
#' measures. The test is based on a mixed-effects model fitted by the `lmer`
#' function from the `lmerTest` package. The function accepts a data object as
#' input and returns a list of tests, one for each alpha diversity index.
#'
#' The mixed-effects model includes the time variable, group variable, and any
#' additional adjustment variables as fixed effects, and the subject variable as
#' a random effect.
#'
#' The output is a list of coefficient tables, one for each alpha diversity
#' index. Each table includes the term, estimate, standard error, t value, and
#' p-value for each fixed effect in the model.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name A character vector with the names of alpha diversity
#' indices to compute. Options include: "shannon", "simpson",
#' "observed_species", "chao1", "ace", and "pielou".
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param group.var A string representing the group variable's name in the
#' metadata.
#' @param time.var A string representing the time variable's name in the
#' metadata. The default is NULL.
#' @param subject.var A string specifying the subject variable column in the metadata.
#' @param adj.vars A character vector with the names of adjustment variables in
#' the metadata.
#' @return A list containing the association tests for each alpha diversity
#' index.
#'
#' @examples
#'
#' data(peerj32.obj)
#' alpha_test_results <- generate_alpha_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' time.var = "time",
#' alpha.name = c("shannon"),
#' subject.var = "subject",
#' group.var = "group",
#' adj.vars = "sex"
#' )
#'
#' @export
generate_alpha_test_pair <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = NULL,
           depth = NULL,
           time.var,
           subject.var,
           group.var,
           adj.vars) {

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
        subject.var, group.var, time.var, adj.vars
      )))

    # Convert the alpha.obj list to a data frame
    alpha_df <-
      dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"),
                        by = c("sample"))

    if (!is.null(group.var)){
      group.levels <- alpha_df %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels() %>% length()
    }

    test.list <- lapply(alpha.name, function(index) {
      try_complex_model <- function(alpha_df, formula_str) {
        tryCatch({
          # 尝试估计复杂模型
          lme.model <- lmerTest::lmer(formula_str, data = alpha_df)
          return(lme.model)
        },
        error = function(e) {
          # 如果发生错误，尝试简化的模型
          message("Complex model failed. Trying a simpler model...")

          correct_formula <-
            function(index,
                     group.var,
                     time.var,
                     subject.var,
                     adj.vars) {
              if (!is.null(group.var)) {
                formula_part <-
                  paste(index,
                        "~",
                        group.var,
                        "*",
                        time.var,
                        " + (1 ",
                        "|",
                        subject.var,
                        ")")
              } else {
                formula_part <-
                  paste(index, "~", time.var, " + (1|", subject.var, ")")
              }

              if (!is.null(adj.vars)) {
                formula_str <- paste(formula_part, "+", adj.vars)
              } else {
                formula_str <- formula_part
              }
              return(as.formula(formula_str))
            }

          new_formula_str <-
            correct_formula(index, group.var, time.var, subject.var, adj.vars)

          lme.model_simple <-
            lmerTest::lmer(new_formula_str, data = alpha_df)
          return(lme.model_simple)
        })
      }

      formula_str <-
        construct_formula(index, group.var, time.var, subject.var, adj.vars)

      lme.model <- try_complex_model(alpha_df, formula_str)

      coef.tab <- extract_coef(lme.model)

    if(!is.null(group.var)){
      # Run ANOVA on the model if group.var is multi-categorical
      if (group.levels > 2) {
        anova.tab <- anova(lme.model) %>%
          as.data.frame() %>%
          rownames_to_column("Term") %>%
          dplyr::select(Term,
                        Statistic = `F value`,
                        P.Value = `Pr(>F)`) %>%
          as_tibble() %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          dplyr::select(
            Term,
            Estimate,
            Std.Error,
            Statistic,
            P.Value
          ) %>%
          dplyr::filter(Term %in% c(group.var, paste0(group.var, ":", time.var)))

        coef.tab <-
          rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
      }
    }

      return(as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list
    names(test.list) <- alpha.name

    return(test.list)
  }
