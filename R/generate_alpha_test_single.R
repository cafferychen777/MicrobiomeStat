#' Perform alpha diversity association tests using linear models
#'
#' This function conducts association tests for multiple alpha diversity indices using linear model analyses.
#' It can be applied to cross-sectional data, a single time point from longitudinal or paired data.
#' It takes a data object as input and returns a list of association tests for each specified alpha diversity index.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name character vector containing the names of alpha diversity indices to calculate.
#'                   Possible values are: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param time.var Character string specifying the column name in metadata containing time variable.
#'                Used to subset to a single time point if provided. Default NULL does not subset data.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var Character string specifying the column name in metadata containing the grouping
#'                 variable. Used as a predictor in the linear models.
#' @param adj.vars Character vector specifying column names in metadata containing variables to adjust
#'                for in the linear models.
#'
#' @return A list containing the association tests for each alpha diversity index.
#'         Each element in the list corresponds to a different alpha diversity index,
#'         and contains a dataframe with the linear model's coefficients, standard errors, t values, and p values.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' generate_alpha_test_single(data.obj = subset_T2D.obj,
#'                            time.var = "visit_number",
#'                            t.level = "   4",
#'                            alpha.name = c("shannon", "simpson"),
#'                            group.var = "subject_race",
#'                            adj.vars = "subject_gender")
#' }
#' @export
generate_alpha_test_single <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = NULL,
           depth = NULL,
           time.var = NULL,
           t.level = NULL,
           group.var,
           adj.vars) {
    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    if (is.null(alpha.obj)) {
      if (!is_rarefied(data.obj)) {
        message(
          "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj, depth = depth)
      }
      otu_tab <- load_data_obj_count(data.obj)
      alpha.obj <- mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
    }

    # Generate tests
    test.list <- lapply(seq_along(alpha.obj), function(i) {
      df <- alpha.obj[[i]]
      # Join the alpha diversity index with metadata
      merged_df <-
        dplyr::left_join(
          df %>% rownames_to_column("sample"),
          data.obj$meta.dat %>% rownames_to_column("sample"),
          by = "sample"
        )

      # Create a formula for lm
      formula <-
        as.formula(paste0(names(merged_df)[2], "~", paste(c(
          adj.vars, group.var
        ), collapse = "+")))

      # Run lm and create a coefficient table
      lm.model <- lm(formula, data = merged_df)
      coef.tab <- broom::tidy(summary(lm.model))

      # Rearrange the table
      coef.tab <-
        coef.tab %>% select(
          Term = term,
          Estimate = estimate,
          Std.Error = std.error,
          Statistic = statistic,
          P.Value = p.value
        )

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(unique(merged_df[[group.var]])) > 2) {
        anova.tab <- broom::tidy(anova(lm.model))

        # Rearrange the table and add missing columns
        anova.tab <- anova.tab %>%
          select(
            term = term,
            Statistic = statistic,
            df = df,
            P.Value = p.value
          ) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA)

        # Reorder the columns to match coef.tab
        anova.tab <- anova.tab %>%
          select(
            Term = term,
            Estimate = Estimate,
            Std.Error = Std.Error,
            Statistic = Statistic,
            P.Value = P.Value
          )

        coef.tab <-
          rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
      }

      return(coef.tab)
    })

    # Assign names to the elements of test.list
    names(test.list) <- alpha.name

    return(test.list)
  }
