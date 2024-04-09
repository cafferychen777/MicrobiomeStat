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
#' # Example where alpha diversity indices are calculated beforehand
#' alpha.obj <- mStat_calculate_alpha_diversity(subset_T2D.obj$feature.tab,
#'                                              c("shannon", "observed_species"))
#' generate_alpha_test_single(data.obj = subset_T2D.obj,
#'                            alpha.obj = alpha.obj,
#'                            alpha.name = c("shannon", "observed_species", "ace"),
#'                            time.var = "visit_number",
#'                            t.level = NULL,
#'                            group.var = "subject_race",
#'                            adj.vars = "subject_gender")
#'
#' # Example where alpha diversity indices are calculated within the function
#' generate_alpha_test_single(data.obj = subset_T2D.obj,
#'                            time.var = "visit_number",
#'                            t.level = "4",
#'                            alpha.name = c("shannon", "observed_species"),
#'                            group.var = "subject_race",
#'                            adj.vars = "subject_gender")
#'
#' data("peerj32.obj")
#' generate_alpha_test_single(data.obj = peerj32.obj,
#'                            time.var = "time",
#'                            t.level = "1",
#'                            alpha.name = c("shannon", "observed_species"),
#'                            group.var = "group",
#'                            adj.vars = "sex")
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

    if (is.null(alpha.name)){
      return()
    }

    if (!is.null(time.var) & !is.null(t.level)) {
      subset.ids <- rownames(data.obj$meta.dat %>%
                               filter(!!sym(time.var) %in% c(t.level)))

      subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)
    }

    if (is.null(alpha.obj)) {
      if (!is_rarefied(data.obj)) {
        message(
          "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj, depth = depth)
      }
      otu_tab <- data.obj$feature.tab
      alpha.obj <- mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
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

    # Generate tests
    test.list <- lapply(alpha.name, function(alpha.name) {

      df <- alpha.obj[[alpha.name]]
      # Join the alpha diversity index with metadata
      merged_df <-
        dplyr::inner_join(
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

      summary <- summary(lm.model)

      coef.tab <- summary$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("Term") %>%
        as_tibble()

      # Rearrange the table
      coef.tab <-
        coef.tab %>% dplyr::select(
          Term,
          Estimate,
          Std.Error = `Std. Error`,
          Statistic = `t value`,
          P.Value = `Pr(>|t|)`
        )

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(na.omit(unique(merged_df[[group.var]]))) > 2) {
        anova <- anova(lm.model)
        anova.tab <- anova %>%
          as.data.frame() %>%
          rownames_to_column("Term") %>%
          dplyr::select(
            Term,
            Statistic = `F value`,
            P.Value = `Pr(>F)`
          ) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA)  %>%
          dplyr::select(
            Term,
            Estimate,
            Std.Error,
            Statistic,
            P.Value
          ) %>%
          dplyr::filter(
            Term == group.var
          ) %>%
          as_tibble()

        coef.tab <-
          rbind(coef.tab, anova.tab)
      }

      return(coef.tab)
    })

    # Assign names to the elements of test.list
    names(test.list) <- alpha.name

    return(test.list)
  }
