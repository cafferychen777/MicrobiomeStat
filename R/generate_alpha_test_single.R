#' Perform alpha diversity association tests using linear models
#'
#' This function conducts association tests for multiple alpha diversity indices using linear model analyses.
#' It can be applied to cross-sectional data, a single time point from longitudinal or paired data.
#' It takes a data object as input and returns a list of association tests for each specified alpha diversity index.
#'
#' @param data.obj A list object containing feature.tab (OTU table with taxa in rows and samples in columns)
#'                 and meta.dat (Metadata table with samples in rows and variables in columns).
#' @param alpha.obj A list object containing pre-calculated alpha diversity indices;
#'                  if not provided, the function will calculate the indices using the 'alpha.name' parameter. Default is NULL.
#' @param time.var character; the name of the time variable in the metadata. Default is NULL.
#' @param t.level character; the level of the time variable to subset the data. Default is NULL.
#' @param alpha.name character vector containing the names of alpha diversity indices to calculate.
#'                   Possible values are: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param group.var character; the name of the group variable in the metadata.
#' @param adj.vars character vector; the names of the adjustment variables in the metadata.
#'
#' @return A list containing the association tests for each alpha diversity index.
#'         Each element in the list corresponds to a different alpha diversity index,
#'         and contains a dataframe with the linear model's coefficients, standard errors, t values, and p values.
#'
#' @examples
#' data("subset_T2D.obj")
#' generate_alpha_test_single(data.obj = subset_T2D.obj,
#'                            time.var = "visit_number",
#'                            t.level = "   4",
#'                            alpha.name = c("shannon", "simpson"),
#'                            group.var = "subject_race",
#'                            adj.vars = "subject_gender")
#'
#' @export
generate_alpha_test_single <-
  function(data.obj,
           alpha.obj = NULL,
           time.var = NULL,
           t.level = NULL,
           alpha.name,
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
        data.obj <- mStat_rarefy_data(data.obj)
      }
      otu_tab <- load_data_obj_count(data.obj)
      alpha.obj <- mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
    }

    # Generate tests
    test.list <- lapply(seq_along(alpha.obj), function(i) {
      df <- alpha.obj[[i]]
      # Join the alpha diversity index with metadata
      merged_df <-
        left_join(
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
          term = term,
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
          mutate(Estimate = NA, Std.Error = NA)

        # Reorder the columns to match coef.tab
        anova.tab <- anova.tab %>%
          select(
            term = term,
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
