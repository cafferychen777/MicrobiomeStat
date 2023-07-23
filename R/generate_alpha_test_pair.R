#' Check if the data has been rarefied
#'
#' This function checks if the data has been rarefied by inspecting if the sum
#' of each column (which represents each sample in the OTU table) is equal.
#'
#' @param data.obj A list object that includes feature.tab (an OTU table with
#'   taxa as rows and samples as columns) and meta.dat (a metadata table with
#'   samples as rows and variables as columns).
#'
#' @return A boolean value indicating whether the data is rarefied. It returns
#'   TRUE if the data is rarefied, FALSE otherwise.
#'
#' @examples
#' # Assuming peerj32.obj is a data object with OTU and metadata tables
#' is_rarefied(peerj32.obj)
#'
#' @export
is_rarefied <- function(data.obj) {
  unique_colsums <- unique(round(colSums(load_data_obj_count(data.obj)), 5))
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
#' @param data.obj A list object that includes feature.tab (an OTU table with
#' taxa as rows and samples as columns) and meta.dat (a metadata table with
#' samples as rows and variables as columns).
#' @param alpha.obj A pre-calculated alpha diversity object. If NULL, the
#' function will calculate alpha diversity based on data.obj.
#' @param time.var A string representing the time variable's name in the
#' metadata. The default is NULL.
#' @param alpha.name A character vector with the names of alpha diversity
#' indices to compute. Options include: "shannon", "simpson",
#' "observed_species", "chao1", "ace", and "pielou".
#' @param group.var A string representing the group variable's name in the
#' metadata.
#' @param subject.var A string specifying the subject variable column in the metadata.
#' @param adj.vars A character vector with the names of adjustment variables in
#' the metadata.
#' @return A list containing the association tests for each alpha diversity
#' index.
#'
#' @examples
#' alpha_test_results <- generate_alpha_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' time.var = "time",
#' alpha.name = c("shannon", "simpson"),
#' subject.var = "subject",
#' group.var = "group",
#' adj.vars = "sex"
#' )
#'
#' @export
generate_alpha_test_pair <-
  function(data.obj,
           alpha.obj = NULL,
           time.var,
           alpha.name,
           subject.var,
           group.var,
           adj.vars) {
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
