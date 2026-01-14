#' Check if Data has been Rarefied
#'
#' Checks if the data has been rarefied by inspecting if all samples have
#' equal total counts.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @return Logical. TRUE if data is rarefied (all samples have equal counts),
#'   FALSE otherwise.
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

#' @title Alpha Diversity Association Test (Paired)
#'
#' @description Tests associations between alpha diversity and time/group variables
#'   using mixed-effects models for paired/longitudinal designs.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @param change.base A value indicating the base level for the time variable.
#'   If provided, the specified level will be used as the reference category in
#'   the model. Default is NULL (first level used).
#'
#' @return A list containing the association tests for each alpha diversity index.
#'
#' @examples
#' data(peerj32.obj)
#' generate_alpha_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' alpha.name = c("shannon", "simpson", "ace"),
#' subject.var = "subject",
#' time.var = "time",
#' group.var = NULL
#' )
#'
#' generate_alpha_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' alpha.name = c("shannon", "simpson", "ace"),
#' subject.var = "subject",
#' time.var = "time",
#' group.var = NULL,
#' change.base = "2"
#' )
#'
#' generate_alpha_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' alpha.name = c("shannon", "simpson", "ace"),
#' subject.var = "subject",
#' time.var = "time",
#' group.var = "group"
#' )
#'
#' generate_alpha_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' alpha.name = c("shannon", "simpson", "ace"),
#' subject.var = "subject",
#' time.var = "time",
#' group.var = "group",
#' adj.vars = "sex"
#' )
#'
#' data("subset_pairs.obj")
#' generate_alpha_test_pair(
#' data.obj = subset_pairs.obj,
#' alpha.obj = NULL,
#' alpha.name = c("shannon", "simpson", "ace"),
#' subject.var = "MouseID",
#' time.var = "Antibiotic",
#' group.var = "Sex"
#' )
#' generate_alpha_test_pair(
#' data.obj = subset_pairs.obj,
#' alpha.obj = NULL,
#' alpha.name = c("shannon", "simpson", "ace"),
#' subject.var = "MouseID",
#' time.var = "Antibiotic",
#' group.var = "Sex",
#' change.base = "Week 2"
#' )
#' @export
generate_alpha_test_pair <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = NULL,
           depth = NULL,
           subject.var,
           time.var,
           group.var,
           adj.vars = NULL,
           change.base = NULL) {

    # Exit the function if no alpha diversity indices are specified
    if (is.null(alpha.name)){
      return()
    }

    # Calculate alpha diversity if not provided
    # This ensures we have the necessary diversity metrics for the analysis
    if (is.null(alpha.obj)) {
      # Perform rarefaction if depth is specified
      # Rarefaction standardizes sampling effort across all samples
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
      
      alpha.obj <-
        mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name, tree = tree)
    } else {
      # Verify that all requested alpha diversity indices are available
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

    # Extract relevant metadata for the analysis
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
        subject.var, group.var, time.var, adj.vars
      )))


    # Change the base level for time.var if change.base is specified
    # This allows for flexibility in choosing the reference time point
    if (!is.null(change.base) && !is.null(time.var)) {
      if (change.base %in% meta_tab[[time.var]]) {
        meta_tab[[time.var]] <- relevel(as.factor(meta_tab[[time.var]]), ref = change.base)
      } else {
        stop("Specified change.base is not a level in the time.var column.", call. = FALSE)
      }
    }

    # Combine alpha diversity data with metadata
    # This creates a comprehensive dataset for our analysis
    alpha_df <-
      dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"),
                        by = c("sample"))

    # Determine the number of levels in the group variable if it exists
    # This information is used later to decide between t-test and ANOVA
    if (!is.null(group.var)){
      group.levels <- alpha_df %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels() %>% length()
    }

    # Perform statistical tests for each alpha diversity index
    test.list <- lapply(alpha.name, function(index) {

      # Function to try fitting a complex mixed-effects model
      # If the complex model fails, it falls back to a simpler model
      try_complex_model <- function(alpha_df, formula_str) {
        tryCatch({
          # Attempt to fit a complex mixed-effects model
          # This model accounts for repeated measures and potential interactions
          lme.model <- lmerTest::lmer(formula_str, data = alpha_df)
          return(lme.model)
        },
        error = function(e) {
          # If the complex model fails, attempt a simpler model
          message("Complex model failed. Trying a simpler model...")

          # Function to construct a simpler formula
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

          # Fit the simpler mixed-effects model
          lme.model_simple <-
            lmerTest::lmer(new_formula_str, data = alpha_df)
          return(lme.model_simple)
        })
      }

      # Construct the formula for the mixed-effects model
      formula_str <-
        construct_formula(index, group.var, time.var, subject.var, adj.vars)

      # Attempt to fit the model, falling back to a simpler model if necessary
      lme.model <- try_complex_model(alpha_df, formula_str)

      # Extract coefficients from the fitted model
      coef.tab <- extract_coef(lme.model)

    # Perform additional analysis if a grouping variable is present
    if(!is.null(group.var)){
      # Run ANOVA if the grouping variable has more than two levels
      # This tests for overall differences among groups, rather than pairwise comparisons
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

        # Combine the coefficient table with the ANOVA results
        coef.tab <-
          rbind(coef.tab, anova.tab)
      }
    }

      return(as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list based on the alpha diversity indices
    names(test.list) <- alpha.name

    return(test.list)
  }