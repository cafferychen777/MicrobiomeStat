#' @title Alpha Diversity Association Test (Single Time Point)
#'
#' @description Performs association tests for alpha diversity indices using linear
#'   models for cross-sectional data or a single time point.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @param t.level Character string specifying the time level/value to subset data to,
#'   if a time variable is provided. Default NULL does not subset data.
#'
#' @return A list containing the association tests for each alpha diversity index.
#'   Each element in the list corresponds to a different alpha diversity index,
#'   and contains a dataframe with the linear model's coefficients, standard errors, t values, and p values.
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
#' generate_alpha_test_single(data.obj = peerj32.obj,
#'                            time.var = "time",
#'                            t.level = "1",
#'                            alpha.name = c("shannon", "observed_species"),
#'                            group.var = "group",
#'                            adj.vars = NULL)
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
           adj.vars = NULL) {

    if (is.null(alpha.name)){
      return()
    }

    prepared <- mStat_prepare_alpha_inputs(
      data.obj = data.obj,
      alpha.obj = alpha.obj,
      alpha.name = alpha.name,
      depth = depth,
      time.var = time.var,
      t.level = t.level
    )
    data.obj <- prepared$data.obj
    alpha.obj <- prepared$alpha.obj

    # Generate tests
    test.list <- lapply(alpha.name, function(alpha.name) {

      df <- alpha.obj[[alpha.name]]
      # Join the alpha diversity index with metadata
      merged_df <-
        dplyr::inner_join(
          df %>% tibble::rownames_to_column("sample"),
          data.obj$meta.dat %>% tibble::rownames_to_column("sample"),
          by = "sample"
        )

      # Create a formula for lm
      formula_vars <- group.var
      if (!is.null(adj.vars)) {
        formula_vars <- c(adj.vars, group.var)
      }
      formula <-
        as.formula(paste0(names(merged_df)[2], "~", paste(formula_vars, collapse = "+")))

      # Run lm and create a coefficient table
      lm.model <- lm(formula, data = merged_df)

      summary <- summary(lm.model)

      coef.tab <- summary$coefficients %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Term") %>%
        tibble::as_tibble()

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
          tibble::rownames_to_column("Term") %>%
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
          tibble::as_tibble()

        coef.tab <-
          rbind(coef.tab, anova.tab)
      }

      return(coef.tab)
    })

    # Assign names to the elements of test.list
    names(test.list) <- alpha.name

    return(test.list)
  }
