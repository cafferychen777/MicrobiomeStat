#' @title Alpha Diversity Change Test (Paired)
#'
#' @description Performs paired tests comparing alpha diversity changes between
#'   two time points, using linear models with optional group comparisons and
#'   covariate adjustment.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @param change.base The baseline time point for calculating changes. If NULL,
#'   the first unique time point in the data will be used.
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. Options include 'log fold change', 'absolute change',
#'   or a custom function taking two arguments (t, t0).
#'
#' @return A list of tables, one for each alpha diversity metric, summarizing the results of the statistical tests.
#' Each table contains the following columns: Term (the name of the variable in the model), Estimate (the estimated coefficient),
#' Std.Error (the standard error of the coefficient), Statistic (the t or F statistic), P.Value (the p-value of the test).
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Example 1: Using both group.var and adj.vars
#' generate_alpha_change_test_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "sex",
#'   adj.vars = NULL,
#'   change.base = "2",
#'   alpha.change.func = "log fold change"
#' )
#'
#' # Rename the time variable in peerj32.obj's metadata
#' peerj32.obj$meta.dat <- peerj32.obj$meta.dat %>%
#'   dplyr::rename(Day = time)
#'
#' # Example 2: Using group.var and adj.vars with a renamed time variable
#' generate_alpha_change_test_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon"),
#'   subject.var = "subject",
#'   time.var = "Day",
#'   group.var = "sex",
#'   adj.vars = c("group"),
#'   change.base = "2",
#'   alpha.change.func = "log fold change"
#' )
#'
#' data("subset_pairs.obj")
#'
#' # Example 3: With group.var and without adj.vars
#' generate_alpha_change_test_pair(
#'   data.obj = subset_pairs.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon"),
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   alpha.change.func = "log fold change"
#' )
#'
#' }
#' @export
generate_alpha_change_test_pair <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = NULL,
           depth = NULL,
           subject.var,
           time.var,
           group.var,
           adj.vars = NULL,
           change.base,
           alpha.change.func = "log fold change") {

    # Check if alpha diversity indices are specified. If not, exit the function.
    if (is.null(alpha.name)){
      return()
    }

    prepared <- mStat_prepare_alpha_inputs(
      data.obj = data.obj,
      alpha.obj = alpha.obj,
      alpha.name = alpha.name,
      depth = depth
    )
    data.obj <- prepared$data.obj
    alpha.obj <- prepared$alpha.obj

    meta_tab <- mStat_prepare_alpha_meta_tab(
      data.obj = data.obj,
      vars = c(subject.var, group.var, time.var, adj.vars)
    )

    # Initialize variable to store information about time-varying covariates.
    time_varying_info <- list(
      time_varying_vars = character(),
      non_time_varying_vars = character()
    )

    # Adjust alpha diversity for non-time-varying covariates if specified.
    # This is an important step to control for potential confounding factors.
    if (!is.null(adj.vars)){
      # Identify which variables are time-varying and which are not.
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # Adjust alpha diversity for non-time-varying covariates.
      # This helps to isolate the effect of the variables of interest.
      if (length(time_varying_info$non_time_varying_vars) > 0){
        alpha.obj <- mStat_calculate_adjusted_alpha_diversity(alpha.obj, meta.dat = meta_tab, time_varying_info$non_time_varying_vars)
      }
    }

    mStat_validate_group_var_contract(
      meta.dat = meta_tab,
      group.var = group.var,
      subject.var = subject.var,
      context = "alpha change testing"
    )

    alpha_df <- mStat_prepare_alpha_data(
      alpha.obj = alpha.obj,
      meta.dat = meta_tab,
      sample_col = "sample",
      join = "inner"
    )

    pair_change <- mStat_prepare_alpha_pair_test_data(
      alpha.df = alpha_df,
      alpha.name = alpha.name,
      subject.var = subject.var,
      time.var = time.var,
      group.var = group.var,
      time_varying_vars = time_varying_info$time_varying_vars,
      change.base = change.base,
      change.func = alpha.change.func,
      context = "alpha change testing"
    )
    combined_alpha <- pair_change$combined_alpha
    change.base <- pair_change$change.base
    change.after <- pair_change$change.after

    if (length(time_varying_info$time_varying_vars) > 0) {
      names_map <- setNames(paste0(time_varying_info$time_varying_vars, "_time_2"), time_varying_info$time_varying_vars)
      drop_names <- paste0(time_varying_info$time_varying_vars, "_time_1")
      combined_alpha <- combined_alpha %>%
        dplyr::select(-any_of(drop_names)) %>%
        rename(!!!names_map)
    }

    test.list <- lapply(alpha.name, function(index) {

      # Create a formula for the linear model.
      # This formula includes the change in alpha diversity as the response variable,
      # and the group variable and time-varying covariates as predictors.
      formula <-
        as.formula(paste0(paste0(index, "_diff"), "~", paste(c(
          time_varying_info$time_varying_vars, group.var
        ), collapse = "+")))

      # Fit the linear model.
      # This model tests for differences in alpha diversity change between groups,
      # while adjusting for time-varying covariates.
      lm.model <- lm(formula, data = combined_alpha)

      # Extract the model summary.
      summary <- summary(lm.model)

      # Create a coefficient table from the model summary.
      # This table includes estimates, standard errors, test statistics, and p-values for each term in the model.
      coef.tab <- summary$coefficients %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Term") %>%
        dplyr::select(
          Term,
          Std.Error = `Std. Error`,
          Statistic = `t value`,
          P.Value = `Pr(>|t|)`,
          Estimate
        ) %>% tibble::as_tibble()

      # Perform ANOVA if the group variable has more than two levels.
      # This tests for overall differences among groups, rather than pairwise comparisons.
      if (length(unique(stats::na.omit(combined_alpha[[group.var]]))) > 2) {
        anova <- anova(lm.model)
        anova.tab <- as.data.frame(anova) %>%
          tibble::rownames_to_column("Term") %>%
          dplyr::select(Term,
                        Statistic = `F value`,
                        P.Value = `Pr(>F)`) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          tibble::as_tibble()

        # Combine the coefficient table and ANOVA table.
        coef.tab <-
          rbind(coef.tab, anova.tab)
      }

      return(coef.tab)
    })

    # Name the elements of the test list according to the alpha diversity indices.
    names(test.list) <- alpha.name

    return(test.list)
  }
