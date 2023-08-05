#' Generate Alpha Diversity Change Test Pair
#'
#' This function generates a pair of statistical tests for comparing alpha diversity metrics
#' at two different time points for a given dataset. The function supports various options
#' for adjusting the tests and calculating the difference in alpha diversity.
#'
#' @param data.obj A data object containing both the count (OTU table) and the metadata.
#' @param alpha.obj A data object containing alpha diversity metrics. If NULL, it will be calculated from the data.obj.
#' @param time.var The name of the variable in the data that represents time points.
#' @param alpha.name The name of the alpha diversity metric to be tested.
#' @param subject.var The name of the variable in the data that represents the subject IDs.
#' @param group.var The name of the variable in the data that represents the group IDs.
#' @param adj.vars A vector of names of the variables in the data that should be used as covariates in the model.
#' @param change.base The name of the baseline time point for calculating changes in alpha diversity. If NULL, the first unique time point in the data will be used.
#' @param change.func The function to be used for calculating changes in alpha diversity. If "lfc", the change will be calculated as the log-fold change. If a function, it will be applied directly.
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
#' generate_alpha_change_test_pair(
#' data.obj = peerj32.obj,
#' alpha.obj = NULL,
#' time.var = "time",
#' alpha.name = c("shannon"),
#' subject.var = "subject",
#' group.var = "group",
#' adj.vars = "sex",
#' change.base = "1"
#' )
#' }
#' @export
generate_alpha_change_test_pair <-
  function(data.obj,
           alpha.obj = NULL,
           time.var,
           alpha.name,
           subject.var,
           group.var,
           adj.vars,
           change.base,
           change.func = "lfc") {

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
      dplyr::bind_cols(alpha.obj) %>% dplyr::bind_cols(tibble("sample" = colnames(otu_tab))) %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"),
                 by = c("sample"))

    if (is.null(change.base)){
      change.base <- unique(alpha_df %>% select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'alpha_df' data frame: ", change.base)
    }

    change.after <-
      unique(alpha_df %>% select(all_of(c(time.var))))[unique(alpha_df %>% select(all_of(c(time.var)))) != change.base]

    alpha_grouped <- alpha_df %>% dplyr::group_by(time)
    alpha_split <- split(alpha_df, f = alpha_grouped$time)

    alpha_time_1 <- alpha_split[[change.base]]
    alpha_time_2 <- alpha_split[[change.after]]

    combined_alpha <- alpha_time_1 %>%
      dplyr::inner_join(
        alpha_time_2,
        by = c(subject.var, group.var),
        suffix = c("_time_1", "_time_2")
      )

    diff_columns <- lapply(alpha.name, function(index) {

      diff_col_name <- paste0(index, "_diff")

      if (is.function(change.func)) {

        combined_alpha <- combined_alpha %>%
          dplyr::mutate(!!diff_col_name := change.func(!!sym(paste0(
            index, "_time_2"
          )), !!sym(paste0(
            index, "_time_1"
          )))) %>%
          select(all_of(diff_col_name))
      } else {

        if (change.func == "lfc") {
          combined_alpha <- combined_alpha %>%
            dplyr::mutate(!!diff_col_name := log(!!sym(paste0(
              index, "_time_2"
            )) / !!sym(paste0(
              index, "_time_1"
            )))) %>%
            select(all_of(diff_col_name))
        } else {
          combined_alpha <- combined_alpha %>%
            dplyr::mutate(!!diff_col_name := !!sym(paste0(index, "_time_2")) -!!sym(paste0(index, "_time_1"))) %>%
            select(all_of(diff_col_name))
        }
      }
    })

    combined_alpha <- dplyr::bind_cols(combined_alpha, diff_columns)

    if (!is.null(adj.vars)) {
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_time_1 %>% select(all_of(c(
          subject.var, adj.vars
        )))
        , by = c(subject.var))
    }

    # Generate tests
    test.list <- lapply(alpha.name, function(index) {

      # Create a formula for lm
      formula <-
        as.formula(paste0(paste0(index, "_diff"), "~", paste(c(
          adj.vars, group.var
        ), collapse = "+")))

      # Run lm and create a coefficient table
      lm.model <- lm(formula, data = combined_alpha)
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
      if (length(unique(combined_alpha[[group.var]])) > 2) {
        anova.tab <- broom::tidy(anova(lm.model))

        # Rearrange the table and add missing columns
        anova.tab <- anova.tab %>%
          select(
            Term = term,
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
