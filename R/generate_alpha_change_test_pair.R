#' Compare alpha diversity between time points
#'
#' This function performs paired tests to compare alpha diversity metrics
#' between two time points in a microbiome dataset.
#'
#' For each alpha diversity metric, the difference between the two time points
#' is calculated for each subject. The difference is used as the response
#' variable in a linear model, with grouping and adjustment variables as
#' predictors.
#'
#' The linear model coefficients are extracted into a results table. If the
#' grouping variable has multiple levels, ANOVA is performed to test for overall
#' significance.
#'
#' Options are provided to customize the alpha diversity difference calculation
#' (e.g. log-fold change) and adjust for covariates.
#'
#' In summary, this provides a flexible paired test for analyzing changes in
#' alpha diversity over time, accounting for groups and covariates.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou". Previously named as `alpha.index`.
#' @param depth An integer specifying the sequencing depth for the "Rarefy" and "Rarefy-TSS" methods.
#' If NULL, no rarefaction is performed.
#' @param time.var Character string specifying the column name in metadata containing
#'                time values for each sample. Required to identify pairs of time
#'                points to calculate changes between.
#' @param subject.var Character string specifying the column name in metadata containing
#'                    unique subject IDs. Required to pair samples from the same subject
#'                    across time points.
#' @param group.var Character string specifying the column name in metadata containing
#'                 grouping categories. Used as a predictor in the models to test for
#'                 differences in changes between groups. Optional, can be NULL.
#' @param adj.vars Character vector specifying column names in metadata containing
#'                covariates to adjust for in the linear models. Optional, can be
#'                left NULL if no adjustment is needed.
#' @param change.base The name of the baseline time point for calculating changes in alpha diversity. If NULL, the first unique time point in the data will be used.
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. This allows flexible options to quantify change:
#'
#'   - If a function is provided: The function will be applied to compare alpha diversity
#'     at timepoint t vs baseline t0. The function should take two arguments
#'     representing the alpha diversity values at t and t0. For instance, a custom function to
#'     calculate the percentage change might look like:
#'     \preformatted{
#'       percentage_change <- function(t, t0) {
#'         return ((t - t0) / t0) * 100
#'       }
#'     }
#'     You can then pass this function as the value for `alpha.change.func`.
#'
#'   - If a string is provided, the following options are supported:
#'     - 'log fold change': Calculates the log2 fold change of alpha diversity at t compared to t0.
#'     - 'absolute change': Calculates the absolute difference in alpha diversity at t compared to t0.
#'     - Any other value: A warning will be given that the provided method is not recognized,
#'       and the default method ('absolute change') will be used.
#'
#'   - Default behavior (if no recognized string or function is provided) is to compute the absolute difference between t and t0.
#' @return A list of tables, one for each alpha diversity metric, summarizing the results of the statistical tests.
#' Each table contains the following columns: Term (the name of the variable in the model), Estimate (the estimated coefficient),
#' Std.Error (the standard error of the coefficient), Statistic (the t or F statistic), P.Value (the p-value of the test).
#'
#' @details
#' This function performs the following statistical tests:
#'
#' - For each alpha diversity metric, a linear model is fitted with the difference in alpha diversity between two time points as the response, and the grouping variable and any adjustment variables as predictors.
#' The linear model coefficients are extracted and formatted into a results table.
#'
#' - If the grouping variable has more than 2 levels, ANOVA is performed on the linear model to obtain the overall significance of the grouping variable.
#' The ANOVA results are also formatted into the results table.
#'
#' So in summary, the function provides a paired test for the change in alpha diversity between two time points,
#' with options to adjust for covariates and test for differences between groups. The output contains coefficient tables
#' summarizing the results of the linear models and ANOVA tests.
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
#'   time.var = "time",
#'   alpha.name = c("shannon"),
#'   subject.var = "subject",
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
#'   time.var = "Day",
#'   alpha.name = c("shannon"),
#'   subject.var = "subject",
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
#'   time.var = "Antibiotic",
#'   alpha.name = c("shannon"),
#'   subject.var = "MouseID",
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
           time.var,
           subject.var,
           group.var,
           adj.vars = NULL,
           change.base,
           alpha.change.func = "log fold change") {

    if (is.null(alpha.name)){
      return()
    }

    if (is.null(alpha.obj)) {
      if (!is.null(depth)) {
        message(
          "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
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

    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
        subject.var, group.var, time.var, adj.vars
      )))

    time_varying_info <- NULL

    if (!is.null(adj.vars)){
      # Use the modified mStat_identify_time_varying_vars function
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      if (length(time_varying_info$non_time_varying_vars) > 0){
        alpha.obj <- mStat_calculate_adjusted_alpha_diversity(alpha.obj, meta.dat = meta_tab, time_varying_info$non_time_varying_vars)
      }
    }

    # Convert the alpha.obj list to a data frame
    alpha_df <-
      dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"),
                 by = c("sample"))

    if (is.null(change.base)){
      change.base <- unique(alpha_df %>% dplyr::select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'alpha_df' data frame: ", change.base)
    }

    change.after <-
      unique(alpha_df %>% dplyr::select(all_of(c(time.var))))[unique(alpha_df %>% dplyr::select(all_of(c(time.var)))) != change.base]

    alpha_grouped <- alpha_df %>% dplyr::group_by(!!sym(time.var))
    alpha_split <- split(alpha_df, f = alpha_grouped[[time.var]])

    alpha_time_1 <- alpha_split[[change.base]]
    alpha_time_2 <- alpha_split[[change.after]]

    combined_alpha <- alpha_time_1 %>%
      dplyr::inner_join(
        alpha_time_2,
        by = c(subject.var, group.var),
        suffix = c("_time_1", "_time_2"),
        relationship = "many-to-many"
      )

    diff_columns <- lapply(alpha.name, function(index) {
      diff_col_name <- paste0(index, "_diff")

      if (is.function(alpha.change.func)) {
        combined_alpha <- combined_alpha %>%
          dplyr::mutate(!!diff_col_name := alpha.change.func(!!sym(paste0(
            index, "_time_2"
          )),!!sym(paste0(
            index, "_time_1"
          )))) %>%
          dplyr::select(all_of(diff_col_name))
      } else {
        if (alpha.change.func == "log fold change") {
          combined_alpha <- combined_alpha %>%
            dplyr::mutate(!!sym(diff_col_name) := log2(!!sym(paste0(
              index, "_time_2"
            )) / !!sym(paste0(
              index, "_time_1"
            )))) %>%
            dplyr::select(all_of(c(diff_col_name)))
        } else
          if (alpha.change.func == "absolute change") {
            combined_alpha <- combined_alpha %>%
              dplyr::mutate(!!diff_col_name := !!sym(paste0(index, "_time_2")) -!!sym(paste0(index, "_time_1"))) %>%
              dplyr::select(all_of(diff_col_name))
          } else {
            message(paste("No valid alpha.change.func provided for", index, ". Defaulting to 'absolute change'."))
            combined_alpha <- combined_alpha %>%
              dplyr::mutate(!!diff_col_name := !!sym(paste0(index, "_time_2")) -!!sym(paste0(index, "_time_1"))) %>%
              dplyr::select(all_of(diff_col_name))
          }
      }
    })

    combined_alpha <- dplyr::bind_cols(combined_alpha, diff_columns)

      if (length(time_varying_info$time_varying_vars) > 0) {
        names_map <- setNames(paste0(time_varying_info$time_varying_vars, "_time_2"), time_varying_info$time_varying_vars)
        combined_alpha <- combined_alpha %>%
          rename(!!!names_map)
      }

    # Generate tests
    test.list <- lapply(alpha.name, function(index) {

      # Create a formula for lm
      formula <-
        as.formula(paste0(paste0(index, "_diff"), "~", paste(c(
          time_varying_info$time_varying_vars, group.var
        ), collapse = "+")))

      # Run lm and create a coefficient table
      lm.model <- lm(formula, data = combined_alpha)

      summary <- summary(lm.model)

      coef.tab <- summary$coefficients %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Term") %>%
        dplyr::select(
          Term,
          Std.Error = `Std. Error`,
          Statistic = `t value`,
          P.Value = `Pr(>|t|)`,
          Estimate
        ) %>% as_tibble()

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(unique(combined_alpha[[group.var]])) > 2) {
        anova <- anova(lm.model)
        anova.tab <- as.data.frame(anova) %>%
          rownames_to_column("Term") %>%
          dplyr::select(Term,
                        Statistic = `F value`,
                        P.Value = `Pr(>F)`) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          as_tibble()

        coef.tab <-
          rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
      }

      return(coef.tab)
    })

    # Assign names to the elements of test.list
    names(test.list) <- alpha.name

    return(test.list)
  }
