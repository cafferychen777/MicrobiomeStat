#' Generate Beta Diversity Volatility Test for Longitudinal Data
#'
#' This function computes a volatility test for longitudinal data on beta diversity.
#' It tests the association between beta diversity volatility and the specified group variable.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param subject.var A string indicating the column name in the metadata table that represents the subject.
#' @param time.var A string indicating the column name in the metadata table that represents the time.
#'   Ensure that it is coded as numeric.
#' @param group.var (Optional) A string indicating the column name in the metadata table that represents the grouping factor.
#' @param adj.vars (Optional) A character vector specifying the column names in the metadata table used for adjustment.
#' @param dist.name A character vector specifying which beta diversity indices to calculate.
#'   Default is "BC" (Bray-Curtis).
#' @param ... (Optional) Additional arguments to pass to internal functions.
#'
#' @return A list containing the result of the volatility test for each specified beta diversity index.
#'   Each element in the list is a tibble with the coefficients extracted from the linear model fitted
#'   for each distance, and an ANOVA table if the group variable is multi-categorical.
#'
#' @details
#' The function starts by validating the input data, processing the time variable, and calculating the beta diversity if necessary.
#' Adjustments are made based on the provided adjusting variables. The volatility of the beta diversity is
#' computed for each subject, and linear models are fitted to test the association between volatility and
#' the specified group variable. The coefficients and ANOVA results are extracted and returned for each
#' beta diversity index specified.
#'
#' @note
#' A warning message will be displayed to ensure that the time variable is coded as numeric.
#' Non-numeric coding may lead to issues in the volatility test computation.
#'
#' @seealso
#' mStat_calculate_beta_diversity, mStat_calculate_adjusted_distance
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_beta_volatility_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   dist.name = c("BC", "Jaccard")
#' )
#'
#' data(subset_T2D.obj)
#' generate_beta_volatility_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   dist.name = c("BC", "Jaccard")
#' )
#' }
#' @export
generate_beta_volatility_test_long <-
  function(data.obj,
           dist.obj = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           adj.vars = NULL,
           dist.name = c("BC"),
           ...) {

    if (is.null(dist.name)){
      return()
    }

    mStat_validate_data(data.obj)

    message(
      "The volatility test in 'generate_beta_volatility_test_long' relies on a numeric time variable.\n",
      "Please ensure that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
      "The time variable will be processed within the function if needed."
    )

    if (is.null(dist.obj)) {
      meta_tab <- data.obj$meta.dat %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      }
    }

    if (nrow(as.matrix(dist.obj[[dist.name[1]]])) > nrow(meta_tab)){
      samIDs <- rownames(meta_tab)
      dist.obj <- mStat_subset_dist(dist.obj = dist.obj, samIDs = samIDs)
    }

    test.list <- lapply(dist.name,function(dist.name){

      dist.df <- as.matrix(dist.obj[[dist.name]])

      dist.df <- dist.df %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      meta_tab <- meta_tab %>% rownames_to_column("sample")

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        dplyr::mutate(min_time_level = min(!!sym(paste0(time.var, ".subject")))[1]) %>%
        dplyr::arrange(!!sym(paste0(time.var, ".sample"))) %>%
        dplyr::mutate(prev_time_level = dplyr::lag(!!sym(paste0(time.var, ".subject")))) %>%
        filter(!!sym(paste0(time.var, ".sample")) == !!sym("prev_time_level")) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        filter(!!sym(paste0(time.var, ".subject")) != !!sym("min_time_level")) %>%
        dplyr::ungroup() %>%
        select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), !!sym(paste0(time.var, ".sample")) ,distance) %>%
        dplyr::rename(
          !!sym(subject.var) := !!sym(paste0(subject.var, ".subject")),
          !!sym(time.var) := !!sym(paste0(time.var, ".subject")),
          !!sym(paste0(time.var, ".before")) := !!sym(paste0(time.var, ".sample"))
        )

      long.df <- long.df %>%
        dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)),
               !!sym(paste0(time.var, ".before")) := as.numeric(!!sym(paste0(time.var, ".before"))))

      volatility_df <- long.df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::arrange(!!sym(time.var)) %>%
        dplyr::mutate(
          time_diff = !!sym(time.var) - !!sym(paste0(time.var, ".before"))
        ) %>%
        dplyr::filter(!is.na(distance), !is.na(time_diff)) %>%
        dplyr::filter(time_diff != 0) %>%
        dplyr::summarise(
          volatility = mean(distance / time_diff, na.rm = TRUE)
        )

      test_df <- volatility_df %>%
        dplyr::left_join(meta_tab %>%
                           select(all_of(c(subject.var, group.var))) %>%
                           dplyr::distinct(), by = subject.var, relationship = "many-to-many")

      # Test the association between the volatility and the grp.var
      formula <- as.formula(paste("volatility ~", group.var))
      test_result <- lm(formula, data = test_df)

      coef.tab <- extract_coef(test_result)

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(unique(test_df[[group.var]])) > 1) {
        anova <- anova(test_result)
        anova.tab <- anova %>% as.data.frame() %>%
          rownames_to_column("Term") %>%
          select(
            Term,
            Statistic = `F value`,
            P.Value = `Pr(>F)`
          ) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          as_tibble() %>%
          select(
            Term,
            Estimate,
            Std.Error,
            Statistic,
            P.Value
          )

        coef.tab <-
          rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
      }

      return(as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list
    names(test.list) <- dist.name

    return(test.list)
  }
