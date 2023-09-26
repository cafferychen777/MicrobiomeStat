#' Compute and analyze changes in beta diversity
#'
#' This function calculates beta diversity measures between time points,
#' performs linear modeling, and generates results.
#'
#' For each subject, it calculates the change in beta diversity between
#' baseline time point and later time point. This within-subject change
#' is used as the outcome in linear models.
#'
#' Adjustment for covariates is supported. Both ANOVA and coefficient
#' tables are provided if a grouping variable is specified.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param time.var The name of the column in metadata containing the time variable.
#'                This should be a column with time points for each sample. Required
#'                to identify pairs of samples for the same subject across time.
#' @param subject.var The name of the column in metadata containing the subject IDs.
#'                    This should uniquely identify each subject in the study. Required
#'                    to identify samples that belong to the same subject.
#' @param group.var The name of the column in metadata containing the grouping variable
#'                 to use in linear modeling. This grouping variable will be used as a
#'                 predictor in the linear models for beta diversity change. Optional.
#' @param adj.vars Vector of names of additional variables in metadata to be included
#'                as covariates in the linear models. Can be empty.
#' @param change.base The baseline time point value in the time variable to be used
#'                   as the reference for calculating beta diversity change. Required
#'                   if time.var contains multiple time points. If NULL, the first
#'                   time point will be used as change.base automatically.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#'
#' @examples
#' \dontrun{
#'
#' # Load packages
#' library(vegan)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # Generate beta diversity test
#' generate_beta_change_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   subject.var = "subject",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   change.base = "1",
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' }
#' @return A named list containing linear modeling results for each beta diversity metric.

#' Each list element corresponds to one of the distance metrics specified in \code{dist.name}.

#' It contains a coefficient table from fitting a linear model with the beta diversity change as
#' response and the \code{group_var} and \code{adj_vars} as predictors.

#' If \code{group_var} has multiple levels, ANOVA results are also included after the coefficients.

#' Column names include:
#' Term, Estimate, Std.Error, Statistic, P.Value
#' @export
#' @name generate_beta_change_test_pair
generate_beta_change_test_pair <-
  function(data.obj,
           dist.obj = NULL,
           time.var = NULL,
           subject.var,
           group.var,
           adj.vars,
           change.base = NULL,
           dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

    if (is.null(dist.obj)&!is.null(data.obj)) {
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      metadata <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var,group.var,time.var, adj.vars))) %>% rownames_to_column("sample")
    } else {
      if (!is.null(data.obj)){
        metadata <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var,group.var,time.var, adj.vars))) %>% rownames_to_column("sample")
      } else {
        metadata <- attr(dist.obj[[dist.name[1]]], "labels")  %>% dplyr::select(all_of(c(subject.var,group.var,time.var,adj.vars))) %>% rownames_to_column("sample")
      }
    }

    if (is.null(change.base)){
      change.base <- unique(metadata %>% dplyr::select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'meta.dat' data frame: ", change.base)
    }

    change.after <-
      unique(metadata %>% dplyr::select(all_of(c(time.var))))[unique(metadata %>% dplyr::select(all_of(c(time.var)))) != change.base]

    test.list <- lapply(dist.name, function(dist.name){
      dist.df <- as.matrix(dist.obj[[dist.name]]) %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(metadata, by = "sample") %>%
        dplyr::left_join(metadata, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      long.df <- long.df %>% dplyr::left_join(metadata %>% dplyr::select(-any_of(c(time.var))) %>% dplyr::distinct(), by = subject.var)

      # Create a formula for lm
      formula <-
        as.formula(paste0("distance", "~", paste(c(
          adj.vars, group.var
        ), collapse = "+")))

      # Run lm and create a coefficient table
      lm.model <- lm(formula, data = long.df)
      coef.tab <- broom::tidy(summary(lm.model))

      # Rearrange the table
      coef.tab <-
        coef.tab %>% dplyr::select(
          Term = term,
          Estimate = estimate,
          Std.Error = std.error,
          Statistic = statistic,
          P.Value = p.value
        )

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(unique(metadata[[group.var]])) > 2) {
        anova.tab <- broom::tidy(anova(lm.model))

        # Rearrange the table and add missing columns
        anova.tab <- anova.tab %>%
          dplyr::select(
            Term = term,
            Statistic = statistic,
            df = df,
            P.Value = p.value
          ) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA)

        # Reorder the columns to match coef.tab
        anova.tab <- anova.tab %>%
          dplyr::select(
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
    names(test.list) <- dist.name

    return(test.list)

  }
