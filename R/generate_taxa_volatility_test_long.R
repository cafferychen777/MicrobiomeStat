#' Longitudinal Taxa Abundance Volatility Test
#'
#' Calculates and tests taxa abundance volatility (variability over time) between groups.
#' Volatility is the mean absolute difference between consecutive time points,
#' normalized by time difference.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @param transform Character; transformation method before volatility calculation.
#'   "CLR" applies CLR transform (default). Any other value skips CLR and uses
#'   filtered abundance values directly.
#' @param ... Additional arguments passed to other methods.
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested comparisons between groups
#'         (e.g., "Level vs Reference (Reference)")
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Effect size for volatility differences between groups
#'           \item \code{SE}: Standard error of the coefficient from the linear model
#'           \item \code{P.Value}: Raw p-value from standard linear model (lm)
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value (Benjamini-Hochberg)
#'           \item \code{Mean.Abundance}: Mean abundance across all samples
#'           \item \code{Prevalence}: Proportion of samples where feature is present (non-zero)
#'         }
#' }
#'
#' This function analyzes VOLATILITY (variability over time) using standard linear models,
#' NOT LinDA. Volatility is calculated as mean absolute differences between consecutive time points.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' test.list <- generate_taxa_volatility_test_long(
#' data.obj = subset_T2D.obj,
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "sample_body_site",
#' prev.filter = 0.1,
#' abund.filter = 0.0001,
#' feature.level = c("Genus"),
#' feature.dat.type = "count",
#' transform = "CLR"
#' )
#' plot.list <- generate_taxa_volatility_volcano_long(data.obj = subset_T2D.obj,
#'                                                    group.var = "subject_race",
#'                                                    test.list = test.list,
#'                                                    feature.sig.level = 0.1,
#'                                                    feature.mt.method = "none")
#'
#' data("ecam.obj")
#' test.list <- generate_taxa_volatility_test_long(
#'   data.obj = ecam.obj,
#'   time.var = "month_num",
#'   subject.var = "subject.id",
#'   group.var = "antiexposedall",
#'   adj.vars = "delivery",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.level = c("Order", "Family", "Genus"),
#'   feature.dat.type = "proportion",
#'   transform = "CLR"
#' )
#' plot.list <- generate_taxa_volatility_volcano_long(
#'   data.obj = ecam.obj,
#'   group.var = "antiexposedall",
#'   test.list = test.list,
#'   feature.sig.level = 0.2,
#'   feature.mt.method = "none"
#' )
#' }
#' @export
generate_taxa_volatility_test_long <- function(data.obj,
                                               time.var,
                                               subject.var,
                                               group.var,
                                               adj.vars = NULL,
                                               prev.filter = 0,
                                               abund.filter = 0,
                                               feature.level,
                                               feature.dat.type = c("count", "proportion", "other"),
                                               transform = "CLR",
                                               ...) {
  # Validate the input data object
  data.obj <- mStat_validate_data(data.obj)

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Inform the user about the importance of numeric time variable
  message(
    "The volatility calculation relies on a numeric time variable.\n",
    "Please check that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
    "You can ensure the time variable is numeric by mutating it in the metadata."
  )

  # Convert time variable to numeric
  data.obj$meta.dat[[time.var]] <- mStat_coerce_time_to_numeric(
    data.obj$meta.dat[[time.var]],
    time.var = time.var,
    context = "taxa volatility analysis"
  )

  # Extract relevant variables from metadata
  meta_tab <- mStat_meta_with_sample(
    data.obj$meta.dat %>%
      dplyr::select(all_of(c(
        time.var, subject.var, group.var, adj.vars
      )))
  )

  # Extract group levels and set reference level
  group_level <- meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels
  reference_level <- group_level[1]

  # CLR-transformed data contains negative values, so abundance filtering must
  # be disabled explicitly.
  if (transform == "CLR"){
    abund.filter <- -Inf
  }

  # Create a formula for the linear model
  formula_str <- paste("volatility ~", group.var)
  if (!is.null(adj.vars)) {
    formula_str <- paste(formula_str, "+", paste(adj.vars, collapse = " + "))
  }
  formula <- as.formula(formula_str)

  # Perform analysis for each taxonomic level
  test.list <- lapply(feature.level, function(feature.level) {

    # Normalize count data if necessary
    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    # Aggregate data to the specified taxonomic level if necessary
    otu_tax_agg <- get_taxa_data(data.obj, feature.level, feature.col = FALSE)

    # Calculate average abundance and prevalence for each feature
    prop_prev_data <- mStat_summarize_taxa_features(
      feature.dat = otu_tax_agg,
      feature.level = feature.level
    )

    if (identical(transform, "CLR")) {
      otu_tax_long <- mStat_prepare_taxa_clr_long_data(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        prev.filter = prev.filter,
        abund.filter = abund.filter,
        value_col = "value"
      )
    } else {
      otu_tax_agg_filtered <- mStat_filter(
        otu_tax_agg,
        prev.filter = prev.filter,
        abund.filter = abund.filter
      )
      otu_tax_long <- mStat_prepare_taxa_long_data(
        feature.dat = otu_tax_agg_filtered,
        feature.level = feature.level,
        value_col = "value",
        meta.dat = NULL,
        feature_in_column = FALSE,
        join = "left"
      )
    }

    taxa.levels <- otu_tax_long %>% select(all_of(feature.level)) %>% pull() %>% unique()

    # Perform volatility analysis for each taxon
    sub_test.list <-
      lapply(taxa.levels, function(taxon) {

        taxa_df <- otu_tax_long %>%
          dplyr::filter(!!sym(feature.level) == taxon) %>%
          dplyr::left_join(meta_tab, by = "sample")

        # Calculate volatility for each subject
        volatility_df <- taxa_df %>%
          dplyr::arrange(!!sym(subject.var),!!sym(time.var)) %>%
          dplyr::group_by(!!sym(subject.var)) %>%
          dplyr::mutate(
            diff_value = abs(value - dplyr::lag(value)),
            diff_time = !!sym(time.var) - dplyr::lag(!!sym(time.var))
          ) %>%
          dplyr::filter(!is.na(diff_value),!is.na(diff_time)) %>%
          dplyr::filter(diff_time != 0) %>%
          dplyr::summarize(volatility = mean(diff_value / diff_time),
                    .groups = 'drop')

        test_df <- volatility_df %>%
          dplyr::left_join(meta_tab %>%
                             select(all_of(c(subject.var, group.var, adj.vars))) %>%
                             dplyr::distinct(),
                           by = subject.var)

        if (nrow(test_df) <= 1) {
          return(NULL)
        }

        model_terms <- c(group.var, adj.vars)
        model_terms <- model_terms[!is.null(model_terms)]

        valid_terms <- model_terms[vapply(model_terms, function(term) {
          term_values <- test_df[[term]]
          if (is.numeric(term_values)) {
            return(stats::var(term_values, na.rm = TRUE) > 0)
          }
          length(unique(stats::na.omit(term_values))) > 1
        }, logical(1))]

        model_formula <- if (length(valid_terms) > 0) {
          stats::as.formula(paste("volatility ~", paste(valid_terms, collapse = " + ")))
        } else {
          stats::as.formula("volatility ~ 1")
        }

        # Fit linear model
        test_result <- lm(model_formula, data = test_df)

        coef.tab <- extract_coef(test_result)

        # Perform ANOVA if group variable has more than two levels
        if (length(unique(taxa_df[[group.var]])) > 2) {
          anova <- anova(test_result)
          anova.tab <- anova %>% as.data.frame() %>%
            tibble::rownames_to_column("Term") %>%
            dplyr::select(
                   Term,
                   Statistic = `F value`,
                   P.Value = `Pr(>F)`) %>%
            tibble::as_tibble() %>%
            dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
            select(
              Term,
              Estimate,
              Std.Error,
              Statistic,
              P.Value
            )

          coef.tab <- rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
        }
        return(tibble::as_tibble(coef.tab))
      })

    # Assign names to the elements of test.list
    names(sub_test.list) <- otu_tax_long %>% select(all_of(feature.level)) %>% pull() %>% unique()

    # Find all unique terms related to the group variable
    unique_terms <- grep(paste0("^", group.var, "$|^", group.var, ".*"), unique(unlist(lapply(sub_test.list, function(df) unique(df$Term)))), value = TRUE)

    # Extract and process results for each term
    result_list <- lapply(unique_terms, function(term) {
      do.call(rbind, lapply(sub_test.list, function(df) {
        df %>% dplyr::filter(Term == term)
      })) %>%
        dplyr::mutate(!!sym(feature.level) := names(sub_test.list)) %>%
        dplyr::left_join(prop_prev_data, by = feature.level) %>%
        dplyr::mutate(Adjusted.P.Value = p.adjust(P.Value, method = "fdr")) %>%
        dplyr::select(all_of(c(feature.level, "Estimate", "Std.Error", "P.Value", "Adjusted.P.Value", "avg_abundance", "prevalence"))) %>%
        dplyr::rename(Coefficient = Estimate,
                      SE = Std.Error,
                      Variable = feature.level,
                      Mean.Abundance = avg_abundance,
                      Prevalence = prevalence)
    })

    # Name the result list
    names(result_list) <- unique_terms

    # Rename the result list elements to include reference level information
    new_names <- sapply(names(result_list), function(name) {
      if (grepl(paste0("^", group.var), name) && !grepl(paste0("^", group.var, "$"), name)) {
        sub_name <- sub(paste0(group.var), "", name)
        return(paste(sub_name, "vs", reference_level, "(Reference)"))
      }
      return(name)
    })

    names(result_list) <- new_names

    return(result_list)
  })

  names(test.list) <- feature.level

  return(test.list)
}
