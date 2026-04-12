#' Longitudinal Per-Time-Point Differential Abundance Test
#'
#' Performs differential abundance testing at each time point separately in
#' longitudinal microbiome data using per-time LinDA fixed-effects models.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @param ... Additional arguments passed to other methods.
#' @details
#' The function integrates various data manipulations, normalization procedures, and statistical tests to assess the significance of taxa changes over time or between groups. It allows for the adjustment of covariates and handles both count and proportion data types.
#'
#' The function constructs a fixed-effects model formula based on the provided variables and performs filtering based on prevalence and abundance thresholds, with optional normalization for count data.
#'
#' Importantly, the function conducts differential abundance analysis separately for each time point in the longitudinal data. This approach allows for the identification of taxa that show significant changes at specific time points, providing insights into the dynamics of the microbiome over time.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by time points (\code{time.var} levels)
#'   \item Second level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Third level: Named by tested comparisons between groups
#'         (e.g., "Level vs Reference (Reference)" for categorical variables,
#'         or variable name for continuous variables)
#'   \item Each final element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Log2 fold change (categorical) or slope (continuous)
#'           \item \code{SE}: Standard error of the coefficient
#'           \item \code{P.Value}: Raw p-value from LinDA's statistical test
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value (Benjamini-Hochberg)
#'           \item \code{Mean.Abundance}: Mean abundance across all samples at that time point
#'           \item \code{Prevalence}: Proportion of samples where feature is present (non-zero)
#'         }
#' }
#'
#' Analysis is performed separately for each time point using LinDA fixed-effects models.
#'
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset
#' data("ecam.obj")
#'
#' # Analyzing the impact of delivery method on microbial composition over months
#' result1 <- generate_taxa_per_time_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   adj.vars = "diet",
#'   feature.level = c("Phylum", "Class"),
#'   feature.dat.type = "proportion"
#' )
#'
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_taxa_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = result1,
#'   group.var = "delivery",
#'   time.var = "month_num",
#'   feature.level = c("Phylum", "Class")
#' )
#'
#' # Example 2: Analyzing the Type 2 Diabetes dataset
#' data("subset_T2D.obj")
#'
#' # Longitudinal analysis of microbial changes in different racial groups
#' result2 <- generate_taxa_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus", "Family"),
#'   feature.dat.type = "count"
#' )
#'
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_taxa_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result2,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   feature.level = c("Genus", "Family")
#' )
#' }
#' @export
generate_taxa_per_time_test_long <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           ...) {
    # Validate the input data object to ensure it meets the required format
    data.obj <- mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Generate the fixed-effects formula for each per-time analysis
    formula <- paste(attr(stats::terms(mStat_build_formula(".response", c(group.var, adj.vars))), "term.labels"), collapse = " + ")
    if (!nzchar(formula)) {
      formula <- "1"
    }

    mStat_validate_time_var_contract(
      meta.dat = data.obj$meta.dat,
      time.var = time.var,
      context = "taxa per-time testing"
    )

    # Extract ordered time levels from the data.
    time.levels <- mStat_order_time_labels(data.obj$meta.dat[[time.var]])

    # Perform analysis for each time point.
    mStat_run_per_time_analysis(
      time.levels = time.levels,
      context = "taxa per-time testing",
      analysis_fn = function(t.level) {
        subset_data.obj <- mStat_subset_analysis_inputs_by_meta_values(
          data.obj = data.obj,
          var = time.var,
          values = t.level,
          prune.features = TRUE
        )$data.obj

        meta_tab <- subset_data.obj$meta.dat %>%
          dplyr::select(all_of(c(time.var, group.var, adj.vars, subject.var)))

        mStat_validate_group_var_contract(
          meta.dat = meta_tab,
          group.var = group.var,
          context = paste0("taxa per-time testing at ", t.level)
        )

        analysis_data.obj <- mStat_normalize_count_data_if_needed(subset_data.obj, feature.dat.type)
        analysis_feature.dat.type <- if (feature.dat.type == "count") "proportion" else feature.dat.type

        test.list <- lapply(feature.level, function(feature.level) {
          otu_tax_agg_filter <- get_taxa_data(
            analysis_data.obj,
            feature.level,
            prev.filter,
            abund.filter,
            feature.col = FALSE
          )

          linda.obj <- .mStat_run_linda(
            feature.dat = otu_tax_agg_filter,
            meta.dat = meta_tab,
            formula = formula,
            feature.dat.type = analysis_feature.dat.type,
            prev.filter = prev.filter,
            mean.abund.filter = abund.filter,
            extra_args = list(...)
          )

          reference_level <- .mStat_get_group_reference_level(
            meta.dat = meta_tab,
            group.var = group.var
          )
          group_is_categorical <- !is.null(reference_level)

          prop_prev_data <- mStat_summarize_taxa_features(
            feature.dat = otu_tax_agg_filter,
            feature.level = feature.level
          )

          sub_test.list <- .mStat_extract_pair_linda_outputs(
            linda_output = linda.obj$output,
            group_var = group.var,
            time_var = NULL,
            reference_level = if (group_is_categorical) reference_level else NULL
          )

          if (length(sub_test.list) == 0) {
            return(list())
          }

          lapply(sub_test.list, function(df) {
            mStat_format_linda_feature_results(
              result.df = df,
              feature.level = feature.level,
              feature.stats = prop_prev_data
            )
          })
        })

        names(test.list) <- feature.level
        test.list
      }
    )
  }
