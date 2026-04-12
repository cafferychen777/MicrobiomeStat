#' Longitudinal Taxa Association Test
#'
#' Tests associations between taxa abundances and a grouping variable in
#' longitudinal microbiome data using LinDA mixed-effects models.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @details
#' Based on whether adj.vars is NULL, the formula tests:
#'
#' - When adj.vars is NOT NULL:
#'   - Tests group.var and adj.vars main effects.
#'   - Adjusted for adj.vars.
#'
#' - When adj.vars is NULL:
#'   - Tests group.var main effects only.
#'   - Unadjusted analysis.
#'
#' Subject variability is accounted for through random effects.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested comparisons between groups
#'         (e.g., "Level vs Reference (Reference)" for categorical variables,
#'         or variable name for continuous variables)
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Log2 fold change (categorical) or slope (continuous)
#'           \item \code{SE}: Standard error of the coefficient
#'           \item \code{P.Value}: Raw p-value from LinDA's statistical test
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value (Benjamini-Hochberg)
#'           \item \code{Mean.Abundance}: Mean abundance across all samples
#'           \item \code{Prevalence}: Proportion of samples where feature is present (non-zero)
#'         }
#' }
#'
#' Analysis uses LinDA mixed-effects models to test associations between taxa abundances
#' and the grouping variable, accounting for subject-level random effects.
#'
#' @examples
#' \dontrun{
#' # Example 1: Generate taxa association tests and volcano plots for the ecam dataset
#' data("ecam.obj")
#' test.list <- generate_taxa_association_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   group.var = "delivery",
#'   feature.level = c("Phylum", "Class"),
#'   feature.dat.type = c("count")
#' )
#'
#' volcano_plots_ecam <- generate_taxa_volcano_single(
#'   data.obj = ecam.obj,
#'   group.var = "delivery",
#'   test.list = test.list,
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "fdr"
#' )
#'
#'
#' # Example 2: Generate taxa association tests and volcano plots for a subset of the T2D dataset
#' data("subset_T2D.obj")
#' test.list_T2D <- generate_taxa_association_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   feature.level = "Genus",
#'   group.var = "subject_race",
#'   feature.dat.type = c("count"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.001
#' )
#'
#' volcano_plots_T2D <- generate_taxa_volcano_single(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   test.list = test.list_T2D,
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "none"
#' )
#' }
#' @export
generate_taxa_association_test_long <-
  function(data.obj,
           subject.var,
           group.var,
           adj.vars = NULL,
           prev.filter = 0,
           abund.filter = 0,
           feature.level,
           feature.dat.type = c("count", "proportion", "other"),
           ...) {
    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Extract relevant metadata and ensure it remains a data frame
    meta_tab <- data.obj$meta.dat %>% 
        as.data.frame() %>%
        select(all_of(c(group.var, adj.vars, subject.var)))

    formula <- .mStat_build_taxa_linda_formula(
      group.var = group.var,
      adj.vars = adj.vars,
      time.var = NULL,
      subject.var = subject.var,
      random_slopes = FALSE
    )

    analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)
    analysis_feature.dat.type <- if (feature.dat.type == "count") "proportion" else feature.dat.type

    # Perform analysis for each feature level
    test.list <- lapply(feature.level, function(feature.level) {

      # Aggregate features by taxonomy if not already done
      otu_tax_agg_filter <- get_taxa_data(
        analysis_data.obj,
        feature.level,
        prev.filter,
        abund.filter,
        feature.col = FALSE
      )

      # Convert filtered data back to matrix
      otu_tax_agg_filter <- as.matrix(otu_tax_agg_filter)
      meta_tab_level <- meta_tab

      if (any(colSums(otu_tax_agg_filter) == 0)) {
        pruned_inputs <- .mStat_prune_zero_total_samples(
          feature.dat = otu_tax_agg_filter,
          meta.dat = meta_tab_level
        )
        otu_tax_agg_filter <- pruned_inputs$feature.dat
        meta_tab_level <- pruned_inputs$meta.dat
      }

      linda.obj <- .mStat_run_linda(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab_level,
        formula = formula,
        feature.dat.type = analysis_feature.dat.type,
        prev.filter = 0,
        mean.abund.filter = 0,
        extra_args = list(...)
      )

      reference_level <- .mStat_get_group_reference_level(
        meta.dat = meta_tab_level,
        group.var = group.var
      )
      group_is_categorical <- !is.null(reference_level)

      # Calculate average abundance and prevalence for each feature
      prop_prev_data <- mStat_summarize_taxa_features(
        feature.dat = otu_tax_agg_filter,
        feature.level = feature.level
      )

      # Extract data frames from LInDA output
      sub_test.list <- .mStat_extract_pair_linda_outputs(
        linda_output = linda.obj$output,
        group_var = group.var,
        time_var = NULL,
        reference_level = if (group_is_categorical) reference_level else NULL
      )

      # Process and format the extracted data frames
      sub_test.list <- lapply(sub_test.list, function(df){
        mStat_format_linda_feature_results(
          result.df = df,
          feature.level = feature.level,
          feature.stats = prop_prev_data
        )
      })

      return(sub_test.list)

    })

    # Assign names to the elements of test.list based on feature levels
    names(test.list) <- feature.level

    return(test.list)
  }
