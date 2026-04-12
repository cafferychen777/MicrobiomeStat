#' @title Paired/Longitudinal Taxa Differential Abundance Test
#' @description Performs differential abundance analysis for paired or longitudinal data using
#'   linear mixed-effects models via LinDA, accounting for subject-level correlations.
#' @inheritParams mStat_data_obj_doc
#' @param change.base Value indicating the base/reference level for the time variable.
#'   If NULL, the first level is used.
#' @param ref.level Character string specifying the reference level for categorical group.var.
#'   If NULL, the first level alphabetically is used. Ignored for continuous variables.
#' @param feature.mt.method Character string specifying the multiple-testing correction method.
#'   One of "fdr", "bonferroni", or "none".
#' @param feature.sig.level Numeric significance threshold used by the testing procedure.
#' @param ... Additional parameters passed to the linda function.
#' @return A nested list of differential abundance test results by feature level and comparison.
#' @name generate_taxa_test_pair
#' @export
generate_taxa_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           change.base = NULL,
           group.var,
           ref.level = NULL,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           feature.mt.method = c("fdr", "bonferroni", "none"),
           feature.sig.level = 0.05,
           ...) {
    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)
    feature.mt.method <- match.arg(feature.mt.method)

    p.adj.method <- switch(
      feature.mt.method,
      fdr = "BH",
      bonferroni = "bonferroni",
      none = "none"
    )

    extra_args <- list(...)
    extra_args[c("p.adj.method", "alpha")] <- NULL

    # Extract relevant metadata
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    mStat_validate_group_var_contract(
      meta.dat = meta_tab,
      group.var = group.var,
      subject.var = subject.var,
      context = "taxa paired testing"
    )

    # Set reference level for group variable if it is categorical
    if (!is.null(group.var)) {
      if (is.factor(meta_tab[[group.var]]) || is.character(meta_tab[[group.var]])) {
        # Convert to factor if character
        meta_tab[[group.var]] <- as.factor(meta_tab[[group.var]])

        # Get available levels
        available_levels <- levels(meta_tab[[group.var]])

        # Validate and set reference level
        if (!is.null(ref.level)) {
          if (!(ref.level %in% available_levels)) {
            stop(
              "ref.level '", ref.level, "' not found in group.var '", group.var, "'. ",
              "Available levels: ", paste(available_levels, collapse = ", ")
            )
          }
          meta_tab[[group.var]] <- relevel(meta_tab[[group.var]], ref = ref.level)
          message("Reference level for '", group.var, "': ", ref.level)
        } else {
          message(
            "Reference level for '", group.var, "': ", available_levels[1],
            " (alphabetically first)"
          )
        }
      } else if (!is.null(ref.level)) {
        # Warn if ref.level is specified for non-categorical variable
        warning(
          "ref.level is ignored because group.var '", group.var,
          "' is not categorical (factor or character)."
        )
      }
    }

    # Build primary and fallback formulas (fallback removes random slopes)
    formula <- .mStat_build_taxa_linda_formula(
      group.var = group.var,
      adj.vars = adj.vars,
      time.var = time.var,
      subject.var = subject.var,
      random_slopes = TRUE
    )

    formula_corrected <- .mStat_build_taxa_linda_formula(
      group.var = group.var,
      adj.vars = adj.vars,
      time.var = time.var,
      subject.var = subject.var,
      random_slopes = FALSE
    )

    analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)
    analysis_feature.dat.type <- if (feature.dat.type == "count") "proportion" else feature.dat.type

    # Perform analysis for each taxonomic level
    test.list <- lapply(feature.level, function(feature.level) {
      # Aggregate data by taxonomy if necessary
      otu_tax_agg_filter <- get_taxa_data(
        analysis_data.obj,
        feature.level,
        prev.filter,
        abund.filter,
        feature.col = FALSE
      )
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
        fallback_formula = formula_corrected,
        fallback_message = "Due to the above error, a simpler model will be used for fitting.",
        feature.dat.type = analysis_feature.dat.type,
        extra_args = c(
          list(
            p.adj.method = p.adj.method,
            alpha = feature.sig.level
          ),
          extra_args
        )
      )

      reference_level <- .mStat_get_group_reference_level(
        meta.dat = meta_tab_level,
        group.var = group.var
      )

      # Calculate average abundance and prevalence for each feature
      prop_prev_data <- mStat_summarize_taxa_features(
        feature.dat = otu_tax_agg_filter,
        feature.level = feature.level
      )

      # Extract and process results from LInDA output
      sub_test.list <- .mStat_extract_pair_linda_outputs(
        linda_output = linda.obj$output,
        group_var = group.var,
        time_var = time.var,
        reference_level = if (!is.null(group.var)) reference_level else NULL
      )

      # Format and annotate results
      sub_test.list <- lapply(sub_test.list, function(df){
        mStat_format_linda_feature_results(
          result.df = df,
          feature.level = feature.level,
          feature.stats = prop_prev_data,
          include_significant = TRUE
        )
      })

      return(sub_test.list)
    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)
  }
