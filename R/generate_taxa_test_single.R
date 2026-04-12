#' Helper function to perform linear model analysis for "other" data type
#'
#' @param feature.dat Feature data matrix
#' @param meta.dat Metadata data frame
#' @param formula Formula string for the model
#' @param group.var Group variable name
#' @return List with analysis results
#' @noRd
perform_lm_analysis <- function(feature.dat,
                                meta.dat,
                                formula,
                                group.var,
                                p.adj.method = "BH",
                                alpha = 0.05) {
  # Check if group.var is categorical (factor or character)
  is_categorical <- is.factor(meta.dat[, group.var]) || is.character(meta.dat[, group.var])

  # Get reference level only if categorical
  if (is_categorical) {
    reference_level <- .mStat_get_group_reference_level(
      meta.dat = meta.dat,
      group.var = group.var
    )
    # Get group levels (excluding reference)
    group_levels <- levels(as.factor(meta.dat[, group.var]))
    comparison_levels <- group_levels[group_levels != reference_level]
  } else {
    # For continuous variables, no reference level needed
    reference_level <- NULL
    comparison_levels <- NULL
  }

  # Initialize results list
  results <- list()

  # Perform linear regression for each feature
  feature_results <- list()

  for (i in 1:nrow(feature.dat)) {
    feature_name <- rownames(feature.dat)[i]
    feature_values <- as.numeric(feature.dat[i, ])

    # Create data frame for regression
    reg_data <- data.frame(
      y = feature_values,
      meta.dat
    )

    # Fit linear model
    lm_formula <- as.formula(paste("y ~", formula))
    tryCatch({
      lm_fit <- lm(lm_formula, data = reg_data)
      lm_summary <- summary(lm_fit)

      # Extract coefficients for group variable
      coef_table <- lm_summary$coefficients
      group_prefix_pattern <- paste0("^", mStat_escape_regex(group.var))
      group_coefs <- coef_table[grep(group_prefix_pattern, rownames(coef_table)), , drop = FALSE]

      if (nrow(group_coefs) > 0) {
        for (j in 1:nrow(group_coefs)) {
          coef_name <- rownames(group_coefs)[j]

          if (is_categorical) {
            # For categorical variables, extract group value from coefficient name
            group_value <- gsub(group_prefix_pattern, "", coef_name)
          } else {
            # For continuous variables, use the variable name itself
            group_value <- group.var
          }

          if (!group_value %in% names(feature_results)) {
            feature_results[[group_value]] <- data.frame(
              feature = character(),
              log2FoldChange = numeric(),
              lfcSE = numeric(),
              pvalue = numeric(),
              stringsAsFactors = FALSE
            )
          }

          feature_results[[group_value]] <- rbind(
            feature_results[[group_value]],
            data.frame(
              feature = feature_name,
              log2FoldChange = group_coefs[j, "Estimate"],
              lfcSE = group_coefs[j, "Std. Error"],
              pvalue = group_coefs[j, "Pr(>|t|)"],
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }, error = function(e) {
      # Skip features that cause errors
      warning(paste("Error fitting model for feature", feature_name, ":", e$message))
    })
  }

  # Format results to match linda output structure
  output <- list()
  for (group_level in names(feature_results)) {
    df <- feature_results[[group_level]]
    if (nrow(df) > 0) {
      # Add multiple testing correction
      df$padj <- p.adjust(df$pvalue, method = p.adj.method)
      df$reject <- df$padj <= alpha
      df$baseMean <- NA  # Not applicable for linear models
      df$stat <- df$log2FoldChange / df$lfcSE
      df$df <- nrow(meta.dat) - length(all.vars(as.formula(paste("~", formula)))) - 1

      # Set rownames
      rownames(df) <- df$feature
      df$feature <- NULL

      # Reorder columns to match linda output
      df <- df[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "reject", "df")]

      # Format output name based on variable type
      if (is_categorical) {
        output[[paste0(group.var, group_level, ":")]] <- df
      } else {
        # For continuous variables, use variable name only
        output[[paste0(group.var, ":")]] <- df
      }
    }
  }

  return(list(output = output))
}

#' Differential Abundance Testing for Single Time Point
#'
#' Performs differential abundance analysis using LinDA (for count/proportion data)
#' or standard linear models (for other data types). Automatically detects whether
#' the predictor is categorical or continuous.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @param t.level Character string specifying the time level/value to subset data to.
#'   Default NULL does not subset data.
#' @param ref.level Character string specifying the reference level for categorical group.var.
#'   If NULL, the first level alphabetically is used. Ignored for continuous variables.
#' @param feature.mt.method Character string specifying the multiple-testing correction method.
#'   One of "fdr", "bonferroni", or "none".
#' @param feature.sig.level Numeric significance threshold used by the testing procedure.
#' @param ... Additional arguments passed to the linda function.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by comparisons
#'         \itemize{
#'           \item For categorical \code{group.var}: One element per non-reference level,
#'                 named as "Level vs Reference (Reference)" (e.g., "Treatment vs Control (Reference)")
#'           \item For continuous \code{group.var}: One element named by the variable name (e.g., "age")
#'         }
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Log2 fold change (for categorical) or slope coefficient (for continuous)
#'           \item \code{SE}: Standard error of the coefficient
#'           \item \code{P.Value}: Raw p-value from the statistical test
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value using Benjamini-Hochberg method
#'           \item \code{Mean.Abundance}: Mean abundance of the feature across all samples
#'           \item \code{Prevalence}: Proportion of samples where the feature is present (non-zero)
#'         }
#' }
#'
#' @details
#' ## Statistical Methods
#'
#' The function automatically selects the appropriate statistical method based on \code{feature.dat.type}:
#'
#' **For "count" or "proportion" data:**
#' Uses LinDA (Linear models for Differential Abundance), which:
#' \itemize{
#'   \item Handles compositional nature of microbiome data
#'   \item Performs internal normalization (no pre-normalization needed)
#'   \item Handles zero-inflation via pseudo-count or imputation
#'   \item Corrects bias from centered log-ratio transformation
#' }
#'
#' **For "other" data:**
#' Uses standard linear regression models, suitable for:
#' \itemize{
#'   \item Pre-transformed data (e.g., log-transformed, CLR-transformed)
#'   \item Data where compositional adjustments are not needed
#' }
#'
#' ## Variable Type Handling
#'
#' The function automatically detects the type of \code{group.var}:
#'
#' **Categorical variables (factor or character):**
#' \itemize{
#'   \item Uses dummy coding with k-1 coefficients for k levels
#'   \item First level (alphabetically) becomes the reference
#'   \item Each coefficient compares a level to the reference
#'   \item Example: For treatment groups A, B, C, output includes "B vs A (Reference)" and "C vs A (Reference)"
#' }
#'
#' **Continuous variables (numeric or integer):**
#' \itemize{
#'   \item Uses single slope coefficient
#'   \item Tests for linear association with abundance
#'   \item Coefficient represents effect per unit increase
#'   \item Example: For age, coefficient shows abundance change per year
#' }
#'
#' ## Analysis Workflow
#'
#' 1. **Data Subsetting**: If \code{time.var} and \code{t.level} provided, subsets to that timepoint
#'
#' 2. **Feature Aggregation**: Aggregates features to specified taxonomic level if not "original"
#'
#' 3. **Filtering**: Applies prevalence and abundance filters
#'
#' 4. **Statistical Testing**: Runs LinDA or linear models depending on \code{feature.dat.type}
#'
#' 5. **Results Compilation**: Extracts coefficients, standard errors, p-values, and calculates
#'    FDR-adjusted p-values using Benjamini-Hochberg method
#'
#' 6. **Metadata Addition**: Adds mean abundance and prevalence for each feature
#'
#' The function supports covariate adjustment via \code{adj.vars} and allows taxonomic
#' aggregation at multiple levels for customized analyses.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' test.list <- generate_taxa_test_single(
#'     data.obj = peerj32.obj,
#'     time.var = "time",
#'     t.level = "2",
#'     group.var = "group",
#'     adj.vars = "sex",
#'     feature.dat.type = "count",
#'     feature.level = c("Phylum","Genus","Family"),
#'     prev.filter = 0.1,
#'     abund.filter = 0.0001,
#' )
#' plot.list <- generate_taxa_volcano_single(
#'     data.obj = peerj32.obj,
#'     group.var = "group",
#'     test.list = test.list,
#'     feature.sig.level = 0.1,
#'     feature.mt.method = "none"
#' )
#' }
#' @export
generate_taxa_test_single <- function(data.obj,
                                      time.var = NULL,
                                      t.level = NULL,
                                      group.var,
                                      ref.level = NULL,
                                      adj.vars = NULL,
                                      prev.filter = 0,
                                      abund.filter = 0,
                                      feature.level,
                                      feature.dat.type = c("count", "proportion", "other"),
                                      feature.mt.method = c("fdr", "bonferroni", "none"),
                                      feature.sig.level = 0.05,
                                      ...) {
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

  # Validate the input data object
  data.obj <- mStat_validate_data(data.obj)

  context <- mStat_prepare_taxa_single_context(
    data.obj = data.obj,
    time.var = time.var,
    t.level = t.level,
    group.var = group.var
  )
  data.obj <- context$data.obj

  # Extract relevant variables from the metadata
  meta_tab <-
    select_meta_vars(data.obj$meta.dat, time.var, group.var, adj.vars)

  # Set reference level for group variable if it is categorical
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

  # Construct the formula for the linear model
  formula <- group.var

  # Add adjustment variables to the formula if specified
  if (!is.null(adj.vars)) {
    adj.vars_string <- paste(adj.vars, collapse = " + ")
    formula <- paste(formula, "+", adj.vars_string)
  }

  # Set filters to 0 if the feature data type is "other"
  if (feature.dat.type == "other") {
    prev.filter <- 0
    abund.filter <- 0
  }

  abund.filter <- mStat_adjust_other_abundance_filter(
    data.obj = data.obj,
    feature.dat.type = feature.dat.type,
    abund.filter = abund.filter
  )

  # Note: Normalization is handled internally by the LinDA function
  # We skip pre-normalization to preserve pre-computed feature.agg.list
  if (feature.dat.type == "count") {
    message(
      "Your data is in raw count format. Normalization will be handled ",
      "internally by the LinDA analysis function."
    )
  }

  # Perform differential abundance testing for each specified taxonomic level
  test.list <- lapply(feature.level, function(feature.level) {
    # Aggregate data to the specified taxonomic level if necessary.
    otu_tax_agg_filter <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter, feature.col = FALSE)
    meta_tab_level <- meta_tab

    # Add a check before running linda.
    if (nrow(otu_tax_agg_filter) == 0 || ncol(otu_tax_agg_filter) == 0) {
      warning("No features remain after filtering. Consider using less stringent filter thresholds.")
      return(list())
    }

    # Remove samples with zero total abundance for this feature level only.
    if (any(colSums(otu_tax_agg_filter) == 0)) {
      pruned_inputs <- .mStat_prune_zero_total_samples(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab_level
      )
      if (is.null(pruned_inputs)) {
        warning("No samples remain after filtering.")
        return(list())
      }
      otu_tax_agg_filter <- pruned_inputs$feature.dat
      meta_tab_level <- pruned_inputs$meta.dat
    }

    # Choose analysis method based on feature data type.
    if (feature.dat.type == "other") {
      # Use linear models for "other" data type (e.g., log-transformed data).
      lm.results <- perform_lm_analysis(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab_level,
        formula = formula,
        group.var = group.var,
        p.adj.method = p.adj.method,
        alpha = feature.sig.level
      )
    } else {
      # Perform LinDA (Linear models for Differential Abundance) analysis.
      # Muffle linda's own "all filtered" warning; we emit our own below.
      linda.obj <- .mStat_run_linda(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab_level,
        formula = formula,
        feature.dat.type = feature.dat.type,
        prev.filter = prev.filter,
        mean.abund.filter = abund.filter,
        extra_args = c(
          list(
            p.adj.method = p.adj.method,
            alpha = feature.sig.level
          ),
          extra_args
        ),
        muffle_all_filtered_warning = TRUE
      )
    }

    # Check if linda returned empty results (all features filtered internally).
    if (feature.dat.type != "other" && length(linda.obj$output) == 0) {
      warning("No features remain after filtering. Consider using less stringent filter thresholds.")
      return(list())
    }

    reference_level <- .mStat_get_group_reference_level(
      meta.dat = meta_tab_level,
      group.var = group.var
    )
    group_is_categorical <- !is.null(reference_level)

    # Set the analysis object based on the method used.
    if (feature.dat.type == "other") {
      analysis.obj <- lm.results
    } else {
      analysis.obj <- linda.obj
    }

    # Calculate mean abundance and prevalence for each feature.
    prop_prev_data <- mStat_summarize_taxa_features(
      feature.dat = otu_tax_agg_filter,
      feature.level = feature.level
    )

    # Extract relevant data frames from analysis output.
    sub_test.list <- .mStat_extract_pair_linda_outputs(
      linda_output = analysis.obj$output,
      group_var = group.var,
      time_var = NULL,
      reference_level = if (group_is_categorical) reference_level else NULL
    )

    # Process each data frame in the list
    sub_test.list <- lapply(sub_test.list, function(df) {
      mStat_format_linda_feature_results(
        result.df = df,
        feature.level = feature.level,
        feature.stats = prop_prev_data,
        include_significant = TRUE
      )
    })

    return(sub_test.list)
  })

  # Name the elements of the test list with the feature levels
  names(test.list) <- feature.level

  return(test.list)
}
