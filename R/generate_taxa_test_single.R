#' Helper function to perform linear model analysis for "other" data type
#'
#' @param feature.dat Feature data matrix
#' @param meta.dat Metadata data frame
#' @param formula Formula string for the model
#' @param group.var Group variable name
#' @return List with analysis results
#' @noRd
perform_lm_analysis <- function(feature.dat, meta.dat, formula, group.var) {
  # Check if group.var is categorical (factor or character)
  is_categorical <- is.factor(meta.dat[, group.var]) || is.character(meta.dat[, group.var])

  # Get reference level only if categorical
  if (is_categorical) {
    reference_level <- levels(as.factor(meta.dat[, group.var]))[1]
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
      group_coefs <- coef_table[grep(paste0("^", group.var), rownames(coef_table)), , drop = FALSE]

      if (nrow(group_coefs) > 0) {
        for (j in 1:nrow(group_coefs)) {
          coef_name <- rownames(group_coefs)[j]

          if (is_categorical) {
            # For categorical variables, extract group value from coefficient name
            group_value <- gsub(paste0("^", group.var), "", coef_name)
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
      df$padj <- p.adjust(df$pvalue, method = "BH")
      df$reject <- df$padj <= 0.05
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
                                      ...) {
  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Validate the input data object
  mStat_validate_data(data.obj)

  # Subset the data if a specific time point is specified
  if (!is.null(time.var)) {
    if (!is.null(t.level)) {
      # Create a condition string for subsetting
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      # Subset the data object based on the condition
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
  }

  # Extract relevant variables from the metadata
  meta_tab <-
    data.obj$meta.dat %>% select(all_of(c(time.var, group.var, adj.vars)))

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

  # For "other" data type, check if data contains negative values
  # If so, adjust abundance filter to handle negative values appropriately
  if (feature.dat.type == "other") {
    # Check if any feature table contains negative values
    has_negative <- FALSE
    if (!is.null(data.obj$feature.tab)) {
      has_negative <- any(data.obj$feature.tab < 0, na.rm = TRUE)
    }
    if (!has_negative && !is.null(data.obj$feature.agg.list)) {
      for (agg_table in data.obj$feature.agg.list) {
        if (any(agg_table < 0, na.rm = TRUE)) {
          has_negative <- TRUE
          break
        }
      }
    }

    if (has_negative) {
      message("Note: Negative values detected in 'other' data type. Abundance filtering is disabled to preserve all features.")
      abund.filter <- -Inf  # Set to negative infinity to include all features regardless of abundance
    }
  }

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
    # Aggregate data to the specified taxonomic level if necessary
    if (is.null(data.obj$feature.agg.list[[feature.level]]) &
        feature.level != "original") {
      data.obj <-
        mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    # Extract the appropriate feature table
    if (feature.level != "original") {
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- data.obj$feature.tab
    }

    # Apply prevalence and abundance filters
    otu_tax_agg_filter <-  otu_tax_agg %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter)

    # Add a check before running linda
    if (nrow(otu_tax_agg_filter) == 0 || ncol(otu_tax_agg_filter) == 0) {
      warning("No features remain after filtering. Consider using less stringent filter thresholds.")
      return(list())
    }

    # Add this check before linda analysis
    if (any(colSums(otu_tax_agg_filter) == 0)) {
      keep_samples <- colSums(otu_tax_agg_filter) > 0
      if (sum(keep_samples) == 0) {
        warning("No samples remain after filtering.")
        return(list())
      }
      otu_tax_agg_filter <- otu_tax_agg_filter[, keep_samples]
      allvars <- names(meta_tab)
      meta_tab <- as.data.frame(meta_tab[keep_samples, ])
      names(meta_tab) <- allvars
    }

    # Choose analysis method based on feature data type
    if (feature.dat.type == "other") {
      # Use linear models for "other" data type (e.g., log-transformed data)
      lm.results <- perform_lm_analysis(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab,
        formula = formula,
        group.var = group.var
      )
    } else {
      # Perform LinDA (Linear models for Differential Abundance) analysis
      # Muffle linda's own "all filtered" warning; we emit our own below
      linda.obj <- withCallingHandlers(
        linda(
          feature.dat = otu_tax_agg_filter,
          meta.dat = meta_tab,
          formula = paste("~", formula),
          feature.dat.type = feature.dat.type,
          prev.filter = prev.filter,
          mean.abund.filter = abund.filter,
          ...
        ),
        warning = function(w) {
          if (grepl("All features were filtered out", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
    }

    # Check if linda returned empty results (all features filtered internally)
    if (feature.dat.type != "other" && length(linda.obj$output) == 0) {
      warning("No features remain after filtering. Consider using less stringent filter thresholds.")
      return(list())
    }

    # Determine the reference level for the group variable (only for categorical variables)
    if (!is.null(group.var)) {
      if (is.factor(meta_tab[, group.var]) || is.character(meta_tab[, group.var])) {
        # Only get reference level for categorical variables
        reference_level <- levels(as.factor(meta_tab[, group.var]))[1]
      } else {
        # For continuous variables, no reference level
        reference_level <- NULL
      }
    }

    # Set the analysis object based on the method used
    if (feature.dat.type == "other") {
      analysis.obj <- lm.results
    } else {
      analysis.obj <- linda.obj
    }

    # Calculate mean abundance and prevalence for each feature
    prop_prev_data <-
      otu_tax_agg %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::group_by(Var1) %>%
      dplyr::summarise(avg_abundance = mean(Freq),
                       prevalence = sum(Freq != 0) / dplyr::n(),
                       .groups = "drop") %>%
      column_to_rownames("Var1") %>%
      rownames_to_column(feature.level)

    # Function to extract relevant data frames from LinDA output
    extract_data_frames <-
      function(linda_object, group_var = NULL, reference_level = NULL) {
        result_list <- list()

        # Find data frames related to the group variable
        matching_dfs <-
          grep(paste0(group_var), names(linda_object$output), value = TRUE)

        for (df_name in matching_dfs) {
          group_prefix <- paste0(group_var)

          # Extract the group value from the data frame name
          group_value <- unlist(strsplit(df_name, split = ":"))[1]
          group_value <-
            gsub(pattern = group_prefix,
                 replacement = "",
                 x = group_value)

          # Store the data frame with a descriptive name
          if (!is.null(reference_level) && reference_level != "") {
            # For categorical variables, show comparison vs reference
            result_list[[paste0(group_value, " vs ", reference_level, " (Reference)")]] <-
              linda_object$output[[df_name]]
          } else {
            # For continuous variables, just use the variable name
            result_list[[group_var]] <- linda_object$output[[df_name]]
          }
        }

        return(result_list)
      }

    # Extract relevant data frames from analysis output
    sub_test.list <-
      extract_data_frames(linda_object = analysis.obj, group_var = group.var, reference_level = reference_level)

    # Process each data frame in the list
    sub_test.list <- lapply(sub_test.list, function(df) {
      df <- df %>%
        # Add feature level as a column
        rownames_to_column(feature.level) %>%
        # Join with prevalence and abundance data
        dplyr::left_join(prop_prev_data, by = feature.level) %>%
        # Select relevant columns
        dplyr::select(all_of(all_of(
          c(
            feature.level,
            "log2FoldChange",
            "lfcSE",
            "pvalue",
            "padj",
            "avg_abundance",
            "prevalence"
          )
        ))) %>%
        # Rename columns for clarity
        dplyr::rename(
          Variable = feature.level,
          Coefficient = log2FoldChange,
          SE = lfcSE,
          P.Value = pvalue,
          Adjusted.P.Value = padj,
          Mean.Abundance = avg_abundance,
          Prevalence = prevalence
        )

      return(df)
    })

    return(sub_test.list)
  })

  # Name the elements of the test list with the feature levels
  names(test.list) <- feature.level

  return(test.list)
}
