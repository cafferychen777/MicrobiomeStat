#' Longitudinal Taxa Trend Test
#'
#' Conducts longitudinal trend tests to analyze how microbial taxa abundance changes
#' over time and across groups using linear mixed-effects models.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @param ref.level Character string specifying the reference level for categorical group.var.
#'   If NULL, the first level alphabetically is used. Ignored for continuous variables.
#' @param ... Additional arguments passed to downstream functions.
#' @details
#' This function requires \\code{time.var} and always evaluates temporal trends.
#' The fixed-effects structure is:
#' \\itemize{
#'   \\item \\code{time.var} when \\code{group.var = NULL}
#'   \\item \\code{group.var * time.var} when \\code{group.var} is provided
#'   \\item optional additive adjustment terms from \\code{adj.vars}
#' }
#' Random effects are modeled as \\code{(1 + time.var | subject.var)}.
#'
#' Output includes time main effect (when available), group main effects, and
#' group-time interaction terms when \\code{group.var} is specified.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested effects (e.g., "time", "group:time interaction")
#'         \itemize{
#'           \item For categorical \code{group.var}: Elements named as
#'                 "Level vs Reference (Reference) [Main Effect]" and
#'                 "Level vs Reference (Reference) [Interaction]"
#'           \item For continuous \code{group.var}: Elements named by model terms
#'                 (e.g., \code{group.var}, \code{group.var:time.var})
#'           \item \code{time.var} main effect is included when present in model output
#'         }
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Effect size. For time effects, represents change per unit time.
#'                 For group effects, represents log2 fold change. For interactions, represents
#'                 difference in trends between groups
#'           \item \code{SE}: Standard error of the coefficient
#'           \item \code{P.Value}: Raw p-value from the statistical test (mixed-effects model or LinDA)
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value using Benjamini-Hochberg method
#'           \item \code{Mean.Abundance}: Mean abundance of the feature across all samples
#'           \item \code{Prevalence}: Proportion of samples where the feature is present (non-zero)
#'         }
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1
#' data("ecam.obj")
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   adj.vars = "diet",
#'   feature.level = c("Phylum","Class"),
#'   feature.dat.type = c("proportion")
#' )
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   feature.level = c("Phylum","Class"),
#'   feature.dat.type = c("proportion")
#' )
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = NULL,
#'   feature.level = c("Phylum","Class"),
#'   feature.dat.type = c("proportion")
#' )
#'
#' # Example 2
#' data("subset_T2D.obj")
#' test.list <- generate_taxa_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus","Family"),
#'   feature.dat.type = c("count")
#' )
#'
#' plot.list <- generate_taxa_trend_volcano_long(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   test.list = test.list,
#'   feature.mt.method = "none")
#'
#' }
#' @noRd
.mStat_extract_trend_linda_outputs <- function(linda_output,
                                               group_var = NULL,
                                               time_var = NULL,
                                               reference_level = NULL) {
  result_list <- list()
  output_names <- names(linda_output)

  if (is.null(group_var)) {
    if (!is.null(time_var) && time_var %in% output_names) {
      result_list[[time_var]] <- linda_output[[time_var]]
    }
    return(result_list)
  }

  if (!is.null(time_var) && time_var %in% output_names) {
    result_list[[time_var]] <- linda_output[[time_var]]
  }

  matching_dfs <- grep(paste0("^", group_var), output_names, value = TRUE)
  for (df_name in matching_dfs) {
    term_head <- strsplit(df_name, split = ":", fixed = TRUE)[[1]][1]
    group_value <- sub(paste0("^", group_var), "", term_head)
    is_interaction <- !is.null(time_var) && grepl(paste0(":", time_var), df_name, fixed = TRUE)

    if (nzchar(group_value)) {
      base_label <- paste0(group_value, " vs ", reference_level, " (Reference)")
      label <- if (is_interaction) {
        paste0(base_label, " [Interaction]")
      } else {
        paste0(base_label, " [Main Effect]")
      }
    } else {
      label <- if (is_interaction) {
        paste0(group_var, ":", time_var)
      } else {
        group_var
      }
    }

    result_list[[label]] <- linda_output[[df_name]]
  }

  result_list
}

#' @export
generate_taxa_trend_test_long <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var = NULL,
           ref.level = NULL,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           ...) {
    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    if (is.null(time.var)) {
      stop("`time.var` is required for generate_taxa_trend_test_long.", call. = FALSE)
    }

    # Inform the user about the importance of numeric time variable
    message(
      "The trend test calculation relies on a numeric time variable.\n",
      "Please check that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the test results.\n",
      "You can ensure the time variable is numeric by mutating it in the metadata."
    )

    # Convert time variable to numeric
    data.obj$meta.dat[[time.var]] <- mStat_coerce_time_to_numeric(
      data.obj$meta.dat[[time.var]],
      time.var = time.var,
      context = "taxa trend analysis"
    )

    # Extract relevant variables from metadata
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

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

    # Build longitudinal formula
    formula <- .mStat_build_taxa_linda_formula(
      group.var = group.var,
      adj.vars = adj.vars,
      time.var = time.var,
      subject.var = subject.var,
      random_slopes = TRUE
    )

    analysis_data.obj <- data.obj
    analysis_feature.dat.type <- feature.dat.type
    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      analysis_data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      analysis_feature.dat.type <- "proportion"
    } else if (feature.dat.type == "other") {
      message(
        "CAUTION: You have selected 'other' as feature.dat.type. This assumes your data has been appropriately pre-processed.\n",
        "Please ensure your data:\n",
        "1. Does not contain zero values (required for CLR transformation)\n",
        "2. Has been properly normalized or transformed for compositional data analysis\n",
        "3. Is suitable for trend analysis across time points\n",
        "Non-standard preprocessing may affect result interpretation and comparability."
      )
    }

    # Perform analysis for each taxonomic level
    test.list <- lapply(feature.level, function(feature.level) {
      # Aggregate data to the specified taxonomic level if necessary
      otu_tax_agg_filter <- get_taxa_data(
        analysis_data.obj,
        feature.level,
        prev.filter,
        abund.filter,
        feature.col = FALSE
      )
      meta_tab_level <- meta_tab

      # Handle zero values in the data if using 'other' mode
      if (analysis_feature.dat.type == "other" && any(otu_tax_agg_filter == 0)) {
        warning(
          "Zero values detected in your data with feature.dat.type='other'.\n",
          "Applying half-minimum imputation for zeros to allow CLR transformation.\n",
          "For full control over zero handling, please pre-process your data before analysis."
        )
        
        # Apply half-minimum approach for zero values (similar to proportion handling in LinDA)
        otu_tax_agg_filter <- t(apply(otu_tax_agg_filter, 1, function(x) {
          if(any(x != 0)) { # Only process if row has some non-zero values
            x[x == 0] <- 0.5 * min(x[x != 0])
          } else {
            # If all values are zero (unlikely after filtering), add small pseudo-count
            x[x == 0] <- 1e-8
          }
          return(x)
        }))
      }

      # Add this check before linda analysis
      if (any(colSums(otu_tax_agg_filter) == 0)) {
        keep_samples <- colSums(otu_tax_agg_filter) > 0
        otu_tax_agg_filter <- otu_tax_agg_filter[, keep_samples]
        meta_tab_level <- meta_tab_level[keep_samples, , drop = FALSE]
      }

      # Perform LinDA (Linear models for Differential Abundance) analysis
      linda.obj <- linda(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab_level,
        formula = paste("~", formula),
        feature.dat.type = analysis_feature.dat.type,
        prev.filter = prev.filter,
        mean.abund.filter = abund.filter,
        ...
      )

      # Determine the reference level for the group variable
      if (!is.null(group.var)){
        reference_level <- levels(as.factor(meta_tab_level[,group.var]))[1]
      }

      # Calculate mean abundance and prevalence for each feature
      prop_prev_data <- mStat_summarize_taxa_features(
        feature.dat = otu_tax_agg_filter,
        feature.level = feature.level
      )

      # Extract relevant data frames from LinDA output
      sub_test.list <- .mStat_extract_trend_linda_outputs(
        linda_output = linda.obj$output,
        group_var = group.var,
        time_var = time.var,
        reference_level = if (!is.null(group.var)) reference_level else NULL
      )

      if (length(sub_test.list) == 0 && !is.null(time.var)) {
        warning(paste("No data frame found for time variable:", time.var))
      }

      # Process each data frame in the list
      sub_test.list <- lapply(sub_test.list, function(df){
        mStat_format_linda_feature_results(
          result.df = df,
          feature.level = feature.level,
          feature.stats = prop_prev_data
        )
      })

      return(sub_test.list)
    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)
  }
