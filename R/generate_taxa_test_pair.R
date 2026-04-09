#' Paired/Longitudinal Taxa Differential Abundance Test
#'
#' Performs differential abundance analysis for paired or longitudinal data using
#' linear mixed-effects models via LinDA, accounting for subject-level correlations.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @param change.base Value indicating the base/reference level for the time variable.
#'   If NULL, the first level is used.
#' @param ref.level Character string specifying the reference level for categorical group.var.
#'   If NULL, the first level alphabetically is used. Ignored for continuous variables.
#' @param ... Additional parameters passed to the linda function.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' test.list <- generate_taxa_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.dat.type = "count"
#' )
#' plot.list <-
#' generate_taxa_volcano_single(
#'  data.obj = peerj32.obj,
#'  group.var = "group",
#'  test.list = test.list,
#'  feature.sig.level = 0.1,
#'  feature.mt.method = "none"
#')
#'
#' data("subset_pairs.obj")
#' test.list <- generate_taxa_test_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.dat.type = "count"
#' )
#' plot.list <-
#' generate_taxa_volcano_single(
#'  data.obj = subset_pairs.obj,
#'  group.var = "Sex",
#'  test.list = test.list,
#'  feature.sig.level = 0.1,
#'  feature.mt.method = "none"
#')
#' }
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested comparisons and effects
#'         \itemize{
#'           \item For categorical \code{group.var}: Elements named as
#'                 "Level vs Reference (Reference) [Main Effect]" or "[Interaction]"
#'           \item For continuous \code{group.var}: Element named by the variable
#'           \item When \code{time.var} is provided: Both main effects and interaction effects
#'         }
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Log2 fold change for the comparison or time effect
#'           \item \code{SE}: Standard error of the coefficient
#'           \item \code{P.Value}: Raw p-value from LinDA's statistical test
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value using Benjamini-Hochberg method
#'           \item \code{Mean.Abundance}: Mean abundance of the feature across all samples
#'           \item \code{Prevalence}: Proportion of samples where the feature is present (non-zero)
#'         }
#' }
#'
#' @details This function performs paired/longitudinal differential abundance analysis.
#' Each list element corresponds to a taxonomic level specified in \code{feature.level}.
#' The function uses linear mixed-effects models via LinDA to account for subject-level
#' correlations in paired or longitudinal data.
#' @noRd
.mStat_build_taxa_linda_formula <- function(group.var = NULL,
                                            adj.vars = NULL,
                                            time.var = NULL,
                                            subject.var = NULL,
                                            random_slopes = TRUE) {
  adj.vars_str <- if (!is.null(adj.vars)) paste(adj.vars, collapse = " + ") else NULL

  if (is.null(time.var)) {
    if (is.null(group.var)) {
      fixed_effects <- if (!is.null(adj.vars_str)) adj.vars_str else "1"
    } else if (!is.null(adj.vars_str)) {
      fixed_effects <- paste(adj.vars_str, "+", group.var)
    } else {
      fixed_effects <- group.var
    }
    random_effects <- paste("(1 |", subject.var, ")")
  } else {
    if (is.null(group.var)) {
      fixed_effects <- if (!is.null(adj.vars_str)) paste(adj.vars_str, "+", time.var) else time.var
    } else if (!is.null(adj.vars_str)) {
      fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
    } else {
      fixed_effects <- paste(group.var, "*", time.var)
    }

    random_effects <- if (random_slopes) {
      paste("(1 +", time.var, "|", subject.var, ")")
    } else {
      paste("(1 |", subject.var, ")")
    }
  }

  paste(fixed_effects, random_effects, sep = " + ")
}

#' @noRd
.mStat_extract_pair_linda_outputs <- function(linda_output,
                                              group_var = NULL,
                                              time_var = NULL,
                                              reference_level = NULL) {
  result_list <- list()
  output_names <- names(linda_output)

  if (is.null(group_var)) {
    for (df_name in output_names) {
      result_list[[df_name]] <- linda_output[[df_name]]
    }
    return(result_list)
  }

  matching_dfs <- grep(paste0("^", group_var), output_names, value = TRUE)
  for (df_name in matching_dfs) {
    group_value <- strsplit(df_name, split = ":", fixed = TRUE)[[1]][1]
    group_value <- sub(paste0("^", group_var), "", group_value)

    base_label <- if (nzchar(group_value)) {
      paste0(group_value, " vs ", reference_level, " (Reference)")
    } else {
      group_var
    }

    is_interaction <- !is.null(time_var) && grepl(paste0(":", time_var), df_name, fixed = TRUE)
    label <- if (is_interaction) {
      paste0(base_label, " [Interaction]")
    } else {
      paste0(base_label, " [Main Effect]")
    }

    result_list[[label]] <- linda_output[[df_name]]
  }

  result_list
}

#' @export
#' @name generate_taxa_test_pair
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

    analysis_data.obj <- data.obj
    analysis_feature.dat.type <- feature.dat.type
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      analysis_data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      analysis_feature.dat.type <- "proportion"
    }

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

      # Add this check before linda analysis
      if (any(colSums(otu_tax_agg_filter) == 0)) {
        keep_samples <- colSums(otu_tax_agg_filter) > 0
        otu_tax_agg_filter <- otu_tax_agg_filter[, keep_samples]
        meta_tab_level <- meta_tab_level[keep_samples, , drop = FALSE]
      }

      # Perform linear mixed model analysis using LInDA
      # If the original model fails, a simpler model is used as a fallback
      linda.obj <- tryCatch({
        do.call(
          linda,
          c(
            list(
              feature.dat = otu_tax_agg_filter,
              meta.dat = meta_tab_level,
              formula = paste("~", formula),
              feature.dat.type = analysis_feature.dat.type,
              p.adj.method = p.adj.method,
              alpha = feature.sig.level
            ),
            extra_args
          )
        )
      }, error = function(e) {
        message("Error in linda: ", e)
        message("Due to the above error, a simpler model will be used for fitting.")
        do.call(
          linda,
          c(
            list(
              feature.dat = otu_tax_agg_filter,
              meta.dat = meta_tab_level,
              formula = paste("~", formula_corrected),
              feature.dat.type = analysis_feature.dat.type,
              p.adj.method = p.adj.method,
              alpha = feature.sig.level
            ),
            extra_args
          )
        )
      })

      # Determine reference level for group comparisons
      if (!is.null(group.var)){
        reference_level <- levels(as.factor(meta_tab_level[,group.var]))[1]
      }

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
