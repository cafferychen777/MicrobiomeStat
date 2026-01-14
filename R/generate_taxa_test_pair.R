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
#' @export
#' @name generate_taxa_test_pair
generate_taxa_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           change.base = NULL,
           group.var = NULL,
           ref.level = NULL,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           ...) {
    # Validate the input data object
    mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Extract relevant metadata
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

    # Function to generate the formula for the linear mixed model
    # This function constructs the fixed and random effects based on the provided variables
    generate_formula <- function(group.var = NULL, adj.vars = NULL, time.var = NULL, subject.var = NULL) {
      # Initialize variables for fixed and random effects
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine adjustment variables into a single string if they exist
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Construct the formula based on the presence or absence of time variable
      if (is.null(time.var)) {
        # For cross-sectional analysis (no time variable)
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model if no predictors
          }
        } else {
          # Include group variable and adjustment variables in fixed effects
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var)
          } else {
            fixed_effects <- group.var
          }
        }
        # Add random intercept for subject
        random_effects <- paste("(1 |", subject.var, ")")
      } else {
        # For longitudinal analysis (time variable present)
        if (is.null(group.var)) {
          # Time effect and adjustment variables
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          # Interaction between group and time, plus adjustment variables
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        # Random intercept and slope for subject
        random_effects <- paste("(1 +", time.var, "|", subject.var, ")")
      }

      # Combine fixed and random effects into full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    # Generate the formula for the current analysis
    formula <- generate_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = time.var,
                                subject.var = subject.var)

    # Function to generate a corrected formula with simplified random effects
    # This is used as a fallback if the original model fails to converge
    correct_formula <- function(group.var = NULL, adj.vars = NULL, time.var = NULL, subject.var = NULL) {
      # Initialize variables for fixed and random effects
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine adjustment variables into a single string if they exist
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Construct the formula based on the presence or absence of time variable
      if (is.null(time.var)) {
        # For cross-sectional analysis (no time variable)
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model if no predictors
          }
        } else {
          # Include group variable and adjustment variables in fixed effects
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var)
          } else {
            fixed_effects <- group.var
          }
        }
        # Add random intercept for subject
        random_effects <- paste("(1 |", subject.var, ")")
      } else {
        # For longitudinal analysis (time variable present)
        if (is.null(group.var)) {
          # Time effect and adjustment variables
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          # Interaction between group and time, plus adjustment variables
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        # Random intercept for subject
        random_effects <- paste("(1 |", subject.var, ")")
      }

      # Combine fixed and random effects into full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    # Generate the corrected formula
    formula_corrected <- correct_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = time.var,
                                subject.var = subject.var)

    # Perform analysis for each taxonomic level
    test.list <- lapply(feature.level, function(feature.level) {

      # Normalize count data if necessary
      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      }

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Get the appropriate feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Apply abundance and prevalence filters
      otu_tax_agg_filter <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter)

      if (feature.dat.type == "count"){
        feature.dat.type <- "proportion"
      }

      # Add this check before linda analysis
      if (any(colSums(otu_tax_agg_filter) == 0)) {
        keep_samples <- colSums(otu_tax_agg_filter) > 0
        otu_tax_agg_filter <- otu_tax_agg_filter[, keep_samples]
        meta_tab <- meta_tab[keep_samples, ]
      }

      # Perform linear mixed model analysis using LInDA
      # If the original model fails, a simpler model is used as a fallback
      linda.obj <- tryCatch({
        linda(feature.dat = otu_tax_agg_filter,
              meta.dat = meta_tab,
              formula = paste("~", formula),
              feature.dat.type = feature.dat.type,
              ...)
      }, error = function(e) {
        message("Error in linda: ", e)
        message("Due to the above error, a simpler model will be used for fitting.")
        linda(feature.dat = otu_tax_agg_filter,
              meta.dat = meta_tab,
              formula = paste("~", formula_corrected),
              feature.dat.type = feature.dat.type)
      })

      # Determine reference level for group comparisons
      if (!is.null(group.var)){
        reference_level <- levels(as.factor(meta_tab[,group.var]))[1]
      }

      # Calculate average abundance and prevalence for each feature
      prop_prev_data <-
        otu_tax_agg %>%
        as.matrix() %>%
        as.table() %>%
        as.data.frame() %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(
          avg_abundance = mean(Freq),
          prevalence = sum(Freq > 0) / dplyr::n()
        ) %>% column_to_rownames("Var1") %>%
        rownames_to_column(feature.level)

      # Function to extract relevant data frames from LInDA output
      extract_data_frames <- function(linda_object, group_var = NULL) {
        # Initialize an empty list to store the extracted dataframes
        result_list <- list()

        # Get all matching dataframe names
        matching_dfs <- grep(paste0(group_var), names(linda_object$output), value = TRUE)

        # Iteratively traverse all matching dataframe names and extract them
        for (df_name in matching_dfs) {
          # Extract group values from data frame names
          group_prefix <- paste0(group_var)

          # Extract the content after "group_prefix" and stop before the ":"
          group_value <- unlist(strsplit(df_name, split = ":"))[1]
          group_value <- gsub(pattern = group_prefix, replacement = "", x = group_value)

          if (grepl(pattern = time.var, x = df_name)){
            result_list[[paste0(group_value," vs ", reference_level, " (Reference) [", "Interaction", "]")]] <- linda_object$output[[df_name]]
          } else {
            # Add the data frame to the result list
            result_list[[paste0(group_value," vs ", reference_level, " (Reference) [", "Main Effect", "]")]] <- linda_object$output[[df_name]]
          }
        }

        return(result_list)
      }

      # Extract and process results from LInDA output
      sub_test.list <- extract_data_frames(linda_object = linda.obj, group_var = group.var)

      # Format and annotate results
      sub_test.list <- lapply(sub_test.list, function(df){
        df <- df %>%
          rownames_to_column(feature.level) %>%
          dplyr::left_join(prop_prev_data, by = feature.level) %>%
          dplyr::select(all_of(all_of(c(feature.level,"log2FoldChange","lfcSE","pvalue","padj","avg_abundance","prevalence")))) %>%
          dplyr::rename(Variable = feature.level,
                        Coefficient = log2FoldChange,
                        SE = lfcSE,
                        P.Value = pvalue,
                        Adjusted.P.Value = padj,
                        Mean.Abundance = avg_abundance,
                        Prevalence = prevalence)

        return(df)
      })

      return(sub_test.list)
    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)
  }
