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
#' Based on whether group.var, adj.vars, and time.var are NULL, the formula tests:
#'
#' - When time.var is NULL:
#'   - When group.var is NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars main effects only.
#'     - Adjusted for adj.vars but not for group.var or time.var.
#'   - When group.var is NOT NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars and group.var main effects.
#'     - Adjusted for adj.vars but not for time.var.
#'   - When group.var is NOT NULL and adj.vars is NULL:
#'     - Tests group.var main effects only.
#'     - Unadjusted for adj.vars but adjusted for group.var.
#'   - When both group.var and adj.vars are NULL:
#'     - Tests the intercept only.
#'     - Unadjusted analysis.
#'
#' - When time.var is NOT NULL:
#'   - When group.var is NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars and time.var main effects.
#'     - Adjusted for adj.vars but not for group.var.
#'   - When group.var is NOT NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars main effects.
#'     - Tests group.var and time.var main effects.
#'     - Tests group.var x time.var interaction.
#'     - Adjusted for adj.vars.
#'   - When group.var is NOT NULL and adj.vars is NULL:
#'     - Tests group.var and time.var main effects.
#'     - Tests group.var x time.var interaction.
#'     - Unadjusted analysis.
#'   - When both group.var and adj.vars are NULL:
#'     - Tests time.var main effect only.
#'     - Unadjusted analysis.
#'
#' The formula combines the appropriate terms based on which variables are NULL.
#' Subject variability is accounted for through random effects.
#'
#' When group.var = NULL and adj.vars = NULL and time.var is NOT NULL,
#' the slope of time.var is tested without adjusting for any additional covariates.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested effects (e.g., "time", "group:time interaction")
#'         \itemize{
#'           \item For categorical \code{group.var}: Separate elements for each non-reference level
#'                 (e.g., "Treatment vs Control (Reference)")
#'           \item For continuous \code{group.var}: Single element named by the variable
#'           \item Time effects and interactions are also included depending on the model
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
    mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Inform the user about the importance of numeric time variable
    message(
      "The trend test calculation relies on a numeric time variable.\n",
      "Please check that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the test results.\n",
      "You can ensure the time variable is numeric by mutating it in the metadata."
    )

    # Convert time variable to numeric if it exists
    if (!is.null(time.var)){
      data.obj$meta.dat <-
        data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))
    }

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

    # Function to generate the formula for linear mixed-effects model
    generate_formula <- function(group.var = NULL, adj.vars = NULL, time.var = NULL, subject.var = NULL) {
      # Initialize fixed and random effects
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine adjustment variables
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Generate formula based on available variables
      if (is.null(time.var)) {
        # Case: No time variable
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var)
          } else {
            fixed_effects <- group.var
          }
        }
        random_effects <- paste("(1 |", subject.var, ")")
      } else {
        # Case: Time variable present
        if (is.null(group.var)) {
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        random_effects <- paste("(1 +", time.var, "|", subject.var, ")")
      }

      # Combine fixed and random effects into full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    # Generate the formula for the linear mixed-effects model
    formula <- generate_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = time.var,
                                subject.var = subject.var)

    # Perform analysis for each taxonomic level
    test.list <- lapply(feature.level, function(feature.level) {
      # Handle data based on feature.dat.type
      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
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

      # Aggregate data to the specified taxonomic level if necessary
      otu_tax_agg_filter <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter, feature.col = FALSE)

      # Set feature data type to proportion if it was originally count
      if (feature.dat.type == "count"){
        feature.dat.type = "proportion"
      }
      
      # Handle zero values in the data if using 'other' mode
      if (feature.dat.type == "other" && any(otu_tax_agg_filter == 0)) {
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
        meta_tab <- meta_tab[keep_samples, ]
      }

      # Perform LinDA (Linear models for Differential Abundance) analysis
      linda.obj <- linda(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab,
        formula = paste("~", formula),
        feature.dat.type = feature.dat.type,
        prev.filter = prev.filter,
        mean.abund.filter = abund.filter,
        ...
      )

      # Determine the reference level for the group variable
      if (!is.null(group.var)){
        reference_level <- levels(as.factor(meta_tab[,group.var]))[1]
      }

      # Calculate mean abundance and prevalence for each feature
      prop_prev_data <-
        otu_tax_agg_filter %>%
        as.matrix() %>%
        as.table() %>%
        as.data.frame() %>%
        dplyr::group_by(Var1) %>%  # Var1 represents taxa
        dplyr::summarise(
          avg_abundance = mean(Freq),
          prevalence = sum(Freq > 0) / dplyr::n()
        ) %>% column_to_rownames("Var1") %>%
        rownames_to_column(feature.level)

      # Function to extract relevant data frames from LinDA output
      extract_data_frames <- function(linda_object, group_var = NULL, time_var) {
        result_list <- list()

        if (!is.null(group_var)) {
          # Extract data frames for group comparisons
          matching_dfs <- grep(paste0(group_var, ".*:", time_var), names(linda_object$output), value = TRUE)

          for (df_name in matching_dfs) {
            group_prefix <- paste0(group_var)
            group_value <- unlist(strsplit(df_name, split = ":"))[1]
            group_value <- gsub(pattern = group_prefix, replacement = "", x = group_value)

            result_list[[paste0(group_value," vs ", reference_level, " (Reference)")]] <- linda_object$output[[df_name]]
          }
        } else {
          # Extract data frame for time effect when no group variable is present
          df <- linda_object$output[[time_var]]
          if (!is.null(df)) {
            result_list[[time_var]] <- df
          } else {
            warning(paste("No data frame found with the name:", time_var))
          }
        }

        return(result_list)
      }

      # Extract relevant data frames from LinDA output
      sub_test.list <- extract_data_frames(linda_object = linda.obj, group_var = group.var, time_var = time.var)

      # Process each data frame in the list
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
