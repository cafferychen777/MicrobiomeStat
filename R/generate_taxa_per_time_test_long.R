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

    build_per_time_formula <- function(group.var = NULL, adj.vars = NULL) {
      terms <- c(group.var, adj.vars)
      terms <- terms[!is.null(terms)]

      if (length(terms) == 0) {
        return("1")
      }

      paste(terms, collapse = " + ")
    }

    # Generate the fixed-effects formula for each per-time analysis
    formula <- build_per_time_formula(group.var = group.var, adj.vars = adj.vars)

    mStat_validate_time_var_contract(
      meta.dat = data.obj$meta.dat,
      time.var = time.var,
      context = "taxa per-time testing"
    )

    # Extract unique time levels from the data
    time.levels <- data.obj$meta.dat %>% select(all_of(c(time.var))) %>% unique() %>% pull()

    # Perform analysis for each time point
    test.list <- lapply(time.levels, function(t.level){
      tryCatch({
        # Subset data for the current time point
        subset.ids <- get_sample_ids(data.obj, time.var, t.level)
        subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)
        
        # Extract relevant metadata for the current subset
        meta_tab <-
          subset_data.obj$meta.dat %>% select(all_of(c(
            time.var, group.var, adj.vars, subject.var
          )))

        mStat_validate_group_var_contract(
          meta.dat = meta_tab,
          group.var = group.var,
          context = paste0("taxa per-time testing at ", t.level)
        )

        analysis_data.obj <- subset_data.obj
        analysis_feature.dat.type <- feature.dat.type
        if (feature.dat.type == "count"){
          message(
            "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
          )
          analysis_data.obj <- mStat_normalize_data(subset_data.obj, method = "TSS")$data.obj.norm
          analysis_feature.dat.type <- "proportion"
        }

        # Perform analysis for each taxonomic level
        test.list <- lapply(feature.level, function(feature.level) {

        # Aggregate, extract, and filter feature data
        otu_tax_agg_filter <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter, feature.col = FALSE)

        # Perform per-time fixed-effects analysis using LinDA
        linda.obj <- linda(
          feature.dat = otu_tax_agg_filter,
          meta.dat = meta_tab,
          formula = paste("~", formula),
          feature.dat.type = analysis_feature.dat.type,
          prev.filter = prev.filter,
          mean.abund.filter = abund.filter,
          ...
        )

        group_is_categorical <- is.factor(meta_tab[[group.var]]) || is.character(meta_tab[[group.var]])
        reference_level <- if (group_is_categorical) {
          levels(as.factor(meta_tab[[group.var]]))[1]
        } else {
          NULL
        }

        # Calculate average abundance and prevalence for each feature
        prop_prev_data <- mStat_summarize_taxa_features(
          feature.dat = otu_tax_agg_filter,
          feature.level = feature.level
        )

        # Function to extract relevant data frames from LInDA output
        extract_data_frames <- function(linda_object,
                                        group_var = NULL,
                                        reference_level = NULL,
                                        is_categorical = TRUE) {
          result_list <- list()

          if (is.null(group_var)) {
            return(result_list)
          }

          output_names <- names(linda_object$output)
          matching_dfs <- grep(paste0("^", group_var), output_names, value = TRUE)

          for (df_name in matching_dfs) {
            group_value <- strsplit(df_name, split = ":", fixed = TRUE)[[1]][1]
            group_value <- sub(paste0("^", group_var), "", group_value)

            label <- if (is_categorical && nzchar(group_value)) {
              paste0(group_value, " vs ", reference_level, " (Reference)")
            } else {
              group_var
            }

            result_list[[label]] <- linda_object$output[[df_name]]
          }

          result_list
        }

        # Extract and process results from LInDA output
        sub_test.list <- extract_data_frames(
          linda_object = linda.obj,
          group_var = group.var,
          reference_level = reference_level,
          is_categorical = group_is_categorical
        )

        if (length(sub_test.list) == 0) {
          return(list())
        }

        # Format and annotate results
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
      }, error = function(e) {
        warning(
          "Error analyzing time point ", t.level, ": ", conditionMessage(e), "\n",
          "Skipping this time point and continuing with others."
        )
        return(NULL) # Return NULL for this time point if an error occurs
      })
    })

    # Keep time points that produced at least one non-empty feature-level result
    has_results <- vapply(test.list, function(tp_result) {
      if (is.null(tp_result)) {
        return(FALSE)
      }
      any(vapply(tp_result, length, integer(1)) > 0)
    }, logical(1))

    valid_time_indices <- which(has_results)

    # Remove empty/failed time points from the list
    test.list <- test.list[has_results]
    
    # If all time points were skipped, return a message
    if (length(test.list) == 0) {
      stop("No time points could be analyzed. Check if you have enough observations per subject at each time point.")
    }
    
    # FIXED: Use the indices recorded from the original list to assign correct time point names
    # This ensures that the names correspond to the actual time points that were analyzed
    if (length(valid_time_indices) > 0) {
      names(test.list) <- time.levels[valid_time_indices]
    }

    return(test.list)
  }
