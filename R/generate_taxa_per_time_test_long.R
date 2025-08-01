#' Longitudinal Differential Abundance Test in Microbiome Data
#'
#' This function performs differential abundance testing across multiple time points in longitudinal microbiome data. It is tailored to analyze how the abundance of microbial taxa varies over time, within different groups or under various conditions.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param subject.var A string specifying the column name in meta.dat that uniquely identifies each subject.
#' @param time.var Optional; a string representing the time variable in the meta.dat. If provided, enables longitudinal analysis.
#' @param group.var Optional; a string specifying the group variable in meta.dat for between-group comparisons.
#' @param adj.vars Optional; a vector of strings representing covariates in meta.dat for adjustment in the analysis.
#' @param feature.level A string or vector of strings indicating the taxonomic level(s) for analysis (e.g., "Phylum", "Class").
#' @param prev.filter Numeric; a minimum prevalence threshold for taxa inclusion in the analysis.
#' @param abund.filter Numeric; a minimum abundance threshold for taxa inclusion in the analysis.
#' @param feature.dat.type Character; "count" or "proportion", indicating the type of feature data.
#' @param ... Additional arguments passed to other methods.
#' @details
#' The function integrates various data manipulations, normalization procedures, and statistical tests to assess the significance of taxa changes over time or between groups. It allows for the adjustment of covariates and handles both count and proportion data types.
#'
#' The function constructs a mixed-effects model formula based on the provided variables, handling fixed and random effects to account for repeated measures in subjects. It performs filtering based on prevalence and abundance thresholds and applies normalization and aggregation procedures as necessary.
#'
#' Importantly, the function conducts differential abundance analysis separately for each time point in the longitudinal data. This approach allows for the identification of taxa that show significant changes at specific time points, providing insights into the dynamics of the microbiome over time.
#'
#' @return
#' A nested list structure. The top level of the list corresponds to different time points, and each element contains a list of dataframes for each taxonomic level. Each dataframe provides statistical analysis results for taxa at that level and time point.
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
           group.var = NULL,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion"),
           ...) {
    # Validate the input data object to ensure it meets the required format
    mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Define a function to generate the formula for the linear mixed model
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
                                time.var = NULL,
                                subject.var = subject.var)

    # Extract unique time levels from the data
    time.levels <- data.obj$meta.dat %>% select(all_of(c(time.var))) %>% unique() %>% pull()

    # Perform analysis for each time point
    test.list <- lapply(time.levels, function(t.level){
      tryCatch({
        # Subset data for the current time point
        subset.ids <- rownames(data.obj$meta.dat %>%
                                 filter(!!sym(time.var) %in% c(t.level)))
        
        subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)
        
        # Extract relevant metadata for the current subset
        meta_tab <-
          subset_data.obj$meta.dat %>% select(all_of(c(
            time.var, group.var, adj.vars, subject.var
          )))
        
        # Check if we have enough data points for mixed effects model
        # Count unique subjects and observations
        n_subjects <- length(unique(meta_tab[[subject.var]]))
        n_observations <- nrow(meta_tab)
        
        if (n_subjects >= n_observations) {
          warning(
            "At time point ", t.level, ", the number of unique subjects (", n_subjects, ") ",
            "is greater than or equal to the number of observations (", n_observations, "). ",
            "This makes it impossible to fit a mixed effects model. ",
            "Skipping this time point."
          )
          return(NULL) # Return NULL for this time point
        }
        
        # Perform analysis for each taxonomic level
        test.list <- lapply(feature.level, function(feature.level) {

        # Normalize count data if necessary
        if (feature.dat.type == "count"){
          message(
            "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
          )
          subset_data.obj <- mStat_normalize_data(subset_data.obj, method = "TSS")$data.obj.norm
        }

        # Aggregate data by taxonomy if necessary
        if (is.null(subset_data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
          subset_data.obj <- mStat_aggregate_by_taxonomy(data.obj = subset_data.obj, feature.level = feature.level)
        }

        # Select appropriate feature table
        if (feature.level != "original"){
          otu_tax_agg <- subset_data.obj$feature.agg.list[[feature.level]]
        } else {
          otu_tax_agg <- subset_data.obj$feature.tab
        }

        # Apply prevalence and abundance filters
        otu_tax_agg_filter <-  otu_tax_agg %>%
          as.data.frame() %>%
          mStat_filter(prev.filter = prev.filter,
                       abund.filter = abund.filter)

        if (feature.dat.type == "count"){
          feature.dat.type = "proportion"
        }

        # Perform linear mixed model analysis using LInDA
        linda.obj <- linda(
          feature.dat = otu_tax_agg_filter,
          meta.dat = meta_tab,
          formula = paste("~", formula),
          feature.dat.type = feature.dat.type,
          prev.filter = prev.filter,
          mean.abund.filter = abund.filter,
          ...
        )

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

          # Get all matching data frame names
          matching_dfs <- grep(group_var, names(linda_object$output), value = TRUE)

          # Iterate through all matching dataframe names and extract them
          for (df_name in matching_dfs) {
            # Extract group values from the data frame name
            group_prefix <- paste0(group_var)

            # Extract the content after "group_prefix", and stop before the ":"
            group_value <- unlist(strsplit(df_name, split = ":"))[1]
            group_value <- gsub(pattern = group_prefix, replacement = "", x = group_value)

            # Add the data frame to the result list
            result_list[[paste0(group_value," vs ", reference_level, " (Reference)")]] <- linda_object$output[[df_name]]
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
      }, error = function(e) {
        warning(
          "Error analyzing time point ", t.level, ": ", conditionMessage(e), "\n",
          "Skipping this time point and continuing with others."
        )
        return(NULL) # Return NULL for this time point if an error occurs
      })
    })

    # FIXED: Record which time points have valid results BEFORE filtering
    # This preserves the correct mapping to time.levels indices
    valid_time_indices <- which(!sapply(test.list, is.null))
    
    # Remove NULL entries from the list (time points that were skipped)
    test.list <- Filter(Negate(is.null), test.list)
    
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