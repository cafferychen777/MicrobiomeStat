#' Generate Taxa Test Pair
#'
#' This function, `generate_taxa_test_pair`, is designed for analysis in microbiome studies.
#' It takes a MicrobiomeStat data object as input and performs several key steps in the analysis of microbiome data.
#' The function filters taxa based on prevalence and abundance thresholds, aggregates taxon abundances by sample,
#' and applies the linda method to fit linear mixed effects models. These models are used to identify significant
#' taxon changes across different groups over time, taking into account various covariates.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A string that specifies the name of the subject variable column in the metadata.
#' @param time.var A string that specifies the name of the time variable column in the metadata. If not provided, it's NULL by default.
#' @param change.base A value indicating the base level for the time variable.
#' If provided, the specified level will be used as the reference category in
#' the model. Default is NULL, which means the first level of the factor will
#' be used.
#' @param group.var A string that specifies the name of the grouping variable column in the metadata for linear modelling.
#' @param adj.vars A vector of strings that specify the names of additional variables to be used as covariates in the analysis.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param ... Additional parameters to be passed to the linda function.
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
#' @return A named list containing data frames summarizing taxon test results for each taxonomic level.
#' @details Each list element corresponds to a taxonomic level specified in `feature.level`.
#' The data frame contains columns for taxon name, log2 fold change, p-values, adjusted p-values,
#' mean abundance, mean prevalence, and the output element from `linda` where the taxon was found significant.
#' @export
#' @name generate_taxa_test_pair
generate_taxa_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           change.base,
           group.var,
           adj.vars,
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