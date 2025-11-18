#' Longitudinal Taxa Association Test Generation
#'
#' This function performs association testing between taxa abundances and a grouping variable in longitudinal microbiome data. It discerns how taxa abundances differ between experimental or observational groups over time.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string that indicates the column name in the metadata which uniquely identifies each subject or sample.
#' @param group.var A character string specifying the grouping variable column in the metadata.
#'                  Can be either:
#'                  \itemize{
#'                    \item Categorical (factor or character): Creates pairwise comparisons between groups
#'                    \item Continuous (numeric or integer): Tests for linear association
#'                  }
#'                  This variable differentiates between different experimental or observational groups.
#' @param adj.vars A vector of character strings. Each string should denote a column name in the metadata that will serve as a covariate in the analysis. These variables might account for potential confounding influences. Default is NULL.
#' @param feature.level A character string indicating the taxonomic resolution for analysis (e.g., "Phylum", "Class"). This choice will determine the granularity of the analysis.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled.
#' Should be one of:
#' \itemize{
#'   \item "count": Raw count data. This function will first apply TSS (Total Sum Scaling) normalization,
#'         then LinDA performs zero-handling using half-minimum approach for statistical testing
#'   \item "proportion": Pre-normalized proportional data. LinDA performs zero-handling using
#'                       half-minimum approach without additional normalization
#' }
#' Default is "count".
#' @param ... Additional arguments to cater to any specialized requirements. For now, these are placeholder and not used.
#' @details
#' Based on whether group.var and adj.vars are NULL, the formula tests:
#'
#' - When group.var is NULL and adj.vars is NOT NULL:
#'   - Tests adj.vars main effects only.
#'   - Adjusted for adj.vars but not group.var.
#'
#' - When group.var is NOT NULL and adj.vars is NOT NULL:
#'   - Tests adj.vars and group.var main effects.
#'   - Adjusted for adj.vars.
#'
#' - When group.var is NOT NULL and adj.vars is NULL:
#'   - Tests group.var main effects only.
#'   - Unadjusted analysis.
#'
#' - When both group.var and adj.vars are NULL:
#'   - Tests the intercept only.
#'   - Unadjusted analysis.
#'
#' The formula combines the appropriate terms based on which variables are NULL.
#' Subject variability is accounted for through random effects.
#'
#' When group.var and adj.vars are NULL, the intercept is tested without adjusting for any covariates.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested comparisons between groups
#'         (e.g., "Level vs Reference (Reference)" for categorical variables,
#'         or variable name for continuous variables)
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Log2 fold change (categorical) or slope (continuous)
#'           \item \code{SE}: Standard error of the coefficient
#'           \item \code{P.Value}: Raw p-value from LinDA's statistical test
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value (Benjamini-Hochberg)
#'           \item \code{Mean.Abundance}: Mean abundance across all samples
#'           \item \code{Prevalence}: Proportion of samples where feature is present (non-zero)
#'         }
#' }
#'
#' Analysis uses LinDA mixed-effects models to test associations between taxa abundances
#' and the grouping variable, accounting for subject-level random effects.
#'
#' @examples
#' \dontrun{
#' # Example 1: Generate taxa association tests and volcano plots for the ecam dataset
#' data("ecam.obj")
#' test.list <- generate_taxa_association_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   group.var = "delivery",
#'   feature.level = c("Phylum", "Class"),
#'   feature.dat.type = c("count")
#' )
#'
#' volcano_plots_ecam <- generate_taxa_volcano_single(
#'   data.obj = ecam.obj,
#'   group.var = "delivery",
#'   test.list = test.list,
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "fdr"
#' )
#'
#'
#' # Example 2: Generate taxa association tests and volcano plots for a subset of the T2D dataset
#' data("subset_T2D.obj")
#' test.list_T2D <- generate_taxa_association_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   feature.level = "Genus",
#'   group.var = "subject_race",
#'   feature.dat.type = c("count"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.001
#' )
#'
#' volcano_plots_T2D <- generate_taxa_volcano_single(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   test.list = test.list_T2D,
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "none"
#' )
#' }
#' @export
generate_taxa_association_test_long <-
  function(data.obj,
           subject.var,
           group.var = NULL,
           adj.vars = NULL,
           prev.filter = 0,
           abund.filter = 0,
           feature.level,
           feature.dat.type = c("count", "proportion"),
           ...) {
    # Validate the input data object
    mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Extract relevant metadata and ensure it remains a data frame
    meta_tab <- data.obj$meta.dat %>% 
        as.data.frame() %>%
        select(all_of(c(group.var, adj.vars, subject.var)))

    # Define a function to generate the formula for statistical modeling
    # This function creates a formula string based on the provided variables
    generate_formula <- function(group.var=NULL, adj.vars=NULL, time.var=NULL, subject.var="Subject") {

      # Initialize variables for fixed and random effects
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine multiple adjustment variables into a single string
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Generate the formula based on the presence or absence of time variable
      if (is.null(time.var)) {
        # For cross-sectional data (no time variable)
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model if no variables are provided
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
        # For longitudinal data (time variable present)
        if (is.null(group.var)) {
          # Include time and adjustment variables
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          # Include interaction between group and time, along with adjustment variables
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        # Add random slope and intercept for subject
        random_effects <- paste("(1 +", time.var, "|", subject.var, ")")
      }

      # Combine fixed and random effects into a complete formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    # Generate the formula for the current analysis
    formula <- generate_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = NULL,
                                subject.var = subject.var)

    # Perform analysis for each feature level
    test.list <- lapply(feature.level, function(feature.level) {

      # Normalize count data if necessary
      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      }

      # Aggregate features by taxonomy if not already done
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Select the appropriate feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Filter features based on prevalence and abundance
      otu_tax_agg_filter <- otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter)

      # Convert filtered data back to matrix
      otu_tax_agg_filter <- as.matrix(otu_tax_agg_filter)

      # Add this check before linda analysis
      if (any(colSums(otu_tax_agg_filter) == 0)) {
        keep_samples <- colSums(otu_tax_agg_filter) > 0
        otu_tax_agg_filter <- otu_tax_agg_filter[, keep_samples]
        meta_tab <- meta_tab[keep_samples, , drop = FALSE]
      }

      # Perform linear mixed model analysis using LInDA
      linda.obj <- linda(
        feature.dat = otu_tax_agg_filter,
        meta.dat = meta_tab,
        formula = paste("~", formula),
        feature.dat.type = "proportion",
        prev.filter = 0,
        mean.abund.filter = 0,
        ...
      )

      # Identify the reference level for the group variable
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
        matching_dfs <- grep(group_var, names(linda_object$output), value = TRUE)

        # Iterate through all matching dataframe names and extract them
        for (df_name in matching_dfs) {
          # Extract group values from the data frame name
          group_prefix <- paste0(group_var)

          # Extract the content after the "group_prefix" and stop before the ":".
          group_value <- unlist(strsplit(df_name, split = ":"))[1]
          group_value <- gsub(pattern = group_prefix, replacement = "", x = group_value)

          # Add the data frame to the result list with appropriate naming
          result_list[[paste0(group_value," vs ", reference_level, " (Reference)")]] <- linda_object$output[[df_name]]
        }

        return(result_list)
      }

      # Extract data frames from LInDA output
      sub_test.list <- extract_data_frames(linda_object = linda.obj, group_var = group.var)

      # Process and format the extracted data frames
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

    # Assign names to the elements of test.list based on feature levels
    names(test.list) <- feature.level

    return(test.list)
  }