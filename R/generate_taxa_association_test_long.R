#' Longitudinal Taxa Association Test Generation
#'
#' This function performs association testing between taxa abundances and a grouping variable in longitudinal microbiome data. It discerns how taxa abundances differ between experimental or observational groups over time.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string that indicates the column name in the metadata which uniquely identifies each subject or sample.
#' @param group.var A character string specifying the grouping variable column in the metadata. This variable differentiates between different experimental or observational groups.
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
#' @param feature.dat.type A character string, either "count" or "proportion", indicating the nature of the data in the `data.obj`. This helps the function to determine if normalization is required. Default is "count".
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
#' @return
#' A list of dataframes, with each dataframe representing a specific taxonomic level (as specified in `feature.level`). These dataframes contain essential statistics, including taxa changes, p-values, and other metrics derived from the linear model.
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
    # Extract data
    mStat_validate_data(data.obj)

    feature.dat.type <- match.arg(feature.dat.type)

    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        group.var, adj.vars, subject.var
      )))

    # Function to generate the formula
    generate_formula <- function(group.var=NULL, adj.vars=NULL, time.var=NULL, subject.var="Subject") {

      # Initialize the fixed_effects and random_effects variables
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine multiple adj.vars into a single string if they exist
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Case where time.var is NULL
      if (is.null(time.var)) {
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

        # Case where time.var is NOT NULL
      } else {
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

      # Generate the full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    formula <- generate_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = NULL,
                                subject.var = subject.var)

    test.list <- lapply(feature.level, function(feature.level) {

      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      }

      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter)

      linda.obj <- linda(
        feature.dat = otu_tax_agg,
        meta.dat = meta_tab,
        formula = paste("~", formula),
        feature.dat.type = "proportion",
        prev.filter = prev.filter,
        mean.abund.filter = abund.filter,
        ...
      )

      if (!is.null(group.var)){
        reference_level <- levels(as.factor(meta_tab[,group.var]))[1]
      }

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

          # Add the data frame to the result list.
          result_list[[paste0(group_value," vs ", reference_level, " (Reference)")]] <- linda_object$output[[df_name]]
        }

        return(result_list)
      }

      # Extract data frame using function
      sub_test.list <- extract_data_frames(linda_object = linda.obj, group_var = group.var)

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
