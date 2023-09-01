#' Generate Taxa Test Pair
#'
#' This function takes a MicrobiomeStat data object as input, filters taxa based on
#' prevalence and abundance thresholds, aggregates taxon abundances by sample,
#' fits linear models using the linda method to identify significant taxon changes
#' across groups over time, accounting for specified covariates. It returns data frames
#' summarizing the results for each taxonomic level, including statistics from the
#' linear models.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A string that specifies the name of the subject variable column in the metadata.
#' @param time.var A string that specifies the name of the time variable column in the metadata. If not provided, it's NULL by default.
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
#' library(vegan)
#'
#' data(peerj32.obj)
#' generate_taxa_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   feature.level = "Family",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.dat.type = "count"
#' )
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
           group.var,
           adj.vars,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion"),
           ...) {
    # Extract data
    mStat_validate_data(data.obj)

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    # Create the formula
    fixed_effects <- paste(paste(adj.vars, collapse = " + "), group.var, time.var, sep = " + ")
    random_effects <- paste("(1|", subject.var, ")", sep = "")
    formula <- paste(fixed_effects, random_effects, sep = " + ")

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
        otu_tax_agg <- load_data_obj_count(data.obj)
      }

      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter)

      linda.obj <- linda(feature.dat = otu_tax_agg,
                         meta.dat = meta_tab,
                         formula = paste("~",formula),
                         feature.dat.type = "proportion",
                         ...)

      # Extract relevant information
      significant_taxa <- rownames(linda.obj$feature.dat.use)

      # 计算每个分组的平均丰度
      prop_prev_data <-
        linda.obj$feature.dat.use %>% as.data.frame() %>% rownames_to_column(feature.level) %>%
        tidyr::gather(-!!sym(feature.level),
               key = "sample",
               value = "count") %>%
        dplyr::inner_join(linda.obj$meta.dat.use %>% rownames_to_column("sample"), by = "sample", relationship = "many-to-many") %>%
        dplyr::group_by(!!sym(group.var), !!sym(feature.level)) %>%
        dplyr::summarise(mean_abundance = mean(count),
                  prevalence = sum(count > 0) / dplyr::n())

      # Initialize results table
      results <- data.frame()

      # Iterate over each element in linda.obj$output
      for (output_element in names(linda.obj$output)) {
        current_output <- linda.obj$output[[output_element]]

        for (taxa in significant_taxa) {
          baseMean <- current_output[taxa, "baseMean"]
          log2FoldChange <- current_output[taxa, "log2FoldChange"]
          lfcSE <- current_output[taxa, "lfcSE"]
          stat <- current_output[taxa, "stat"]
          pvalue <- current_output[taxa, "pvalue"]
          padj <- current_output[taxa, "padj"]

          # 检查group.var是不是因子
          if (is_categorical(linda.obj$meta.dat.use[[group.var]])) {
            for (group in unique(linda.obj$meta.dat.use[[group.var]])) {
              group_data <-
                prop_prev_data[which(prop_prev_data[[feature.level]] == taxa &
                                       prop_prev_data[[group.var]] == group), ]
              mean_prop <- group_data$mean_abundance
              mean_prev <- group_data$prevalence

              results <- rbind(
                results,
                data.frame(
                  Variable = taxa,
                  Group = group,
                  Base.Mean = baseMean,
                  Log2.Fold.Change = log2FoldChange,
                  LFC.SE = lfcSE,
                  Stat = stat,
                  P.Value = pvalue,
                  Adjusted.P.Value = padj,
                  Mean.Abundance = mean_prop,
                  Mean.Prevalence = mean_prev,
                  Output.Element = output_element
                )
              )
            }
          } else {
            total_data <-
              prop_prev_data[which(prop_prev_data[[feature.level]] == taxa), ]
            mean_prop <- total_data$mean_abundance
            mean_prev <- total_data$prevalence

            results <- rbind(
              results,
              data.frame(
                Variable = taxa,
                Base.Mean = baseMean,
                Log2.Fold.Change = log2FoldChange,
                LFC.SE = lfcSE,
                Stat = stat,
                P.Value = pvalue,
                Adjusted.P.Value = padj,
                Mean.Abundance = mean_prop,
                Mean.Prevalence = mean_prev,
                Output.Element = output_element
              )
            )
          }
        }
      }

      return(results)

    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)

  }
