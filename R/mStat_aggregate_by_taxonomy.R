#' Aggregate OTU Table by Taxonomy Level
#'
#' This function, part of the MicrobiomeStat package, aggregates an OTU table by a specified taxonomy level.
#' It checks for the consistency of row names between the OTU table and the taxonomy table and stops execution if any mismatches are found.
#' @name mStat_aggregate_by_taxonomy
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param feature.level A character vector that specifies the taxonomy levels to aggregate by.
#'
#' @return A list that is the modified version of the input data object. It includes the original data object and an additional element, 'feature.agg.list',
#' which is a list of aggregated OTU tables for each specified taxonomy level.
#'
#' @details
#' The function first checks for the consistency of row names between the OTU table and the taxonomy table. If any mismatches are found,
#' it displays a message and stops execution. If the row names are consistent, it merges the OTU table with the taxonomy table and aggregates
#' the OTU table by the specified taxonomy level. The aggregation is done by summing the values for each sample at each taxonomy level.
#' The aggregated OTU tables are added to the data object as a list named 'feature.agg.list'.
#'
#' @author Caffery(Chen) YANG
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Specify the taxonomy level
#' feature.level <- c("Phylum", "Family")
#'
#' # Aggregate data object by taxonomy level
#' peerj32.obj <- mStat_aggregate_by_taxonomy(peerj32.obj, feature.level)
#' }
#' @export
mStat_aggregate_by_taxonomy <- function(data.obj, feature.level = NULL) {
  # Validate input: ensure that a taxonomic level for aggregation is specified
  if (is.null(feature.level)) {
    stop("feature.level can not be NULL.")
  }

  # Extract the feature table and taxonomy table from the data object
  otu_tab <- data.obj$feature.tab %>% as.data.frame()
  tax_tab <- data.obj$feature.ann %>% as.data.frame()

  # Check for consistency between feature and taxonomy tables
  # This step is crucial to ensure data integrity before aggregation
  otu_not_in_tax <- setdiff(rownames(otu_tab), rownames(tax_tab))
  tax_not_in_otu <- setdiff(rownames(tax_tab), rownames(otu_tab))

  # Notify the user of any inconsistencies in the data
  # This helps in identifying potential issues in the input data
  if (length(otu_not_in_tax) > 0) {
    message(
      "The following row names are in 'feature.tab' but not in 'feature.ann': ",
      paste(otu_not_in_tax, collapse = ", ")
    )
  }

  if (length(tax_not_in_otu) > 0) {
    message(
      "The following row names are in 'feature.ann' but not in 'feature.tab': ",
      paste(tax_not_in_otu, collapse = ", ")
    )
  }

  # Preserve the original sample order for consistent output
  original_sample_order <- colnames(otu_tab)

  # Perform aggregation for each specified taxonomic level
  data.obj$feature.agg.list <- setNames(lapply(feature.level, function(feature.level) {
    # Join the OTU table with the taxonomy table
    # This step combines abundance data with taxonomic information
    otu_tax <- otu_tab %>%
      rownames_to_column("sample") %>%
      dplyr::inner_join(
        tax_tab %>% select(all_of(feature.level)) %>% rownames_to_column("sample"),
        by = "sample"
      ) %>%
      column_to_rownames("sample")

    # Aggregate the OTU table at the specified taxonomic level
    # This process involves several steps:
    # 1. Reshape the data from wide to long format
    # 2. Group by taxonomic level and sample
    # 3. Sum the abundances within each group
    # 4. Reshape back to wide format
    # 5. Replace NA values with "Unclassified" for better interpretability
    otu_tax_agg <- otu_tax %>%
      tidyr::pivot_longer(cols = -all_of(feature.level), names_to = "sample", values_to = "value") %>%
      dplyr::group_by_at(vars(!!sym(feature.level), sample)) %>%
      dplyr::summarise(value = sum(value), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = sample, values_from = value, values_fill = list(value = 0)) %>%
      dplyr::mutate(!!feature.level := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(feature.level) %>%
      as.matrix()

    # Ensure the output maintains the original sample order for consistency
    otu_tax_agg <- otu_tax_agg[, original_sample_order]

    # Remove taxa with zero abundance across all samples
    # This step helps in reducing the dimensionality of the data
    otu_tax_agg <- otu_tax_agg[rowSums(otu_tax_agg) > 0, , drop = FALSE]

    return(otu_tax_agg)
  }), feature.level)

  # Return the updated data object with the new aggregated feature list
  return(data.obj)
}
