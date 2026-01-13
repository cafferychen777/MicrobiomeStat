#' Aggregate OTU Table by Taxonomy Level
#'
#' This function aggregates an OTU table by a specified taxonomy level.
#' It checks for consistency between OTU and taxonomy tables before aggregation.
#'
#' @param feature.tab OTU table as matrix or data.frame
#' @param feature.ann Taxonomy table as matrix or data.frame
#' @param feature.level Taxonomy level to aggregate by as character vector.
#'
#' @return Aggregated OTU table matrix.
#'
#'
#' @details
#' The function first checks that the row names of OTU and taxonomy tables match.
#' If any mismatches are found, a message is displayed but the function continues.
#' The OTU and taxonomy tables are joined by row names.
#' The OTU table is then aggregated by summing abundance values for each sample at the specified taxonomy level.
#' Rows with all zero values after aggregation are removed.
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Specify the taxonomy level
#' feature.level <- c("Family")
#'
#' # Aggregate data object by taxonomy level
#' mStat_aggregate_by_taxonomy2(peerj32.obj$feature.tab, peerj32.obj$feature.ann, feature.level)
#' }
#' @export
mStat_aggregate_by_taxonomy2 <-
  function (feature.tab,
            feature.ann,
            feature.level = NULL) {
    # Validate input: ensure that a taxonomic level for aggregation is specified
    if (is.null(feature.level)) {
      stop("feature.level can not be NULL.")
    }

    # Ensure only one feature level is provided for aggregation
    # This function is designed to aggregate at a single taxonomic level at a time
    if (length(feature.level) != 1) {
      stop(
        "The function 'mStat_aggregate_by_taxonomy2' only supports aggregation by a single feature level at a time. Please provide exactly one feature level for aggregation."
      )
    }

    # Convert input tables to data frames for consistent processing
    otu_tab <- feature.tab %>% as.data.frame()
    tax_tab <- feature.ann %>% as.data.frame()

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

    # Merge the OTU table with the taxonomy table
    # This step combines abundance data with taxonomic information
    # An inner join is used, which means only features present in both tables are retained
    message("Performing an inner join on two tables. Rows with non-matching names will be excluded.")

    otu_tax <-
      otu_tab %>% rownames_to_column("sample") %>%
      dplyr::inner_join(tax_tab %>%
                          select(all_of(feature.level)) %>%
                          rownames_to_column("sample"),
                        by = "sample") %>%
      column_to_rownames("sample")

    # Aggregate the OTU table at the specified taxonomic level
    # This process involves several steps:
    # 1. Reshape the data from wide to long format
    # 2. Group by sample and the specified taxonomic level
    # 3. Sum the abundances within each group
    # 4. Reshape back to wide format
    # 5. Replace NA values with "Unclassified" for better interpretability
    otu_tax_agg <- otu_tax %>%
      tidyr::gather(key = "sample", value = "value", -all_of(feature.level)) %>%
      dplyr::group_by_at(vars(sample,!!sym(feature.level))) %>%
      dplyr::summarise(value = sum(value)) %>%
      tidyr::spread(key = "sample", value = "value") %>%
      dplyr::mutate(!!feature.level := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(feature.level) %>%
      as.matrix()

    # Remove taxa with zero abundance across all samples
    # This step helps in reducing the dimensionality of the data and focuses on relevant taxa
    otu_tax_agg <- otu_tax_agg[rowSums(otu_tax_agg) > 0,]

    # Return the aggregated OTU table
    return(otu_tax_agg)
  }
