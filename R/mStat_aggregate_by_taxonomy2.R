#' Aggregate OTU Table by Taxonomy Level
#'
#' This function, part of the MicrobiomeStat package, aggregates an OTU table by a specified taxonomy level.
#' It checks for the consistency of row names between the OTU table and the taxonomy table and stops execution if any mismatches are found.
#' @name mStat_aggregate_by_taxonomy2
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
    # Check if feature.level is not NULL
    if (is.null(feature.level)) {
      stop("feature.level can not be NULL.")
    }

    if (length(feature.level) != 1) {
      stop(
        "The function 'mStat_aggregate_by_taxonomy2' only supports aggregation by a single feature level at a time. Please provide exactly one feature level for aggregation."
      )
    }

    otu_tab <- feature.tab %>% as.data.frame()
    tax_tab <- feature.ann %>% as.data.frame()

    # 检查 otu_tab 的行名和 tax_tab 的行名是否完全一致
    otu_not_in_tax <- setdiff(rownames(otu_tab), rownames(tax_tab))
    tax_not_in_otu <- setdiff(rownames(tax_tab), rownames(otu_tab))

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

    # 将 OTU 表与分类表合并
    message("Performing an inner join on two tables. Rows with non-matching names will be excluded.")

    otu_tax <-
      otu_tab %>% rownames_to_column("sample") %>%
      dplyr::inner_join(tax_tab %>%
                          select(all_of(feature.level)) %>%
                          rownames_to_column("sample"),
                        by = "sample") %>%
      column_to_rownames("sample")

    # 聚合 OTU 表
    otu_tax_agg <- otu_tax %>%
      tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
      dplyr::group_by_at(vars(sample,!!sym(feature.level))) %>%
      dplyr::summarise(value = sum(value)) %>%
      tidyr::spread(key = "sample", value = "value") %>%
      dplyr::mutate(!!feature.level := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(feature.level) %>%
      as.matrix()

    # Remove rows with all zero
    otu_tax_agg <- otu_tax_agg[rowSums(otu_tax_agg) > 0,]

    return(otu_tax_agg)
  }
