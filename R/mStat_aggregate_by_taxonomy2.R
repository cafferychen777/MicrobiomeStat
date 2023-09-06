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
