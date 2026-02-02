#' @title Taxa Data Utilities
#' @description Helper functions for retrieving, aggregating, and filtering
#'   feature (OTU/ASV/taxa) tables. Centralizes the aggregation + extraction +
#'   filtering pipeline that was previously duplicated across 30+ files.
#' @name taxa_utils
#' @keywords internal
NULL


#' Retrieve a filtered feature table for a given taxonomic level
#'
#' Performs the three-step pipeline that appears in virtually every
#' \code{generate_taxa_*} function:
#' \enumerate{
#'   \item Aggregate features to \code{feature.level} if not already cached.
#'   \item Extract the aggregated table (or the raw table for \code{"original"}).
#'   \item Optionally filter by prevalence / abundance, then convert rownames
#'         to a column named \code{feature.level}.
#' }
#'
#' @param data.obj A MicrobiomeStat data object containing at least
#'   \code{feature.tab} and optionally \code{feature.agg.list}.
#' @param feature.level Character string specifying the taxonomic level to
#'   extract. Use \code{"original"} to retrieve \code{data.obj$feature.tab}
#'   without aggregation.
#' @param prev.filter Numeric prevalence threshold passed to
#'   \code{\link{mStat_filter}}. Set to 0 to skip filtering. Default 0.
#' @param abund.filter Numeric abundance threshold passed to
#'   \code{\link{mStat_filter}}. Set to 0 to skip filtering. Default 0.
#'
#' @return A data.frame with a leading column named \code{feature.level}
#'   (former rownames) followed by sample columns.
#'
#' @keywords internal
get_taxa_data <- function(data.obj,
                          feature.level,
                          prev.filter = 0,
                          abund.filter = 0) {

  # Step 1: Aggregate if needed (idempotent, cached in data.obj$feature.agg.list)
  if (is.null(data.obj$feature.agg.list[[feature.level]]) &
      feature.level != "original") {
    data.obj <-
      mStat_aggregate_by_taxonomy(data.obj = data.obj,
                                  feature.level = feature.level)
  }

  # Step 2: Extract the appropriate table
  if (feature.level != "original") {
    otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
  } else {
    otu_tax_agg <- data.obj$feature.tab
  }

  # Step 3: Filter and convert rownames to column
  otu_tax_agg %>%
    as.data.frame() %>%
    mStat_filter(prev.filter = prev.filter,
                 abund.filter = abund.filter) %>%
    tibble::rownames_to_column(feature.level)
}
