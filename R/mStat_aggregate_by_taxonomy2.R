#' Aggregate Feature Table by Taxonomy Level (Low-Level Interface)
#'
#' A lower-level interface for taxonomy aggregation that works directly with
#' feature and annotation matrices, without requiring a MicrobiomeStat data object.
#'
#' @param feature.tab Feature table (matrix or data.frame) with taxa as rows
#'   and samples as columns.
#' @param feature.ann Taxonomy annotation table (matrix or data.frame) with taxa
#'   as rows and taxonomy levels as columns.
#' @param feature.level Character string specifying the taxonomy level to aggregate
#'   by. Must match a column name in feature.ann. Supports only a single level;
#'   for multiple levels, call this function multiple times or use
#'   \code{\link{mStat_aggregate_by_taxonomy}}.
#'
#' @return Aggregated feature table as a matrix with taxonomy groups as rows
#'   and samples as columns. The original sample order is preserved.
#'
#' @details
#' This function provides direct access to the taxonomy aggregation algorithm
#' without the data object wrapper. It is useful when:
#' \itemize{
#'   \item You have separate feature and annotation tables
#'   \item You want to aggregate a single level without modifying a data object
#'   \item You are building custom analysis pipelines
#' }
#'
#' The function performs an inner join between feature.tab and feature.ann,
#' so only features present in both tables will be included in the output.
#' NA taxonomy values are replaced with "Unclassified".
#'
#' @author Caffery(Chen) YANG
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Aggregate to Family level
#' family_table <- mStat_aggregate_by_taxonomy2(
#'   feature.tab = peerj32.obj$feature.tab,
#'   feature.ann = peerj32.obj$feature.ann,
#'   feature.level = "Family"
#' )
#'
#' # Aggregate to multiple levels separately
#' phylum_table <- mStat_aggregate_by_taxonomy2(
#'   peerj32.obj$feature.tab, peerj32.obj$feature.ann, "Phylum"
#' )
#' genus_table <- mStat_aggregate_by_taxonomy2(
#'   peerj32.obj$feature.tab, peerj32.obj$feature.ann, "Genus"
#' )
#' }
#'
#' @seealso \code{\link{mStat_aggregate_by_taxonomy}} for the higher-level interface
#'   that works with MicrobiomeStat data objects and supports multiple levels.
#'
#' @export
mStat_aggregate_by_taxonomy2 <- function(feature.tab,
                                          feature.ann,
                                          feature.level = NULL) {
  # Validate inputs
  if (is.null(feature.level)) {
    stop("feature.level cannot be NULL. Please specify a taxonomy level.")
  }

  if (length(feature.level) != 1) {
    stop(
      "mStat_aggregate_by_taxonomy2 supports only a single feature level. ",
      "For multiple levels, use mStat_aggregate_by_taxonomy() or call this ",
      "function multiple times."
    )
  }

  if (is.null(feature.tab)) {
    stop("feature.tab cannot be NULL.")
  }

  if (is.null(feature.ann)) {
    stop("feature.ann cannot be NULL.")
  }

  # Delegate to core function (defined in mStat_aggregate_by_taxonomy.R)
  .aggregate_single_level(
    feature.tab = feature.tab,
    feature.ann = feature.ann,
    feature.level = feature.level
  )
}
