#' @title Core Taxonomy Aggregation Function
#'
#' @description Internal function that performs the actual aggregation of a feature
#' table by a single taxonomy level. This is the computational core shared by all
#' public aggregation interfaces.
#'
#' @param feature.tab Feature table (matrix or data.frame) with taxa as rows and samples as columns.
#' @param feature.ann Taxonomy annotation table (matrix or data.frame) with taxa as rows.
#' @param feature.level Single character string specifying the taxonomy level to aggregate by.
#'
#' @return Aggregated feature table as a matrix with taxonomy groups as rows and samples as columns.
#'
#' @details
#' The aggregation process:
#' 1. Validates consistency between feature and taxonomy tables
#' 2. Joins tables by row names (inner join - only matching features retained)
#' 3. Reshapes to long format, groups by taxonomy level, sums abundances
#' 4. Reshapes back to wide format
#' 5. Replaces NA taxonomy values with "Unclassified"
#' 6. Removes zero-sum rows
#'
#' @keywords internal
.aggregate_single_level <- function(feature.tab, feature.ann, feature.level) {
  # Convert to data frames for consistent processing

  otu_tab <- as.data.frame(feature.tab)
  tax_tab <- as.data.frame(feature.ann)

  # Preserve original sample order for consistent output

  original_sample_order <- colnames(otu_tab)

  # Validate that feature.level exists in taxonomy table

if (!feature.level %in% colnames(tax_tab)) {
    stop(
      "feature.level '", feature.level, "' not found in feature.ann. ",
      "Available levels: ", paste(colnames(tax_tab), collapse = ", ")
    )
  }

  # Check for row name consistency and report mismatches
  otu_not_in_tax <- setdiff(rownames(otu_tab), rownames(tax_tab))
  tax_not_in_otu <- setdiff(rownames(tax_tab), rownames(otu_tab))

  if (length(otu_not_in_tax) > 0) {
    message(
      "Note: ", length(otu_not_in_tax), " features in feature.tab not found in feature.ann ",
      "(will be excluded from aggregation)"
    )
  }

  if (length(tax_not_in_otu) > 0) {
    message(
      "Note: ", length(tax_not_in_otu), " features in feature.ann not found in feature.tab ",
      "(ignored)"
    )
  }

  # Join feature table with taxonomy annotation
  # Using inner_join: only features present in both tables are retained
  # Convert taxonomy column to character to handle factor types consistently
  tax_subset <- tax_tab %>%
    dplyr::select(dplyr::all_of(feature.level)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    tibble::rownames_to_column("feature_id")

  otu_tax <- otu_tab %>%
    tibble::rownames_to_column("feature_id") %>%
    dplyr::inner_join(tax_subset, by = "feature_id") %>%
    tibble::column_to_rownames("feature_id")

  # Aggregate by taxonomy level:
  # 1. Reshape wide -> long
  # 2. Group by taxonomy level and sample
  # 3. Sum abundances within each group
  # 4. Reshape long -> wide
  # 5. Handle NA taxonomy values
  otu_tax_agg <- otu_tax %>%
    tidyr::pivot_longer(
      cols = -dplyr::all_of(feature.level),
      names_to = "sample",
      values_to = "value"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(feature.level, "sample")))) %>%
    dplyr::summarise(value = sum(value), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = "sample",
      values_from = "value",
      values_fill = list(value = 0)
    ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(feature.level),
        ~ tidyr::replace_na(., "Unclassified")
      )
    ) %>%
    tibble::column_to_rownames(feature.level) %>%
    as.matrix()

  # Restore original sample order
  otu_tax_agg <- otu_tax_agg[, original_sample_order, drop = FALSE]

  # Remove zero-abundance taxa
  otu_tax_agg <- otu_tax_agg[rowSums(otu_tax_agg) > 0, , drop = FALSE]

  return(otu_tax_agg)
}


#' Aggregate Feature Table by Taxonomy Level
#'
#' Aggregates a feature table within a MicrobiomeStat data object by one or more
#' specified taxonomy levels. The aggregated tables are stored in the
#' `feature.agg.list` component of the data object.
#'
#' @param data.obj A MicrobiomeStat data object (list) containing at minimum:
#'   \itemize{
#'     \item \code{feature.tab}: Feature table (matrix) with taxa as rows
#'     \item \code{feature.ann}: Taxonomy annotation (matrix/data.frame)
#'   }
#' @param feature.level Character vector specifying one or more taxonomy levels
#'   to aggregate by (e.g., c("Phylum", "Family", "Genus")).
#'
#' @return The input data.obj with an added/updated `feature.agg.list` component
#'   containing aggregated feature tables for each specified taxonomy level.
#'
#' @details
#' This function is the primary interface for taxonomy aggregation in MicrobiomeStat.
#' It supports aggregating to multiple taxonomy levels in a single call, storing
#' results in a named list for easy access.
#'
#' The aggregation uses an inner join between feature.tab and feature.ann,
#' meaning only features present in both tables will be included in the output.
#' A message is displayed if any features are excluded due to mismatches.
#'
#' @author Caffery(Chen) YANG
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Aggregate by multiple taxonomy levels
#' peerj32.obj <- mStat_aggregate_by_taxonomy(
#'   peerj32.obj,
#'   feature.level = c("Phylum", "Family", "Genus")
#' )
#'
#' # Access aggregated tables
#' phylum_table <- peerj32.obj$feature.agg.list$Phylum
#' family_table <- peerj32.obj$feature.agg.list$Family
#' }
#'
#' @seealso \code{\link{mStat_aggregate_by_taxonomy2}} for a lower-level interface
#'   that works directly with matrices.
#'
#' @export
mStat_aggregate_by_taxonomy <- function(data.obj, feature.level = NULL) {
  # Validate input
  if (is.null(feature.level)) {
    stop("feature.level cannot be NULL. Please specify at least one taxonomy level.")
  }

  if (is.null(data.obj$feature.tab)) {
    stop("data.obj must contain a 'feature.tab' component.")
  }

  if (is.null(data.obj$feature.ann)) {
    stop("data.obj must contain a 'feature.ann' component.")
  }

  # Aggregate each level using the core function
  data.obj$feature.agg.list <- stats::setNames(
    lapply(feature.level, function(level) {
      .aggregate_single_level(
        feature.tab = data.obj$feature.tab,
        feature.ann = data.obj$feature.ann,
        feature.level = level
      )
    }),
    feature.level
  )

  return(data.obj)
}
