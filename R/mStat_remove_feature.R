#' Remove Specific Features from Data Object
#'
#' Removes specified features from the feature table and annotations.
#'
#' @inheritParams mStat_data_obj_doc
#' @param featureIDs Character vector of feature IDs to remove.
#' @param feature.level Taxonomic level to match feature IDs. Use "original" for row names.
#'
#' @return A MicrobiomeStat data object with specified features removed.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' featureIDs <- rownames(peerj32.obj$feature.tab)[1:3]
#' peerj32.obj <- mStat_remove_feature(peerj32.obj, featureIDs, feature.level = "original")
#' }
#'
#' @export
#' @importFrom dplyr filter select
mStat_remove_feature <- function (data.obj, featureIDs, feature.level = NULL) {

  # Check if feature IDs are provided in a character vector
  if (is.numeric(featureIDs)) {
    warning('featureIDs should be characters!\n')
  }

  # If feature.level is 'original', then filter by rownames
  if ("original" %in% feature.level) {
    # Save the initial number of features
    initial_features <- nrow(data.obj$feature.tab)

    data.obj$feature.tab <- data.obj$feature.tab[!rownames(data.obj$feature.tab) %in% featureIDs, ]
    data.obj$feature.ann <- data.obj$feature.ann[!rownames(data.obj$feature.ann) %in% featureIDs, ]

    # Compute the number of removed features
    removed_features <- initial_features - nrow(data.obj$feature.tab)

    # Print the message
    message(paste(removed_features, "features were removed from the", feature.level, "level. The remaining number of features is", nrow(data.obj$feature.tab), "."))

    return(data.obj)
  }

  # Get the feature levels from the feature annotation dataframe
  feature.levels <- colnames(data.obj$feature.ann)
  found.levels <- feature.levels[sapply(feature.levels, function(x) any(data.obj$feature.ann[, x] %in% featureIDs))]

  if (is.null(feature.level)) {
    message(paste("The feature IDs were found in the following levels:", paste(found.levels, collapse = ", "),
               ". Please specify the 'feature.level' parameter for the level you want to remove."))
    return()
  } else if (!(feature.level %in% found.levels)) {
    message(paste("The specified 'feature.level' does not contain any of the feature IDs. The feature IDs were found in the following levels:",
               paste(found.levels, collapse = ", "), "."))
    return()
  }

  # Save the initial number of features before filtering
  initial_features <- nrow(data.obj$feature.tab)

  # Use the tidyverse function 'filter' to subset the data
  data.obj$feature.tab <- dplyr::filter(data.obj$feature.tab %>% as.data.frame(), !(data.obj$feature.ann[, feature.level] %in% featureIDs))
  data.obj$feature.ann <- dplyr::filter(data.obj$feature.ann %>% as.data.frame(), !(data.obj$feature.ann[, feature.level] %in% featureIDs))

  # Message about removed features (correct calculation)
  removed_features <- initial_features - nrow(data.obj$feature.tab)
  message(paste(removed_features, "features were removed from the", feature.level, "level. The remaining number of features is", nrow(data.obj$feature.tab), "."))

  # Check if the feature.agg.list exists in data.obj
  if (!is.null(data.obj$feature.agg.list)) {
    # Retrieve the names of all elements in feature.agg.list
    feature_levels_in_agg_list <- names(data.obj$feature.agg.list)

    # Inform the user that feature.agg.list is being recalculated
    message("Recalculating feature.agg.list using mStat_aggregate_by_taxonomy function...")

    data.obj <- mStat_aggregate_by_taxonomy(data.obj, feature_levels_in_agg_list)

    # Inform the user that the recalculation is complete
    message("Recalculation of feature.agg.list is complete. The updated feature.agg.list contains", length(feature_levels_in_agg_list), "elements.")
  }

  return(data.obj)
}
