#' Remove Specific Features from a MicrobiomeStat Data Object
#'
#' This function removes specific features from a data object in the MicrobiomeStat package. It can handle different feature levels and adjust accordingly. If 'original' is specified as the feature level, the function will remove features based on the original row names of 'feature.tab' and 'feature.ann'. In all other cases, the function will filter based on the specified feature level.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param featureIDs A character vector of feature IDs to be removed. The function will issue a warning if numeric IDs are provided.
#' @param feature.level An optional character string that specifies the feature level to remove features from. This can be 'original' for removing based on the original row names or any other feature level that corresponds to a column in 'feature.ann'. If no level is specified and multiple levels contain the provided feature IDs, the function will stop and ask for a specific level.
#'
#' @return A list that is the modified version of the input data object. It includes the subsetted metadata, feature table, feature names, full feature name list, and abundance list, with the specified features removed. If 'original' was specified as the feature level, the returned object will also have adjusted 'feature.tab' and 'feature.ann'.
#'
#' @examples
#' \dontrun{
#' # Load the necessary libraries
#'
#' data(peerj32.obj)
#'
#' # Use the remove_feature function
#' # Here we take the first 3 features as an example
#' featureIDs <- rownames(peerj32.obj$feature.tab)[1:3]
#' peerj32.obj <- mStat_remove_feature(peerj32.obj, featureIDs, feature.level = "original")
#' }
#' @details
#' The function first checks if the specified feature level is 'original'. If so, it removes the specified features from the 'feature.tab' and 'feature.ann' based on the original row names. If the feature level is specified and is not 'original', it removes the specified features based on this level. After removal, it recalculates 'feature.agg.list' if it exists in the data object.
#'
#' @author Jun Chen, Chen Yang
#' @seealso \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}
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

  # Use the tidyverse function 'filter' to subset the data
  data.obj$feature.tab <- dplyr::filter(data.obj$feature.tab %>% as.data.frame(), !(data.obj$feature.ann[, feature.level] %in% featureIDs))
  data.obj$feature.ann <- dplyr::filter(data.obj$feature.ann %>% as.data.frame(), !(data.obj$feature.ann[, feature.level] %in% featureIDs))

  # Message about removed features
  removed_features <- length(featureIDs) - nrow(data.obj$feature.tab)
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
