#' Load Metadata from a MicrobiomeStat Data Object
#'
#' This internal function loads metadata from a MicrobiomeStat data object.
#' It checks for the presence of various variable names in the data object
#' and returns the data associated with the first match it finds.
#' @name load_data_obj_metadata
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#'
#' @return A data frame that contains the metadata. If no metadata is found,
#' the function returns NULL and displays a message.
#'
#' @details
#' The function first checks for the presence of a set of primary variable names in the MicrobiomeStat data object.
#' If it finds a match, it returns the data associated with the first match.
#' If it doesn't find a match, it checks for the presence of a set of possible confusing variable names.
#' If it finds multiple matches, it displays a message and returns the data associated with the first match.
#' If it finds a single match, it returns the associated data.
#' If it doesn't find a match, it performs a fuzzy search using a set of search words.
#' If it finds a match, it returns the associated data. If it doesn't find a match, it displays a message and returns NULL.
#' @keywords internal
#' @noRd
load_data_obj_metadata <- function(data.obj) {
  # Define a list of primary variable names to check first.
  # These names are commonly used in various omics studies for metadata.
  # The list includes variations to account for different naming conventions.
  primary_names <-
    c(
      "meta.dat",
      "meta_dat",
      "metadata",
      "meta_tab",
      "meta.data",
      "sample.data",
      "sample_information",
      "phenodata",
      "sample.metadata",
      "experiment.data",
      "experiment_information",
      "samplemeta",
      "experimentmeta",
      "sample_meta",
      "experiment_meta",
      "phenotype_data"
    )

  # Define a list of potentially confusing variable names.
  # These are alternative or shortened names that might be used for metadata in different contexts.
  # Including these allows for more flexible metadata detection.
  possible_confusing_names <-
    c(
      "meta",
      "meta_info",
      "meta_data",
      "sample_data",
      "sample_info",
      "exp_info",
      "exp_data",
      "phenotype",
      "samplemetadata",
      "experimentmetadata",
      "sample_metainfo",
      "experiment_metainfo",
      "phenotype_info"
    )

  # Define words for fuzzy searching.
  # These words are commonly associated with metadata and can be used for partial matching.
  # This provides a last resort for finding metadata when naming is non-standard.
  fuzzy_search_words <- c("sample", "meta")

  # Extract the names of all variables in the input data object.
  keys_in_data_obj <- names(data.obj)

  # Find the intersection between primary names and the names in the data object.
  # This identifies which of the preferred variable names are present in the data.
  existing_primary_names <-
    intersect(primary_names, keys_in_data_obj)

  # If any of the primary names are found, return the corresponding data.
  # This prioritizes the use of standard, expected variable names for metadata.
  if (length(existing_primary_names) > 0) {
    message("Using table '", existing_primary_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_primary_names[1]]])
  }

  # If no primary names are found, check for potentially confusing names.
  # This allows for flexibility in metadata naming conventions across different studies or data formats.
  existing_confusing_names <-
    intersect(possible_confusing_names, keys_in_data_obj)

  # If multiple potentially confusing names are found, warn the user and use the first one.
  # This helps manage ambiguity when multiple possible metadata tables are present in the data object.
  if (length(existing_confusing_names) > 1) {
    message(
      "Multiple potential metadata tables detected: ",
      paste(existing_confusing_names, collapse = ", "),
      ". Returning the first one."
    )
  }

  # If any potentially confusing names are found, return the corresponding data.
  # This allows the function to work with less standardized naming conventions.
  if (length(existing_confusing_names) > 0) {
    message("Using table '",
            existing_confusing_names[1],
            "' in 'data.obj'.")
    return(data.obj[[existing_confusing_names[1]]])
  }

  # If no exact matches are found, attempt fuzzy matching.
  # This provides a last resort for finding metadata when naming is highly non-standard.
  fuzzy_matches <-
    grep(paste(fuzzy_search_words, collapse = "|"), keys_in_data_obj)

  # If any fuzzy matches are found, return the first matching data.
  # This allows the function to work even with unconventional naming schemes.
  if (length(fuzzy_matches) > 0) {
    message("Using table '", names(data.obj)[fuzzy_matches[1]], "' in 'data.obj'.")
    return(data.obj[[fuzzy_matches[1]]])
  }

  # If no matches are found at all, inform the user and return NULL.
  # This indicates that no suitable metadata could be identified in the input object.
  message("No potential metadata tables detected.")
  return(NULL)
}
