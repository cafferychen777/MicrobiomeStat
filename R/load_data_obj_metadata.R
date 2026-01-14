#' Load Metadata from a MicrobiomeStat Data Object
#'
#' Internal function to extract metadata from a MicrobiomeStat data object.
#' Searches for metadata using various common naming conventions.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @return A data frame containing the metadata, or NULL if not found.
#'
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
