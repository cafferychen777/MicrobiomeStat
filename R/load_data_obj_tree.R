#' Load Tree Data from a MicrobiomeStat Data Object
#'
#' Internal function to extract phylogenetic tree from a MicrobiomeStat data object.
#' Searches for tree data using various common naming conventions.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @return A phylogenetic tree object, or NULL if not found.
#'
#' @keywords internal
#' @noRd
load_data_obj_tree <- function(data.obj) {

  # Define a list of primary variable names to check first.
  # These names are commonly used in phylogenetic studies for tree data.
  # The list includes variations to account for different naming conventions across studies.
  primary_names <- c("tree", "tree.data", "tree_info", "tree_table", "phylo_tree",
                     "phylo.data", "phylo_info", "phylo.tab", "phylogenetic_tree",
                     "phylogeny.data", "phylogeny_info", "phylogeny.tab")

  # Define a list of potentially confusing variable names.
  # These are alternative or shortened names that might be used for tree data in different contexts.
  # Including these allows for more flexible tree data detection.
  possible_confusing_names <- c("tree_data", "tree_info", "phylo", "phylogeny",
                                "phylo_data", "phylogeny_data")

  # Define words for fuzzy searching.
  # These words are commonly associated with phylogenetic tree data and can be used for partial matching.
  # This provides a last resort for finding tree data when naming is non-standard.
  fuzzy_search_words <- c("tree", "phylo", "phylogeny")

  # Extract the names of all variables in the input data object.
  keys_in_data_obj <- names(data.obj)

  # Find the intersection between primary names and the names in the data object.
  # This identifies which of the preferred variable names are present in the data.
  existing_primary_names <- intersect(primary_names, keys_in_data_obj)

  # If any of the primary names are found, return the corresponding data.
  # This prioritizes the use of standard, expected variable names for tree data.
  if (length(existing_primary_names) > 0) {
    return(data.obj[[existing_primary_names[1]]])
  }

  # If no primary names are found, check for potentially confusing names.
  # This allows for flexibility in tree data naming conventions across different studies or data formats.
  existing_confusing_names <- intersect(possible_confusing_names, keys_in_data_obj)

  # If multiple potentially confusing names are found, warn the user and use the first one.
  # This helps manage ambiguity when multiple possible tree objects are present in the data object.
  if (length(existing_confusing_names) > 1) {
    message("Multiple potential tree objects detected: ", paste(existing_confusing_names, collapse=", "), ". Returning the first one.")
  }

  # If any potentially confusing names are found, return the corresponding data.
  # This allows the function to work with less standardized naming conventions.
  if (length(existing_confusing_names) > 0) {
    return(data.obj[[existing_confusing_names[1]]])
  }

  # If no exact matches are found, attempt fuzzy matching.
  # This provides a last resort for finding tree data when naming is highly non-standard.
  fuzzy_matches <- grep(paste(fuzzy_search_words, collapse = "|"), keys_in_data_obj)

  # If any fuzzy matches are found, return the first matching data.
  # This allows the function to work even with unconventional naming schemes.
  if (length(fuzzy_matches) > 0) {
    return(data.obj[[fuzzy_matches[1]]])
  }

  # If no matches are found at all, inform the user and return NULL.
  # This indicates that no suitable tree data could be identified in the input object.
  message("No potential tree objects detected.")
  return(NULL)
}
