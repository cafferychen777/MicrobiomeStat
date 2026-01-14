#' Load Count Data from a MicrobiomeStat Data Object
#'
#' Internal function to extract the count table from a MicrobiomeStat data object.
#' Searches for count data using various common naming conventions.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @return The count table (matrix) if found, NULL otherwise.
#'
#' @keywords internal
#' @noRd
load_data_obj_count <- function(data.obj) {
  # Define a list of primary variable names to check first.
  # These names are commonly used in microbiome and omics studies for count data.
  primary_names <-
    c(
      "feature.tab",
      "otu.tab",
      "otu_tab",
      "expression.tab",
      "exp.tab",
      "gene.counts",
      "metabolite.tab",
      "metab.tab",
      "metabolite.counts",
      "protein.tab",
      "prot.tab",
      "protein.counts",
      "transcript.tab",
      "transcript.counts",
      "rna.tab",
      "rna.counts",
      "peptide.tab",
      "peptide.counts"
    )

  # Define a list of potentially confusing variable names.
  # These are alternative names that might be used for count data in different contexts.
  possible_confusing_names <-
    c(
      "otu.table",
      "otu_matrix",
      "otu_mat",
      "otu_data",
      "otu_list",
      "otu_counts",
      "otu_abundance",
      "exp_matrix",
      "exp_mat",
      "exp_data",
      "metab_matrix",
      "metab_mat",
      "metab_data",
      "prot_matrix",
      "prot_mat",
      "prot_data",
      "transcript_matrix",
      "transcript_mat",
      "transcript_data",
      "rna_matrix",
      "rna_mat",
      "rna_data",
      "peptide_matrix",
      "peptide_mat",
      "peptide_data"
    )

  # Define words for fuzzy searching.
  # These words are commonly associated with count data and can be used for partial matching.
  fuzzy_search_words <- c("count", "abundance")

  # Extract the names of all variables in the input data object.
  keys_in_data_obj <- names(data.obj)

  # Find the intersection between primary names and the names in the data object.
  # This identifies which of the preferred variable names are present in the data.
  existing_primary_names <-
    intersect(primary_names, keys_in_data_obj)

  # If any of the primary names are found, return the corresponding data.
  # This prioritizes the use of standard, expected variable names.
  if (length(existing_primary_names) > 0) {
    message("Using table '", existing_primary_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_primary_names[1]]])
  }

  # If no primary names are found, check for potentially confusing names.
  # This allows for flexibility in variable naming conventions.
  existing_confusing_names <-
    intersect(possible_confusing_names, keys_in_data_obj)

  # If multiple potentially confusing names are found, warn the user and use the first one.
  # This helps manage ambiguity when multiple possible count tables are present.
  if (length(existing_confusing_names) > 1) {
    message(
      "Multiple potential count tables detected: ",
      paste(existing_confusing_names, collapse = ", "),
      ". Returning the first one."
    )
  }

  # If any potentially confusing names are found, return the corresponding data.
  if (length(existing_confusing_names) > 0) {
    message("Using table '", existing_confusing_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_confusing_names[1]]])
  }

  # If no exact matches are found, attempt fuzzy matching.
  # This provides a last resort for finding count data when naming is non-standard.
  fuzzy_matches <-
    grep(paste(fuzzy_search_words, collapse = "|"), keys_in_data_obj)

  # If any fuzzy matches are found, return the first matching data.
  if (length(fuzzy_matches) > 0) {
    message("Using table '", names(data.obj)[fuzzy_matches[1]], "' in 'data.obj'.")
    return(data.obj[[fuzzy_matches[1]]])
  }

  # If no matches are found at all, inform the user and return NULL.
  # This indicates that no suitable count data could be identified in the input object.
  message("No potential count tables detected.")
  return(NULL)
}
