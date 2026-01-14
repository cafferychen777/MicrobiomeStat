#' Update Count Table in Data Object
#'
#' Internal function to replace the count table in a MicrobiomeStat data object.
#'
#' @inheritParams mStat_data_obj_doc
#' @param new_count_table Matrix/data.frame with new count data (must have row and column names).
#'
#' @return A MicrobiomeStat data object with updated count table.
#'
#' @keywords internal
#' @noRd
update_data_obj_count <- function(data.obj, new_count_table) {
  # Check if the new count table is valid
  if (is.null(new_count_table)) {
    stop("The provided count table is NULL.")
  }
  # Check if the new count table has row and column names
  if (is.null(rownames(new_count_table)) || is.null(colnames(new_count_table))) {
    stop("The provided count table must have row and column names.")
  }
  # Define a list of priority variable names to check
  primary_names <- c("feature.tab", "otu.tab", "otu_tab", "expression.tab", "exp.tab", "gene.counts",
                     "metabolite.tab", "metab.tab", "metabolite.counts", "protein.tab", "prot.tab", "protein.counts",
                     "count", "counts", "abundance")
  # Define a list of other variable names that might cause confusion
  possible_confusing_names <- c("otu.table", "otu_matrix", "otu_mat", "otu_data", "otu_list", "otu_counts",
                                "otu_abundance", "exp_matrix", "exp_mat", "exp_data", "metab_matrix",
                                "metab_mat", "metab_data", "prot_matrix", "prot_mat", "prot_data")
  # Check if any priority variable names exist in data.obj
  keys_in_data_obj <- names(data.obj)
  existing_primary_names <- intersect(primary_names, keys_in_data_obj)
  # If priority variable names exist, update the corresponding data
  if (length(existing_primary_names) > 0) {
    data.obj[[existing_primary_names[1]]] <- new_count_table
    return(data.obj)
  }
  # Otherwise, check if any potentially confusing variable names exist in data.obj
  existing_confusing_names <- intersect(possible_confusing_names, keys_in_data_obj)
  # If multiple conflicting names exist, send a warning message and update the data corresponding to the first conflicting name
  if (length(existing_confusing_names) > 1) {
    message("Multiple potential count tables detected: ", paste(existing_confusing_names, collapse=", "), ". Updating the first one.")
  }
  if (length(existing_confusing_names) > 0) {
    # Update the data corresponding to the first existing conflicting name
    data.obj[[existing_confusing_names[1]]] <- new_count_table
    return(data.obj)
  }
  # If no conflicting names exist, send a warning message and return NULL
  message("No potential count tables detected. The count table was not updated.")
  return(data.obj)
}
