#' @title Update the count table in a MicrobiomeStat data object
#'
#' @description This function updates the count table within a MicrobiomeStat data object based on the provided new count table. It checks for the presence of a count table within the data object using a list of pre-defined primary and possible names. If a count table is identified, it is updated with the provided new count table. If multiple potential count tables are detected, a warning message is generated and the first one is updated.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param new_count_table A matrix or data.frame representing the new count table, which will replace the existing count table in the data.obj. It must have row and column names.
#'
#' @return A MicrobiomeStat data object with the updated count table. If no count table is detected within the data object, the function returns the original data object without any modification, and a message is printed.
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
