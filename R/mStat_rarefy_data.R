#' Rarefy data to a specified sampling depth
#'
#' This function takes an input object containing an OTU table and optionally a specified rarefaction depth.
#' If no depth is specified, the smallest column sum is chosen as the depth.
#' It rarefies the data to the specified depth and returns the object with the rarefied OTU table.
#' @name mStat_rarefy_data
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param depth Depth at which to rarefy; if not provided, defaults to the smallest column sum.
#'
#' @return The input object with the rarefied OTU table added.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and data
#' library(vegan)
#' data(peerj32.obj)
#' # Perform aggregation to create feature.agg.list with Phylum and Family
#' peerj32.obj <- mStat_aggregate_by_taxonomy(peerj32.obj, c("Phylum", "Family"))
#' # Perform rarefication, remove all-zero rows and update feature.agg.list
#' rarefy_peerj32.obj <- mStat_rarefy_data(data.obj = peerj32.obj)
#' }
#' @export
mStat_rarefy_data <- function(data.obj, depth = NULL) {

  # Check if data.obj is the correct type
  if (!is.list(data.obj)) {
    stop("data.obj should be a list.")
  }

  otu_tab <- as.data.frame(load_data_obj_count(data.obj))

  # If no depth is specified, use the smallest column sum as the depth
  if (is.null(depth)) {
    depth = min(colSums(otu_tab))
    message(paste("No depth specified, using the smallest column sum: ", depth))
  } else if (!is.numeric(depth) || depth < 0) {
    stop("Depth should be a non-negative number.")
  }

  # Check if all samples have a column sum greater than the specified depth
  if (any(colSums(otu_tab) < depth)) {
    stop("Not all samples have enough reads for the specified depth.")
  }

  # Rarefy the data
  rarefied_otu_tab <- t(vegan::rrarefy(t(otu_tab), sample = depth))

  # Identify all-zero rows
  zero_row_indices <- which(rowSums(rarefied_otu_tab) == 0)

  # Update the data.obj
  if(length(zero_row_indices) > 0) {
    rarefied_otu_tab <- rarefied_otu_tab[-zero_row_indices, ]
    data.obj <- update_data_obj_count(data.obj, rarefied_otu_tab)
    if(!is.null(data.obj$feature.ann)) {
      data.obj$feature.ann <- data.obj$feature.ann[-zero_row_indices, ]
    }
    message(paste("Removed ", length(zero_row_indices), " rows from the OTU table and feature annotations due to all-zero counts after rarefication. The removed rows are: ", paste0(zero_row_indices, collapse = ", "), "."))
  } else {
    data.obj <- update_data_obj_count(data.obj, rarefied_otu_tab)
  }

  non_zero_col_names <- colnames(rarefied_otu_tab)[which(colSums(rarefied_otu_tab) != 0)]

  if (length(non_zero_col_names) > 0){
    data.obj <- mStat_subset_data(data.obj = data.obj, samIDs = non_zero_col_names)
  }

  # Normalize feature.agg.list if it exists
  if ('feature.agg.list' %in% names(data.obj)) {
    data.obj <- mStat_aggregate_by_taxonomy(data.obj, names(data.obj$feature.agg.list))
    message("feature.agg.list has been updated.")
  }

  return(data.obj)
}
