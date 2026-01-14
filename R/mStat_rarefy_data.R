#' Rarefy data to a specified sampling depth
#'
#' This function takes an input object containing a feature table and optionally a specified rarefaction depth.
#' If no depth is specified, the smallest column sum is chosen as the depth.
#' It rarefies the data to the specified depth and returns the object with the rarefied feature table.
#' @name mStat_rarefy_data
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param depth Depth at which to rarefy; if not provided, defaults to the smallest column sum.
#'
#' @return The input object with the rarefied feature table added.
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

  # Validate input data structure
  # Ensuring the input is a list is crucial for maintaining the expected data format
  if (!is.list(data.obj)) {
    stop("data.obj should be a list.")
  }

  # Extract the feature table from the input data object
  # The feature table is the core data structure in microbiome analysis, representing taxon abundances across samples
  otu_tab <- as.data.frame(data.obj$feature.tab)

  # Determine the rarefaction depth
  # Rarefaction is a technique used to standardize sampling effort across different samples
  if (is.null(depth)) {
    # If no depth is specified, use the smallest sample size (column sum) as the rarefaction depth
    # This approach ensures that all samples can be rarefied without introducing NA values
    depth = min(colSums(otu_tab))
    message(paste("No depth specified, using the smallest column sum: ", depth))
  } else if (!is.numeric(depth) || depth < 0) {
    # Validate that the specified depth is a non-negative number
    stop("Depth should be a non-negative number.")
  }

  # Check if all samples have enough reads for the specified rarefaction depth
  # This step prevents potential issues with samples that have fewer reads than the specified depth
  if (any(colSums(otu_tab) < depth)) {
    stop("Not all samples have enough reads for the specified depth.")
  }

  # Perform rarefaction
  # Rarefaction helps to account for uneven sequencing depth across samples
  if (all(round(colSums(otu_tab),5) == 1)){
    # If the data is already in relative abundance format, no rarefaction is needed
    rarefied_otu_tab <- as.matrix(otu_tab)
  } else {
    # Use the vegan package to perform rarefaction
    # This function randomly subsamples the feature table to the specified depth
    rarefied_otu_tab <- t(vegan::rrarefy(t(otu_tab), sample = depth))
  }

  # Identify features (rows) that become all-zero after rarefaction
  # These features are typically rare taxa that are lost due to the subsampling process
  zero_row_indices <- which(rowSums(rarefied_otu_tab) == 0)

  # Update the data object with the rarefied feature table
  if(length(zero_row_indices) > 0) {
    # Remove all-zero rows from the rarefied feature table
    rarefied_otu_tab <- rarefied_otu_tab[-zero_row_indices, ]
    # Update the data object with the filtered rarefied feature table
    data.obj <- update_data_obj_count(data.obj, rarefied_otu_tab)
    # Remove corresponding rows from the feature annotations if they exist
    if(!is.null(data.obj$feature.ann)) {
      data.obj$feature.ann <- data.obj$feature.ann[-zero_row_indices, ]
    }
    # Inform the user about the removed features
    message(paste("Removed ", length(zero_row_indices), " rows from the feature table and feature annotations due to all-zero counts after rarefaction. The removed rows are: ", paste0(zero_row_indices, collapse = ", "), "."))
  } else {
    # If no all-zero rows were created, simply update the data object with the rarefied feature table
    data.obj <- update_data_obj_count(data.obj, rarefied_otu_tab)
  }

  # Identify samples that have non-zero total counts after rarefaction
  # This step ensures that we don't keep samples that might have become all-zero due to rarefaction
  non_zero_col_names <- colnames(rarefied_otu_tab)[which(colSums(rarefied_otu_tab) != 0)]

  # Subset the data object to include only non-zero samples
  if (length(non_zero_col_names) > 0){
    data.obj <- mStat_subset_data(data.obj = data.obj, samIDs = non_zero_col_names)
  }

  # Update aggregated feature lists if they exist
  # This step ensures that any pre-existing aggregations are updated with the rarefied data
  if ('feature.agg.list' %in% names(data.obj)) {
    data.obj <- mStat_aggregate_by_taxonomy(data.obj, names(data.obj$feature.agg.list))
    message("feature.agg.list has been updated.")
  }

  # Return the updated data object with rarefied feature table
  return(data.obj)
}
