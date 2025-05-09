#' Update Metadata in a MicrobiomeStat Data Object
#'
#' This function updates the metadata in a MicrobiomeStat data object. It either
#' reads the metadata from a file (CSV or tab-delimited) or directly uses an
#' input dataframe.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param map.file a character string or dataframe. The path to the metadata file (CSV or TSV) or a dataframe containing the metadata.
#' @param meta.sep a character string. The field separator character in the metadata file. Default is a tab ("\\t").
#' @param quote a character string. The set of quoting characters for the metadata file. Default is a double quote ('"').
#' @param comment a character string. The comment character for the metadata file. Lines beginning with this character are ignored.
#' @param ... further arguments to be passed to the read.csv or read.table function.
#'
#' @return a list. An updated MicrobiomeStat data object that includes the new meta.dat.
#'
#' @examples
#' \dontrun{
#'   # Load required libraries
#'   library(vegan)
#'   data(peerj32.obj)
#'
#'   # Update metadata using a CSV file
#'   # peerj32.obj <- mStat_update_meta_data(peerj32.obj, "metadata.csv")
#'
#'   # Update metadata using a dataframe
#'   metadata <- data.frame(Treatment = sample(c("Control", "Treatment"),
#'                          length(colnames(peerj32.obj$feature.tab)), replace = TRUE))
#'   rownames(metadata) <- rownames(peerj32.obj$meta.dat)
#'   peerj32.obj <- mStat_update_meta_data(peerj32.obj, metadata)
#' }
#' @export
mStat_update_meta_data <- function (data.obj, map.file, meta.sep='\t', quote="\"", comment="", ...) {
  # Inform the user that the metadata file is being loaded.
  # This provides feedback on the current operation, which is useful for tracking progress.
  cat("Load meta file...\n")

  # Determine the type of input for metadata and process accordingly.
  # This allows flexibility in how users can provide metadata information.
  if (is.character(map.file)) {
    # If the input is a file path, read the file based on its type.
    # This accommodates different file formats commonly used for metadata.
    if (grepl("csv$", map.file)) {
      # For CSV files, use read.csv function.
      # Parameters are set to ensure proper reading of the file structure.
      meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote, ...)
    } else {
      # For other file types (e.g., tab-delimited), use read.table function.
      # This allows for more flexible file format handling.
      meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote, ...)
    }
  } else {
    # If the input is already a data frame, use it directly.
    # This allows users to pass pre-loaded metadata directly to the function.
    meta.dat <- map.file
  }

  # Update the metadata in the data object.
  # This step replaces any existing metadata with the newly loaded information.
  data.obj$meta.dat <- meta.dat

  # Find the intersection of sample names between the metadata and the feature table.
  # This step is crucial for ensuring that the metadata and abundance data correspond to the same samples.
  samIDs <- intersect(rownames(meta.dat), colnames(data.obj$feature.tab))

  # Check if there are any common sample names between metadata and feature table.
  # If no common samples are found, it indicates a mismatch between metadata and abundance data.
  if (length(samIDs) == 0)  stop('Sample names in the meta file and biom file differ?\n')

  # Subset the data object to include only the samples present in both metadata and feature table.
  # This ensures consistency across all components of the data object.
  data.obj <- mStat_subset_data(data.obj, samIDs)

  # Return the updated data object.
  # The returned object now has updated metadata and is consistent across all its components.
  return(data.obj)
}
