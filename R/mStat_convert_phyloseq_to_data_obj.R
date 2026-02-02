#' Convert a Phyloseq Object to a MicrobiomeStat's Data Object
#'
#' This function is a part of the MicrobiomeStat package. It takes a phyloseq object, extracts relevant information and repackages it in the MicrobiomeStat's data object format. This enables easy use of phyloseq data with MicrobiomeStat's analysis functions.
#'
#' @name mStat_convert_phyloseq_to_data_obj
#' @param phylo.obj A phyloseq object to be converted. This should contain an OTU (operational taxonomic unit) table, sample data, taxonomy table, and a phylogenetic tree.
#'
#' @return A MicrobiomeStat data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix of the feature table. Rows with a sum of zero are removed, so only the features present in the samples are included.
#'   \item meta.dat: A data frame of the sample data. This contains the metadata for each of the samples.
#'   \item feature.ann: A matrix of the feature annotation table (taxonomy table). Only the rows that exist in the feature table are included.
#'   \item tree: A phylogenetic tree. The tree is rooted by midpointing if it is not already rooted. Tips not present in the feature table are dropped.
#' }
#'
#' @examples
#' \dontrun{
#'   #library(microbiome)
#'   # Load phangorn for midpoint rooting if tree is present
#'   #library(phangorn)
#'   #data(peerj32)
#'   #peerj32.phy <- peerj32$phyloseq
#'   #data.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#' }
#'
#' @details
#' This function checks each component (feature table, sample data, taxonomy table, and phylogenetic tree) of the phyloseq object for null values. If a component is not null, it is converted to the appropriate format and added to the MicrobiomeStat data object. The feature and taxonomy tables are converted to matrices, while the sample data is converted to a data frame.
#'
#' The function automatically handles the feature table orientation by checking the \code{taxa_are_rows} slot of the phyloseq object. If \code{taxa_are_rows} is FALSE (i.e., samples are rows and taxa are columns in the phyloseq object), the feature table will be automatically transposed to ensure that the resulting feature table has taxa as rows and samples as columns, which is the required format for MicrobiomeStat functions.
#'
#' The phylogenetic tree is checked if it is rooted, and if not, it is rooted by midpointing. Tips not present in the feature table are dropped from the tree. This ensures the output data object is consistent and ready for further microbiome statistical analysis.
#'
#' @author Jun Chen
#' @references McMurdie PJ, Holmes S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 2013;8(4):e61217.
#'
#' @export
mStat_convert_phyloseq_to_data_obj <- function (phylo.obj) {

  # Initialize an empty list to store the converted data
  # This list will contain various components of the phyloseq object in a format suitable for MicrobiomeStat
  data.obj <- list()

  # Process the OTU (Operational Taxonomic Unit) table if it exists
  if (!is.null(phylo.obj@otu_table)) {
    # Convert the feature table to a matrix format
    # This step ensures compatibility with downstream analyses and improves computational efficiency
    otu_matrix <- phylo.obj@otu_table %>%
      as.data.frame() %>%
      as.matrix()

    # Check if taxa are rows in the phyloseq object
    # If taxa_are_rows is FALSE (samples as rows), transpose the matrix
    # This ensures the feature table always has taxa as rows and samples as columns
    if (!phylo.obj@otu_table@taxa_are_rows) {
      message("Note: taxa_are_rows is FALSE in the phyloseq object. Transposing feature table to ensure taxa as rows...")
      otu_matrix <- t(otu_matrix)
    }

    data.obj$feature.tab <- otu_matrix

    # Remove features (rows) with zero counts across all samples
    # This step is crucial for reducing sparsity in the data, which can improve statistical power and reduce computational burden in subsequent analyses
    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  # Process the sample data if it exists
  if (!is.null(phylo.obj@sam_data)) {
    # Convert the sample data to a data frame format for easier manipulation.
    # Use data.frame() directly to preserve factor levels and column types;
    # the old as.matrix() pipeline coerced everything to character.
    data.obj$meta.dat <- data.frame(phylo.obj@sam_data, stringsAsFactors = FALSE)

    # Remove the "sample" column if it exists
    # This step prevents redundancy, as sample information is typically contained in the row names
    if ("sample" %in% colnames(data.obj$meta.dat)) {
      data.obj$meta.dat <- data.obj$meta.dat %>% select(-sample)
    }
  }

  # Process the taxonomy table if it exists
  if (!is.null(phylo.obj@tax_table)) {
    # Convert the taxonomy table to a matrix format
    data.obj$feature.ann <- phylo.obj@tax_table %>%
      as.data.frame() %>%
      as.matrix()

    # Ensure that the taxonomy table only includes features present in the feature table
    # This step maintains consistency between the feature table and feature annotations
    if (exists("feature.tab", data.obj)) {
      data.obj$feature.ann <- data.obj$feature.ann[rownames(data.obj$feature.ann) %in% rownames(data.obj$feature.tab), ]
    }
  }

  # Process the phylogenetic tree if it exists
  if (!is.null(phylo.obj@phy_tree)) {
    data.obj$tree <- phylo.obj@phy_tree
    
    # Check if the tree is rooted, and if not, root it by midpoint
    # Rooting the tree is important for many phylogenetic analyses and ensures consistency across different studies
    if (!ape::is.rooted(data.obj$tree)) {
      message('Root the tree by midpointing ...')
      data.obj$tree <- midpoint(data.obj$tree)
    }
    
    # Remove tree tips (leaves) that are not present in the feature table
    # This step ensures that the tree and feature table are consistent, which is crucial for phylogenetic diversity analyses
    if (exists("feature.tab", data.obj)) {
      absent <- data.obj$tree$tip.label[!(data.obj$tree$tip.label %in% rownames(data.obj$feature.tab))]
      if (length(absent) != 0) {
        message('Drop features not in the feature table ...')
        data.obj$tree <- ape::drop.tip(data.obj$tree, absent)
      }
    }
  }

  # Return the processed data object
  # This object contains all the components of the phyloseq object, reformatted for use with MicrobiomeStat functions
  return(data.obj)
}
