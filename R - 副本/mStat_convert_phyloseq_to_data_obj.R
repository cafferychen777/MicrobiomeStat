#' Convert a Phyloseq Object to a MicrobiomeStat's Data Object
#'
#' This function is a part of the MicrobiomeStat package. It takes a phyloseq object, extracts relevant information and repackages it in the MicrobiomeStat's data object format. This enables easy use of phyloseq data with MicrobiomeStat's analysis functions.
#'
#' @name mStat_convert_phyloseq_to_data_obj
#' @param phylo.obj A phyloseq object to be converted. This should contain an OTU (operational taxonomic unit) table, sample data, taxonomy table, and a phylogenetic tree.
#'
#' @return A MicrobiomeStat data object (a list) containing the following elements:
#' \itemize{
#'   \item feature.tab: A matrix of the feature table (OTU table). Rows with a sum of zero are removed, so only the features present in the samples are included.
#'   \item meta.dat: A data frame of the sample data. This contains the metadata for each of the samples.
#'   \item feature.ann: A matrix of the feature annotation table (taxonomy table). Only the rows that exist in the feature table are included.
#'   \item tree: A phylogenetic tree. The tree is rooted by midpointing if it is not already rooted. Tips not present in the OTU table are dropped.
#' }
#'
#' @examples
#' \dontrun{
#'   library(microbiome)
#'   data(peerj32)
#'   peerj32.phy <- peerj32$phyloseq
#'   data.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#' }
#'
#' @details
#' This function checks each component (OTU table, sample data, taxonomy table, and phylogenetic tree) of the phyloseq object for null values. If a component is not null, it is converted to the appropriate format and added to the MicrobiomeStat data object. The OTU and taxonomy tables are converted to matrices, while the sample data is converted to a data frame. The phylogenetic tree is checked if it is rooted, and if not, it is rooted by midpointing. Tips not present in the OTU table are dropped from the tree. This ensures the output data object is consistent and ready for further microbiome statistical analysis.
#'
#' @author Jun Chen
#' @seealso \code{\link[phyloseq]{phyloseq-class}}
#' @references McMurdie PJ, Holmes S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLoS ONE. 2013;8(4):e61217.
#'
#' @export
mStat_convert_phyloseq_to_data_obj <- function (phylo.obj) {

  data.obj <- list()

  if (!is.null(phylo.obj@otu_table)) {
    data.obj$feature.tab <- phylo.obj@otu_table %>%
      as.data.frame() %>%
      as.matrix()

    data.obj$feature.tab <- data.obj$feature.tab[rowSums(data.obj$feature.tab) > 0, ]
  }

  if (!is.null(phylo.obj@sam_data)) {
    data.obj$meta.dat <- phylo.obj@sam_data %>% as.matrix() %>%
      as.data.frame()
    # 检查是否存在"sample"列，如果存在，就删除
    if ("sample" %in% colnames(data.obj$meta.dat)) {
      data.obj$meta.dat <- data.obj$meta.dat %>% select(-sample)
    }
  }

  if (!is.null(phylo.obj@tax_table)) {
    data.obj$feature.ann <- phylo.obj@tax_table %>%
      as.data.frame() %>%
      as.matrix()

    if (exists("feature.tab", data.obj)) {
      data.obj$feature.ann <- data.obj$feature.ann[rownames(data.obj$feature.ann) %in% rownames(data.obj$feature.tab), ]
    }
  }

  if (!is.null(phylo.obj@phy_tree)) {
    data.obj$tree <- phylo.obj@phy_tree
    if (!is.rooted(data.obj$tree)) {
      message('Root the tree by midpointing ...')
      data.obj$tree <- midpoint(data.obj$tree)
    }
    if (exists("feature.tab", data.obj)) {
      absent <- data.obj$tree$tip.label[!(data.obj$tree$tip.label %in% rownames(data.obj$feature.tab))]
      if (length(absent) != 0) {
        message('Drop features not in the feature table ...')
        data.obj$tree <- drop.tip(data.obj$tree, absent)
      }
    }
  }

  return(data.obj)
}
