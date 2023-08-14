#' Calculate various beta diversity indices
#'
#' This function calculates a variety of beta diversity indices based on the input data, such as
#' Bray-Curtis (BC), Jaccard, unweighted UniFrac (UniFrac), generalized UniFrac (GUniFrac),
#' weighted UniFrac (WUniFrac), and Jensen-Shannon divergence (JS).
#' @name mStat_calculate_beta_diversity
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#'   "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac),
#'   "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence).
#'
#' @return A list containing the calculated beta diversity indices. The indices are named with the abbreviation.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(GUniFrac)
#' library(ape)
#' library(philentropy)
#'
#' # Load example data
#' data(GlobalPatterns)
#'
#' GlobalPatterns.obj <- mStat_convert_phyloseq_to_data_obj(GlobalPatterns)
#'
#' # Calculate various beta diversity indices
#' dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('BC'))
#' dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c("Jaccard"))
#' dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('UniFrac'))
#' dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('WUniFrac'))
#' dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('JS'))
#' }
#'
#' @export
mStat_calculate_beta_diversity <- function(data.obj, dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

  otu_tab <- load_data_obj_count(data.obj)
  tax_tab <- load_data_obj_taxonomy(data.obj)
  meta_tab <- load_data_obj_metadata(data.obj)

  message("Initializing distance objects...")
  dist.obj <- list()

  if ('BC' %in% dist.name) {
    message("Calculating Bray-Curtis dissimilarity...")
    dist.obj$BC <- vegan::vegdist(t(otu_tab), method = 'bray')
  }

  if ('Jaccard' %in% dist.name) {
    message("Calculating Jaccard dissimilarity...")
    dist.obj$Jaccard <- vegan::vegdist(t(otu_tab), method = 'jaccard')
  }

  if ('UniFrac' %in% dist.name || 'GUniFrac' %in% dist.name || 'WUniFrac' %in% dist.name) {
    message("Calculating UniFrac dissimilarities...")
    phy_tree <- load_data_obj_tree(data.obj)
    if (is.null(phy_tree)) {
      stop("Phylogenetic tree is required for UniFrac, GUniFrac and WUniFrac calculations.")
    }
    alpha_values <- c()
    if ('UniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 1)
    }
    if ('GUniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 0)
    }
    if ('WUniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 0.5)
    }
    message("Performing GUniFrac calculations...")
    unifracs <- GUniFrac(t(otu_tab), phy_tree, alpha = unique(alpha_values))$unifracs
  }

  if ('UniFrac' %in% dist.name) {
    message("Assigning UniFrac results...")
    dist.obj$UniFrac <- unifracs[, , "d_1"]
  }

  if ('GUniFrac' %in% dist.name) {
    message("Assigning GUniFrac results...")
    dist.obj$GUniFrac <- unifracs[, , "d_0"]
  }

  if ('WUniFrac' %in% dist.name) {
    message("Assigning WUniFrac results...")
    dist.obj$WUniFrac <- unifracs[, , "d_0.5"]
  }

  if ('JS' %in% dist.name) {
    message("Calculating Jensen-Shannon divergence...")
    jsd <- philentropy::JSD(as.matrix(t(otu_tab)))
    rownames(jsd) <- colnames(otu_tab)
    colnames(jsd) <- colnames(otu_tab)
    dist.obj$JS <- as.dist(jsd)
  }

  message("All calculations complete.")
  return(dist.obj)
}



