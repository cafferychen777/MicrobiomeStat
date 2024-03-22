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
#' library(vegan) # Used for community ecology analyses
#' library(GUniFrac) # For generalized UniFrac distances
#' library(philentropy) # For distance measures like Jensen-Shannon divergence
#' # library(phyloseq) # For handling and analyzing phylogenetic sequencing data
#' # Load example data
#' # data(GlobalPatterns) # An example dataset from the phyloseq package
#' # Convert the phyloseq object to a MicrobiomeStat data object
#' # This step is crucial for making the dataset compatible with MicrobiomeStat functions
#' # GlobalPatterns.obj <- mStat_convert_phyloseq_to_data_obj(GlobalPatterns)
#' # Calculate various beta diversity indices
#' # Beta diversity measures the difference in microbial communities across samples
#' # Bray-Curtis dissimilarity (BC)
#' # A commonly used measure of dissimilarity based on counts
#' # dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('BC'))
#' # Jaccard index
#' # A measure based on presence/absence, useful for binary data
#' # dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c("Jaccard"))
#' # UniFrac distance
#' # A phylogenetic measure of community dissimilarity
#' # Requires a phylogenetic tree as part of the input
#' # dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('UniFrac'))
#' # Weighted UniFrac distance (WUniFrac)
#' # A variation of UniFrac that accounts for relative abundance
#' # Also requires a phylogenetic tree
#' # dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('WUniFrac'))
#' # Jensen-Shannon divergence (JS)
#' # A symmetric and smoothed version of the Kullback-Leibler divergence
#' # Useful for comparing probability distributions
#' # dist.obj <- mStat_calculate_beta_diversity(GlobalPatterns.obj, dist.name = c('JS'))
#' }
#' @export
mStat_calculate_beta_diversity <- function(data.obj,
                                           dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

  if (is.null(dist.name)){
    return()
  }

  otu_tab <- data.obj$feature.tab
  tax_tab <- data.obj$feature.ann
  meta_tab <- data.obj$meta.dat

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
    phy_tree <- data.obj$tree
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

    # Normalize OTU table to relative abundances
    otu_tab_norm <- sweep(otu_tab, 2, colSums(otu_tab), FUN = "/")

    # Define KLD and JSD functions
    KLD <- function(p, q) {
      non_zero <- p > 0
      sum(p[non_zero] * log(p[non_zero] / q[non_zero]))
    }

    JSD <- function(p, q) {
      m <- (p + q) / 2
      (KLD(p, m) + KLD(q, m)) / 2
    }

    # Calculate JSD for each pair of samples
    num_samples <- ncol(otu_tab_norm)
    jsd_matrix <- matrix(0, nrow = num_samples, ncol = num_samples)

    for (i in 1:(num_samples - 1)) {
      for (j in (i + 1):num_samples) {
        jsd_matrix[i, j] <- JSD(otu_tab_norm[, i], otu_tab_norm[, j])
        jsd_matrix[j, i] <- jsd_matrix[i, j]  # JSD is symmetric
      }
    }

    rownames(jsd_matrix) <- colnames(otu_tab)
    colnames(jsd_matrix) <- colnames(otu_tab)

    # Convert to a dist object
    jsd_dist <- as.dist(jsd_matrix)
    dist.obj$JS <- jsd_dist
  }

  message("All calculations complete.")
  return(dist.obj)
}
