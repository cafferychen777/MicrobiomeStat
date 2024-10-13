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

  # Check if dist.name is NULL and return early if so
  if (is.null(dist.name)){
    return()
  }

  # Extract the OTU table from the data object
  otu_tab <- data.obj$feature.tab

  # Check if the data has been rarefied by comparing the column sums
  # Rarefaction ensures equal sampling depth across all samples, which is important for many beta diversity measures
  if (length(unique(colSums(otu_tab))) != 1) {
    warning("It appears the data may not have been rarefied. Please verify.")
  }

  # Initialize a list to store the calculated distance matrices
  message("Initializing distance objects...")
  dist.obj <- list()

  # Calculate Bray-Curtis dissimilarity if requested
  # Bray-Curtis quantifies the compositional dissimilarity between two different sites, based on counts at each site
  if ('BC' %in% dist.name) {
    message("Calculating Bray-Curtis dissimilarity...")
    dist.obj$BC <- vegan::vegdist(t(otu_tab), method = 'bray')
  }

  # Calculate Jaccard dissimilarity if requested
  # Jaccard index measures dissimilarity between sample sets, and is defined as the size of the intersection divided by the size of the union of the sample sets
  if ('Jaccard' %in% dist.name) {
    message("Calculating Jaccard dissimilarity...")
    dist.obj$Jaccard <- vegan::vegdist(t(otu_tab), method = 'jaccard')
  }

  # Calculate UniFrac distances if any of UniFrac, GUniFrac, or WUniFrac are requested
  # UniFrac incorporates phylogenetic distances between observed organisms in the computation
  if ('UniFrac' %in% dist.name || 'GUniFrac' %in% dist.name || 'WUniFrac' %in% dist.name) {
    message("Calculating UniFrac dissimilarities...")
    phy_tree <- data.obj$tree
    if (is.null(phy_tree)) {
      stop("Phylogenetic tree is required for UniFrac, GUniFrac and WUniFrac calculations.")
    }
    
    # Determine which alpha values to use based on the requested UniFrac variants
    alpha_values <- c()
    if ('UniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 1)  # Unweighted UniFrac
    }
    if ('GUniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 0)  # Generalized UniFrac
    }
    if ('WUniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 0.5)  # Weighted UniFrac
    }
    
    # Perform GUniFrac calculations for all requested alpha values
    message("Performing GUniFrac calculations...")
    unifracs <- GUniFrac(t(otu_tab), phy_tree, alpha = unique(alpha_values))$unifracs
  }

  # Assign UniFrac results to the distance object if requested
  if ('UniFrac' %in% dist.name) {
    message("Assigning UniFrac results...")
    dist.obj$UniFrac <- unifracs[, , "d_1"]
  }

  # Assign GUniFrac results to the distance object if requested
  if ('GUniFrac' %in% dist.name) {
    message("Assigning GUniFrac results...")
    dist.obj$GUniFrac <- unifracs[, , "d_0"]
  }

  # Assign WUniFrac results to the distance object if requested
  if ('WUniFrac' %in% dist.name) {
    message("Assigning WUniFrac results...")
    dist.obj$WUniFrac <- unifracs[, , "d_0.5"]
  }

  # Calculate Jensen-Shannon divergence if requested
  # Jensen-Shannon divergence is a method of measuring the similarity between two probability distributions
  if ('JS' %in% dist.name) {
    message("Calculating Jensen-Shannon divergence...")

    # Normalize OTU table to relative abundances
    otu_tab_norm <- sweep(otu_tab, 2, colSums(otu_tab), FUN = "/")

    # Define Kullback-Leibler divergence (KLD) function
    # KLD measures how one probability distribution diverges from a second, expected probability distribution
    KLD <- function(p, q) {
      non_zero <- p > 0
      sum(p[non_zero] * log(p[non_zero] / q[non_zero]))
    }

    # Define Jensen-Shannon divergence (JSD) function
    # JSD is a symmetrized and smoothed version of the KLD
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

    # Assign row and column names to the JSD matrix
    rownames(jsd_matrix) <- colnames(otu_tab)
    colnames(jsd_matrix) <- colnames(otu_tab)

    # Convert the JSD matrix to a dist object
    jsd_dist <- as.dist(jsd_matrix)
    dist.obj$JS <- jsd_dist
  }

  # Inform the user that all calculations are complete
  message("All calculations complete.")
  
  # Return the list of calculated beta diversity indices
  return(dist.obj)
}
