#' Calculate Beta Diversity Indices
#'
#' Calculates various beta diversity distance matrices from feature abundance data.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#'
#' @return A named list of distance matrices (class "dist").
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
  
  # Define valid distance metrics with their correct capitalization
  valid_metrics <- c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')
  
  # Create a case-insensitive mapping for validation
  valid_metrics_lower <- tolower(valid_metrics)
  names(valid_metrics_lower) <- valid_metrics
  
  # Check for invalid or incorrectly capitalized distance metrics
  invalid_metrics <- c()
  corrected_dist_name <- c()
  
  for (metric in dist.name) {
    metric_lower <- tolower(metric)
    if (metric_lower %in% tolower(valid_metrics)) {
      # Find the correct capitalization
      correct_metric <- valid_metrics[which(tolower(valid_metrics) == metric_lower)]
      
      # If the metric is not correctly capitalized, add it to the list for warning
      if (metric != correct_metric) {
        invalid_metrics <- c(invalid_metrics, metric)
        message(paste0("Warning: '", metric, "' is not correctly capitalized. Using '", correct_metric, "' instead."))
      }
      
      # Add the correctly capitalized metric to the new list
      corrected_dist_name <- c(corrected_dist_name, correct_metric)
    } else {
      # If the metric is not valid, warn and skip it
      warning(paste0("Invalid distance metric: '", metric, "'. Valid options are: ", paste(valid_metrics, collapse=", ")))
    }
  }
  
  # If no valid metrics remain after correction, stop execution
  if (length(corrected_dist_name) == 0) {
    stop("No valid distance metrics provided. Please use one or more of: ", paste(valid_metrics, collapse=", "))
  }
  
  # Replace the original dist.name with the corrected version
  dist.name <- corrected_dist_name
  
  # Extract the feature table from the data object
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
    dist.obj$BC <- mStat_attach_dist_metadata(
      vegan::vegdist(t(otu_tab), method = 'bray'),
      data.obj$meta.dat
    )
  }

  # Calculate Jaccard dissimilarity if requested
  # Jaccard index measures dissimilarity between sample sets, and is defined as the size of the intersection divided by the size of the union of the sample sets
  if ('Jaccard' %in% dist.name) {
    message("Calculating Jaccard dissimilarity...")
    dist.obj$Jaccard <- mStat_attach_dist_metadata(
      vegan::vegdist(t(otu_tab), method = 'jaccard', binary = TRUE),
      data.obj$meta.dat
    )
  }

  if ('UniFrac' %in% dist.name || 'GUniFrac' %in% dist.name || 'WUniFrac' %in% dist.name) {
    message("Calculating UniFrac dissimilarities...")
    phy_tree <- data.obj$tree
    if (is.null(phy_tree)) {
      stop("Phylogenetic tree is required for UniFrac, GUniFrac and WUniFrac calculations.")
    }
    
    # Check if the tree is binary and convert if necessary
    if (!ape::is.binary(phy_tree)) {
      message("Converting multifurcating tree to binary tree for UniFrac calculations...")
      phy_tree <- ape::multi2di(phy_tree)
      # Update the tree in the data object
      data.obj$tree <- phy_tree
    }
    
    alpha_values <- c()
    if ('GUniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 0)
    }
    if ('WUniFrac' %in% dist.name) {
      alpha_values <- c(alpha_values, 1)
    }
    if (length(alpha_values) == 0) {
      alpha_values <- 0
    }

    message("Performing GUniFrac calculations...")
    unifracs <- GUniFrac::GUniFrac(t(otu_tab), phy_tree, alpha = unique(alpha_values))$unifracs
  }

  if ('UniFrac' %in% dist.name) {
    message("Assigning UniFrac results...")
    dist.obj$UniFrac <- mStat_attach_dist_metadata(
      stats::as.dist(unifracs[, , "d_UW"]),
      data.obj$meta.dat
    )
  }

  if ('GUniFrac' %in% dist.name) {
    message("Assigning GUniFrac results...")
    dist.obj$GUniFrac <- mStat_attach_dist_metadata(
      stats::as.dist(unifracs[, , "d_0"]),
      data.obj$meta.dat
    )
  }

  if ('WUniFrac' %in% dist.name) {
    message("Assigning WUniFrac results...")
    dist.obj$WUniFrac <- mStat_attach_dist_metadata(
      stats::as.dist(unifracs[, , "d_1"]),
      data.obj$meta.dat
    )
  }

  # Calculate Jensen-Shannon divergence if requested
  # Jensen-Shannon divergence is a method of measuring the similarity between two probability distributions
  if ('JS' %in% dist.name) {
    message("Calculating Jensen-Shannon divergence...")

    sample_sums <- colSums(otu_tab, na.rm = TRUE)
    zero_samples <- which(sample_sums == 0)

    if (length(zero_samples) > 0) {
      warning(paste("Found", length(zero_samples), "samples with zero total counts.",
                    "These samples will be handled by adding a small pseudocount (1e-10)."))
      otu_tab[, zero_samples] <- otu_tab[, zero_samples] + 1e-10
      sample_sums[zero_samples] <- colSums(otu_tab[, zero_samples, drop = FALSE])
    }

    otu_tab_norm <- sweep(otu_tab, 2, sample_sums, FUN = "/")

    if (any(is.nan(otu_tab_norm)) || any(is.infinite(otu_tab_norm))) {
      warning("NaN or Inf values detected after normalization. Setting to 0.")
      otu_tab_norm[is.nan(otu_tab_norm) | is.infinite(otu_tab_norm)] <- 0
    }

    KLD <- function(p, q) {
      valid_idx <- (p > 0) & (q > 0)
      valid_idx[is.na(valid_idx)] <- FALSE

      if (!any(valid_idx)) {
        return(0)
      }

      sum(p[valid_idx] * log(p[valid_idx] / q[valid_idx]))
    }

    JSD <- function(p, q) {
      m <- (p + q) / 2
      (KLD(p, m) + KLD(q, m)) / 2
    }

    num_samples <- ncol(otu_tab_norm)
    jsd_matrix <- matrix(0, nrow = num_samples, ncol = num_samples)

    if (num_samples > 1) {
      for (i in seq_len(num_samples - 1)) {
        for (j in seq.int(i + 1, num_samples)) {
          jsd_value <- sqrt(JSD(otu_tab_norm[, i], otu_tab_norm[, j]))
          if (!is.finite(jsd_value)) {
            warning(paste("Invalid JSD value between samples", i, "and", j, ". Setting to 0."))
            jsd_value <- 0
          }
          jsd_matrix[i, j] <- jsd_value
          jsd_matrix[j, i] <- jsd_value
        }
      }
    }

    rownames(jsd_matrix) <- colnames(otu_tab)
    colnames(jsd_matrix) <- colnames(otu_tab)

    if (any(is.na(jsd_matrix)) || any(is.infinite(jsd_matrix))) {
      warning("Jensen-Shannon divergence matrix contains NA or infinite values. This may cause issues in downstream analyses.")
    }

    jsd_dist <- stats::as.dist(jsd_matrix)
    dist.obj$JS <- mStat_attach_dist_metadata(jsd_dist, data.obj$meta.dat)
  }

  # Inform the user that all calculations are complete
  message("All calculations complete.")
  
  # Return the list of calculated beta diversity indices
  return(dist.obj)
}
