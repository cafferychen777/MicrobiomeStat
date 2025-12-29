#' Check if the Data Matrix Contains Count Data
#'
#' This internal function determines whether the given data matrix contains
#' non-negative integer values. It's used to verify if the data represents counts,
#' e.g., read counts in a sequencing context.
#'
#' @param data_mat A matrix containing the data to check.
#'
#' @return A logical value. Returns `TRUE` if the data matrix contains only non-negative
#'         integer values, and `FALSE` otherwise.
#'
#' @noRd
is_count_data <- function(data_mat) {
  # Check if all values in the matrix are non-negative integers
  all(data_mat == floor(data_mat) & data_mat >= 0)
}

#' Normalize a MicrobiomeStat Data Object
#'
#' This function is part of the MicrobiomeStat package. It normalizes a data object based on the chosen method.
#' @name mStat_normalize_data
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param method A string. The normalization method to be applied. It must be one of the following: "Rarefy-TSS", "Rarefy", "TSS", "GMPR", "CSS", "DESeq", "TMM". The default is "Rarefy-TSS".
#' - "Rarefy-TSS": Rarefaction followed by Total Sum Scaling normalization.
#' - "Rarefy": Rarefaction normalization only.
#' - "TSS": Total Sum Scaling normalization only.
#' - "GMPR": Geometric Mean of Pairwise Ratios normalization method.
#' - "CSS": Cumulative Sum Scaling normalization method (Paulson et al. 2013). Requires metagenomeSeq package.
#'   CSS normalizes by the cumulative sum of counts up to a data-driven quantile threshold,
#'   making it robust to high-abundance taxa and varying library sizes.
#' - "DESeq": Normalization using the DESeq method for RNA-seq data.
#' - "TMM": Normalization using the Trimmed Mean of M-values (TMM) method from the edgeR package.
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#'
#' @return A list. The normalized data object and the scale factor used for normalization.
#'
#' @examples
#' \dontrun{
#' # Load example data object
#' data(peerj32.obj)
#'
#' # Applying Total Sum Scaling (TSS) normalization
#' norm_result_tss <- mStat_normalize_data(data.obj = peerj32.obj, method = "TSS")
#' print(norm_result_tss$data.obj.norm)  # Display normalized data object
#'
#' # Applying Rarefaction followed by Total Sum Scaling (Rarefy-TSS) with a specified depth
#' norm_result_rarefy_tss <- mStat_normalize_data(data.obj = peerj32.obj,
#' method = "Rarefy-TSS", depth = 5000)
#' print(norm_result_rarefy_tss$data.obj.norm)  # Display normalized data object
#'
#' # Normalization using Geometric Mean of Pairwise Ratios (GMPR)
#' norm_result_gmpr <- mStat_normalize_data(data.obj = peerj32.obj, method = "GMPR")
#' print(norm_result_gmpr$data.obj.norm)  # Display normalized data object
#'
#' # Utilizing the DESeq normalization method
#' # This is particularly useful for RNA-seq data from microbiome studies
#' norm_result_deseq <- mStat_normalize_data(data.obj = peerj32.obj, method = "DESeq")
#' print(norm_result_deseq$data.obj.norm)  # Display normalized data object
#'
#' # Example of error handling when an incorrect depth is specified for the "Rarefy" method
#' tryCatch({
#'   norm_result_error <- mStat_normalize_data(data.obj = peerj32.obj,
#'   method = "Rarefy", depth = 10000000)
#'   print(norm_result_error$data.obj.norm)
#' }, error = function(e) {
#'   print(e$message)  # Print the error message if depth is not feasible
#' })
#'
#' }
#' @details
#' The function first checks if 'data.obj' is a list. It then retrieves the OTU table and estimates the normalization/scale factor based on the chosen method. The data object is then updated with the normalized OTU table and the chosen method is added as 'norm.status'. The function returns the normalized data object and the scale factor.
#'
#' @export
mStat_normalize_data <-
  function(data.obj,
           method = c("Rarefy-TSS", "Rarefy", "TSS", "GMPR", "CSS", "DESeq", "TMM"),
           depth = NULL) {
    # Validate input data structure
    # Ensuring the input is a list is crucial for maintaining the expected data format
    if (!is.list(data.obj)) {
      stop("data.obj should be a list.")
    }

    # Extract the OTU (Operational Taxonomic Unit) table from the input data object
    # The OTU table is the core data structure in microbiome analysis, representing taxon abundances across samples
    otu_tab <- as.data.frame(data.obj$feature.tab)

    # Validate and select the normalization method
    # This step ensures that only supported methods are used
    method <- match.arg(method)

    # Normalization process
    # Different normalization methods are applied based on the selected method
    # Each method addresses different aspects of microbiome data variability

    if (method == "Rarefy-TSS") {
      # Rarefaction followed by Total Sum Scaling
      # This method first standardizes sampling depth, then converts to relative abundance
      if (is.null(depth)) {
        depth <- min(colSums(otu_tab))
      } else if (depth > min(colSums(otu_tab))) {
        stop("Depth is greater than the smallest total count across samples.")
      }
      rarefy_depth <- ifelse(is.null(depth), min(colSums(otu_tab)), depth)
      
      # Check if data is already in relative abundance format
      if (all(round(colSums(otu_tab),5) == 1)){
        rarefied_otu_tab <- as.matrix(otu_tab)
      } else {
        # Perform rarefaction using vegan package
        rarefied_otu_tab <- t(vegan::rrarefy(t(otu_tab), sample = rarefy_depth))
      }
      # Convert to relative abundance
      rarefied_otu_tab <- rarefied_otu_tab / rarefy_depth
      scale_factor <- rarefy_depth
    } else if (method == "Rarefy") {
      # Rarefaction only
      # This method standardizes sampling depth across all samples
      if (is.null(depth)) {
        depth <- min(colSums(otu_tab))
      } else if (depth > min(colSums(otu_tab))) {
        stop("Depth is greater than the smallest total count across samples.")
      }
      rarefy_depth <- ifelse(is.null(depth), min(colSums(otu_tab)), depth)
      
      # Check if data is already in relative abundance format
      if (all(round(colSums(otu_tab),5) == 1)){
        rarefied_otu_tab <- as.matrix(otu_tab)
      } else {
        # Perform rarefaction using vegan package
        rarefied_otu_tab <- t(vegan::rrarefy(t(otu_tab), sample = rarefy_depth))
      }
      scale_factor <- rarefy_depth
    } else if (method == "TSS") {
      # Total Sum Scaling
      # This method converts counts to relative abundances
      scale_factor <- colSums(otu_tab, na.rm = TRUE)
      
      # Check for zero-sum samples and handle them gracefully
      zero_samples <- which(scale_factor == 0)
      if (length(zero_samples) > 0) {
        warning(paste("Found", length(zero_samples), "samples with zero total counts.",
                      "Setting scale factor to 1 for these samples to avoid division by zero."))
        scale_factor[zero_samples] <- 1
      }
    } else if (method == "GMPR") {
      # Geometric Mean of Pairwise Ratios
      # This method is robust to compositional effects and uneven sequencing depth
      scale_factor <- GUniFrac::GMPR(otu_tab)
    } else if (method == "CSS") {
      # Cumulative Sum Scaling (Paulson et al. 2013, Nature Methods)
      # CSS normalizes by the cumulative sum up to a data-driven quantile
      # This approach is robust to high-abundance taxa and varying library sizes

      # Check for metagenomeSeq package
      if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
        stop(
          "Package 'metagenomeSeq' required for CSS normalization.\n",
          "Install with: BiocManager::install('metagenomeSeq')"
        )
      }

      # Convert to MRexperiment object (required format for metagenomeSeq)
      mr_obj <- metagenomeSeq::newMRexperiment(counts = as.matrix(otu_tab))

      # Calculate optimal quantile threshold using data-driven approach
      # cumNormStatFast determines the quantile where features stabilize
      p <- metagenomeSeq::cumNormStatFast(mr_obj)

      # Calculate CSS normalization factors
      # Returns cumulative sum up to quantile p for each sample as a data.frame
      css_factors_df <- metagenomeSeq::calcNormFactors(mr_obj, p = p)

      # Extract normalization factors as a vector
      scale_factor <- css_factors_df$normFactors
      names(scale_factor) <- rownames(css_factors_df)

      message(paste0("CSS normalization using quantile threshold p = ", round(p, 3)))
    } else if (method == "DESeq") {
      # DESeq normalization
      # This method is particularly useful for RNA-seq data from microbiome studies
      scale_factor <- apply(otu_tab, 2, function(x) {
        positive_vals <- x[x > 0]
        # Handle edge case: if no positive values, return NA
        if (length(positive_vals) == 0) {
          return(NA_real_)
        }
        sum(x) / exp(mean(log(positive_vals)))
      })
      # Warn about samples with NA scale factors (all zeros)
      if (any(is.na(scale_factor))) {
        n_na <- sum(is.na(scale_factor))
        warning(
          n_na, " sample(s) have all zero counts and cannot be normalized with DESeq method. ",
          "These samples will have NA scale factors."
        )
      }
    } else if (method == "TMM") {
      # TMM normalization (Trimmed Mean of M-values)
      # This method is robust to compositional effects and uneven sequencing depth

      # Check for edgeR package
      if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop(
          "Package 'edgeR' required for TMM normalization.\n",
          "Install with: BiocManager::install('edgeR')"
        )
      }

      # Calculate TMM normalization factors
      scale_factor <- edgeR::calcNormFactors(otu_tab, method = "TMM")
    } else {
      stop("Invalid normalization method.")
    }

    if (method %in% c("TSS", "GMPR", "CSS", "DESeq", "TMM")) {
      # Normalize the data
      normalized_tab <- as.matrix(sweep(otu_tab, 2, scale_factor, "/"))
      
      # Check for and handle NaN/Inf values that may result from division
      if (any(is.nan(normalized_tab)) || any(is.infinite(normalized_tab))) {
        warning("NaN or Inf values detected after normalization. Setting these values to 0.")
        normalized_tab[is.nan(normalized_tab) | is.infinite(normalized_tab)] <- 0
      }
      
      data.obj.norm <- update_data_obj_count(data.obj, normalized_tab)
    } else {
      data.obj.norm <-
        update_data_obj_count(data.obj, as.matrix(rarefied_otu_tab))
    }

    # Normalize feature.agg.list if it exists
    if ('feature.agg.list' %in% names(data.obj.norm)) {
        data.obj.norm <-
          mStat_aggregate_by_taxonomy(data.obj.norm, names(data.obj.norm$feature.agg.list))
        message("The feature.agg.list has been re-aggregated based on the normalized data.")
    }

    data.obj <- data.obj.norm

    message(paste0(
      "Data has been successfully normalized using ",
      method,
      " method."
    ))

    # Return the normalized data and the scale factor
    return(list(data.obj.norm = data.obj, scale_factor = scale_factor))
  }
