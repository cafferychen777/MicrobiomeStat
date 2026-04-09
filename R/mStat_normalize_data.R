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

mStat_sanitize_scale_factor <- function(scale_factor, otu_tab, method) {
  scale_factor <- as.numeric(scale_factor)
  names(scale_factor) <- colnames(otu_tab)

  zero_samples <- colSums(otu_tab, na.rm = TRUE) == 0
  if (any(zero_samples)) {
    warning(
      sum(zero_samples),
      " sample(s) have zero total counts and cannot be meaningfully normalized with ",
      method,
      ". Using scale factor 1 to preserve zero profiles."
    )
    scale_factor[zero_samples] <- 1
  }

  invalid <- !is.finite(scale_factor) | scale_factor <= 0
  invalid[zero_samples] <- FALSE
  if (any(invalid)) {
    stop(
      "Failed to compute positive finite scale factors for ",
      method,
      " normalization."
    )
  }

  scale_factor
}

mStat_apply_scale_factor <- function(otu_tab, scale_factor) {
  normalized_tab <- as.matrix(sweep(otu_tab, 2, scale_factor, "/"))

  non_finite <- !is.finite(normalized_tab)
  if (any(non_finite)) {
    warning("NaN or Inf values detected after normalization. Setting these values to 0.")
    normalized_tab[non_finite] <- 0
  }

  normalized_tab
}

mStat_compute_deseq_size_factors <- function(otu_tab) {
  otu_mat <- as.matrix(otu_tab)
  positive_everywhere <- rowSums(otu_mat <= 0) == 0

  if (any(positive_everywhere)) {
    geo_means <- exp(rowMeans(log(otu_mat[positive_everywhere, , drop = FALSE])))
    ratio_mat <- sweep(otu_mat[positive_everywhere, , drop = FALSE], 1, geo_means, "/")
    size_factor <- apply(ratio_mat, 2, function(x) {
      stats::median(x[is.finite(x) & x > 0], na.rm = TRUE)
    })
  } else {
    warning(
      "Standard DESeq median-ratio normalization requires features that are positive in every sample. ",
      "Falling back to a positive-count geometric mean variant for sparse data."
    )

    positive_mask <- otu_mat > 0
    positive_per_feature <- rowSums(positive_mask)
    valid_features <- positive_per_feature > 0

    if (!any(valid_features)) {
      size_factor <- rep(1, ncol(otu_mat))
      names(size_factor) <- colnames(otu_mat)
      return(size_factor)
    }

    log_counts <- matrix(
      0,
      nrow = nrow(otu_mat),
      ncol = ncol(otu_mat),
      dimnames = dimnames(otu_mat)
    )
    log_counts[positive_mask] <- log(otu_mat[positive_mask])

    geo_means <- exp(rowSums(log_counts) / positive_per_feature)
    ratio_mat <- sweep(otu_mat[valid_features, , drop = FALSE], 1, geo_means[valid_features], "/")
    ratio_mat[ratio_mat <= 0] <- NA_real_

    size_factor <- apply(ratio_mat, 2, function(x) {
      stats::median(x[is.finite(x) & x > 0], na.rm = TRUE)
    })
  }

  names(size_factor) <- colnames(otu_mat)

  finite_positive <- is.finite(size_factor) & size_factor > 0
  if (any(finite_positive)) {
    size_factor[finite_positive] <- size_factor[finite_positive] /
      exp(mean(log(size_factor[finite_positive])))
  }

  size_factor
}

mStat_compute_tmm_scale_factors <- function(otu_tab) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(
      "Package 'edgeR' required for TMM normalization.\n",
      "Install with: BiocManager::install('edgeR')"
    )
  }

  dge <- edgeR::DGEList(counts = as.matrix(otu_tab))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  scale_factor <- dge$samples$lib.size * dge$samples$norm.factors
  names(scale_factor) <- colnames(otu_tab)
  scale_factor
}

#' Normalize a MicrobiomeStat Data Object
#'
#' Normalizes feature abundance data using various methods.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param method Normalization method. One of:
#'   \itemize{
#'     \item "Rarefy-TSS": Rarefaction + Total Sum Scaling (default)
#'     \item "Rarefy": Rarefaction only
#'     \item "TSS": Total Sum Scaling only
#'     \item "GMPR": Geometric Mean of Pairwise Ratios
#'     \item "CSS": Cumulative Sum Scaling (requires metagenomeSeq)
#'     \item "DESeq": DESeq normalization
#'     \item "TMM": Trimmed Mean of M-values (requires edgeR)
#'   }
#'
#' @return A list with normalized data object and scale factors.
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
#' The function first checks if 'data.obj' is a list. It then retrieves the feature table and estimates the normalization/scale factor based on the chosen method. The data object is then updated with the normalized feature table and the chosen method is added as 'norm.status'. The function returns the normalized data object and the scale factor.
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
    # The feature table is the core data structure in microbiome analysis, representing taxon abundances across samples
    otu_tab <- as.matrix(data.obj$feature.tab)

    # Validate and select the normalization method
    # This step ensures that only supported methods are used
    method <- match.arg(method)

    # Normalization process
    # Different normalization methods are applied based on the selected method
    # Each method addresses different aspects of microbiome data variability

    if (method %in% c("Rarefy-TSS", "Rarefy")) {
      # Rarefaction-based normalization methods
      # - Rarefy: Standardizes sampling depth across all samples
      # - Rarefy-TSS: Rarefaction followed by Total Sum Scaling (converts to relative abundance)

      # Validate and determine rarefaction depth
      min_depth <- min(colSums(otu_tab))

      if (is.null(depth)) {
        depth <- min_depth
      } else if (!is.numeric(depth) || length(depth) != 1 || !is.finite(depth) || depth < 0) {
        stop("Depth should be a single non-negative finite number.")
      } else if (depth > min_depth) {
        stop("Depth is greater than the smallest total count across samples.")
      }

      # Check if data is already in relative abundance format
      if (all(round(colSums(otu_tab), 5) == 1)) {
        rarefied_otu_tab <- as.matrix(otu_tab)
      } else {
        # Perform rarefaction using vegan package
        rarefied_otu_tab <- t(vegan::rrarefy(t(otu_tab), sample = depth))
      }

      # Apply TSS normalization (convert to relative abundance) for Rarefy-TSS only
      if (method == "Rarefy-TSS") {
        rarefied_otu_tab <- rarefied_otu_tab / depth
      }

      scale_factor <- depth
    } else if (method == "TSS") {
      # Total Sum Scaling
      # This method converts counts to relative abundances
      scale_factor <- mStat_sanitize_scale_factor(
        colSums(otu_tab, na.rm = TRUE),
        otu_tab,
        method
      )
    } else if (method == "GMPR") {
      # Geometric Mean of Pairwise Ratios
      # This method is robust to compositional effects and uneven sequencing depth
      scale_factor <- mStat_sanitize_scale_factor(
        GUniFrac::GMPR(otu_tab),
        otu_tab,
        method
      )
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
      scale_factor <- mStat_sanitize_scale_factor(
        setNames(css_factors_df$normFactors, rownames(css_factors_df)),
        otu_tab,
        method
      )

      message(paste0("CSS normalization using quantile threshold p = ", round(p, 3)))
    } else if (method == "DESeq") {
      # DESeq median-ratio normalization.
      # Use a sparse-data fallback when no feature is positive in every sample.
      scale_factor <- mStat_sanitize_scale_factor(
        mStat_compute_deseq_size_factors(otu_tab),
        otu_tab,
        method
      )
    } else if (method == "TMM") {
      # TMM normalization uses edgeR composition factors together with library sizes.
      scale_factor <- mStat_sanitize_scale_factor(
        mStat_compute_tmm_scale_factors(otu_tab),
        otu_tab,
        method
      )
    } else {
      stop("Invalid normalization method.")
    }

    if (method %in% c("TSS", "GMPR", "CSS", "DESeq", "TMM")) {
      # Normalize the data
      normalized_tab <- mStat_apply_scale_factor(otu_tab, scale_factor)

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
