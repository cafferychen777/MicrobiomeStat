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
#' - "CSS": Cumulative Sum Scaling normalization method.
#' - "DESeq": Normalization using the DESeq method for RNA-seq data.
#' - "TMM": Normalization using the Trimmed Mean of M-values (TMM) method from the edgeR package.
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count across samples is used as the rarefaction depth.
#'
#' @return A list. The normalized data object and the scale factor used for normalization.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#' peerj32.obj <- mStat_aggregate_by_taxonomy(peerj32.obj, feature.level = c("Phylum","Family"))
#' }
#'
#' @details
#' The function first checks if 'data.obj' is a list. It then retrieves the OTU table and estimates the normalization/scale factor based on the chosen method. The data object is then updated with the normalized OTU table and the chosen method is added as 'norm.status'. The function returns the normalized data object and the scale factor.
#'
#' @export
mStat_normalize_data <-
  function(data.obj,
           method = c("Rarefy-TSS", "Rarefy", "TSS", "GMPR", "CSS", "DESeq", "TMM"),
           depth = NULL) {
    # Check if data.obj is the correct type
    if (!is.list(data.obj)) {
      stop("data.obj should be a list.")
    }

    # Get the OTU table
    otu_tab <- as.data.frame(load_data_obj_count(data.obj))

    # Check method
    method <- match.arg(method)

    # Estimate the normalization/scale factor using otu_tab
    if (method == "Rarefy-TSS") {
      if (is.null(depth)) {
        depth <- min(colSums(otu_tab))
      } else if (depth > min(colSums(otu_tab))) {
        stop("Depth is greater than the smallest total count across samples.")
      }
      rarefy_depth <-
        ifelse(is.null(depth), min(colSums(otu_tab)), depth)
      rarefied_otu_tab <-
        t(rrarefy(t(otu_tab), sample = rarefy_depth))
      rarefied_otu_tab <- rarefied_otu_tab / rarefy_depth
      scale_factor <- rarefy_depth
    } else if (method == "Rarefy") {
      if (is.null(depth)) {
        depth <- min(colSums(otu_tab))
      } else if (depth > min(colSums(otu_tab))) {
        stop("Depth is greater than the smallest total count across samples.")
      }
      rarefy_depth <-
        ifelse(is.null(depth), min(colSums(otu_tab)), depth)
      rarefied_otu_tab <-
        t(rrarefy(t(otu_tab), sample = rarefy_depth))
      scale_factor <- rarefy_depth
    } else if (method == "TSS") {
      scale_factor <- colSums(otu_tab)
    } else if (method == "GMPR") {
      scale_factor <- GMPR(otu_tab)
    } else if (method == "CSS") {
      scale_factor <- apply(otu_tab, 2, function(x) {
        sum(x) / median(x[x > 0])
      })
    } else if (method == "DESeq") {
      scale_factor <- apply(otu_tab, 2, function(x) {
        sum(x) / exp(mean(log(x[x > 0])))
      })
    } else if (method == "TMM") {
      scale_factor <- edgeR::calcNormFactors(otu_tab, method = "TMM")
    } else {
      stop("Invalid normalization method.")
    }

    if (method %in% c("TSS", "GMPR", "CSS", "DESeq", "TMM")) {
      # Normalize the data
      data.obj.norm <-
        update_data_obj_count(data.obj, sweep(otu_tab, 2, scale_factor, "/"))
    } else {
      data.obj.norm <-
        update_data_obj_count(data.obj, rarefied_otu_tab)
    }

    # Normalize feature.agg.list if it exists
    if ('feature.agg.list' %in% names(data.obj.norm)) {
        data.obj.norm <-
          mStat_aggregate_by_taxonomy(data.obj.norm, names(data.obj.norm$feature.agg.list))
      message("feature.agg.list has been normalized.")
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
