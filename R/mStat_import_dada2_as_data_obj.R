#' Convert DADA2 Data into MicrobiomeStat's Data Object
#'
#' The `mStat_import_dada2_as_data_obj` function acts as a bridge between DADA2 output and MicrobiomeStat, effectively enabling smooth data transition between the two. It accepts essential DADA2 outputs and converts them into a MicrobiomeStat compatible format for subsequent microbiome data analysis.
#' @name mStat_import_dada2_as_data_obj
#' @param seq_tab A sequence table from DADA2, where rows represent samples, and columns represent sequences (features). This table is crucial for subsequent MicrobiomeStat analyses.
#' @param tax_tab (Optional) A taxonomy assignment table derived from DADA2. Each row corresponds to a feature (sequence), and each column corresponds to a different taxonomic rank. This table enriches the feature information.
#' @param sam_tab (Optional) A sample metadata table. Each row corresponds to a sample, and each column provides different metadata attributes. This table enriches the sample information.
#' @param phy_tree (Optional) A phylogenetic tree of class 'phylo'. This tree represents the phylogenetic relationships between features (sequences) and can be incorporated into certain analyses.
#'
#' @return A MicrobiomeStat data object, containing the following elements:
#'   \itemize{
#'     \item feature.tab: A matrix of the feature (sequence) table.
#'     \item meta.dat: A data frame of the sample metadata.
#'     \item feature.ann: A matrix of the feature (sequence) taxonomy annotations.
#'     \item tree: A phylogenetic tree, if provided.
#'   }
#' This object is readily accepted by MicrobiomeStat's downstream analysis functions.
#'
#'
#'
#' @examples
#' \dontrun{
#'  # library(Biostrings)
#'  # seq_tab <- readRDS(
#'  #   system.file(
#'  #       "extdata", "dada2_seqtab.rds",
#'  #       package = "microbiomeMarker"
#'  #   )
#'  # )
#'  # tax_tab <- readRDS(
#'  #   system.file(
#'  #       "extdata", "dada2_taxtab.rds",
#'  #       package = "microbiomeMarker"
#'  #   )
#'  # )
#'  # sam_tab <- read.table(
#'  #   system.file(
#'  #       "extdata", "dada2_samdata.txt",
#'  #       package = "microbiomeMarker"
#'  #   ),
#'  #   sep = "\t",
#'  #   header = TRUE,
#'  #   row.names = 1
#'  # )
#'  # data_obj <- mStat_import_dada2_as_data_obj(seq_tab = seq_tab,
#'  # tax_tab = tax_tab, sam_tab = sam_tab)
#'
#' }
#'
#' @export
mStat_import_dada2_as_data_obj <- function(seq_tab,
                                           tax_tab = NULL,
                                           sam_tab = NULL,
                                           phy_tree = NULL) {

  # refseq
  refseq <- colnames(seq_tab)
  # set refseq and taxa names to ASV_1, ASV_2,...
  refseq_nm <- paste0("ASV", seq_along(refseq))
  colnames(seq_tab) <- refseq_nm
  names(refseq) <- refseq_nm

  if (!is.null(tax_tab)) {
    if (!identical(refseq_nm, row.names(tax_tab))) {
      tax_tab <- tax_tab[match(refseq, row.names(tax_tab)), ,
                         drop = FALSE]
    }
    row.names(tax_tab) <- refseq_nm
  }

  # refseq to XStringSet
  refseq <- DNAStringSet(refseq)

  if (!is.null(phy_tree) && inherits(phy_tree, "character")) {
    phy_tree <- read_tree(phy_tree)
  }

  asv_tab <- seq_tab

  data.obj <- list()
  data.obj$feature.tab <- t(asv_tab)
  data.obj$meta.dat <- sam_tab
  data.obj$feature.ann <- tax_tab
  data.obj$tree <- phy_tree

  return(data.obj)
}
