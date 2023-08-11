#' Import a BIOM formatted data file into an mStat data object
#'
#' The `mStat_import_biom_as_data_obj` function facilitates the conversion of a BIOM (Biological Observation Matrix) file to an mStat data object. This includes importing and parsing the relevant metadata, taxonomic information, and optionally a phylogenetic tree.
#'
#' @param BIOMfilename A character string specifying the path to the BIOM file, or a 'biom-class' object. This file contains the essential biological observation matrix and related metadata.
#' @param treefilename (Optional) A character string specifying the path to the phylogenetic tree file, or a 'phylo' object. This file represents the phylogenetic tree associated with the biological observation. If not provided, the function continues without incorporating a phylogenetic tree.
#' @param parseFunction A user-specified function for parsing the taxonomy from the metadata in the BIOM file. Default function is 'parse_taxonomy_default'.
#' @param ... Additional arguments passed to 'read_tree' function when reading the phylogenetic tree file.
#'
#' @return A list representing an mStat data object. The list includes the following elements:
#'   \itemize{
#'     \item feature.tab: A matrix that contains the biological observations (OTU table).
#'     \item feature.ann: A matrix that contains parsed taxonomic annotations for each feature.
#'     \item meta.dat: A data frame that contains metadata for each sample.
#'     \item tree: A phylogenetic tree, if provided.
#'   }
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   #rich_dense_biom <- system.file("extdata", "rich_dense_otu_table.biom",  package="phyloseq")
#'   #treefilename <- system.file("extdata", "biom-tree.phy",  package="phyloseq")
#'   #refseqfilename <- system.file("extdata", "biom-refseq.fasta",  package="phyloseq")
#'   #data.obj <- mStat_import_biom_as_data_obj(rich_dense_biom,
#'   #treefilename, refseqfilename, parseFunction=parse_taxonomy_greengenes)
#' }
#'
#' @details
#' This function reads a BIOM file, which typically contains a biological observation matrix and its related metadata. It also optionally accepts a phylogenetic tree file. It creates a list object that follows the mStat data object structure for subsequent analysis. If any of the required elements are missing in the BIOM file, warnings will be issued, but the function will continue to run. Please ensure the BIOM file and optional phylogenetic tree file are properly formatted for successful conversion.
#'
#' @author Jun Chen
#' @references The Biological Observation Matrix (BIOM) format or: how I learned to stop worrying and love the ome-ome. Daniel McDonald et al. GigaScience 2012.
mStat_import_biom_as_data_obj <-
  function (BIOMfilename,
            treefilename = NULL,
            parseFunction = parse_taxonomy_default,
            ...)
  {
    data.obj <- list()
    if (inherits(BIOMfilename, "character")) {
      x = read_biom(biom_file = BIOMfilename)
    }
    else if (inherits(BIOMfilename, "biom")) {
      x = BIOMfilename
    }
    else {
      stop("import_biom requires a 'character' string to a biom file or a 'biom-class' object")
    }

    data.obj$feature.tab = as(biom_data(x), "matrix")

    if (all(sapply(sapply(x$rows, function(i) {
      i$metadata
    }), is.null))) {
      taxtab <- NULL
    }
    else {
      taxlist = lapply(x$rows, function(i) {
        parseFunction(i$metadata$taxonomy)
      })
      names(taxlist) = sapply(x$rows, function(i) {
        i$id
      })
      taxtab = build_mStat_tax_table(taxlist)
    }
    data.obj$feature.ann <- taxtab

    if (is.null(sample_metadata(x))) {
      samdata <- NULL
    }
    else {
      samdata = sample_metadata(x)
    }

    data.obj$meta.dat <- samdata

    if (!is.null(treefilename)) {
      if (inherits(treefilename, "phylo")) {
        tree <- treefilename
      }
      else {
        tree <- read_tree(treefilename, ...)
      }
      if (is.null(tree)) {
        warning("treefilename failed import. It not included.")
      }
      else {
        data.obj$tree <- tree
      }
    }

    return(data.obj)
  }

#' Build mStat tax table
#'
#' This function takes a list of taxonomy information and creates a taxonomic table compatible with mStat.
#'
#' @param taxlist A list where each element is a named vector representing the taxonomy of a feature.
#'
#' @return A matrix representing the taxonomic table.
#' @export
#' @examples
#' # Assume `taxlist` is a list of taxonomy information
#' # tax_table <- build_mStat_tax_table(taxlist)
build_mStat_tax_table <- function(taxlist) {
  columns = unique(unlist(lapply(taxlist, names)))
  taxmat <-
    matrix(NA_character_,
           nrow = length(taxlist),
           ncol = length(columns))
  colnames(taxmat) = columns
  for (i in 1:length(taxlist)) {
    if (length(taxlist[[i]]) > 0) {
      taxmat[i, names(taxlist[[i]])] <- taxlist[[i]]
    }
  }
  taxmat[taxmat == ""] <- NA_character_
  taxmat <- as(taxmat, "matrix")
  rownames(taxmat) = names(taxlist)
  return(taxmat)
}
