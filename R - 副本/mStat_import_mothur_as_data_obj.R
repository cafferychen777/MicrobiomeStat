#' Convert Mothur Data into MicrobiomeStat's Data Object
#'
#' This function serves as a bridge between Mothur and MicrobiomeStat, allowing seamless integration of Mothur data with the powerful MicrobiomeStat analysis toolkit. It accepts various Mothur microbiome study-related files and converts them into a MicrobiomeStat data object.
#'
#' MicrobiomeStat is a comprehensive and user-friendly tool designed for robust microbiome data analysis. The resulting data object from this function can be directly fed into a series of MicrobiomeStat's statistical tools for further analysis.
#'
#' @param mothur_list_file (Optional) A list file from Mothur.
#' @param mothur_group_file (Optional) A group file from Mothur.
#' @param mothur_tree_file (Optional) A tree file from Mothur.
#' @param cutoff (Optional) A cutoff value for OTU.
#' @param mothur_shared_file (Optional) A shared file from Mothur.
#' @param mothur_constaxonomy_file (Optional) A constaxonomy file from Mothur.
#' @param parseFunction (Optional) A function for parsing taxonomy. Defaults to parse_taxonomy_default.
#'
#' @return A MicrobiomeStat data object that contains the feature table, metadata, feature annotations, phylogenetic tree, and normalization status. This object is ready to be used in MicrobiomeStat's downstream analysis functions.
#'
#' @examples
#' \dontrun{
#' # library(plyr)
#' # path_to_list_file <- "/Users/apple/Downloads/esophagus/Esophagus/esophagus.fn.list"
#' # path_to_group_file <- "/Users/apple/Downloads/esophagus/Esophagus/esophagus.good.groups"
#' # path_to_tree_file <- "/Users/apple/Downloads/esophagus/Esophagus/esophagus.tree"
#' # path_to_shared_file <- "/Users/apple/Downloads/esophagus/Esophagus/esophagus.fn.shared"
#'
#' # data_obj <- mStat_import_mothur_as_data_obj(mothur_list_file = path_to_list_file,
#'                                            # mothur_group_file = path_to_group_file,
#'                                            # mothur_tree_file = path_to_tree_file,
#'                                            # mothur_shared_file = path_to_shared_file
#'                                            # )
#' }
#' @seealso MicrobiomeStat documentation for further details on how to use the returned data object.
#'
#' @export
#'
#' @author Caffery Yang
mStat_import_mothur_as_data_obj <- function(mothur_list_file = NULL, mothur_group_file = NULL,
                                            mothur_tree_file = NULL, cutoff = NULL, mothur_shared_file = NULL,
                                            mothur_constaxonomy_file = NULL, parseFunction = parse_taxonomy_default)
{
  data.obj <- list()
  if (!is.null(mothur_group_file) & !is.null(mothur_list_file)) {
    groupOTU = import_mothur_otu_table(mothur_list_file,
                                       mothur_group_file, cutoff)
    data.obj$feature.tab <- groupOTU
  } else if (!is.null(mothur_shared_file)) {
    OTUshared <- import_mothur_shared(mothur_shared_file)
    data.obj$feature.tab <- OTUshared
  }

  if (!is.null(mothur_tree_file)) {
    tree <- read_tree(mothur_tree_file)
    data.obj$tree <- tree
  }

  if (!is.null(mothur_constaxonomy_file)) {
    tax <- import_mothur_constaxonomy(mothur_constaxonomy_file,
                                      parseFunction)
    data.obj$feature.ann <- tax
  }
  return(data.obj)
}


#' Build Taxonomy Table (Internal)
#'
#' This function accepts a list of taxonomic assignments for each feature and transforms it into a matrix format. It ensures the correct matching of taxonomic ranks across different features. This function is particularly useful for creating a unified taxonomy table in microbiome studies.
#'
#' @param taxlist A list of taxonomy vectors. Each element of the list corresponds to a single feature (OTU or ASV). The names of the list should correspond to the identifiers of the features.
#'
#' @return A character matrix representing the taxonomic assignments of all features. Each row of the matrix corresponds to a feature, and each column corresponds to a specific taxonomic rank.
#'
#' @examples
#' \dontrun{
#'   # Sample list of taxonomy vectors
#'   tax_list <- list(
#'     "OTU1" = c("k__Bacteria", "p__Firmicutes", "c__Bacilli"),
#'     "OTU2" = c("k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria")
#'   )
#'
#'   # Apply function
#'   tax_table <- build_mStat_tax_table(tax_list)
#'   print(tax_table)
#' }
#'
#' @details
#' This function forms an essential part of microbiome data preprocessing in the MicrobiomeStat toolkit. It facilitates the organization of taxonomy information for downstream analyses, such as diversity and compositionality analyses, taxon significance testing, etc.
#'
#' @keywords internal
#' @noRd
build_mStat_tax_table <- function(taxlist) {
  columns = unique(unlist(lapply(taxlist, names)))
  taxmat <- matrix(NA_character_, nrow = length(taxlist), ncol = length(columns))
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

#' Import Mothur Consensus Taxonomy (Internal)
#'
#' This internal function reads in a Mothur consensus taxonomy file and applies the user-specified or default taxonomy parsing function to it. It transforms the raw taxonomy strings into a format that can be utilized by other functions in the MicrobiomeStat toolkit.
#'
#' @param mothur_constaxonomy_file A character string specifying the path to the Mothur consensus taxonomy file.
#' @param parseFunction A function for parsing the taxonomy. Default is 'parse_taxonomy_default'.
#'
#' @return A character matrix representing the taxonomic assignments of all features. Each row of the matrix corresponds to a feature, and each column corresponds to a specific taxonomic rank.
#'
#' @examples
#' \dontrun{
#'   # Replace 'path_to_your_file.txt' with your actual file path
#'   tax_table <- import_mothur_constaxonomy('path_to_your_file.txt')
#'   print(tax_table)
#' }
#'
#' @details
#' This function is an integral part of microbiome data preprocessing in the MicrobiomeStat toolkit. It allows users to import taxonomy information directly from Mothur software output, further enhancing the compatibility and ease of use of the toolkit.
#'
#' @keywords internal
#' @noRd
import_mothur_constaxonomy <- function(mothur_constaxonomy_file, parseFunction=parse_taxonomy_default) {
  read.table(mothur_constaxonomy_file)
  rawtab = read.table(mothur_constaxonomy_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)[, "Taxonomy", drop=FALSE]
  if( identical(parseFunction, parse_taxonomy_default) ){
    # Proceed with default parsing stuff.
    # Remove the confidence strings inside the parentheses, if present
    rawtab[, "Taxonomy"] = gsub("\\([[:digit:]]+\\)", "", rawtab[, "Taxonomy"])
    # Remove the quotation marks, if present
    rawtab[, "Taxonomy"] = gsub("\"", "", rawtab[, "Taxonomy"])
    # Remove trailing semicolon
    rawtab[, "Taxonomy"] = gsub(";$", "", rawtab[, "Taxonomy"])
    # Split on semicolon
    taxlist = strsplit(rawtab[, "Taxonomy"], ";", fixed=TRUE)
    taxlist = lapply(taxlist, parseFunction)
  } else {
    taxlist = lapply(rawtab[, "Taxonomy"], parseFunction)
  }
  names(taxlist) <- rownames(rawtab)
  return(build_mStat_tax_table(taxlist)) # Use new function name
}


#' Import Mothur OTU List (Internal)
#'
#' This function imports an OTU list file created by Mothur, turning it into a list object in R. The list object returned by this function requires further processing before it can be used by other phyloseq functions, hence this function is mostly intended for troubleshooting and inspection purposes rather than as a direct data import tool.
#'
#' @param mothur_list_file A character string specifying the file name or path of the Mothur OTU list file.
#' @param cutoff A character string specifying the cutoff value to be used for OTU clustering. If multiple cutoff values are available in the file, the largest value will be used by default. If only one cutoff value is present, this argument is not required.
#'
#' @return A list where each element is a character vector of one or more sequence identifiers. These identifiers show how each sequence has been grouped into OTUs by Mothur. The choice of cutoff can have a significant impact on the resulting OTU groups.
#'
#' @examples
#' \dontrun{
#'   # Please replace 'path_to_your_file.list' with your actual file path
#'   otu_list <- import_mothur_otulist('path_to_your_file.list')
#'   print(otu_list)
#' }
#'
#' @seealso
#' \code{show_mothur_cutoffs()}, \code{import_mothur()}
#'
#' @keywords internal
#' @noRd
import_mothur_otulist <- function(mothur_list_file, cutoff=NULL){
  # mothur_list_file = system.file("extdata", "esophagus.fn.list.gz", package="phyloseq")
  # cutoff = 0.04
  cutoffs = show_mothur_cutoffs(mothur_list_file)
  cutoff = select_mothur_cutoff(cutoff, cutoffs)
  # Read only the line corresponding to that cutoff
  inputline = which(cutoffs == cutoff)
  rawlines = scan(mothur_list_file, "character", sep="\t", skip=(inputline-1), nlines=1, na.strings="", quiet=TRUE)
  rawlines = rawlines[!is.na(rawlines)]
  # The first two elements are the cutoff and the number of OTUs. skip, and read to first comma for OTUnames
  OTUnames = scan(text=rawlines, what="character", comment.char=",", quiet=TRUE)[3:as.integer(rawlines[2])]
  # split each element on commas
  OTUs <- strsplit(rawlines[3:as.integer(rawlines[2])], ",", fixed=TRUE)
  # Name each OTU (currently as the first seq name in each cluster), and return the list
  names(OTUs) <- OTUnames
  # return as-is
  return(OTUs)
}

#' Import Mothur Shared File (Internal)
#'
#' This internal function reads in a shared file created by Mothur. The output is a processed OTU table that can be used for further analysis in phyloseq.
#'
#' @param mothur_shared_file A character string specifying the file name or path of the Mothur shared file.
#' @param cutoff A character string specifying the cutoff value to be used. If multiple cutoff values are available in the file, the largest value will be used by default. If only one cutoff value is present, this argument is not required.
#'
#' @return An OTU table (as an instance of the otu_table class) with rows representing OTUs and columns representing samples. This table is ready for further analysis in phyloseq.
#'
#' @examples
#' \dontrun{
#'   # Please replace 'path_to_your_file.shared' with your actual file path
#'   otu_table <- import_mothur_shared('path_to_your_file.shared')
#'   print(otu_table)
#' }
#'
#' @seealso
#' \code{show_mothur_cutoffs()}, \code{select_mothur_cutoff()}
#'
#' @keywords internal
#' @noRd
import_mothur_shared = function(mothur_shared_file, cutoff=NULL){
  #mothur_shared_file = "~/github/phyloseq/inst/extdata/esophagus.fn.shared.gz"
  # Check that cutoff is in cutoffs, or select a cutoff if none given.
  cutoffs = show_mothur_cutoffs(mothur_shared_file)
  cutoffs = cutoffs[!cutoffs %in% "label"]
  cutoff = select_mothur_cutoff(cutoff, cutoffs)
  x = readLines(mothur_shared_file)
  rawtab = read.table(text=x[grep(paste0("^", cutoff), x)], header=FALSE, row.names=2, stringsAsFactors=FALSE)[, -(1:2)]
  colnames(rawtab) <- strsplit(x[1], "\t")[[1]][4:(ncol(rawtab)+3)]
  return(otu_table(t(as.matrix(rawtab)), taxa_are_rows=TRUE))
}

#' Select or Verify Cutoff Value from Mothur Outputs (Internal)
#'
#' This internal function is used to select a cutoff value if none is provided by the user, or verify the cutoff value provided by the user. This function is used with mothur's output files that provide multiple cutoffs.
#'
#' @param cutoff The cutoff value provided by the user. If not provided, the function will select the largest non-"unique" cutoff value. If "unique" is the only cutoff value in the file, it will be selected.
#' @param cutoffs A character vector of cutoff values available in the file.
#'
#' @return A character string representing the selected or verified cutoff value.
#'
#' @examples
#' \dontrun{
#'   # Assume 'cutoffs_available' contains cutoff values from your file
#'   selected_cutoff <- select_mothur_cutoff(cutoff = NULL, cutoffs = cutoffs_available)
#'   print(selected_cutoff)
#' }
#'
#' @seealso
#' \code{show_mothur_cutoffs()}
#'
#' @keywords internal
#' @noRd
select_mothur_cutoff = function(cutoff, cutoffs){
  if( is.null(cutoff) ){
    # cutoff was NULL, need to select one.
    if( length(cutoffs) > 1 ){
      # Select the largest value, avoiding the "unique" option.
      selectCutoffs <- as(cutoffs[cutoffs != "unique"], "numeric")
      cutoff <- as.character(max(selectCutoffs))
    } else {
      # There is only one cutoff value, so use it.
      # Don't have to specify a cutoff, in this case
      cutoff <- cutoffs
    }
  } else {
    # Provided by user, non-null. Coerce to character for indexing
    cutoff <- as.character(cutoff)
    # Check that it is in set of available cutoffs.
    if( !cutoff %in% cutoffs ){
      stop("The cutoff value you provided is not among those available. Try show_mothur_cutoffs()")
    }
  }
}

#' Import Mothur List and Group Files as an OTU Table (Internal)
#'
#' This internal function imports a list and group file produced by Mothur and returns an otu_table object in R. The table contains information about the sample source of individual sequences and species/taxa abundance.
#'
#' @param mothur_list_file The path and name of the list file produced by Mothur.
#' @param mothur_group_file The path and name of the group file produced by Mothur's make.group() function.
#' @param cutoff An optional character string specifying the cutoff value to be used. If not provided, the function will select the largest value among the cutoff values contained in the list file.
#'
#' @return An otu_table object.
#'
#' @examples
#' \dontrun{
#'   # Assume 'list_file' and 'group_file' point to your Mothur output files
#'   otu_table <- import_mothur_otu_table(list_file, group_file)
#'   print(otu_table)
#' }
#'
#' @seealso
#' \code{show_mothur_cutoffs()}, \code{import_mothur()}
#'
#' @keywords internal
#' @noRd
import_mothur_otu_table <- function(mothur_list_file, mothur_group_file, cutoff=NULL){
  otulist       <- import_mothur_otulist(mothur_list_file, cutoff)
  mothur_groups <- read.table(mothur_group_file, sep="\t", as.is=TRUE, stringsAsFactors=FALSE, colClasses="character", row.names=1)
  # Initialize abundance matrix with zeros for sparse assignment
  samplenames = unique(mothur_groups[, 1])
  mothur_otu_table <- matrix(0, nrow=length(otulist), ncol=length(samplenames))
  colnames(mothur_otu_table) <- samplenames
  rownames(mothur_otu_table) <- names(otulist)

  # Write a sparse versino of the abundance table
  df = ldply(otulist, function(x){data.frame(read=x, stringsAsFactors=FALSE)})
  colnames(df)[1] <- "OTU"
  df = data.frame(df, sample=mothur_groups[df[, "read"], 1], stringsAsFactors=FALSE)
  adf = ddply(df, c("OTU", "sample"), function(x){
    # x = subset(df, OTU=="59_3_17" & sample=="C")
    data.frame(x[1, c("OTU", "sample"), drop=FALSE], abundance=nrow(x))
  })

  # Vectorized for speed using matrix indexing.
  # See help("Extract") for details about matrix indexing. Diff than 2-vec index.
  mothur_otu_table[as(adf[, c("OTU", "sample")], "matrix")] <- adf[, "abundance"]

  # Finally, return the otu_table as a phyloseq otu_table object.
  return(as.matrix(mothur_otu_table))
}

#' Read a Tree File (Internal)
#'
#' This internal function reads a tree file in either the "phylo" or NEXUS format, and returns it as an object of class "phylo". If the tree file cannot be read, the function either returns NULL or throws an error, depending on the value of 'errorIfNULL'.
#'
#' @param treefile The path and name of the tree file to be read. If an object of class "phylo" is provided, the function simply returns it.
#' @param errorIfNULL If TRUE, the function stops and throws an error when it fails to read the tree file. If FALSE (default), it simply returns NULL.
#' @param ... Further arguments passed to either read.nexus or read.tree function.
#'
#' @return An object of class "phylo" representing the tree.
#'
#' @examples
#' \dontrun{
#'   # Assume 'tree_path' points to your tree file
#'   tree <- read_tree(tree_path)
#'   print(tree)
#' }
#'
#' @keywords internal
#' @noRd
read_tree <- function (treefile, errorIfNULL = FALSE, ...)
{
  if (class(treefile)[1] %in% c("phylo")) {
    tree <- treefile
  }
  else {
    tree <- NULL
    try(tree <- read.nexus(treefile, ...), TRUE)
    if (is.null(tree))
      try(tree <- read.tree(treefile, ...), TRUE)
  }
  if (errorIfNULL & is.null(tree)) {
    stop("tree file could not be read.\nPlease retry with valid tree.")
  }

  return(tree)
}

#' Parse Taxonomy Vector (Internal)
#'
#' This function standardizes the taxonomy vector by removing leading and trailing spaces. It assigns names in the form "Rank1", "Rank2", etc., to each element of the taxonomy vector. If the input vector is empty, a warning message will be thrown.
#'
#' @param char.vec A character vector representing a taxonomy assignment for a single feature. Each element in the vector corresponds to a different taxonomic rank (from highest to lowest).
#'
#' @return A character vector representing the taxonomy assignment for a single feature with cleaned strings and assigned rank names. The names of the elements follow the "Rank1", "Rank2", etc., format.
#'
#' @examples
#' \dontrun{
#'   # Sample taxonomy vector
#'   tax_vec <- c(" k__Bacteria", "p__Firmicutes ", "c__Bacilli")
#'
#'   # Apply function
#'   cleaned_tax_vec <- parse_taxonomy_default(tax_vec)
#'   print(cleaned_tax_vec)
#' }
#'
#' @details
#' This function is primarily used to preprocess the taxonomy information obtained from metagenomic sequencing. Proper taxonomy parsing is essential for the subsequent annotation and analysis of features (e.g., Operational Taxonomic Units, OTUs, or Amplicon Sequence Variants, ASVs) in microbiome studies.
#'
#' @keywords internal
#' @noRd
parse_taxonomy_default <- function (char.vec)
{
  char.vec = gsub("^[[:space:]]{1,}", "", char.vec)
  char.vec = gsub("[[:space:]]{1,}$", "", char.vec)
  if (length(char.vec) > 0) {
    names(char.vec) = paste("Rank", 1:length(char.vec), sep = "")
  }
  else {
    warning("Empty taxonomy vector encountered.")
  }
  return(char.vec)
}

