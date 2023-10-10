#' Convert Qiime2 Data into MicrobiomeStat's Data Object
#'
#' This function serves as a bridge between Qiime2 and MicrobiomeStat, allowing seamless integration of Qiime2 data with the powerful MicrobiomeStat analysis toolkit. It accepts various Qiime2 microbiome study-related files and converts them into a MicrobiomeStat data object.
#'
#' MicrobiomeStat is a comprehensive and user-friendly tool designed for robust microbiome data analysis. The resulting data object from this function can be directly fed into a series of MicrobiomeStat's statistical tools for further analysis.
#'
#' @param otu_qza A feature table with observations (OTUs, ASVs, etc.) in .qza format from Qiime2.
#' @param taxa_qza (Optional) A taxonomy assignment table in .qza format from Qiime2.
#' @param sam_tab (Optional) A sample metadata table.
#' @param refseq_qza (Optional) A representative sequences table in .qza format from Qiime2.
#' @param tree_qza (Optional) A phylogenetic tree in .qza format from Qiime2.
#'
#' @return A MicrobiomeStat data object that contains the feature table, metadata, feature annotations, phylogenetic tree, and normalization status. This object is ready to be used in MicrobiomeStat's downstream analysis functions.
#'
#' @examples
#' \dontrun{
#' # library(Biostrings)
#' # library(yaml)
#' data_obj <- mStat_import_qiime2_as_data_obj(otu_qza = "path_to_otu.qza",
#'                                              taxa_qza = "path_to_taxa.qza",
#'                                              sam_tab = "path_to_sample_metadata",
#'                                              refseq_qza = "path_to_refseq.qza",
#'                                              tree_qza = "path_to_tree.qza",
#'                                              norm.status = "Raw")
#'
#' otuqza_file <- system.file("extdata", "table.qza", package = "microbiomeMarker")
#' taxaqza_file <- system.file("extdata", "taxonomy.qza", package = "microbiomeMarker")
#' sample_file <- system.file("extdata", "sample-metadata.tsv", package = "microbiomeMarker")
#' treeqza_file <- system.file("extdata", "tree.qza", package = "microbiomeMarker")
#' data.obj <- mStat_import_qiime2_as_data_obj(
#'     otu_qza = otuqza_file, taxa_qza = taxaqza_file,
#'     sam_tab = sample_file, tree_qza = treeqza_file,
#' )
#' }
#' @seealso MicrobiomeStat documentation for further details on how to use the returned data object.
#'
#' @export
#'
#' @author Caffery Yang
mStat_import_qiime2_as_data_obj <- function(otu_qza,
                                            taxa_qza = NULL,
                                            sam_tab = NULL,
                                            refseq_qza = NULL,
                                            tree_qza = NULL) {

  message("Loading feature table...")
  feature.tab <- read_qza(otu_qza)

  if (!is.null(taxa_qza)) {
    message("Loading feature annotations...")
    feature.ann <- read_qza(file = taxa_qza)
    feature.ann <- feature.ann %>% as.data.frame() %>% filter(Feature.ID %in% rownames(feature.tab))
    feature.ann <- parse_q2taxonomy(feature.ann)
  } else {
    message("No feature annotations provided. Skipping...")
    taxa_tab <- NULL
  }

  if (!is.null(sam_tab)) {
    message("Loading sample metadata...")
    sam_tab <- read_q2sample_meta(sam_tab)
  } else {
    message("No sample metadata provided. Skipping...")
    sam_tab <- NULL
  }

  if (!is.null(tree_qza)) {
    message("Loading phylogenetic tree...")
    tree <- read_qza(tree_qza)
  } else {
    message("No phylogenetic tree provided. Skipping...")
    tree <- NULL
  }

  # Create data.obj list
  message("Creating MicrobiomeStat data object...")
  data.obj <- list(feature.tab = feature.tab,
                   meta.dat = sam_tab,
                   feature.ann = feature.ann,
                   tree = tree)


  message("MicrobiomeStat data object created successfully!")
  return(data.obj)
}



#' Read the qza file output from qiime2 and import into MicrobiomeStat
#'
#' This function imports the Qiime2 artifacts to R, preparing them for use in MicrobiomeStat.
#' @param file A character string indicating the path of the input qza file. Currently supports files in format
#'   of 'BIOMV210DirFmt' (feature table), 'TSVTaxonomyDirectoryFormat' (taxonomic table), 'NewickDirectoryFormat'
#'   (phylogenetic tree), and 'DNASequencesDirectoryFormat' (representative sequences).
#' @param temp A character string specifying a temporary directory where the qza file will be decompressed.
#'   Default is `tempdir()`.
#' @return Depending on the input file format, the function returns different types of output suitable for
#'   further analysis in MicrobiomeStat:
#'   - For a feature table file ('BIOMV210DirFmt'), it returns a feature table.
#'   - For a taxonomic table file ('TSVTaxonomyDirectoryFormat'), it returns a taxonomic table.
#'   - For a phylogenetic tree file ('NewickDirectoryFormat'), it returns a phylogenetic tree.
#'   - For a representative sequences file ('DNASequencesDirectoryFormat'), it returns a DNA sequence set.
#' @noRd
#' @examples
#' # Please replace 'path_to_your_file.qza' with your actual file path
#' # library(yaml)
#' # library(Biostrings)
#' # read_qza('path_to_your_file.qza')
read_qza <- function(file, temp = tempdir()) {
  # Unzip the file
  message("Unzipping the .qza file...")
  unzipped_file_paths <- utils::unzip(file, exdir = temp)

  # Find metadata.yaml file in the unzipped files
  message("Locating the metadata file...")
  metadata_file_path <- grep("metadata.yaml", unzipped_file_paths, value = TRUE)

  # Read the metadata
  message("Reading the metadata...")
  metadata <- yaml::read_yaml(metadata_file_path[1])
  file_format <- metadata$format
  file_uuid <- metadata$uuid

  # Check the format of the file and process accordingly
  message("Processing the file based on its format...")
  if (grepl("BIOMV", file_format)) {
    biom_file_path <- file.path(temp, file_uuid, "data/feature-table.biom")
    processed_file <- read_q2biom(biom_file_path)
  } else if (file_format == "TSVTaxonomyDirectoryFormat") {
    taxonomy_file_path <- file.path(temp, file_uuid, "data/taxonomy.tsv")
    processed_file <- read_q2taxa(taxonomy_file_path)
  } else if (file_format == "NewickDirectoryFormat") {
    tree_file_path <- file.path(temp, file_uuid, "data/tree.nwk")
    processed_file <- read_tree(tree_file_path)
  } else if (file_format == "DNASequencesDirectoryFormat") {
    dna_sequences_file_path <- file.path(temp, file_uuid, "data/dna-sequences.fasta")
    processed_file <- readDNAStringSet(dna_sequences_file_path)
  } else {
    stop(
      "Only files in format of 'BIOMV210DirFmt' ",
      "'TSVTaxonomyDirectoryFormat', ",
      "'NewickDirectoryFormat' and 'DNASequencesDirectoryFormat'",
      " are supported."
    )
  }

  message("File processed successfully!")
  return(processed_file)
}


#' Read qiime2 feature table
#'
#' This function imports a Qiime2 feature table file and converts it into a matrix, which can be used as an input for MicrobiomeStat.
#' @param file A character string specifying the file name or path of the biom file. The file should be in Qiime2 biom format.
#' @return A matrix object containing the feature table. Rows represent features (OTUs or ASVs), and columns represent samples. This output is compatible with MicrobiomeStat functions.
#' @noRd
#' @examples
#' # Please replace 'path_to_your_file.biom' with your actual file path
#' # feature_table <- read_q2biom('path_to_your_file.biom')
read_q2biom <- function(file) {
  biomobj <- read_biom(file)
  feature_tab <- as(biom_data(biomobj), "matrix")

  return(feature_tab)
}

#' Read qiime2 taxa file
#'
#' This function imports a Qiime2 taxa file and converts it into a matrix, which can be used as an input for MicrobiomeStat.
#' @param file A character string specifying the file name or path of the taxa file. The file should be in Qiime2 taxa format.
#' @keywords internal
#' @return A matrix object containing the taxonomic annotations. Rows represent features (OTUs or ASVs), and columns represent taxonomic ranks. This output is compatible with MicrobiomeStat functions.
#' @noRd
#' @examples
#' # Please replace 'path_to_your_file.tsv' with your actual file path
#' # taxa_matrix <- read_q2taxa('path_to_your_file.tsv')
read_q2taxa <- function(file) {
  taxa <- utils::read.table(file, sep = "\t", header = TRUE)
  if ("Confidence" %in% names(taxa)) {
    taxa$Confidence <- NULL
  }
  feature.ann <- as.matrix(taxa)
  return(feature.ann)
}

#' Read qiime2 sample meta data file
#'
#' This function imports a Qiime2 sample metadata file and converts it into a dataframe, which can be used as an input for MicrobiomeStat.
#' @param file A character string specifying the file name or path of the metadata file. The file should be in Qiime2 sample metadata format.
#' @keywords internal
#' @return A dataframe object containing the sample metadata. Rows represent samples, and columns represent different metadata categories. This output is compatible with MicrobiomeStat functions.
#' @noRd
#' @examples
#' # Please replace 'path_to_your_file.tsv' with your actual file path
#' # metadata_df <- read_q2sample_meta('path_to_your_file.tsv')
read_q2sample_meta <- function(file) {
  QiimeMap <- read.table(file = file, header = TRUE,
                         sep = "\t", comment.char = "")
  rownames(QiimeMap) <- as.character(QiimeMap[, 1])
  return(as.data.frame(QiimeMap))
}

#' Parse qiime2 taxa in different taxonomic levels
#'
#' This function processes a taxonomic table from Qiime2, separating the taxonomic ranks into dplyr::distinct columns, and converts it into a matrix compatible with MicrobiomeStat.
#' @param taxa A dataframe or matrix where the rows represent features (OTUs or ASVs) and the 'Taxon' column contains taxonomic strings.
#' @param sep A character string containing a regular expression, acting as the separator between different taxonomic levels. The default is set to "; |;", which is compatible with both GreenGenes and SILVA taxonomies.
#' @param trim_rank_prefix A logical parameter that determines whether to remove the leading characters from taxonomic ranks, such as "k__" or "D_0__". Default is `TRUE`.
#' @keywords internal
#' @return A matrix object containing the parsed taxonomic annotations. Rows represent features (OTUs or ASVs), and columns represent taxonomic ranks (from Kingdom to Species). This output is compatible with MicrobiomeStat functions.
#' @noRd
#' @examples
#' # Please replace 'taxa_df' with your actual taxa dataframe or matrix
#' # parsed_taxa <- parse_q2taxonomy(taxa_df)
parse_q2taxonomy <- function(taxa, sep = "; |;", trim_rank_prefix = TRUE) {
  taxa <- data.frame(taxa)
  if (trim_rank_prefix) {
    # remove leading characters from GG
    taxa$Taxon <- gsub("[kpcofgs]__", "", taxa$Taxon)
    # remove leading characters from SILVA
    taxa$Taxon <- gsub("D_\\d__", "", taxa$Taxon)
  }

  taxa <- tidyr::separate(taxa, .data$Taxon,
                          c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                          sep = sep,
                          fill = "right"
  )
  taxa <- apply(taxa, 2, function(x) ifelse(x == "", NA, x))
  taxa <- as.data.frame(taxa)

  rownames(taxa) <- taxa$Feature.ID
  taxa$Feature.ID <- NULL

  return(as.matrix(taxa))
}


