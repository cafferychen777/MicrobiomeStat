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
#' # This function requires the `Biostrings` and `yaml` packages.
#' # If you encounter issues when running the example code,
#' # please ensure that you have installed and loaded these packages.
#' # library(Biostrings)
#' # library(yaml)
#' # library(biomformat)
#' #data_obj <- mStat_import_qiime2_as_data_obj(otu_qza = "path_to_otu.qza",
#' #                                             taxa_qza = "path_to_taxa.qza",
#' #                                              sam_tab = "path_to_sample_metadata",
#' #                                             refseq_qza = "path_to_refseq.qza",
#' #                                              tree_qza = "path_to_tree.qza",
#' #                                            norm.status = "Raw")
#'
#' #otuqza_file <- system.file("extdata", "table.qza", package = "microbiomeMarker")
#' #taxaqza_file <- system.file("extdata", "taxonomy.qza", package = "microbiomeMarker")
#' #sample_file <- system.file("extdata", "sample-metadata.tsv", package = "microbiomeMarker")
#' #treeqza_file <- system.file("extdata", "tree.qza", package = "microbiomeMarker")
#' #data.obj <- mStat_import_qiime2_as_data_obj(
#' #     otu_qza = otuqza_file, taxa_qza = taxaqza_file,
#' #     sam_tab = sample_file, tree_qza = treeqza_file,
#' #)
#' }
#' @seealso MicrobiomeStat documentation for further details on how to use the returned data object.
#'
#' @export
#'
mStat_import_qiime2_as_data_obj <- function(otu_qza,
                                            taxa_qza = NULL,
                                            sam_tab = NULL,
                                            refseq_qza = NULL,
                                            tree_qza = NULL) {

  message("Loading feature table...")
  feature.tab <- read_qza(otu_qza)

  if (is.null(feature.tab) || nrow(feature.tab) == 0 || ncol(feature.tab) == 0) {
    stop("Failed to read feature table or table is empty")
  }

  if (!is.null(taxa_qza)) {
    message("Loading feature annotations...")
    taxa_data <- read_qza_taxonomy(taxa_qza)
    if (!"Feature.ID" %in% colnames(taxa_data)) {
      stop("Taxonomy file must contain a 'Feature.ID' column. Available columns: ",
           paste(colnames(taxa_data), collapse = ", "))
    }
    if (!"Taxon" %in% colnames(taxa_data)) {
      stop("Taxonomy file must contain a 'Taxon' column")
    }
    
    common_features <- intersect(rownames(feature.tab), taxa_data$Feature.ID)
    if (length(common_features) == 0) {
      stop("No matching features between feature table and taxonomy data")
    }
    
    feature.tab <- feature.tab[common_features, , drop = FALSE]
    taxa_data <- taxa_data[match(common_features, taxa_data$Feature.ID), , drop = FALSE]
  } else {
    message("No feature annotations provided. Skipping...")
    taxa_data <- NULL
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
                   feature.ann = taxa_data,
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
read_qza <- function(file) {
  temp <- tempdir()
  message("Unzipping the .qza file...")
  utils::unzip(file, exdir = temp)
  
  # 列出所有文件
  all_files <- list.files(temp, recursive = TRUE, full.names = TRUE)
  message("Found files: ", paste(all_files, collapse = "\n"))
  
  # 查找 metadata.yaml
  metadata_file <- list.files(temp, pattern = "metadata.yaml$", 
                             full.names = TRUE, recursive = TRUE)[1]
  if (is.na(metadata_file)) {
    stop("Could not find metadata.yaml in qza archive")
  }
  
  metadata <- yaml::read_yaml(metadata_file)
  message("Format from metadata: ", metadata$format)
  
  # 处理 BIOMV210 格式
  if (grepl("BIOMV210", metadata$format)) {
    # 查找任意目录下的 biom 文件
    biom_file <- list.files(temp, pattern = "\\.biom$", 
                           full.names = TRUE, recursive = TRUE)[1]
    
    if (!is.na(biom_file)) {
      message("Found BIOM file: ", biom_file)
      message("Reading BIOM file...")
      return(read_q2biom(biom_file))
    } else {
      stop("Could not find BIOM file in qza archive")
    }
  }
  
  # 处理分类数据
  if (grepl("TSVTaxonomyDirectoryFormat", metadata$format)) {
    # 查找任意目录下的 taxonomy.tsv 文件
    tsv_file <- list.files(temp, pattern = "taxonomy\\.tsv$", 
                          full.names = TRUE, recursive = TRUE)[1]
    
    if (!is.na(tsv_file)) {
      message("Found taxonomy file: ", tsv_file)
      message("Reading taxonomy TSV file...")
      data <- utils::read.table(tsv_file, header = TRUE, sep = "\t", 
                               quote = "", comment.char = "")
      
      # 重命名列
      col_mapping <- c(
        "#OTU ID" = "Feature.ID",
        "OTU ID" = "Feature.ID",
        "ID" = "Feature.ID",
        "Taxonomy" = "Taxon"
      )
      
      for (old_name in names(col_mapping)) {
        if (old_name %in% colnames(data)) {
          data[[col_mapping[old_name]]] <- data[[old_name]]
          data[[old_name]] <- NULL
        }
      }
      
      if (!"Feature.ID" %in% colnames(data)) {
        stop("Could not find Feature.ID column in taxonomy file")
      }
      
      return(data)
    } else {
      stop("Could not find taxonomy.tsv file in qza archive")
    }
  }
  
  # 处理 Newick 树
  if (grepl("NewickDirectoryFormat", metadata$format)) {
    tree_file <- list.files(temp, pattern = "tree\\.nwk$", 
                           full.names = TRUE, recursive = TRUE)[1]
    if (!is.na(tree_file)) {
      message("Found tree file: ", tree_file)
      return(ape::read.tree(tree_file))
    } else {
      stop("Could not find tree.nwk file in qza archive")
    }
  }
  
  stop("Unsupported format or file not found: ", metadata$format)
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
read_q2biom <- function(biom_file) {
  message("\n=== BIOM File Reading Debug Info ===")
  message("BIOM file path: ", biom_file)
  message("File exists: ", file.exists(biom_file))
  message("File size: ", file.size(biom_file), " bytes")
  
  # Check file content type
  file_header <- readBin(biom_file, "raw", n = 10)
  message("File header (hex): ", paste(as.hexmode(file_header), collapse = " "))
  
  # Try biomformat first
  biom_data <- tryCatch({
    message("\nAttempting to read with biomformat...")
    biom <- biomformat::read_biom(biom_file)
    message("BIOM object class: ", class(biom))
    message("BIOM object structure:")
    str(biom)
    
    # Convert to matrix and then data frame
    mat <- as.matrix(biom)
    message("\nMatrix dimensions: ", nrow(mat), " x ", ncol(mat))
    message("Sample of matrix data:")
    print(mat[1:min(5, nrow(mat)), 1:min(5, ncol(mat))])
    
    df <- as.data.frame(mat)
    message("Successfully converted to data frame")
    return(df)
  }, error = function(e) {
    message("biomformat reading failed with error: ", e$message)
    message("Error class: ", class(e))
    message("Full error:")
    print(e)
    return(NULL)
  })
  
  if (!is.null(biom_data)) {
    message("Successfully read BIOM file with biomformat")
    return(biom_data)
  }
  
  # Try HDF5
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    message("\nAttempting to read as HDF5...")
    rhdf5::h5closeAll()
    
    # Check if file is HDF5
    is_hdf5 <- tryCatch({
      rhdf5::H5Fis_hdf5(biom_file)
    }, error = function(e) {
      message("HDF5 format check failed: ", e$message)
      FALSE
    })
    
    message("Is HDF5 format: ", is_hdf5)
    
    if (is_hdf5) {
      # List HDF5 structure
      message("HDF5 file structure:")
      print(rhdf5::h5ls(biom_file))
      
      hdf5_data <- tryCatch({
        mat <- rhdf5::h5read(biom_file, "observation/matrix")
        message("Matrix dimensions: ", nrow(mat), " x ", ncol(mat))
        
        obs_ids <- rhdf5::h5read(biom_file, "observation/ids")
        message("Number of observation IDs: ", length(obs_ids))
        message("Sample observation IDs: ", paste(head(obs_ids), collapse = ", "))
        
        sam_ids <- rhdf5::h5read(biom_file, "sample/ids")
        message("Number of sample IDs: ", length(sam_ids))
        message("Sample sample IDs: ", paste(head(sam_ids), collapse = ", "))
        
        df <- as.data.frame(mat)
        colnames(df) <- sam_ids
        rownames(df) <- obs_ids
        
        message("Successfully created data frame from HDF5")
        return(df)
      }, error = function(e) {
        message("HDF5 reading failed with error: ", e$message)
        message("Error class: ", class(e))
        message("Full error:")
        print(e)
        return(NULL)
      }, finally = {
        rhdf5::h5closeAll()
      })
      
      if (!is.null(hdf5_data)) {
        return(hdf5_data)
      }
    }
  }
  
  # Try reading as JSON text
  message("\nAttempting to read as JSON text...")
  json_data <- tryCatch({
    json_text <- readLines(biom_file, warn = FALSE)
    message("First few lines of file:")
    print(head(json_text))
    
    # Add JSON parsing here if needed
    
  }, error = function(e) {
    message("JSON reading failed: ", e$message)
    return(NULL)
  })
  
  message("\n=== End of BIOM File Reading Debug Info ===\n")
  stop("Failed to read BIOM file in any format")
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

read_qza_taxonomy <- function(taxa_qza) {
  # Extract taxonomy file
  unzip_dir <- tempfile()
  dir.create(unzip_dir)
  unzip(taxa_qza, exdir = unzip_dir)
  
  # Find taxonomy file
  taxa_file <- list.files(unzip_dir, 
                         pattern = "taxonomy.tsv$", 
                         recursive = TRUE, 
                         full.names = TRUE)[1]
  
  if (is.na(taxa_file)) {
    stop("Cannot find taxonomy.tsv file in the qza archive")
  }
  
  message("Found taxonomy file: ", taxa_file)
  
  # Read taxonomy data
  taxa_data <- read.table(taxa_file, 
                         sep = "\t", 
                         header = TRUE, 
                         quote = "", 
                         comment.char = "",
                         stringsAsFactors = FALSE)
  
  message("Found columns: ", paste(colnames(taxa_data), collapse = ", "))
  
  # Ensure Feature.ID column exists
  if (!"Feature.ID" %in% colnames(taxa_data)) {
    # Try common alternative names
    alt_names <- c("Feature ID", "#OTU ID", "ID", "OTU", "ASV")
    for (name in alt_names) {
      if (name %in% colnames(taxa_data)) {
        taxa_data$Feature.ID <- taxa_data[[name]]
        break
      }
    }
    
    # If still no Feature.ID, use first column
    if (!"Feature.ID" %in% colnames(taxa_data)) {
      taxa_data$Feature.ID <- taxa_data[[1]]
    }
  }
  
  # Clean up temporary directory
  unlink(unzip_dir, recursive = TRUE)
  
  return(taxa_data)
}


