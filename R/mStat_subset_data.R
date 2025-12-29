#' @title Subset Data Object by Sample IDs in MicrobiomeStat
#'
#' @description
#' This function is part of the MicrobiomeStat package. It subsets a multi-omics data object by
#' a specified set of sample IDs. It checks for the presence of various data components in the
#' data object and subsets them if they exist.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' table (`feature.tab`), feature names, a full feature name list, and a feature aggregation list
#' (`feature.agg.list`).
#' @param samIDs A vector of sample IDs to subset by. This can be a logical vector, a numeric
#' vector, or a character vector of sample IDs.
#' @param condition A logical expression string describing rows to keep. This allows subsetting samples based on metadata, similarly to dplyr::filter. Only one of 'samIDs' and 'condition' should be provided.
#' @return A list that is the modified version of the input data object. It includes the subsetted
#' metadata, feature table, feature names, full feature name list, and feature aggregation list.
#'
#' @examples
#' \dontrun{
#' # Load the required libraries
#' library(MicrobiomeStat)
#' # Prepare data for the function
#' data(peerj32.obj)
#' peerj32.obj$meta.dat <- peerj32.obj$meta.dat %>%
#' dplyr::select(all_of("subject")) %>% dplyr::distinct() %>%
#' dplyr::mutate(cons = runif(dplyr::n(),0,5)) %>%
#' dplyr::left_join(peerj32.obj$meta.dat %>% rownames_to_column("sample"),by = "subject") %>%
#' tibble::column_to_rownames("sample")
#' # Example 1: Subset data for a specific time point
#' # Subset to include only samples from time point 1
#' subset_time1 <- mStat_subset_data(data.obj = peerj32.obj, condition = "time == '1'")
#' # Example 2: Subset data for a specific group
#' # Subset to include only samples from the 'LGG' group
#' subset_LGG <- mStat_subset_data(data.obj = peerj32.obj, condition = "group == 'LGG'")
#' # Example 3: Subset data for a specific sex
#' # Subset to include only male samples
#' subset_male <- mStat_subset_data(data.obj = peerj32.obj, condition = "sex == 'male'")
#' # Example 4: Complex condition subsetting
#' # Subset based on multiple conditions: male samples from 'Placebo' group at time point 2
#' complex_condition <- "sex == 'male' & group == 'Placebo' & time == '2'"
#' subset_complex <- mStat_subset_data(data.obj = peerj32.obj, condition = complex_condition)
#' # Example 5: Subset data using a combination of sample IDs and condition
#' # First, get a subset of sample IDs (e.g., first 10 samples)
#' subset_ids <- rownames(peerj32.obj$meta.dat)[1:10]
#' # Subset the data object based on these sample IDs
#' subset_by_ids <- mStat_subset_data(data.obj = peerj32.obj, samIDs = subset_ids)
#' # Then, further subset the result to include only those samples from the 'Placebo' group
#' subset_ids_condition <- mStat_subset_data(data.obj = subset_by_ids,
#' condition = "group == 'Placebo'")
#' # Example 6: Subset data for female samples in the 'Placebo' group
#' # Subset to include only female samples from the 'Placebo' group
#' subset_female_placebo <- mStat_subset_data(data.obj = peerj32.obj,
#' condition = "sex == 'female' & group == 'Placebo'")
#' # Example 7: Subset data excluding certain subjects
#' # Subset to exclude subjects S1 and S2
#' subset_exclude_subjects <- mStat_subset_data(data.obj = peerj32.obj,
#' condition = "!subject %in% c('S1', 'S2')")
#' # Example 8: Subset data for a specific range of the 'cons' variable
#' # Subset to include samples with 'cons' values greater than 2
#' subset_cons_gt_2 <- mStat_subset_data(data.obj = peerj32.obj, condition = "cons > 2")
#' # Example 9: Subset data for samples with even-numbered IDs
#' # Subset to include samples with even-numbered IDs
#' even_sample_ids <- rownames(peerj32.obj$meta.dat)[as.integer(gsub('sample-', '',
#' rownames(peerj32.obj$meta.dat))) %% 2 == 0]
#' subset_even_samples <- mStat_subset_data(data.obj = peerj32.obj, samIDs = even_sample_ids)
#' # Example 10: Combine multiple conditions with logical operators
#' # Subset to include male samples from 'LGG' group at time point 1 with 'cons' less than 3
#' complex_multiple_conditions <- "sex == 'male' & group == 'LGG' & time == '1' & cons < 3"
#' subset_complex_multiple <- mStat_subset_data(data.obj = peerj32.obj,
#' condition = complex_multiple_conditions)
#' }
#' @details
#' The function first checks if samIDs is logical or numeric, and if so, converts it to a
#' character vector of sample IDs. It then subsets the metadata by the sample IDs. If a feature
#' table exists, it subsets the feature table and the feature names by the sample IDs. If a full
#' feature name list exists, it subsets the full feature name list by the sample IDs. If a feature
#' aggregation list exists, it subsets each feature aggregation table in the list by the sample IDs.
#' The function returns the subsetted data object.
#'
#' @seealso \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr filter select
mStat_subset_data <- function (data.obj, samIDs = NULL, condition = NULL) {

  # Store the original sample names
  original_samIDs <- rownames(data.obj$meta.dat)

  # Check that not both samIDs and condition are provided
  if (!is.null(samIDs) & !is.null(condition)) {
    stop("Only one of 'samIDs' and 'condition' should be provided.")
  }

  if (!is.null(condition)) {
    # If condition is provided, filter meta.dat by condition
    data.obj$meta.dat <- dplyr::filter(data.obj$meta.dat, !!rlang::parse_expr(condition))
    samIDs <- rownames(data.obj$meta.dat)
    message("Data has been subsetted based on the provided condition.")
  } else {
    # If samIDs is logical or numeric, convert it to character form of sample IDs
    if (is.logical(samIDs) | is.numeric(samIDs)) {
      samIDs <- rownames(data.obj$meta.dat)[samIDs]
    }
    message("Data has been subsetted based on the provided samIDs.")
  }

  # Extract subset of metadata that matches samIDs
  data.obj$meta.dat <- dplyr::filter(data.obj$meta.dat, rownames(data.obj$meta.dat) %in% samIDs)
  message("Updated metadata to match the subsetted data.")

  # Get the sample names that were excluded
  excluded_samIDs <- setdiff(original_samIDs, samIDs)
  message(paste("The following samples were excluded:", paste(excluded_samIDs, collapse = ", ")))

  # If feature table exists, extract subset of feature table that matches samIDs
  if (!is.null(data.obj$feature.tab)) {
    # Extract sample subset from feature table
    data.obj$feature.tab <- data.obj$feature.tab %>%
      as.data.frame() %>%
      dplyr::select(all_of(samIDs)) %>%
      dplyr::filter(rowSums(data.obj$feature.tab) != 0) %>%
      as.matrix()

    message("Updated feature table to match the subsetted data.")

  }

  # Extract corresponding feature annotation subset
  if (!is.null(data.obj$feature.tab) & !is.null(data.obj$feature.ann)) {
    data.obj$feature.ann <- data.obj$feature.ann %>%
      as.data.frame() %>%
      dplyr::filter(rownames(data.obj$feature.ann) %in% rownames(data.obj$feature.tab)) %>%
      as.matrix()

    message("Updated feature annotation to match the subsetted data.")
  }

  # If feature aggregation list exists, perform the same subset extraction operation for each feature aggregation table
  if (!is.null(data.obj$feature.agg.list)) {
    data.obj$feature.agg.list <- lapply(data.obj$feature.agg.list, function(x) {
      x %>% as.data.frame() %>%
        dplyr::select(all_of(samIDs)) %>%
        dplyr::filter(rowSums(x) != 0) %>%
        as.matrix()
    })

    message("Updated feature aggregation list to match the subsetted data.")
  }

  message("Data subsetting complete. Returning updated data object.")
  # Return the processed data object
  return(data.obj)
}
