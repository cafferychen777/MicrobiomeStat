#' mStat_combine_data
#'
#' This function is designed specifically for combining two MicrobiomeStat data objects,
#' assuming that they are both in 'Raw' format. Each MicrobiomeStat data object should be a list that contains
#' 'feature.tab' (a matrix of features),
#' 'meta.dat' (a data frame of metadata), and 'feature.ann' (a matrix of feature annotations).
#'
#' @param data.obj1 The first MicrobiomeStat data object to be combined.
#' @param data.obj2 The second MicrobiomeStat data object to be combined.
#'
#' @return A list with the combined 'feature.tab', 'meta.dat' and 'feature.ann'.
#'
#' @examples
#' # combined_data <- mStat_combine_data(data_obj1, data_obj2)
#'
#' @export
mStat_combine_data <- function(data.obj1, data.obj2) {

  # Inform the user about the expected format of input data objects
  message("Both data objects should be in 'Raw' format...")

  # Define an inner function to combine common data elements
  # This function is used for combining feature tables, feature annotations, and metadata
  combine_data_common <- function(data1, data2) {
    # Convert input data to data frames and add a 'feature' column from row names
    # This step ensures consistent data structure for further processing
    data1 <- data1 %>% as.data.frame() %>% rownames_to_column("feature")
    data2 <- data2 %>% as.data.frame() %>% rownames_to_column("feature")

    # Identify common features and samples between the two datasets
    # This is crucial for maintaining data integrity during combination
    common_features <- intersect(data1$feature, data2$feature)
    common_samples <- intersect(colnames(data1)[-1], colnames(data2)[-1])

    # Check for consistency in common features and samples
    # This ensures that the data for overlapping elements is identical in both datasets
    if (length(common_features) > 0 & length(common_samples) > 0) {
      data1_common <- data1 %>% filter(feature %in% common_features) %>% select(common_samples)
      data2_common <- data2 %>% filter(feature %in% common_features) %>% select(common_samples)
      if (!identical(data1_common, data2_common)) {
        stop("Inconsistent data for common features and samples.")
      } else {
        message("Consistent data for common features and samples.")
      }
    } else {
      message("No common features and samples found.")
    }

    # Identify columns to be replaced with zeros in case of missing data
    cols_to_replace <- setdiff(names(data1), "feature")

    # Combine the datasets using a series of data manipulation steps
    combined_data <- data1 %>%
      # Perform a full join to include all features from both datasets
      dplyr::full_join(data2, by = "feature") %>%
      # Replace NA values with zeros for numeric columns
      tidyr::replace_na(setNames(as.list(rep(0, length(cols_to_replace))), cols_to_replace)) %>%
      # Reshape the data to long format
      tidyr::gather(key = "sample", value = "count", -feature) %>%
      # Reshape back to wide format, filling in missing values
      tidyr::spread(key = "sample", value = "count") %>% 
      # Convert 'feature' column back to row names
      column_to_rownames("feature")

    message("Data combined.")

    # Return the combined data as a matrix
    return(as.matrix(combined_data))
  }

  # Combine feature tables from both data objects
  combined_feature_tab <- combine_data_common(data1 = data.obj1$feature.tab, data2 = data.obj2$feature.tab)
  
  # Combine feature annotations from both data objects
  combined_feature_ann <- combine_data_common(data.obj1$feature.ann, data.obj2$feature.ann)
  
  # Combine metadata from both data objects
  combined_meta_dat <- combine_data_common(data.obj1$meta.dat, data.obj2$meta.dat)

  # Inform the user that the combined data object is being returned
  message("Returning the combined data object.")
  
  # Return a list containing the combined feature table, metadata, and feature annotations
  return(list(feature.tab = combined_feature_tab, meta.dat = combined_meta_dat, feature.ann = combined_feature_ann))
}
