#' mStat_combine_data
#'
#' This function is designed specifically for combining two MicrobiomeStat data objects,
#' assuming that they are both in 'Raw' format. Each MicrobiomeStat data object should be a list that contains
#' 'norm.status' (indicating the status of normalization, should be 'Raw'), 'feature.tab' (a matrix of features),
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
  if (!(data.obj1$norm.status == 'Raw' & data.obj2$norm.status == 'Raw')) {
    stop("This function is only applicable for data objects in 'Raw' format.")
  } else {
    message("Both data objects are in 'Raw' format, proceeding with combination...")
  }

  combine_data_common <- function(data1, data2) {
    data1 <- data1 %>% as.data.frame() %>% rownames_to_column("feature")
    data2 <- data2 %>% as.data.frame()%>% rownames_to_column("feature")

    common_features <- intersect(data1$feature, data2$feature)
    common_samples <- intersect(colnames(data1)[-1], colnames(data2)[-1])

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

    cols_to_replace <- setdiff(names(data1), "feature")

    combined_data <- data1 %>%
      full_join(data2, by = "feature") %>%
      tidyr::replace_na(setNames(as.list(rep(0, length(cols_to_replace))), cols_to_replace)) %>%
      tidyr::gather(key = "sample", value = "count", -feature) %>%
      tidyr::spread(key = "sample", value = "count") %>% column_to_rownames("feature")

    message("Data combined.")

    return(as.matrix(combined_data))
  }

  combined_feature_tab <- combine_data_common(data1 = data.obj1$feature.tab, data2 = data.obj2$feature.tab)
  combined_feature_ann <- combine_data_common(data.obj1$feature.ann, data.obj2$feature.ann)
  combined_meta_dat <- combine_data_common(data.obj1$meta.dat, data.obj2$meta.dat)

  message("Returning the combined data object.")
  return(list(feature.tab = combined_feature_tab, meta.dat = combined_meta_dat, feature.ann = combined_feature_ann))
}
