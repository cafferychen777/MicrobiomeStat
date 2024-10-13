#' @title Aggregate Data Function
#' @description This function aggregates data by a specified subject and optional strata.
#' @name aggregate_data
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character specifying the name of the subject variable.
#' @param strata.var A character specifying the name of the strata variable. If NULL, no stratification is performed.
#'
#' @return A new list object with aggregated data.
#'
#'
#' @examples
#' \dontrun{
#' # Prepare data for the function
#' data(peerj32.obj)
#'
#' # Call the function
#' aggregated_data <- mStat_aggregate_data(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   strata.var = NULL
#' )
#' }
#' @export
mStat_aggregate_data <- function (data.obj, subject.var, strata.var = NULL) {

  message("Aggregate data by a factor ...")

  # Aggregate pre-aggregated feature data if it exists
  # This step handles cases where data has been pre-aggregated at different taxonomic levels
  if (!is.null(data.obj$feature.agg.list)) {
    feature.agg.list <- list()
    for (name in names(data.obj$feature.agg.list)) {
      abund <- data.obj$feature.agg.list[[name]]
      
      # Aggregate data by subject only if no strata variable is provided
      if (is.null(strata.var)) {
        obj <- dplyr::group_by(as.data.frame(t(abund)), data.obj$meta.dat[, subject.var]) %>%
          dplyr::summarise(dplyr::across(everything(), list(~sum(.x, na.rm = TRUE))))
        feature.agg.list[[name]] <- t(as.matrix(obj[, -1]))
      } else {
        # Aggregate data by both subject and strata if a strata variable is provided
        if (!is.null(data.obj$meta.dat[, strata.var])) {
          obj <- dplyr::group_by(as.data.frame(t(abund)), data.obj$meta.dat[, subject.var], data.obj$meta.dat[, strata.var]) %>%
            dplyr::summarise(dplyr::across(everything(), list(~sum(.x, na.rm = TRUE))))
          feature.agg.list[[name]] <- t(as.matrix(obj[, -(1:2)]))
          colnames(feature.agg.list[[name]]) <- paste(obj[, 1], obj[, 2], sep="_")
        } else {
          stop("Input data object does not contain strata variable.")
        }
      }
    }
    data.obj.new <- list(feature.agg.list = feature.agg.list)
  } else {
    data.obj.new <- list()
  }

  # Aggregate the main feature table
  # This step is crucial for combining data from the same subject or subject-strata combination
  if (!("feature.tab" %in% names(data.obj)) || !("subject" %in% names(data.obj$meta.dat))) {
    stop("Input data object does not contain feature.tab or subject variable.")
  }
  abund <- data.obj$feature.tab
  if (is.null(strata.var)) {
    # Aggregate by subject only
    obj <- dplyr::group_by(as.data.frame(t(abund)), data.obj$meta.dat %>% select(all_of(subject.var)) %>% dplyr::pull()) %>%
      dplyr::summarise(dplyr::across(everything(), list(~sum(.x, na.rm = TRUE))))
    feature.tab <- t(as.matrix(obj[, -1]))
    rownames(feature.tab) <- sub("\\._1$", "", rownames(feature.tab))
    colnames(feature.tab) <- obj[, 1] %>% dplyr::pull()
  } else {
    # Aggregate by both subject and strata
    if (!is.null(data.obj$meta.dat[, strata.var])) {
      obj <- dplyr::group_by(as.data.frame(t(abund)), data.obj$meta.dat[, subject.var], data.obj$meta.dat[, strata.var]) %>%
        dplyr::summarise(dplyr::across(everything(), list(~sum(.x, na.rm = TRUE))))
      feature.tab <- t(as.matrix(obj[, -(1:2)]))
      colnames(feature.tab) <- paste(obj[, 1] %>% dplyr::pull(), obj[, 2] %>% dplyr::pull(), sep="_")
      rownames(feature.tab) <- sub("\\_1$", "", rownames(feature.tab))
    } else {
      stop("Input data object does not contain strata variable.")
    }
  }
  data.obj.new$feature.tab <- feature.tab

  # Preserve tree and feature annotation data
  # These elements are not affected by the aggregation process
  data.obj.new$tree <- data.obj$tree
  data.obj.new$feature.ann <- data.obj$feature.ann

  # Aggregate the metadata
  # This step ensures that the metadata matches the aggregated feature data
  if (is.null(strata.var)) {
    # Aggregate metadata by subject only
    meta.dat <- dplyr::group_by(data.obj$meta.dat, data.obj$meta.dat[, subject.var]) %>%
      dplyr::summarise(dplyr::across(everything(), ~.[1])) %>%
      dplyr::ungroup()
    colnames(meta.dat)[1] <- c(subject.var)
  } else {
    # Aggregate metadata by both subject and strata
    if (!is.null(data.obj$meta.dat[, strata.var])) {
      meta.dat <- data.obj$meta.dat %>% dplyr::group_by(data.obj$meta.dat[, subject.var], data.obj$meta.dat[, strata.var]) %>%
        dplyr::summarise(dplyr::across(everything(), ~.[1])) %>%
        dplyr::ungroup() %>% select(-all_of(c(subject.var,strata.var)))
      colnames(meta.dat)[1:2] <- c(subject.var,strata.var)
    } else {
      stop("Input data object does not contain strata variable.")
    }
  }
  # Ensure meta.dat is a data.frame for consistency
  meta.dat <- as.data.frame(meta.dat)

  # Align metadata with aggregated feature data
  # This step ensures that the row names of the metadata match the column names of the feature table
  rownames(meta.dat) <- colnames(feature.tab)
  data.obj.new$meta.dat <- meta.dat

  message("Finished!")
  return(data.obj.new)
}
