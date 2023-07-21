#' @title Update Sample Names in MicrobiomeStat Data Object
#'
#' @description
#' This function is part of the MicrobiomeStat package. It updates the sample names
#' in a given data object. The data object should contain three components:
#' a metadata table (`meta.dat`), a feature aggregation list (`feature.agg.list`),
#' and a feature table (`feature.tab`).
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' `meta.dat`, `feature.agg.list`, and `feature.tab` components.
#' @param new.name A character vector containing the new sample names.
#' The length of new.name should match the number of rows in `meta.dat`, and there should be no duplicates.
#'
#' @return A data object with updated sample names. The object is similar to the input object,
#' but the sample names in `meta.dat`, `feature.agg.list`, and `feature.tab` are updated to `new.name`.
#'
#' @examples
#' \dontrun{
#' # Load the required libraries
#' library(microbiome)
#' library(MicrobiomeStat)
#'
#' # Prepare data for the function
#' data(peerj32)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#'
#' # Update sample names
#' new.names <- paste0("new-", colnames(peerj32.obj$feature.tab))
#' updated_peerj32.obj <- MicrobiomeStat::update_sample_name(data.obj = peerj32.obj,
#' new.name = new.names)
#' }
#'
#' @export
mStat_update_sample_name <- function (data.obj, new.name) {
  # 检查新的样本名的数量是否与元数据表中的样本数量一致
  # 如果不一致，那么就停止执行并给出错误信息
  if (length(new.name) != nrow(data.obj$meta.dat)) stop('The number of sample do not agree!\n')

  # 检查新的样本名是否有重复
  # 如果有重复，那么就停止执行并给出错误信息
  if (length(new.name) != length(unique(new.name))) stop ('The new names have duplicates!\n')

  # 更新元数据表的行名（即样本名）
  rownames(data.obj$meta.dat) <- new.name

  # 如果data.obj中存在feature.agg.list，就进行下面的操作，否则跳过
  if ("feature.agg.list" %in% names(data.obj)) {
    # 更新feature.agg.list中每个元素的列名（即样本名）
    # 这里用到了lapply函数，它会将一个函数应用到列表的每一个元素上
    data.obj$feature.agg.list <- lapply(data.obj$feature.agg.list, function (x) {colnames(x) <- new.name; x})
  }

  # 更新OTU表的列名（即样本名）
  colnames(data.obj$feature.tab) <- new.name

  # 返回更新过后的数据对象
  return(data.obj)
}
