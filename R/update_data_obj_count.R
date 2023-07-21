#' @title Update the count table in a MicrobiomeStat data object
#'
#' @description This function updates the count table within a MicrobiomeStat data object based on the provided new count table. It checks for the presence of a count table within the data object using a list of pre-defined primary and possible names. If a count table is identified, it is updated with the provided new count table. If multiple potential count tables are detected, a warning message is generated and the first one is updated.
#'
#' @param data.obj A list object specific to the MicrobiomeStat package, which contains a count table that needs to be updated.
#' @param new_count_table A matrix or data.frame representing the new count table, which will replace the existing count table in the data.obj. It must have row and column names.
#'
#' @return A MicrobiomeStat data object with the updated count table. If no count table is detected within the data object, the function returns the original data object without any modification, and a message is printed.
#'
#' @keywords internal
#' @noRd
update_data_obj_count <- function(data.obj, new_count_table) {

  # 检查新的计数表是否有效
  if (is.null(new_count_table)) {
    stop("The provided count table is NULL.")
  }

  # 检查新的计数表是否有行和列的名字
  if (is.null(rownames(new_count_table)) || is.null(colnames(new_count_table))) {
    stop("The provided count table must have row and column names.")
  }

  # 定义优先检查的变量名列表
  primary_names <- c("feature.tab", "otu.tab", "otu_tab", "expression.tab", "exp.tab", "gene.counts",
                     "metabolite.tab", "metab.tab", "metabolite.counts", "protein.tab", "prot.tab", "protein.counts",
                     "count", "counts", "abundance")

  # 定义可能导致混淆的其他变量名列表
  possible_confusing_names <- c("otu.table", "otu_matrix", "otu_mat", "otu_data", "otu_list", "otu_counts",
                                "otu_abundance", "exp_matrix", "exp_mat", "exp_data", "metab_matrix",
                                "metab_mat", "metab_data", "prot_matrix", "prot_mat", "prot_data")

  # 查看 data.obj 中是否存在优先检查的变量名
  keys_in_data_obj <- names(data.obj)
  existing_primary_names <- intersect(primary_names, keys_in_data_obj)

  # 如果存在优先检查的变量名，则更新对应的数据
  if (length(existing_primary_names) > 0) {
    data.obj[[existing_primary_names[1]]] <- new_count_table
    return(data.obj)
  }

  # 否则，查看 data.obj 中是否存在其他可能导致混淆的变量名
  existing_confusing_names <- intersect(possible_confusing_names, keys_in_data_obj)

  # 如果存在多个冲突名字，发送警告消息并更新第一个冲突名字对应的数据
  if (length(existing_confusing_names) > 1) {
    message("Multiple potential count tables detected: ", paste(existing_confusing_names, collapse=", "), ". Updating the first one.")
  }

  if (length(existing_confusing_names) > 0) {
    # 更新第一个存在的冲突名字对应的数据
    data.obj[[existing_confusing_names[1]]] <- new_count_table
    return(data.obj)
  }

  # 如果不存在冲突名字，发送警告消息并返回 NULL
  message("No potential count tables detected. The count table was not updated.")
  return(data.obj)
}
