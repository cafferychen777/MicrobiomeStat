#' Update Count Data in a MicrobiomeStat Data Object
#'
#' This internal function updates count data in a MicrobiomeStat data object.
#' It checks for the presence of various variable names in the data object
#' and updates the data associated with the first match it finds.
#' @name load_data_obj_count
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param new_count_table A new count table to replace the existing count data in the MicrobiomeStat data object.
#'
#' @return A MicrobiomeStat data object with updated count data. If no count data is found for updating,
#' the function returns the original data object and displays a message.
#'
#' @details
#' The function first checks for the presence of a set of primary variable names in the MicrobiomeStat data object.
#' If it finds a match, it updates the data associated with the first match.
#' If it doesn't find a match, it checks for the presence of a set of possible confusing variable names.
#' If it finds multiple matches, it displays a message and updates the data associated with the first match.
#' If it finds a single match, it updates the associated data.
#' If it doesn't find a match, it displays a message and returns the original data object.
#'
#' @author Chen Yang
#' @seealso \code{\link[base]{intersect}}
#' @keywords internal
#' @noRd
load_data_obj_count <- function(data.obj) {
  # 定义优先检查的变量名列表
  primary_names <-
    c(
      "feature.tab",
      "otu.tab",
      "otu_tab",
      "expression.tab",
      "exp.tab",
      "gene.counts",
      "metabolite.tab",
      "metab.tab",
      "metabolite.counts",
      "protein.tab",
      "prot.tab",
      "protein.counts",
      "transcript.tab",
      "transcript.counts",
      "rna.tab",
      "rna.counts",
      "peptide.tab",
      "peptide.counts"
    )

  # 定义可能导致混淆的其他变量名列表
  possible_confusing_names <-
    c(
      "otu.table",
      "otu_matrix",
      "otu_mat",
      "otu_data",
      "otu_list",
      "otu_counts",
      "otu_abundance",
      "exp_matrix",
      "exp_mat",
      "exp_data",
      "metab_matrix",
      "metab_mat",
      "metab_data",
      "prot_matrix",
      "prot_mat",
      "prot_data",
      "transcript_matrix",
      "transcript_mat",
      "transcript_data",
      "rna_matrix",
      "rna_mat",
      "rna_data",
      "peptide_matrix",
      "peptide_mat",
      "peptide_data"
    )

  # 定义进行模糊搜索的词
  fuzzy_search_words <- c("count", "abundance")

  # 查看 data.obj 中是否存在优先检查的变量名
  keys_in_data_obj <- names(data.obj)
  existing_primary_names <-
    intersect(primary_names, keys_in_data_obj)

  # 如果存在优先检查的变量名，则返回对应的数据
  if (length(existing_primary_names) > 0) {
    message("Using table '", existing_primary_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_primary_names[1]]])
  }

  # 否则，查看 data.obj 中是否存在其他可能导致混淆的变量名
  existing_confusing_names <-
    intersect(possible_confusing_names, keys_in_data_obj)

  # 如果存在多个冲突名字，发送警告消息并返回第一个冲突名字对应的数据
  if (length(existing_confusing_names) > 1) {
    message(
      "Multiple potential count tables detected: ",
      paste(existing_confusing_names, collapse = ", "),
      ". Returning the first one."
    )
  }

  if (length(existing_confusing_names) > 0) {
    # 返回第一个存在的冲突名字对应的数据
    message("Using table '", existing_confusing_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_confusing_names[1]]])
  }

  # 如果没有精确匹配，尝试模糊搜索
  fuzzy_matches <-
    grep(paste(fuzzy_search_words, collapse = "|"), keys_in_data_obj)

  if (length(fuzzy_matches) > 0) {
    # 返回第一个模糊匹配的数据
    message("Using table '", names(data.obj)[fuzzy_matches[1]], "' in 'data.obj'.")
    return(data.obj[[fuzzy_matches[1]]])
  }

  # 如果不存在冲突名字，发送警告消息并返回 NULL
  message("No potential count tables detected.")
  return(NULL)
}
