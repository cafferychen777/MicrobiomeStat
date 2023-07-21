#' Load Taxonomy from a MicrobiomeStat Data Object
#'
#' This internal function loads taxonomy from a MicrobiomeStat data object.
#' It checks for the presence of various variable names in the data object
#' and returns the data associated with the first match it finds.
#'
#' @param data.obj A MicrobiomeStat data object, which is a list that contains various elements related to microbiome data.
#'
#' @return A data frame that contains the taxonomy. If no taxonomy is found,
#' the function returns NULL and displays a message.
#'
#' @details
#' The function first checks for the presence of a set of primary variable names in the MicrobiomeStat data object.
#' If it finds a match, it returns the data associated with the first match.
#' If it doesn't find a match, it checks for the presence of a set of possible confusing variable names.
#' If it finds multiple matches, it displays a message and returns the data associated with the first match.
#' If it finds a single match, it returns the associated data.
#' If it doesn't find a match, it performs a fuzzy search using a set of search words.
#' If it finds a match, it returns the associated data. If it doesn't find a match, it displays a message and returns NULL.
#'
#' @author Chen Yang
#' @seealso \code{\link[base]{intersect}}, \code{\link[base]{grep}}
#' @keywords internal
#' @noRd
load_data_obj_taxonomy <- function(data.obj) {

  # 定义优先检查的变量名列表
  primary_names <- c("feature.ann","otu.name","otu_name", "taxonomy", "tax.data", "tax_table", "tax_tab","tax.tab", "tax_info",
                     "taxa_table", "taxa.tab", "taxa_info", "species_table", "species.tab",
                     "species_info", "species.data")

  # 定义可能导致混淆的其他变量名列表
  possible_confusing_names <- c("tax", "taxa", "species", "species_data", "taxon", "taxon_data",
                                "taxon_table", "taxon.tab", "taxon_info")

  # 定义进行模糊搜索的词
  fuzzy_search_words <- c("tax", "species")

  # 查看 data.obj 中是否存在优先检查的变量名
  keys_in_data_obj <- names(data.obj)
  existing_primary_names <- intersect(primary_names, keys_in_data_obj)

  # 如果存在优先检查的变量名，则返回对应的数据
  if (length(existing_primary_names) > 0) {
    message("Using table '", existing_primary_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_primary_names[1]]])
  }

  # 否则，查看 data.obj 中是否存在其他可能导致混淆的变量名
  existing_confusing_names <- intersect(possible_confusing_names, keys_in_data_obj)

  # 如果存在多个冲突名字，发送警告消息并返回第一个冲突名字对应的数据
  if (length(existing_confusing_names) > 1) {
    message("Multiple potential taxonomy tables detected: ", paste(existing_confusing_names, collapse=", "), ". Returning the first one.")
  }

  if (length(existing_confusing_names) > 0) {
    # 返回第一个存在的冲突名字对应的数据
    message("Using table '", existing_confusing_names[1], "' in 'data.obj'.")
    return(data.obj[[existing_confusing_names[1]]])
  }

  # 如果没有精确匹配，尝试模糊搜索
  fuzzy_matches <- grep(paste(fuzzy_search_words, collapse = "|"), keys_in_data_obj)

  if (length(fuzzy_matches) > 0) {
    # 返回第一个模糊匹配的数据
    message("Using table '", names(data.obj)[fuzzy_matches[1]], "' in 'data.obj'.")
    return(data.obj[[fuzzy_matches[1]]])
  }

  # 如果不存在冲突名字，发送警告消息并返回 NULL
  message("No potential taxonomy tables detected.")
  return(NULL)
}
