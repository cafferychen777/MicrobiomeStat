#' Filter a Microbiome Data Frame by Prevalence and Abundance
#'
#' This function filters taxa in a microbiome data frame based on a minimum
#' prevalence and average abundance threshold.
#'
#' @param x A matrix or data frame containing taxa (in rows) by sample (in columns)
#' microbial abundance data.
#' @param prev.filter Numeric, the minimum prevalence threshold for a taxon to be
#' retained. Prevalence is calculated as the proportion of samples where the
#' taxon is present.
#' @param abund.filter Numeric, the minimum average abundance threshold for a taxon
#' to be retained.
#'
#' @return A matrix or data frame with filtered taxa based on the specified thresholds.
#'
#' @details
#' The function first reshapes the input data to long format and then groups by taxa
#' to compute the average abundance and prevalence. After this, it filters out taxa
#' based on the provided prevalence and abundance thresholds.
#'
#' @examples
#' \dontrun{
#' data_matrix <- matrix(c(0, 3, 4, 0, 2, 7, 8, 9, 10), ncol=3)
#' colnames(data_matrix) <- c("sample1", "sample2", "sample3")
#' rownames(data_matrix) <- c("taxa1", "taxa2", "taxa3")
#'
#' filtered_data <- mStat_filter(data_matrix, 0.5, 5)
#' print(filtered_data)
#' }
#'
#' @export
mStat_filter <- function(x, prev.filter, abund.filter){
  # 将数据框从宽格式转换为长格式
  x_long <- as.data.frame(as.table(as.matrix(x)))

  # 计算每个taxa的平均abundance和prevalence
  filtered_taxa <- x_long %>%
    dplyr::group_by(Var1) %>%  # Var1是taxa
    dplyr::summarise(
      avg_abundance = mean(Freq),
      prevalence = sum(Freq > 0) / dplyr::n()
    ) %>%
    dplyr::filter(prevalence >= prev.filter, avg_abundance >= abund.filter) %>%
    dplyr::pull("Var1")

  filtered_x <- x[filtered_taxa,]

  return(filtered_x)
}
