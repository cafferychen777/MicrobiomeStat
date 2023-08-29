mStat_filter <- function(x, prev.filter, abund.filter){
  # 将数据框从宽格式转换为长格式
  x_long <- as.data.frame(as.table(as.matrix(x)))

  # 计算每个taxa的平均abundance和prevalence
  filtered_taxa <- x_long %>%
    group_by(Var1) %>%  # Var1是taxa
    summarise(
      avg_abundance = mean(Freq),
      prevalence = sum(Freq > 0) / n()
    ) %>%
    filter(prevalence >= prev.filter, avg_abundance >= abund.filter) %>%
    pull("Var1")

  filtered_x <- x[filtered_taxa,]

  return(filtered_x)
}
