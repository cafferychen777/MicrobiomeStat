mStat_calculate_adjusted_distance <- function (data.obj, dist.obj, adj.vars, dist.name) {
  D.adj <- lapply(dist.name, function(dist.name){
    # 使用经典的多维缩放方法（MDS）对距离矩阵'Distance'进行分析
    # 选择的维度数为'nrow(D) - 1'
    # 但请注意：这里应该使用'Distance'而不是'D'
    D <- dist.obj[[dist.name]]

    obj <- cmdscale(D, k = D %>% as.matrix() %>% nrow() - 1, eig = TRUE)

    # 从MDS结果中提取特征值，并存储它们的符号和绝对值
    eig <- obj$eig
    s <- sign(eig)
    eig <- abs(eig)

    # 从MDS结果中提取点坐标
    xx <- obj$points

    # 使用线性模型(lm)对坐标进行回归分析，以adj.vars为协变量
    res <- residuals(lm(xx ~ data.obj$meta.dat %>% select(all_of(c(adj.vars))) %>% pull()))

    # 使用线性模型的残差计算调整后的距离
    # 正特征值和负特征值被分别处理
    # 这里的正负特征值可能与主坐标分析的解释有关
    D.adj <- dist(t(t(res[s > 0 ,]) * sqrt(eig)[s > 0])) - dist(t(t(res[, s < 0]) * sqrt(eig)[s < 0]))

    return(D.adj)
  })
}
