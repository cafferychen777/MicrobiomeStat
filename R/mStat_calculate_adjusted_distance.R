#' Calculate Adjusted Distances
#'
#' This function calculates the adjusted distances for a given set of variables,
#' using a model matrix dynamically created from the specified adjustment variables.
#'
#' @param data.obj A MicrobiomeStat data object. The heart of the MicrobiomeStat, consisting of several key components:
#'   \itemize{
#'     \item \strong{Feature.tab (Matrix)}: Meeting point for research objects (OTU/ASV/KEGG/Gene, etc.) and samples.
#'     \item \strong{Meta.dat (Data frame)}: Rows correspond to the samples, and columns serve as annotations, describing the samples.
#'     \item \strong{Feature.ann (Matrix)}: Annotations, carrying classification information like Kingdom, Phylum, etc.
#'     \item \strong{Phylogenetic tree (Optional)}: An evolutionary perspective, illuminating the relationships among various research objects.
#'     \item \strong{Feature.agg.list (Optional)}: Aggregated results based on the feature.tab and feature.ann.
#'   }
#' @param dist.obj A list of distance matrices corresponding to the names in `dist.name`.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param dist.name A character vector of names corresponding to distance matrices in `dist.obj`.
#' @return A list of adjusted distance matrices.
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' dist.obj <- mStat_calculate_beta_diversity(subset_T2D.obj, c("BC","Jaccard"))
#' adj.dist.obj <- mStat_calculate_adjusted_distance(
#' data.obj = subset_T2D.obj,
#' dist.obj = dist.obj,
#' adj.vars = c("subject_gender", "subject_race"),
#' dist.name = c("BC"))
#' data("peerj32.obj")
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, c("BC"))
#' adj.dist.obj <- mStat_calculate_adjusted_distance(
#' data.obj = peerj32.obj,
#' dist.obj = dist.obj,
#' adj.vars = c("sex"),
#' dist.name = c("BC"))
#' }
#' @export
mStat_calculate_adjusted_distance <- function (data.obj,
                                               dist.obj,
                                               adj.vars,
                                               dist.name) {
  # Message to inform the user
  message(
    "Calculating adjusted distances using the provided adjustment variables and distance matrices..."
  )

  # Load metadata and select the required adjustment variables
  meta_tab <-
    load_data_obj_metadata(data.obj) %>% select(all_of(c(adj.vars)))

  # Create a formula string from the adjustment variables and convert it into a formula object
  formula_str <- paste("~", paste(adj.vars, collapse = "+"))
  dynamic_formula <- as.formula(formula_str)

  rhs.frame <- model.frame(dynamic_formula, meta_tab, drop.unused.levels = TRUE)

  # Create a model matrix using the dynamic formula
  rhs <- model.matrix(dynamic_formula, rhs.frame)

  # 获取模型矩阵的属性“assign”，它指定了每一列是由哪个公式项生成的
  grps <- attr(rhs, "assign")

  # 对模型矩阵进行QR分解
  qrhs <- qr(rhs)

  # 使用QR分解的结果对模型矩阵进行列排序
  rhs <- rhs[, qrhs$pivot, drop = FALSE]

  # 只保留模型矩阵的秩所对应的列
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]

  # 更新“assign”属性以匹配排序后的模型矩阵
  grps <- grps[qrhs$pivot][1:qrhs$rank]

  # 获取“assign”属性中的唯一值
  u.grps <- unique(grps)

  # 计算公式右侧的唯一项的数量，减1是为了排除截距项
  nterms <- length(u.grps) - 1

  if (nterms < 1)
    stop("right-hand-side of formula has no usable terms")

  # Iterate over the distance names to calculate adjusted distances
  adj.dist.obj <- lapply(dist.name, function(sub_dist.name) {
    # Extract the corresponding distance matrix
    lhs <- dist.obj[[sub_dist.name]]

    dmat <- as.matrix(lhs^2)

    y <- lhs

    X <- rhs

    X<-X[rownames(dmat),]
    X<-as.matrix(X[,-1],nrow=nrow(X))

    # Compute the hat matrix (projection matrix)
    H <- X %*% solve(t(X) %*% X) %*% t(X)

    # Compute the matrix A
    A <- -1 / 2 * as.matrix(y) ^ 2

    # Compute the matrix J
    J <-
      diag(nrow(X)) - matrix(rep(1 / (nrow(X)), length(A)), nrow = nrow(A))

    # Compute the matrix E (adjusted distances)
    E <-
      (diag(nrow(H)) - H) %*% J %*% A %*% J %*% (diag(nrow(H)) - H)

    # Return the adjusted distances as a distance object
    return(as.dist(E))
  })

  names(adj.dist.obj) <- dist.name

  return(adj.dist.obj)
}
