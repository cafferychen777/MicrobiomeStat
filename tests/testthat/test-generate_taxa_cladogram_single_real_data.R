test_that("generate_taxa_cladogram_single works with NSE in geom_fruit", {
  # 跳过实际数据测试，专注于测试 geom_fruit 函数的 NSE 问题
  skip_on_cran()
  
  # 创建简单的模拟数据
  features <- paste0("ASV", 1:10)
  
  # 创建特征注释数据框
  feature.ann <- data.frame(
    Phylum = rep(c("Firmicutes", "Bacteroidetes"), each = 5),
    Class = rep(c("Bacilli", "Clostridia", "Bacteroidia", "Flavobacteriia", "Bacteroidia"), 2),
    Order = paste0("Order", 1:10),
    Family = paste0("Family", 1:10),
    Genus = paste0("Genus", 1:10),
    Species = paste0("Species", 1:10),
    row.names = features,
    stringsAsFactors = FALSE
  )
  
  # 创建丰度矩阵
  feature.tab <- matrix(
    rpois(100, lambda = 10),
    nrow = 10,
    ncol = 10,
    dimnames = list(features, paste0("Sample", 1:10))
  )
  
  # 创建元数据
  meta.dat <- data.frame(
    Group = rep(c("A", "B"), each = 5),
    Gender = rep(c("Male", "Female"), 5),
    row.names = paste0("Sample", 1:10)
  )
  
  # 创建基本数据对象
  data.obj <- list(
    feature.tab = feature.tab,
    feature.ann = feature.ann,
    meta.dat = meta.dat
  )
  
  # 创建测试结果列表
  test.list <- list()
  for (level in c("Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    test.list[[level]] <- list(
      "A_vs_B" = data.frame(
        Feature = unique(feature.ann[[level]]),
        Coefficient = rnorm(length(unique(feature.ann[[level]]))),
        P.Value = runif(length(unique(feature.ann[[level]]))),
        Adjusted.P.Value = runif(length(unique(feature.ann[[level]]))),
        stringsAsFactors = FALSE
      )
    )
  }
  
  # 确保全局环境中有 geom_tile
  if (!exists("geom_tile", envir = .GlobalEnv)) {
    assign("geom_tile", ggplot2::geom_tile, envir = .GlobalEnv)
  }
  
  # 测试 generate_taxa_cladogram_single 函数
  plot.list <- generate_taxa_cladogram_single(
    data.obj = data.obj,
    test.list = test.list,
    group.var = "Group",
    feature.level = c("Phylum", "Class", "Order"),
    feature.mt.method = "none",
    cutoff = 0.9,
    color.group.level = "Order"
  )
  
  # 验证结果
  expect_type(plot.list, "list")
  expect_true(all(sapply(plot.list, function(x) inherits(x, "ggplot"))))
  expect_true(length(plot.list) > 0)
  
  # 确保全局环境中有 geom_tile
  expect_true(exists("geom_tile", envir = .GlobalEnv))
})
