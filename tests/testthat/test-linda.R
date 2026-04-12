test_that("linda normalizes zero.handling choice and keeps output shape", {
  set.seed(42)

  feature.tab <- matrix(
    c(10, 0, 5, 2,
      0, 3, 1, 4,
      8, 1, 0, 2),
    nrow = 3,
    dimnames = list(paste0("Feature", 1:3), paste0("Sample", 1:4))
  )
  meta.dat <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  lower <- linda(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    zero.handling = "imputation",
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )
  mixed <- linda(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    zero.handling = "Imputation",
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )

  expect_equal(lower$variables, mixed$variables)
  expect_equal(lower$output[[1]]$log2FoldChange, mixed$output[[1]]$log2FoldChange)
  expect_equal(lower$output[[1]]$pvalue, mixed$output[[1]]$pvalue)
})

test_that("linda2 normalizes zero.handling choice and keeps output shape", {
  set.seed(42)

  feature.tab <- matrix(
    c(10, 0, 5, 2,
      0, 3, 1, 4,
      8, 1, 0, 2),
    nrow = 3,
    dimnames = list(paste0("Feature", 1:3), paste0("Sample", 1:4))
  )
  meta.dat <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  lower <- linda2(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    zero.handling = "pseudo-count",
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )
  mixed <- linda2(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    zero.handling = "Pseudo-count",
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )

  expect_equal(lower$variables, mixed$variables)
  expect_equal(lower$output[[1]]$log2FoldChange, mixed$output[[1]]$log2FoldChange)
  expect_equal(lower$output[[1]]$pvalue, mixed$output[[1]]$pvalue)
})

test_that("linda and linda2 preserve matrix shape after filtering to one feature", {
  feature.tab <- matrix(
    c(
      5, 5, 5, 5,
      0, 0, 0, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Feature1", "Feature2"), paste0("Sample", 1:4))
  )
  meta.dat <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  linda_res <- linda(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    prev.filter = 0.5,
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )

  linda2_res <- linda2(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    prev.filter = 0.5,
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )

  expect_equal(dim(linda_res$feature.dat.use), c(1, 4))
  expect_equal(dim(linda2_res$feature.dat.use), c(1, 4))
  expect_true(length(linda_res$output) > 0)
  expect_true(length(linda2_res$output) > 0)
})


test_that("linda helpers reuse shared sample-total normalization when filtering features", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  linda_src <- paste(readLines(file.path(package_root, "R/linda.R"), warn = FALSE), collapse = "\n")
  linda2_src <- paste(readLines(file.path(package_root, "R/linda2.R"), warn = FALSE), collapse = "\n")

  expect_match(linda_src, "mStat_normalize_feature_matrix_by_sample_total")
  expect_match(linda2_src, "mStat_normalize_feature_matrix_by_sample_total")
  expect_no_match(linda_src, "temp <- t\\(t\\(Y\\) / colSums\\(Y\\)\\)")
  expect_no_match(linda2_src, "temp <- t\\(t\\(Y\\) / colSums\\(Y\\)\\)")
})

test_that("linda2 returns structured empty result when all features are filtered out", {
  feature.tab <- matrix(
    c(
      1, 0, 0, 1,
      0, 1, 1, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Feature1", "Feature2"), paste0("Sample", 1:4))
  )
  meta.dat <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  expect_warning(
    linda2_res <- linda2(
      feature.dat = feature.tab,
      meta.dat = meta.dat,
      formula = "~group",
      prev.filter = 0.75,
      is.winsor = FALSE,
      adaptive = FALSE,
      verbose = FALSE
    ),
    "All features were filtered out in LinDA2"
  )

  expect_equal(dim(linda2_res$feature.dat.use), c(0, 4))
  expect_equal(linda2_res$variables, character(0))
  expect_equal(linda2_res$bias, numeric(0))
  expect_length(linda2_res$output, 0)
  expect_equal(rownames(linda2_res$meta.dat.use), colnames(feature.tab))
})

test_that("linda preserves sample row names with single-column metadata", {
  feature.tab <- matrix(
    c(
      5, 5, 5, 5,
      2, 2, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("Feature1", "Feature2"), paste0("Sample", 1:4))
  )
  meta.dat <- data.frame(
    group = c("A", "A", "B", "B"),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  linda_res <- linda(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group",
    is.winsor = FALSE,
    adaptive = FALSE,
    verbose = FALSE
  )

  expect_equal(rownames(linda_res$meta.dat.use), colnames(feature.tab))
})

test_that("mStat_has_random_effect_term distinguishes transforms from random effects", {
  expect_false(mStat_has_random_effect_term("~log(age) + group"))
  expect_false(mStat_has_random_effect_term(stats::as.formula("~log(age) + group")))
  expect_true(mStat_has_random_effect_term("~group + (1|subject)"))
  expect_true(mStat_has_random_effect_term(stats::as.formula("~group + (1|subject)")))
})

test_that("linda reference level switching works correctly", {
  # Create sample data
  set.seed(123)

  # Create feature table (simulated microbiome data)
  n_samples <- 40
  n_features <- 10
  feature.tab <- matrix(rpois(n_samples * n_features, lambda = 10),
                        nrow = n_features,
                        ncol = n_samples)
  rownames(feature.tab) <- paste0("Feature", 1:n_features)
  colnames(feature.tab) <- paste0("Sample", 1:n_samples)

  # Create metadata
  meta.dat <- data.frame(
    group = factor(rep(c("A", "B"), each = n_samples/2)),
    sex = factor(rep(c("M", "F"), times = n_samples/2))
  )
  rownames(meta.dat) <- colnames(feature.tab)

  # Create feature annotation matrix
  feature.ann <- matrix(
    c(
      paste0("Phylum", rep(1:2, each = 5)),
      paste0("Genus", 1:n_features)
    ),
    nrow = n_features,
    ncol = 2,
    dimnames = list(
      rownames(feature.tab),
      c("Phylum", "Genus")
    )
  )

  # Build data.obj
  data.obj <- list(
    feature.tab = feature.tab,
    meta.dat = meta.dat,
    feature.ann = feature.ann
  )

  # Test 1: Using original order (A as reference)
  test.list1 <- generate_taxa_test_single(
    data.obj = data.obj,
    group.var = "group",
    adj.vars = "sex",
    feature.level = "Genus",
    feature.dat.type = "count",
    prev.filter = 0.1,
    abund.filter = 0.0001
  )

  # Output detailed results for first test
  cat("\n=== Test 1 Results (A as reference) ===\n")
  cat("Result names:", names(test.list1$Genus), "\n")
  cat("\nFirst few rows of results:\n")
  print(head(test.list1$Genus[[1]]))

  # Change reference level
  data.obj$meta.dat$group <- relevel(data.obj$meta.dat$group, ref = "B")

  # Test 2: Using new order (B as reference)
  test.list2 <- generate_taxa_test_single(
    data.obj = data.obj,
    group.var = "group",
    adj.vars = "sex",
    feature.level = "Genus",
    feature.dat.type = "count",
    prev.filter = 0.1,
    abund.filter = 0.0001
  )

  # Output detailed results for second test
  cat("\n=== Test 2 Results (B as reference) ===\n")
  cat("Result names:", names(test.list2$Genus), "\n")
  cat("\nFirst few rows of results:\n")
  print(head(test.list2$Genus[[1]]))

  # Compare coefficients
  cat("\n=== Coefficient Comparison ===\n")
  comparison_df <- data.frame(
    Variable = test.list1$Genus[[1]]$Variable,
    Coef_Test1 = test.list1$Genus[[1]]$Coefficient,
    Coef_Test2 = test.list2$Genus[[1]]$Coefficient,
    P_Value_Test1 = test.list1$Genus[[1]]$P.Value,
    P_Value_Test2 = test.list2$Genus[[1]]$P.Value
  )
  print(comparison_df)

  # Validation
  expect_equal(names(test.list1$Genus), "B vs A (Reference) [Main Effect]")
  expect_equal(names(test.list2$Genus), "A vs B (Reference) [Main Effect]")
  expect_equal(test.list1$Genus[[1]]$Coefficient,
              -test.list2$Genus[[1]]$Coefficient)
  expect_equal(test.list1$Genus[[1]]$P.Value,
              test.list2$Genus[[1]]$P.Value)
  expect_equal(test.list1$Genus[[1]]$SE,
              test.list2$Genus[[1]]$SE)
})
