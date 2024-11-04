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
  expect_equal(names(test.list1$Genus), "B vs A (Reference)")
  expect_equal(names(test.list2$Genus), "A vs B (Reference)")
  expect_equal(test.list1$Genus[[1]]$Coefficient,
              -test.list2$Genus[[1]]$Coefficient)
  expect_equal(test.list1$Genus[[1]]$P.Value,
              test.list2$Genus[[1]]$P.Value)
  expect_equal(test.list1$Genus[[1]]$SE,
              test.list2$Genus[[1]]$SE)
})
