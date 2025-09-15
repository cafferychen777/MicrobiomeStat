# Test Jensen-Shannon divergence calculation in beta diversity
# This test specifically addresses Issue #45: NA values in JS calculation after TSS normalization

test_that("Jensen-Shannon divergence works correctly with TSS normalized data", {
  # Skip if MicrobiomeStat functions are not available
  skip_if_not_installed("MicrobiomeStat")

  # Load required functions
  library(MicrobiomeStat)
  # Create test data with some zeros that will cause issues in the original implementation
  set.seed(123)
  
  # Create feature table with some samples having zero counts for certain features
  feature_tab <- matrix(c(
    100, 0, 50, 200,    # Feature 1: sample 2 has zero
    50, 100, 0, 150,    # Feature 2: sample 3 has zero  
    200, 150, 100, 0,   # Feature 3: sample 4 has zero
    0, 50, 200, 100     # Feature 4: sample 1 has zero
  ), nrow = 4, ncol = 4)
  
  rownames(feature_tab) <- paste0("Feature_", 1:4)
  colnames(feature_tab) <- paste0("Sample_", 1:4)
  
  # Create metadata
  meta_dat <- data.frame(
    sample_id = paste0("Sample_", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(meta_dat) <- paste0("Sample_", 1:4)
  
  # Create feature annotations
  feature_ann <- matrix(paste0("Taxa_", 1:4), nrow = 4, ncol = 1)
  rownames(feature_ann) <- paste0("Feature_", 1:4)
  colnames(feature_ann) <- "Taxonomy"
  
  # Create data object
  test_data <- list(
    feature.tab = feature_tab,
    meta.dat = meta_dat,
    feature.ann = feature_ann
  )
  
  # Test TSS normalization
  norm_result <- mStat_normalize_data(test_data, method = "TSS")
  norm_data <- norm_result$data.obj.norm
  
  # Verify normalization worked correctly
  expect_equal(as.numeric(colSums(norm_data$feature.tab)), rep(1, 4), tolerance = 1e-10)
  
  # Test Jensen-Shannon divergence calculation
  expect_no_error({
    dist_obj <- mStat_calculate_beta_diversity(norm_data, dist.name = "JS")
  })
  
  # Verify the distance object was created
  expect_true("JS" %in% names(dist_obj))
  expect_s3_class(dist_obj$JS, "dist")
  
  # Convert to matrix for easier testing
  js_matrix <- as.matrix(dist_obj$JS)
  
  # Test that the distance matrix is valid
  expect_false(any(is.na(js_matrix)), "Distance matrix should not contain NA values")
  expect_false(any(is.infinite(js_matrix)), "Distance matrix should not contain infinite values")
  
  # Test that diagonal is zero (distance from sample to itself)
  expect_equal(as.numeric(diag(js_matrix)), rep(0, 4), tolerance = 1e-10)
  
  # Test that matrix is symmetric
  expect_equal(js_matrix, t(js_matrix), tolerance = 1e-10)
  
  # Test that all distances are non-negative
  expect_true(all(js_matrix >= 0), "All distances should be non-negative")
  
  # Test that distances are reasonable (JS divergence should be between 0 and 1)
  expect_true(all(js_matrix <= 1), "Jensen-Shannon divergence should be <= 1")
})

test_that("Jensen-Shannon divergence handles edge cases correctly", {
  # Skip if MicrobiomeStat functions are not available
  skip_if_not_installed("MicrobiomeStat")

  # Load required functions
  library(MicrobiomeStat)
  # Test case 1: All samples identical (should give zero distances)
  identical_data <- matrix(rep(c(100, 50, 200, 150), 3), nrow = 4, ncol = 3)
  rownames(identical_data) <- paste0("Feature_", 1:4)
  colnames(identical_data) <- paste0("Sample_", 1:3)
  
  test_obj_identical <- list(
    feature.tab = identical_data,
    meta.dat = data.frame(sample_id = paste0("Sample_", 1:3), row.names = paste0("Sample_", 1:3)),
    feature.ann = matrix(paste0("Taxa_", 1:4), nrow = 4, ncol = 1, dimnames = list(paste0("Feature_", 1:4), "Taxonomy"))
  )
  
  # Normalize and test
  norm_identical <- mStat_normalize_data(test_obj_identical, method = "TSS")$data.obj.norm
  dist_identical <- mStat_calculate_beta_diversity(norm_identical, dist.name = "JS")
  
  # All distances should be zero (or very close to zero)
  js_matrix_identical <- as.matrix(dist_identical$JS)
  expect_true(all(js_matrix_identical < 1e-10), "Identical samples should have zero JS divergence")
  
  # Test case 2: Samples with completely different features (no overlap)
  no_overlap_data <- matrix(c(
    100, 0, 0,    # Feature 1 only in sample 1
    0, 100, 0,    # Feature 2 only in sample 2  
    0, 0, 100     # Feature 3 only in sample 3
  ), nrow = 3, ncol = 3)
  rownames(no_overlap_data) <- paste0("Feature_", 1:3)
  colnames(no_overlap_data) <- paste0("Sample_", 1:3)
  
  test_obj_no_overlap <- list(
    feature.tab = no_overlap_data,
    meta.dat = data.frame(sample_id = paste0("Sample_", 1:3), row.names = paste0("Sample_", 1:3)),
    feature.ann = matrix(paste0("Taxa_", 1:3), nrow = 3, ncol = 1, dimnames = list(paste0("Feature_", 1:3), "Taxonomy"))
  )
  
  # Normalize and test
  norm_no_overlap <- mStat_normalize_data(test_obj_no_overlap, method = "TSS")$data.obj.norm
  
  expect_no_error({
    dist_no_overlap <- mStat_calculate_beta_diversity(norm_no_overlap, dist.name = "JS")
  })
  
  js_matrix_no_overlap <- as.matrix(dist_no_overlap$JS)
  expect_false(any(is.na(js_matrix_no_overlap)), "No overlap case should not produce NA values")
  expect_false(any(is.infinite(js_matrix_no_overlap)), "No overlap case should not produce infinite values")
})

test_that("Jensen-Shannon divergence integrates correctly with beta ordination pipeline", {
  # Skip if MicrobiomeStat functions are not available
  skip_if_not_installed("MicrobiomeStat")

  # Load required functions
  library(MicrobiomeStat)
  # Create test data
  set.seed(456)
  feature_tab <- matrix(rpois(20, lambda = 50), nrow = 5, ncol = 4)
  feature_tab[1, 1] <- 0  # Add some zeros
  feature_tab[2, 2] <- 0
  
  rownames(feature_tab) <- paste0("Feature_", 1:5)
  colnames(feature_tab) <- paste0("Sample_", 1:4)
  
  meta_dat <- data.frame(
    sample_id = paste0("Sample_", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  rownames(meta_dat) <- paste0("Sample_", 1:4)
  
  feature_ann <- matrix(paste0("Taxa_", 1:5), nrow = 5, ncol = 1)
  rownames(feature_ann) <- paste0("Feature_", 1:5)
  colnames(feature_ann) <- "Taxonomy"
  
  test_data <- list(
    feature.tab = feature_tab,
    meta.dat = meta_dat,
    feature.ann = feature_ann
  )
  
  # Normalize data
  norm_data <- mStat_normalize_data(test_data, method = "TSS")$data.obj.norm
  
  # Test that the full beta ordination pipeline works with JS
  expect_no_error({
    plot_result <- generate_beta_ordination_single(
      data.obj = norm_data,
      subject.var = "sample_id",
      group.var = "group",
      dist.name = "JS",
      pdf = FALSE
    )
  })
  
  # The result should be a list of plots
  expect_type(plot_result, "list")
})
