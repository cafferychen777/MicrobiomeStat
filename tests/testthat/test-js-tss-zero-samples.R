# Test for Jensen-Shannon divergence with TSS normalization and zero-sum samples
# Related to GitHub issue #45

library(MicrobiomeStat)

test_that("JS divergence handles zero-sum samples after TSS normalization", {
  
  # Create test data with a zero-sum sample
  feature_tab <- matrix(
    c(100, 50, 30, 20,    # Sample 1: normal
       80, 40, 20, 10,    # Sample 2: normal
       60, 30, 15,  5,    # Sample 3: low counts
        0,  0,  0,  0),   # Sample 4: all zeros (edge case)
    nrow = 4, ncol = 4
  )
  
  colnames(feature_tab) <- paste0("Sample_", 1:4)
  rownames(feature_tab) <- paste0("Feature_", 1:4)
  
  # Create metadata
  meta_dat <- data.frame(
    sample_id = colnames(feature_tab),
    group = c("A", "A", "B", "B"),
    row.names = colnames(feature_tab)
  )
  
  # Create feature annotation
  feature_ann <- matrix(
    paste0("Taxa_", 1:4),
    nrow = 4, ncol = 1,
    dimnames = list(rownames(feature_tab), "Taxonomy")
  )
  
  # Create data object
  data.obj <- list(
    feature.tab = feature_tab,
    meta.dat = meta_dat,
    feature.ann = feature_ann
  )
  
  # Test TSS normalization with zero-sum sample
  expect_warning(
    norm_result <- mStat_normalize_data(data.obj, method = "TSS"),
    regexp = "zero total counts"
  )
  
  norm_data <- norm_result$data.obj.norm
  
  # Verify no NaN values in normalized data
  expect_false(any(is.nan(norm_data$feature.tab)))
  
  # Verify column sums (should be 1 for non-zero samples, 0 for zero sample)
  col_sums <- colSums(norm_data$feature.tab)
  expect_equal(unname(col_sums[1:3]), c(1, 1, 1), tolerance = 1e-10)
  expect_equal(unname(col_sums[4]), 0)
  
  # Test JS divergence calculation
  expect_warning(
    dist_obj <- mStat_calculate_beta_diversity(norm_data, dist.name = "JS"),
    regexp = "zero total counts|rarefied"
  )
  
  # Verify JS matrix exists and has no NA/NaN/Inf values
  expect_true("JS" %in% names(dist_obj))
  js_matrix <- as.matrix(dist_obj$JS)
  expect_false(any(is.na(js_matrix)))
  expect_false(any(is.nan(js_matrix)))
  expect_false(any(is.infinite(js_matrix)))
  
  # Verify matrix properties
  expect_equal(dim(js_matrix), c(4, 4))
  expect_equal(unname(diag(js_matrix)), rep(0, 4))  # Diagonal should be 0
  expect_true(isSymmetric(js_matrix))                # Should be symmetric
  
  # Test that cmdscale works with the JS matrix
  expect_no_error({
    mds_result <- cmdscale(as.dist(js_matrix), k = 2, eig = TRUE)
  })
})

test_that("JS divergence handles multiple edge cases", {
  
  # Create test data with various edge cases
  feature_tab <- matrix(
    c(100, 50, 30,    # Sample 1: normal
        0,  0,  0,    # Sample 2: all zeros
       10,  5,  2,    # Sample 3: very low counts
        1,  0,  0,    # Sample 4: single count
      0.5, 0.5, 0),   # Sample 5: fractional counts (already normalized)
    nrow = 3, ncol = 5
  )
  
  colnames(feature_tab) <- paste0("Sample_", 1:5)
  rownames(feature_tab) <- paste0("Feature_", 1:3)
  
  data.obj <- list(
    feature.tab = feature_tab,
    meta.dat = data.frame(
      sample_id = colnames(feature_tab),
      row.names = colnames(feature_tab)
    ),
    feature.ann = matrix(
      paste0("Taxa_", 1:3),
      nrow = 3, ncol = 1,
      dimnames = list(rownames(feature_tab), "Taxonomy")
    )
  )
  
  # Normalize and calculate JS
  norm_result <- suppressWarnings(mStat_normalize_data(data.obj, method = "TSS"))
  norm_data <- norm_result$data.obj.norm
  
  dist_obj <- suppressWarnings(
    mStat_calculate_beta_diversity(norm_data, dist.name = "JS")
  )
  
  js_matrix <- as.matrix(dist_obj$JS)
  
  # All checks should pass
  expect_false(any(is.na(js_matrix)))
  expect_false(any(is.nan(js_matrix)))
  expect_false(any(is.infinite(js_matrix)))
  
  # JS divergence should be between 0 and 1
  off_diagonal <- js_matrix[upper.tri(js_matrix)]
  expect_true(all(off_diagonal >= 0))
  expect_true(all(off_diagonal <= 1))
})

test_that("Full pipeline works with JS and TSS normalization", {
  
  # Use the peerj32 dataset if available, otherwise create test data
  if (exists("peerj32.obj")) {
    data.obj <- peerj32.obj
  } else {
    # Create simple test data
    set.seed(123)
    feature_tab <- matrix(
      rpois(100, lambda = 10),
      nrow = 10, ncol = 10
    )
    # Add a zero-sum sample
    feature_tab[, 10] <- 0
    
    colnames(feature_tab) <- paste0("Sample_", 1:10)
    rownames(feature_tab) <- paste0("Feature_", 1:10)
    
    data.obj <- list(
      feature.tab = feature_tab,
      meta.dat = data.frame(
        sample_id = colnames(feature_tab),
        group = rep(c("Control", "Treatment"), each = 5),
        row.names = colnames(feature_tab)
      ),
      feature.ann = matrix(
        paste0("Taxa_", 1:10),
        nrow = 10, ncol = 1,
        dimnames = list(rownames(feature_tab), "Taxonomy")
      )
    )
  }
  
  # Test the full pipeline
  expect_no_error({
    # Normalize
    norm_data <- suppressWarnings(
      mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    )
    
    # Calculate beta diversity
    dist_obj <- suppressWarnings(
      mStat_calculate_beta_diversity(norm_data, dist.name = "JS")
    )
    
    # Calculate PC
    pc_obj <- mStat_calculate_PC(dist_obj, method = "mds", k = 2, dist.name = "JS")
    
    # Verify PC calculation succeeded
    expect_true("JS" %in% names(pc_obj))
    expect_true("points" %in% names(pc_obj$JS))
  })
})