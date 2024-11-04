test_that("linda p-value distribution is approximately uniform under null hypothesis", {
  # Create sample data with no true associations (null hypothesis)
  set.seed(123)

  # Generate random count data
  n_samples <- 200  # Larger sample size for better distribution testing
  n_features <- 100
  feature.tab <- matrix(rpois(n_samples * n_features, lambda = 10),
                       nrow = n_features,
                       ncol = n_samples)
  rownames(feature.tab) <- paste0("Feature", 1:n_features)
  colnames(feature.tab) <- paste0("Sample", 1:n_samples)

  # Create random metadata
  meta.dat <- data.frame(
    group = factor(rep(c("A", "B"), each = n_samples/2)),
    age = rnorm(n_samples),
    sex = factor(sample(c("M", "F"), n_samples, replace = TRUE))
  )
  rownames(meta.dat) <- colnames(feature.tab)

  # Create feature annotations
  feature.ann <- matrix(
    c(paste0("Phylum", rep(1:2, each = n_features/2)),
      paste0("Genus", 1:n_features)),
    nrow = n_features,
    ncol = 2,
    dimnames = list(rownames(feature.tab), c("Phylum", "Genus"))
  )

  # Run linda
  linda_result <- linda(
    feature.dat = feature.tab,
    meta.dat = meta.dat,
    formula = "~group + age + sex",
    feature.dat.type = "count"
  )

  # Extract p-values for each variable
  for(var in linda_result$variables) {
    p_values <- linda_result$output[[var]]$pvalue

    # Perform Kolmogorov-Smirnov test against uniform distribution
    ks_test <- ks.test(p_values, "punif")

    # Print diagnostic information
    cat("\nTesting variable:", var, "\n")
    cat("KS test p-value:", ks_test$p.value, "\n")

    # The test should not reject the null hypothesis (p-value should be > 0.05)
    # We use a smaller alpha (0.01) here because KS test can be sensitive
    expect_gt(ks_test$p.value, 0.01)

    # Check for obvious signs of non-uniformity
    # Calculate the ratio between max and min density in histogram
    hist_data <- hist(p_values, breaks = 10, plot = FALSE)
    density_ratio <- max(hist_data$density) / min(hist_data$density)

    cat("Density ratio:", density_ratio, "\n")

    # The ratio shouldn't be too extreme (e.g., not more than 10)
    expect_lt(density_ratio, 10)

    # Check if there's no strong clustering of p-values
    # Calculate the standard deviation of bin counts
    bin_sd <- sd(hist_data$counts)
    bin_mean <- mean(hist_data$counts)
    cv <- bin_sd / bin_mean  # coefficient of variation

    cat("Coefficient of variation:", cv, "\n")

    # CV shouldn't be too large (e.g., not more than 0.7)
    expect_lt(cv, 0.7)

    # Add p-value distribution histogram
    hist(p_values,
         main = paste("P-value distribution for", var),
         xlab = "P-value",
         breaks = 20)
  }
})
