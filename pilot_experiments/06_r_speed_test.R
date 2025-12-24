#!/usr/bin/env Rscript
#' ==============================================================================
#' R Computational Speed Test for LinDA v2 Permutation Testing
#' ==============================================================================
#'
#' Objective: Measure R performance for permutation testing to determine if
#'            Rcpp optimization is necessary for production deployment.
#'
#' Critical Decision Point:
#'   - If pure R is fast enough (< 0.1 sec/feature), deploy in R
#'   - If too slow (> 0.5 sec/feature), rewrite core loop in C++ (Rcpp)
#'
#' Author: Chen Yang (cafferychen777)
#' Date: 2025-12-24
#' ==============================================================================

cat("================================================================================\n")
cat("R Computational Speed Test for LinDA v2\n")
cat("================================================================================\n\n")

set.seed(789)

# ==============================================================================
# Test 1: Basic Permutation Loop Speed
# ==============================================================================

cat("TEST 1: BASIC PERMUTATION LOOP\n")
cat("----------------------------------------------------------------------\n\n")

# Simulate one feature's data (100 subjects, 5 timepoints)
n_subjects <- 100
n_timepoints <- 5
n_total <- n_subjects * n_timepoints

subject_ids <- rep(1:n_subjects, each = n_timepoints)
timepoints <- rep(1:n_timepoints, times = n_subjects)
groups <- rep(c(0, 1), each = (n_subjects / 2) * n_timepoints)
y_values <- rnorm(n_total)

cat("Dataset: 100 subjects × 5 timepoints = 500 observations\n")
cat("Running 1000 permutations...\n\n")

# Permutation test implementation
start_time <- Sys.time()

# Calculate observed statistic
calc_ssd <- function(y, time, group) {
  # Sum of squared differences across timepoints
  mean_0 <- tapply(y[group == 0], time[group == 0], mean)
  mean_1 <- tapply(y[group == 1], time[group == 1], mean)

  # Handle missing timepoints
  common_times <- intersect(names(mean_0), names(mean_1))
  if (length(common_times) == 0) return(0)

  ssd <- sum((mean_1[common_times] - mean_0[common_times])^2)
  return(ssd)
}

observed_stat <- calc_ssd(y_values, timepoints, groups)

# Subject-level group assignment
subject_groups <- groups[seq(1, n_total, by = n_timepoints)]

# Permutation loop
n_permutations <- 1000
null_stats <- numeric(n_permutations)

for (i in 1:n_permutations) {
  # Shuffle group labels at subject level
  perm_groups <- sample(subject_groups)

  # Expand to observation level
  perm_groups_expanded <- rep(perm_groups, each = n_timepoints)

  # Calculate statistic
  null_stats[i] <- calc_ssd(y_values, timepoints, perm_groups_expanded)
}

end_time <- Sys.time()
duration_basic <- as.numeric(end_time - start_time, units = "secs")

# Calculate p-value
p_value <- mean(null_stats >= observed_stat)

cat("Results:\n")
cat(sprintf("  Time for 1 feature: %.3f seconds\n", duration_basic))
cat(sprintf("  Observed statistic: %.4f\n", observed_stat))
cat(sprintf("  P-value: %.4f\n", p_value))
cat(sprintf("\nProjected time for 1000 features: %.1f minutes\n", duration_basic * 1000 / 60))

# Speed assessment
if (duration_basic < 0.1) {
  cat("\n✅ EXCELLENT: Pure R is fast enough!\n")
  cat("   No Rcpp optimization needed.\n")
  speed_grade_basic <- "EXCELLENT"
  rcpp_needed_basic <- FALSE
} else if (duration_basic < 0.5) {
  cat("\n✓ GOOD: Acceptable speed, but Rcpp would help.\n")
  cat("   Consider optimization for large datasets (>500 features).\n")
  speed_grade_basic <- "GOOD"
  rcpp_needed_basic <- FALSE
} else if (duration_basic < 2.0) {
  cat("\n⚠ SLOW: Usable but inefficient.\n")
  cat("   Rcpp optimization recommended for production.\n")
  speed_grade_basic <- "SLOW"
  rcpp_needed_basic <- TRUE
} else {
  cat("\n🚨 TOO SLOW: Must use Rcpp/C++!\n")
  cat("   Pure R implementation not viable for production.\n")
  speed_grade_basic <- "TOO_SLOW"
  rcpp_needed_basic <- TRUE
}

# ==============================================================================
# Test 2: Optimized R Implementation (vectorized)
# ==============================================================================

cat("\n\n")
cat("TEST 2: OPTIMIZED R IMPLEMENTATION (Vectorized)\n")
cat("----------------------------------------------------------------------\n\n")

cat("Testing optimized version with matrix operations...\n\n")

start_time_opt <- Sys.time()

# Vectorized calculation
calc_ssd_vectorized <- function(y, time, group, n_time = 5) {
  # Pre-allocate matrices
  mean_mat <- matrix(NA, nrow = 2, ncol = n_time)

  for (g in 0:1) {
    for (t in 1:n_time) {
      idx <- which(group == g & time == t)
      if (length(idx) > 0) {
        mean_mat[g + 1, t] <- mean(y[idx])
      }
    }
  }

  # Calculate SSD, handling NAs
  diffs <- mean_mat[2, ] - mean_mat[1, ]
  ssd <- sum(diffs^2, na.rm = TRUE)
  return(ssd)
}

observed_stat_opt <- calc_ssd_vectorized(y_values, timepoints, groups, n_timepoints)

# Optimized permutation loop
null_stats_opt <- numeric(n_permutations)

for (i in 1:n_permutations) {
  perm_groups <- sample(subject_groups)
  perm_groups_expanded <- rep(perm_groups, each = n_timepoints)
  null_stats_opt[i] <- calc_ssd_vectorized(y_values, timepoints, perm_groups_expanded, n_timepoints)
}

end_time_opt <- Sys.time()
duration_opt <- as.numeric(end_time_opt - start_time_opt, units = "secs")

p_value_opt <- mean(null_stats_opt >= observed_stat_opt)

cat("Results:\n")
cat(sprintf("  Time for 1 feature: %.3f seconds\n", duration_opt))
cat(sprintf("  Speedup vs basic: %.1fx\n", duration_basic / duration_opt))
cat(sprintf("  P-value: %.4f\n", p_value_opt))
cat(sprintf("\nProjected time for 1000 features: %.1f minutes\n", duration_opt * 1000 / 60))

if (duration_opt < 0.1) {
  speed_grade_opt <- "EXCELLENT"
  rcpp_needed_opt <- FALSE
  cat("\n✅ EXCELLENT: Optimized R is production-ready!\n")
} else if (duration_opt < 0.5) {
  speed_grade_opt <- "GOOD"
  rcpp_needed_opt <- FALSE
  cat("\n✓ GOOD: Optimized version is acceptable.\n")
} else {
  speed_grade_opt <- "SLOW"
  rcpp_needed_opt <- TRUE
  cat("\n⚠ SLOW: Even optimized R needs Rcpp.\n")
}

# ==============================================================================
# Test 3: Parallel Processing Potential
# ==============================================================================

cat("\n\n")
cat("TEST 3: PARALLEL PROCESSING ESTIMATE\n")
cat("----------------------------------------------------------------------\n\n")

# Check available cores
n_cores <- parallel::detectCores()
cat(sprintf("Available CPU cores: %d\n", n_cores))

# Conservative parallelization (use 75% of cores)
n_workers <- floor(n_cores * 0.75)
cat(sprintf("Recommended workers: %d\n\n", n_workers))

# Estimate speedup (accounting for overhead)
parallel_efficiency <- 0.7  # Typical efficiency for embarrassingly parallel tasks
estimated_speedup <- n_workers * parallel_efficiency

time_parallel <- duration_opt / estimated_speedup
time_1000_features_parallel <- time_parallel * 1000 / 60

cat("Estimated parallel performance:\n")
cat(sprintf("  Time per feature: %.3f seconds\n", time_parallel))
cat(sprintf("  Time for 1000 features: %.1f minutes\n", time_1000_features_parallel))
cat(sprintf("  Speedup: %.1fx (with %d workers)\n", estimated_speedup, n_workers))

if (time_1000_features_parallel < 5) {
  cat("\n✅ With parallelization, pure R is production-ready!\n")
  parallel_sufficient <- TRUE
} else if (time_1000_features_parallel < 15) {
  cat("\n✓ With parallelization, acceptable for most use cases.\n")
  parallel_sufficient <- TRUE
} else {
  cat("\n⚠ Even with parallelization, Rcpp optimization recommended.\n")
  parallel_sufficient <- FALSE
}

# ==============================================================================
# Overall Recommendation
# ==============================================================================

cat("\n\n")
cat("================================================================================\n")
cat("OVERALL ASSESSMENT AND RECOMMENDATIONS\n")
cat("================================================================================\n\n")

cat("Performance Summary:\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("  Basic R implementation:      %.3f sec/feature (Grade: %s)\n",
            duration_basic, speed_grade_basic))
cat(sprintf("  Optimized R implementation:  %.3f sec/feature (Grade: %s)\n",
            duration_opt, speed_grade_opt))
cat(sprintf("  Parallel R (estimated):      %.3f sec/feature\n", time_parallel))

cat("\nProjected Analysis Time (1000 features):\n")
cat("----------------------------------------------------------------------\n")
cat(sprintf("  Basic R:           %.1f minutes\n", duration_basic * 1000 / 60))
cat(sprintf("  Optimized R:       %.1f minutes\n", duration_opt * 1000 / 60))
cat(sprintf("  Parallel R:        %.1f minutes\n", time_1000_features_parallel))

cat("\n")
cat("IMPLEMENTATION RECOMMENDATION:\n")
cat("----------------------------------------------------------------------\n")

if (!rcpp_needed_opt && parallel_sufficient) {
  cat("✅ USE PURE R + PARALLELIZATION\n\n")
  cat("Rationale:\n")
  cat("  • Optimized R performance is acceptable\n")
  cat("  • Parallelization provides sufficient speedup\n")
  cat("  • Avoids complexity of C++ integration\n")
  cat("  • Easier to maintain and debug\n\n")
  cat("Implementation:\n")
  cat("  1. Use optimized vectorized R code\n")
  cat("  2. Implement parallel::mclapply() for feature loop\n")
  cat("  3. Set default workers = min(n_features, n_cores * 0.75)\n")
  cat("  4. Provide progress bar (pbapply package)\n\n")
  cat("Expected production performance:\n")
  cat(sprintf("  • 100 features:  ~%.1f minutes\n", time_1000_features_parallel / 10))
  cat(sprintf("  • 500 features:  ~%.1f minutes\n", time_1000_features_parallel / 2))
  cat(sprintf("  • 1000 features: ~%.1f minutes\n", time_1000_features_parallel))

  final_rec <- "PURE_R_PARALLEL"

} else if (!rcpp_needed_opt) {
  cat("✓ PURE R ACCEPTABLE (Optimization Optional)\n\n")
  cat("Rationale:\n")
  cat("  • R performance is workable\n")
  cat("  • Can ship v2.0 in pure R\n")
  cat("  • Consider Rcpp for v2.1 if users complain\n\n")
  cat("Implementation:\n")
  cat("  1. Ship initial version in pure R\n")
  cat("  2. Gather user feedback on performance\n")
  cat("  3. Add Rcpp optimization in future release if needed\n")

  final_rec <- "PURE_R_OK"

} else {
  cat("🚨 RCPP OPTIMIZATION REQUIRED\n\n")
  cat("Rationale:\n")
  cat("  • Pure R too slow for production use\n")
  cat("  • Even parallelization insufficient\n")
  cat("  • Users will not tolerate long wait times\n\n")
  cat("Implementation:\n")
  cat("  1. Rewrite calc_ssd() in C++ (Rcpp)\n")
  cat("  2. Parallelize permutation loop in C++\n")
  cat("  3. Fallback to pure R if Rcpp unavailable\n\n")
  cat("Estimated Rcpp speedup:\n")
  cat("  • Expect 10-50x faster than pure R\n")
  cat("  • Target: < 0.05 sec/feature\n")
  cat(sprintf("  • 1000 features: < %.1f minutes\n", duration_opt * 1000 / 60 / 20))

  final_rec <- "RCPP_REQUIRED"
}

cat("\n")
cat("Next Steps:\n")
cat("----------------------------------------------------------------------\n")

if (final_rec == "PURE_R_PARALLEL") {
  cat("1. Implement optimized R version with parallelization\n")
  cat("2. Add progress bars for user feedback\n")
  cat("3. Test on real datasets (100-1000 features)\n")
  cat("4. Benchmark against ANCOM-BC2, MaAsLin2\n")
  cat("5. Proceed to manuscript preparation\n")
} else if (final_rec == "PURE_R_OK") {
  cat("1. Ship LinDA v2.0 in pure R (optimized)\n")
  cat("2. Monitor user feedback on speed\n")
  cat("3. If complaints, implement Rcpp in v2.1\n")
  cat("4. Document expected runtime in user guide\n")
} else {
  cat("1. Learn Rcpp basics (RcppArmadillo recommended)\n")
  cat("2. Rewrite core permutation loop in C++\n")
  cat("3. Add comprehensive unit tests\n")
  cat("4. Benchmark Rcpp vs pure R\n")
  cat("5. Ensure backward compatibility\n")
}

# ==============================================================================
# Save Results
# ==============================================================================

results <- list(
  basic_r_time = duration_basic,
  optimized_r_time = duration_opt,
  parallel_r_time_est = time_parallel,
  time_1000_features_parallel_min = time_1000_features_parallel,
  speedup_optimization = duration_basic / duration_opt,
  speedup_parallel_est = estimated_speedup,
  n_cores = n_cores,
  n_workers_recommended = n_workers,
  rcpp_needed = rcpp_needed_opt || !parallel_sufficient,
  final_recommendation = final_rec,
  speed_grade_basic = speed_grade_basic,
  speed_grade_optimized = speed_grade_opt
)

saveRDS(results, "pilot_experiments/results/06_r_speed_test_results.rds")

cat("\n")
cat("================================================================================\n")
cat("Results saved to: pilot_experiments/results/06_r_speed_test_results.rds\n")
cat("================================================================================\n")
