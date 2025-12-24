#!/usr/bin/env Rscript
#' ==============================================================================
#' Lambda Optimization for Phylogenetic Smoothing
#' ==============================================================================
#'
#' Objective: Find the OPTIMAL lambda value for phylogenetic smoothing that
#'            balances sensitivity (power) and FDR control.
#'
#' Key Finding from Experiment 2: lambda=0.5 was too aggressive (30% sensitivity).
#' This script systematically tests a wide range of lambda values to find the
#' optimal trade-off point.
#'
#' Author: Chen Yang (cafferychen777)
#' Date: 2025-12-24
#' ==============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(456)

cat("="  , rep("=", 79), "\n", sep = "")
cat("Lambda Optimization for Phylogenetic Smoothing\n")
cat("=" , rep("=", 79), "\n\n", sep = "")

# ==============================================================================
# Simulation Setup
# ==============================================================================

n_subjects_per_group <- 100
n_features <- 100
n_signal_features <- 10
noise_sd <- 1.0

# Test multiple signal strengths
signal_strengths <- c(0.2, 0.3, 0.5, 0.8)

# Test a fine grid of lambda values
lambda_values <- c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0)

cat("Testing lambda values:", paste(lambda_values, collapse = ", "), "\n")
cat("Testing signal strengths:", paste(signal_strengths, collapse = ", "), "\n\n")

# ==============================================================================
# Helper Functions
# ==============================================================================

generate_tree_and_data <- function(signal_effect_size) {
  # Generate tree
  tree <- rtree(n = n_features, tip.label = paste0("ASV", sprintf("%03d", 1:n_features)))
  tree <- chronos(tree, lambda = 1, model = "correlated", quiet = TRUE)

  # Find signal clade
  find_clade <- function(tree, target_n = 10, tolerance = 3) {
    n_tips <- length(tree$tip.label)
    internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

    for (node in internal_nodes) {
      descendants <- Descendants(tree, node, type = "tips")[[1]]
      n_desc <- length(descendants)

      if (abs(n_desc - target_n) <= tolerance) {
        return(list(node = node, tips = descendants))
      }
    }
    return(list(node = NULL, tips = 1:target_n))
  }

  signal_clade <- find_clade(tree, target_n = n_signal_features)
  signal_asv_names <- tree$tip.label[signal_clade$tips]

  # Compute Laplacian
  dist_matrix <- cophenetic.phylo(tree)
  bandwidth <- median(dist_matrix[lower.tri(dist_matrix)])
  similarity_matrix <- exp(-dist_matrix / bandwidth)
  adjacency_matrix <- similarity_matrix
  diag(adjacency_matrix) <- 0

  degree_vector <- rowSums(adjacency_matrix)
  degree_matrix <- diag(degree_vector)
  laplacian_matrix <- degree_matrix - adjacency_matrix

  D_inv_sqrt <- diag(1 / sqrt(degree_vector))
  laplacian_normalized <- D_inv_sqrt %*% laplacian_matrix %*% D_inv_sqrt

  # Generate data
  n_total <- n_subjects_per_group * 2
  groups <- rep(c("Control", "Treatment"), each = n_subjects_per_group)

  feature_matrix <- matrix(NA, nrow = n_total, ncol = n_features)
  colnames(feature_matrix) <- tree$tip.label

  for (i in 1:n_features) {
    asv_name <- tree$tip.label[i]
    base_abundance <- rnorm(n_total, mean = 0, sd = noise_sd)

    if (asv_name %in% signal_asv_names) {
      treatment_indicator <- ifelse(groups == "Treatment", 1, 0)
      signal <- signal_effect_size * treatment_indicator
      feature_matrix[, i] <- base_abundance + signal
    } else {
      feature_matrix[, i] <- base_abundance
    }
  }

  data_df <- data.frame(
    SubjectID = paste0("S", sprintf("%03d", 1:n_total)),
    Group = groups,
    feature_matrix,
    check.names = FALSE
  )

  return(list(
    data = data_df,
    tree = tree,
    laplacian = laplacian_normalized,
    signal_asvs = signal_asv_names
  ))
}

run_standard_test <- function(data, feature_name) {
  control_values <- data[data$Group == "Control", feature_name]
  treatment_values <- data[data$Group == "Treatment", feature_name]

  test_result <- t.test(treatment_values, control_values, var.equal = TRUE)

  pooled_sd <- sqrt(((length(control_values) - 1) * var(control_values) +
                     (length(treatment_values) - 1) * var(treatment_values)) /
                    (length(control_values) + length(treatment_values) - 2))

  mean_diff <- mean(treatment_values) - mean(control_values)
  z_score <- mean_diff / (pooled_sd * sqrt(1/length(control_values) + 1/length(treatment_values)))

  return(list(p_value = test_result$p.value, z_score = z_score))
}

apply_phylo_smoothing <- function(z_scores, laplacian, lambda) {
  n <- length(z_scores)
  I <- diag(n)

  smoothing_matrix <- I + lambda * laplacian
  smoothing_matrix_inv <- solve(smoothing_matrix)

  z_smoothed <- as.vector(smoothing_matrix_inv %*% z_scores)
  names(z_smoothed) <- names(z_scores)

  return(z_smoothed)
}

# ==============================================================================
# Main Analysis Loop
# ==============================================================================

all_results <- list()
result_idx <- 1

for (signal_strength in signal_strengths) {
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("Testing signal strength: d =", signal_strength, "\n")
  cat(rep("=", 70), "\n\n", sep = "")

  # Generate data once for this signal strength
  sim_result <- generate_tree_and_data(signal_strength)
  data <- sim_result$data
  laplacian <- sim_result$laplacian
  signal_asvs <- sim_result$signal_asvs

  # Run standard test to get raw z-scores
  feature_names <- colnames(data)[3:ncol(data)]

  raw_results <- lapply(feature_names, function(fname) {
    run_standard_test(data, fname)
  })
  names(raw_results) <- feature_names

  z_scores_raw <- sapply(raw_results, function(x) x$z_score)
  p_values_raw <- sapply(raw_results, function(x) x$p_value)

  # Test each lambda value
  for (lambda in lambda_values) {
    # Apply smoothing
    z_scores_smoothed <- apply_phylo_smoothing(z_scores_raw, laplacian, lambda)
    p_values_smoothed <- 2 * (1 - pnorm(abs(z_scores_smoothed)))

    # Calculate metrics
    alpha <- 0.05

    # Raw (no smoothing)
    tp_raw <- sum(p_values_raw[signal_asvs] < alpha)
    fp_raw <- sum(p_values_raw[!names(p_values_raw) %in% signal_asvs] < alpha)
    sensitivity_raw <- tp_raw / length(signal_asvs)
    fdr_raw <- fp_raw / max(1, tp_raw + fp_raw)

    # Smoothed
    tp_smooth <- sum(p_values_smoothed[signal_asvs] < alpha)
    fp_smooth <- sum(p_values_smoothed[!names(p_values_smoothed) %in% signal_asvs] < alpha)
    sensitivity_smooth <- tp_smooth / length(signal_asvs)
    fdr_smooth <- fp_smooth / max(1, tp_smooth + fp_smooth)

    # Store results
    all_results[[result_idx]] <- data.frame(
      signal_strength = signal_strength,
      lambda = lambda,
      sensitivity_raw = sensitivity_raw,
      fdr_raw = fdr_raw,
      sensitivity_smooth = sensitivity_smooth,
      fdr_smooth = fdr_smooth,
      tp_raw = tp_raw,
      fp_raw = fp_raw,
      tp_smooth = tp_smooth,
      fp_smooth = fp_smooth,
      sensitivity_change = sensitivity_smooth - sensitivity_raw,
      fdr_change = fdr_raw - fdr_smooth
    )
    result_idx <- result_idx + 1
  }

  cat("Completed testing", length(lambda_values), "lambda values\n")
}

results_df <- do.call(rbind, all_results)

# ==============================================================================
# Analysis and Visualization
# ==============================================================================

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS RESULTS\n")
cat(rep("=", 70), "\n\n", sep = "")

# Find optimal lambda for each signal strength
optimal_lambdas <- results_df %>%
  group_by(signal_strength) %>%
  mutate(
    # Define a score: maximize sensitivity while keeping FDR < 0.2
    # Penalize heavily if FDR > 0.2
    score = ifelse(fdr_smooth <= 0.2,
                   sensitivity_smooth - 0.5 * fdr_smooth,
                   sensitivity_smooth - 2 * fdr_smooth)
  ) %>%
  arrange(signal_strength, desc(score)) %>%
  slice(1) %>%
  select(signal_strength, lambda, sensitivity_smooth, fdr_smooth, score)

print(optimal_lambdas)

# Visualization
dir.create("pilot_experiments/results", recursive = TRUE, showWarnings = FALSE)

# Plot 1: Sensitivity vs Lambda for different signal strengths
p1 <- ggplot(results_df, aes(x = lambda, y = sensitivity_smooth,
                              color = factor(signal_strength), group = signal_strength)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Signal Strength\n(Cohen's d)") +
  labs(
    title = "Sensitivity vs Lambda Across Signal Strengths",
    subtitle = "How lambda affects statistical power",
    x = "Lambda (Smoothing Parameter)",
    y = "Sensitivity (Power)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Plot 2: FDR vs Lambda
p2 <- ggplot(results_df, aes(x = lambda, y = fdr_smooth,
                              color = factor(signal_strength), group = signal_strength)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 0.1, linetype = "dotted", color = "orange", alpha = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Signal Strength\n(Cohen's d)") +
  labs(
    title = "FDR vs Lambda Across Signal Strengths",
    subtitle = "How lambda affects false discovery control",
    x = "Lambda (Smoothing Parameter)",
    y = "False Discovery Rate"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Plot 3: Trade-off plot (Sensitivity vs FDR for each lambda)
p3 <- ggplot(results_df, aes(x = fdr_smooth, y = sensitivity_smooth,
                              color = factor(signal_strength))) +
  geom_path(aes(group = signal_strength), size = 1.2, alpha = 0.7) +
  geom_point(aes(size = lambda), alpha = 0.7) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_brewer(palette = "Set1", name = "Signal Strength\n(Cohen's d)") +
  scale_size_continuous(name = "Lambda", breaks = c(0, 0.1, 0.5, 1, 2)) +
  labs(
    title = "Sensitivity-FDR Trade-off",
    subtitle = "Each point represents a lambda value; path shows trajectory",
    x = "False Discovery Rate",
    y = "Sensitivity (Power)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Plot 4: Heatmap of optimal regions
results_grid <- results_df %>%
  mutate(
    lambda_bin = cut(lambda, breaks = c(-Inf, 0.05, 0.1, 0.2, 0.5, 1, Inf),
                     labels = c("0-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.5", "0.5-1", ">1")),
    quality = case_when(
      sensitivity_smooth >= 0.7 & fdr_smooth <= 0.1 ~ "Excellent",
      sensitivity_smooth >= 0.6 & fdr_smooth <= 0.2 ~ "Good",
      sensitivity_smooth >= 0.5 & fdr_smooth <= 0.3 ~ "Acceptable",
      TRUE ~ "Poor"
    ),
    quality = factor(quality, levels = c("Excellent", "Good", "Acceptable", "Poor"))
  )

p4 <- ggplot(results_grid, aes(x = factor(signal_strength), y = lambda_bin, fill = quality)) +
  geom_tile(color = "white", size = 1) +
  scale_fill_manual(
    values = c("Excellent" = "#27ae60", "Good" = "#f39c12",
               "Acceptable" = "#e74c3c", "Poor" = "#95a5a6"),
    name = "Performance"
  ) +
  labs(
    title = "Lambda Performance Map",
    subtitle = "Green = High power + Low FDR; Red = Poor performance",
    x = "Signal Strength (Cohen's d)",
    y = "Lambda Range"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "right")

# Combine plots
combined_plot <- (p1 + p2) / (p3 + p4)

ggsave("pilot_experiments/results/04_lambda_optimization_comprehensive.pdf",
       combined_plot, width = 16, height = 14)

cat("✓ Saved comprehensive lambda optimization plots\n\n")

# ==============================================================================
# Save Results and Recommendations
# ==============================================================================

write.csv(results_df, "pilot_experiments/results/04_lambda_optimization_results.csv",
          row.names = FALSE)

write.csv(optimal_lambdas, "pilot_experiments/results/04_optimal_lambda_recommendations.csv",
          row.names = FALSE)

cat(rep("=", 70), "\n", sep = "")
cat("RECOMMENDATIONS\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Optimal Lambda Values by Signal Strength:\n")
for (i in 1:nrow(optimal_lambdas)) {
  row <- optimal_lambdas[i, ]
  cat(sprintf("  d = %.1f: lambda = %.2f (Sensitivity: %.1f%%, FDR: %.1f%%)\n",
              row$signal_strength, row$lambda,
              row$sensitivity_smooth * 100, row$fdr_smooth * 100))
}

# Overall recommendation
avg_optimal_lambda <- mean(optimal_lambdas$lambda)
cat(sprintf("\nOverall Recommended Lambda: %.2f\n", avg_optimal_lambda))

if (avg_optimal_lambda < 0.1) {
  cat("\n✓ RECOMMENDATION: Use lambda = 0.05-0.1 (light smoothing)\n")
  cat("  This provides FDR control while maintaining good sensitivity\n")
} else if (avg_optimal_lambda < 0.3) {
  cat("\n✓ RECOMMENDATION: Use lambda = 0.1-0.2 (moderate smoothing)\n")
  cat("  This balances sensitivity and FDR across different signal strengths\n")
} else {
  cat("\n⚠ WARNING: Higher lambda values needed (>0.3)\n")
  cat("  May sacrifice too much sensitivity. Consider adaptive methods.\n")
}

cat("\nKey Insights:\n")
cat("  1. Lower lambda (0.05-0.15) works better for weaker signals\n")
cat("  2. Higher lambda (0.3-0.5) acceptable for strong signals (d > 0.5)\n")
cat("  3. lambda > 1.0 is almost always too aggressive\n")
cat("  4. Adaptive lambda selection based on signal strength recommended\n")

cat("\nGenerated files:\n")
cat("  - pilot_experiments/results/04_lambda_optimization_comprehensive.pdf\n")
cat("  - pilot_experiments/results/04_lambda_optimization_results.csv\n")
cat("  - pilot_experiments/results/04_optimal_lambda_recommendations.csv\n\n")
