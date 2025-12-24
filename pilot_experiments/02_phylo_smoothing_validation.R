#!/usr/bin/env Rscript
#' ==============================================================================
#' Pilot Experiment 2: Phylogenetic Smoothing Validation
#' ==============================================================================
#'
#' Objective: Prove that phylogenetic tree smoothing can recover weak signals
#'            when related ASVs show coordinated (but individually weak) changes.
#'
#' Design:
#' - Simulate microbiome data with a phylogenetic tree
#' - Plant weak signals in a small clade (10 ASVs)
#' - Compare Method A (No smoothing) vs Method B (Tree-based smoothing)
#' - Measure sensitivity improvement and FDR control
#'
#' Author: Chen Yang (cafferychen777)
#' Date: 2025-12-24
#' ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggtree)
  library(patchwork)
})

set.seed(123)

# ==============================================================================
# 1. Simulation Parameters
# ==============================================================================

n_subjects_per_group <- 100  # Large sample for weak signal detection
n_features <- 100            # Total number of ASVs/features
n_signal_features <- 10      # ASVs in the signal clade
signal_effect_size <- 0.3    # WEAK effect (Cohen's d)
noise_sd <- 1.0              # Residual noise

# Smoothing parameter
lambda <- 0.5  # Regularization strength (0 = no smoothing, 1 = full smoothing)

cat("=== Pilot Experiment 2: Phylogenetic Smoothing Validation ===\n\n")
cat("Simulation Parameters:\n")
cat("  - Subjects per group:", n_subjects_per_group, "\n")
cat("  - Total features (ASVs):", n_features, "\n")
cat("  - Signal clade size:", n_signal_features, "\n")
cat("  - Effect size (Cohen's d):", signal_effect_size, "⚠ WEAK\n")
cat("  - Smoothing lambda:", lambda, "\n\n")

# ==============================================================================
# 2. Generate Phylogenetic Tree
# ==============================================================================

cat("Generating phylogenetic tree...\n")

# Create a random tree
tree <- rtree(n = n_features, tip.label = paste0("ASV", sprintf("%03d", 1:n_features)))

# Make it ultrametric (for easier interpretation)
tree <- chronos(tree, lambda = 1, model = "correlated")

# Identify a clade for signal placement
# Pick a random internal node that has ~10 descendants
find_clade_with_n_tips <- function(tree, target_n = 10, tolerance = 3) {
  n_tips <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)

  for (node in internal_nodes) {
    descendants <- Descendants(tree, node, type = "tips")[[1]]
    n_desc <- length(descendants)

    if (abs(n_desc - target_n) <= tolerance) {
      return(list(node = node, tips = descendants))
    }
  }

  # Fallback: use first 10 tips
  return(list(node = NULL, tips = 1:target_n))
}

signal_clade <- find_clade_with_n_tips(tree, target_n = n_signal_features)
signal_tip_indices <- signal_clade$tips
signal_asv_names <- tree$tip.label[signal_tip_indices]

cat("  Signal clade ASVs:", paste(signal_asv_names, collapse = ", "), "\n")
cat("  Clade node:", signal_clade$node, "\n\n")

# Visualize tree with signal clade highlighted
tree_data <- as_tibble(tree)
tree_data$is_signal <- tree_data$label %in% signal_asv_names

p_tree <- ggtree(tree, layout = "circular") +
  geom_tippoint(aes(color = label %in% signal_asv_names), size = 3) +
  scale_color_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#95a5a6"),
    labels = c("TRUE" = "Signal Clade", "FALSE" = "Background"),
    name = "ASV Type"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ggtitle("Phylogenetic Tree with Signal Clade Highlighted")

dir.create("pilot_experiments/results", recursive = TRUE, showWarnings = FALSE)
ggsave("pilot_experiments/results/02a_phylo_tree_signal_clade.pdf",
       p_tree, width = 10, height = 10)
cat("✓ Saved phylogenetic tree plot\n\n")

# ==============================================================================
# 3. Compute Tree-based Distance and Laplacian Matrix
# ==============================================================================

cat("Computing tree-based matrices...\n")

# Compute patristic distance matrix
dist_matrix <- cophenetic.phylo(tree)

# Convert to similarity (inverse distance)
# Use exponential kernel: S_ij = exp(-d_ij / bandwidth)
bandwidth <- median(dist_matrix[lower.tri(dist_matrix)])
similarity_matrix <- exp(-dist_matrix / bandwidth)

# Normalize to adjacency matrix
adjacency_matrix <- similarity_matrix
diag(adjacency_matrix) <- 0  # Remove self-loops

# Compute degree matrix
degree_vector <- rowSums(adjacency_matrix)
degree_matrix <- diag(degree_vector)

# Compute graph Laplacian: L = D - A
laplacian_matrix <- degree_matrix - adjacency_matrix

# Normalize Laplacian: L_norm = D^(-1/2) * L * D^(-1/2)
D_inv_sqrt <- diag(1 / sqrt(degree_vector))
laplacian_normalized <- D_inv_sqrt %*% laplacian_matrix %*% D_inv_sqrt

cat("  Distance matrix:", nrow(dist_matrix), "x", ncol(dist_matrix), "\n")
cat("  Median patristic distance:", round(bandwidth, 3), "\n")
cat("  Laplacian eigenvalues range:",
    round(range(eigen(laplacian_normalized)$values), 3), "\n\n")

# ==============================================================================
# 4. Simulate Feature Abundance Data
# ==============================================================================

cat("Generating simulated abundance data...\n")

simulate_abundance_data <- function() {

  # Create subject IDs
  n_total <- n_subjects_per_group * 2
  subject_ids <- paste0("S", sprintf("%03d", 1:n_total))
  groups <- rep(c("Control", "Treatment"), each = n_subjects_per_group)

  # Initialize feature matrix
  feature_matrix <- matrix(NA, nrow = n_total, ncol = n_features)
  colnames(feature_matrix) <- tree$tip.label

  # Simulate features
  for (i in 1:n_features) {
    asv_name <- tree$tip.label[i]

    # Base abundance (null)
    base_abundance <- rnorm(n_total, mean = 0, sd = noise_sd)

    # Add weak signal to clade members
    if (asv_name %in% signal_asv_names) {
      treatment_indicator <- ifelse(groups == "Treatment", 1, 0)
      signal <- signal_effect_size * treatment_indicator
      feature_matrix[, i] <- base_abundance + signal
    } else {
      feature_matrix[, i] <- base_abundance
    }
  }

  # Create data frame
  result <- data.frame(
    SubjectID = subject_ids,
    Group = groups,
    feature_matrix,
    check.names = FALSE
  )

  return(result)
}

sim_data <- simulate_abundance_data()
cat("  Data dimensions:", nrow(sim_data), "subjects x",
    (ncol(sim_data) - 2), "features\n\n")

# ==============================================================================
# 5. Method A: Standard t-test (No Smoothing)
# ==============================================================================

cat("Running Method A: Standard t-test (no smoothing)...\n")

run_standard_test <- function(data, feature_name) {
  # Simple t-test
  control_values <- data[data$Group == "Control", feature_name]
  treatment_values <- data[data$Group == "Treatment", feature_name]

  test_result <- t.test(treatment_values, control_values, var.equal = TRUE)

  # Extract z-score (standardized effect)
  pooled_sd <- sqrt(((length(control_values) - 1) * var(control_values) +
                     (length(treatment_values) - 1) * var(treatment_values)) /
                    (length(control_values) + length(treatment_values) - 2))

  mean_diff <- mean(treatment_values) - mean(control_values)
  z_score <- mean_diff / (pooled_sd * sqrt(1/length(control_values) + 1/length(treatment_values)))

  return(list(
    p_value = test_result$p.value,
    z_score = z_score,
    mean_diff = mean_diff
  ))
}

# Test all features
feature_names <- tree$tip.label
results_standard <- lapply(feature_names, function(fname) {
  run_standard_test(sim_data, fname)
})
names(results_standard) <- feature_names

# Extract vectors
p_values_standard <- sapply(results_standard, function(x) x$p_value)
z_scores_standard <- sapply(results_standard, function(x) x$z_score)

cat("  Significant features (P < 0.05):", sum(p_values_standard < 0.05), "\n")
cat("  True positives detected:",
    sum(p_values_standard[signal_asv_names] < 0.05), "/", n_signal_features, "\n")
cat("  Mean |Z| for signal clade:", round(mean(abs(z_scores_standard[signal_asv_names])), 3), "\n\n")

# ==============================================================================
# 6. Method B: Phylogenetic Smoothing
# ==============================================================================

cat("Running Method B: Tree-based smoothing...\n")

apply_phylo_smoothing <- function(z_scores, laplacian, lambda) {
  # Smoothing formula: z_smooth = (I + lambda * L)^(-1) * z_raw
  #
  # Interpretation:
  # - lambda = 0: No smoothing, z_smooth = z_raw
  # - lambda > 0: Borrow strength from phylogenetically close neighbors
  # - Higher lambda = more smoothing

  n <- length(z_scores)
  I <- diag(n)

  smoothing_matrix <- I + lambda * laplacian
  smoothing_matrix_inv <- solve(smoothing_matrix)

  z_smoothed <- as.vector(smoothing_matrix_inv %*% z_scores)
  names(z_smoothed) <- names(z_scores)

  return(z_smoothed)
}

# Apply smoothing
z_scores_smoothed <- apply_phylo_smoothing(z_scores_standard, laplacian_normalized, lambda)

# Convert back to p-values (two-tailed)
p_values_smoothed <- 2 * (1 - pnorm(abs(z_scores_smoothed)))

cat("  Significant features (P < 0.05):", sum(p_values_smoothed < 0.05), "\n")
cat("  True positives detected:",
    sum(p_values_smoothed[signal_asv_names] < 0.05), "/", n_signal_features, "\n")
cat("  Mean |Z| for signal clade (smoothed):",
    round(mean(abs(z_scores_smoothed[signal_asv_names])), 3), "\n\n")

# ==============================================================================
# 7. Performance Comparison
# ==============================================================================

cat("=== Performance Comparison ===\n\n")

# Create results table
results_df <- data.frame(
  ASV = feature_names,
  TrueSignal = feature_names %in% signal_asv_names,
  Z_Raw = z_scores_standard,
  Z_Smoothed = z_scores_smoothed,
  P_Raw = p_values_standard,
  P_Smoothed = p_values_smoothed
) %>%
  mutate(
    Sig_Raw = P_Raw < 0.05,
    Sig_Smoothed = P_Smoothed < 0.05,
    Z_Change = Z_Smoothed - Z_Raw,
    P_Log_Ratio = log10(P_Raw / pmax(P_Smoothed, 1e-100))
  )

# Calculate metrics
alpha <- 0.05

# True Positives
tp_raw <- sum(results_df$TrueSignal & results_df$Sig_Raw)
tp_smoothed <- sum(results_df$TrueSignal & results_df$Sig_Smoothed)

# False Positives
fp_raw <- sum(!results_df$TrueSignal & results_df$Sig_Raw)
fp_smoothed <- sum(!results_df$TrueSignal & results_df$Sig_Smoothed)

# Sensitivity
sensitivity_raw <- tp_raw / n_signal_features
sensitivity_smoothed <- tp_smoothed / n_signal_features

# FDR
fdr_raw <- fp_raw / max(1, sum(results_df$Sig_Raw))
fdr_smoothed <- fp_smoothed / max(1, sum(results_df$Sig_Smoothed))

# Print comparison
comparison_table <- data.frame(
  Metric = c("True Positives", "False Positives", "Sensitivity (Power)",
             "False Discovery Rate", "Total Significant"),
  No_Smoothing = c(
    tp_raw,
    fp_raw,
    sprintf("%.1f%%", sensitivity_raw * 100),
    sprintf("%.1f%%", fdr_raw * 100),
    sum(results_df$Sig_Raw)
  ),
  With_Smoothing = c(
    tp_smoothed,
    fp_smoothed,
    sprintf("%.1f%%", sensitivity_smoothed * 100),
    sprintf("%.1f%%", fdr_smoothed * 100),
    sum(results_df$Sig_Smoothed)
  ),
  stringsAsFactors = FALSE
)

print(comparison_table)
cat("\n")

# Sensitivity improvement
sensitivity_improvement <- (sensitivity_smoothed - sensitivity_raw) / max(0.01, sensitivity_raw)
cat("Sensitivity Improvement:", sprintf("%.1f%%", sensitivity_improvement * 100), "\n\n")

# ==============================================================================
# 8. Visualization: Z-score Boost in Signal Clade
# ==============================================================================

# Plot Z-score changes for signal clade
signal_clade_data <- results_df %>%
  filter(TrueSignal) %>%
  select(ASV, Z_Raw, Z_Smoothed) %>%
  pivot_longer(cols = c(Z_Raw, Z_Smoothed),
               names_to = "Method",
               values_to = "Z_Score") %>%
  mutate(
    Method = factor(Method,
                   levels = c("Z_Raw", "Z_Smoothed"),
                   labels = c("No Smoothing", "With Smoothing"))
  )

p_z_boost <- ggplot(signal_clade_data, aes(x = ASV, y = Z_Score, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_hline(yintercept = -1.96, linetype = "dashed", color = "red", alpha = 0.6) +
  scale_fill_manual(values = c("No Smoothing" = "#3498db",
                               "With Smoothing" = "#e74c3c")) +
  labs(
    title = "Z-score Amplification in Signal Clade",
    subtitle = paste0("Smoothing lambda = ", lambda, " | True signal: weak (d = ", signal_effect_size, ")"),
    x = "ASV",
    y = "Z-score",
    fill = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) +
  annotate("text", x = 1, y = 2.5, label = "Significance threshold (±1.96)",
           color = "red", size = 3.5, hjust = 0)

ggsave("pilot_experiments/results/02b_z_score_boost.pdf",
       p_z_boost, width = 10, height = 6)
cat("✓ Saved Z-score boost plot\n\n")

# ==============================================================================
# 9. Visualization: P-value Comparison Scatter
# ==============================================================================

p_pvalue_scatter <- ggplot(results_df,
                           aes(x = pmax(P_Raw, 1e-10),
                               y = pmax(P_Smoothed, 1e-10),
                               color = TrueSignal)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#95a5a6"),
    labels = c("TRUE" = "Signal Clade", "FALSE" = "Background")
  ) +
  labs(
    title = "P-value Comparison: Phylogenetic Smoothing Effect",
    subtitle = paste0("Signal clade size: ", n_signal_features, " | Effect size: ", signal_effect_size),
    x = "P-value (No Smoothing)",
    y = "P-value (With Smoothing)",
    color = "ASV Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold")
  ) +
  annotate("text", x = 0.1, y = 1e-6, label = "Smoothing Helps",
           color = "#27ae60", fontface = "bold", size = 5)

ggsave("pilot_experiments/results/02c_pvalue_scatter.pdf",
       p_pvalue_scatter, width = 10, height = 8)
cat("✓ Saved p-value scatter plot\n\n")

# ==============================================================================
# 10. Visualization: Smoothing Effect Across the Tree
# ==============================================================================

# Create tree visualization with Z-score changes
tree_plot_data <- data.frame(
  label = tree$tip.label,
  Z_Change = results_df$Z_Change,
  Is_Signal = results_df$TrueSignal
)

p_tree_effect <- ggtree(tree, layout = "rectangular") +
  geom_tippoint(aes(color = tree_plot_data$Z_Change), size = 4) +
  scale_color_gradient2(
    low = "#3498db", mid = "gray90", high = "#e74c3c",
    midpoint = 0,
    name = "Z-score\nChange"
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ggtitle("Phylogenetic Smoothing Effect Across Tree")

ggsave("pilot_experiments/results/02d_tree_smoothing_effect.pdf",
       p_tree_effect, width = 12, height = 10)
cat("✓ Saved tree smoothing effect plot\n\n")

# ==============================================================================
# 11. Lambda Tuning: Sensitivity Analysis
# ==============================================================================

cat("Running lambda tuning analysis...\n")

lambda_values <- c(0, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0)

lambda_tuning_results <- lapply(lambda_values, function(lam) {

  # Apply smoothing with current lambda
  z_smooth <- apply_phylo_smoothing(z_scores_standard, laplacian_normalized, lam)
  p_smooth <- 2 * (1 - pnorm(abs(z_smooth)))

  # Calculate metrics
  tp <- sum(p_smooth[signal_asv_names] < 0.05)
  fp <- sum(p_smooth[!names(p_smooth) %in% signal_asv_names] < 0.05)

  sensitivity <- tp / n_signal_features
  fdr <- fp / max(1, tp + fp)

  data.frame(
    Lambda = lam,
    Sensitivity = sensitivity,
    FDR = fdr,
    True_Positives = tp,
    False_Positives = fp
  )
}) %>% bind_rows()

print(lambda_tuning_results)
cat("\n")

# Plot lambda tuning
p_lambda <- ggplot(lambda_tuning_results, aes(x = Lambda)) +
  geom_line(aes(y = Sensitivity, color = "Sensitivity"), size = 1.5) +
  geom_point(aes(y = Sensitivity, color = "Sensitivity"), size = 4) +
  geom_line(aes(y = FDR, color = "FDR"), size = 1.5) +
  geom_point(aes(y = FDR, color = "FDR"), size = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray50", alpha = 0.5) +
  scale_color_manual(
    values = c("Sensitivity" = "#27ae60", "FDR" = "#e74c3c"),
    name = "Metric"
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Lambda Tuning: Sensitivity vs FDR Trade-off",
    subtitle = "Optimal lambda balances power and false discovery control",
    x = "Smoothing Parameter (λ)",
    y = "Rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave("pilot_experiments/results/02e_lambda_tuning.pdf",
       p_lambda, width = 10, height = 6)
cat("✓ Saved lambda tuning plot\n\n")

# ==============================================================================
# 12. Save Results
# ==============================================================================

# Save detailed results
write.csv(results_df,
          "pilot_experiments/results/02_detailed_results.csv",
          row.names = FALSE)

# Save summary
summary_results <- list(
  simulation_params = list(
    n_subjects_per_group = n_subjects_per_group,
    n_features = n_features,
    n_signal_features = n_signal_features,
    signal_effect_size = signal_effect_size,
    lambda = lambda
  ),
  performance = comparison_table,
  sensitivity_improvement_pct = sensitivity_improvement * 100,
  lambda_tuning = lambda_tuning_results
)

saveRDS(summary_results, "pilot_experiments/results/02_summary.rds")

cat("=== Experiment Complete ===\n\n")
cat("Key Findings:\n")
cat("  - Smoothing detected", tp_smoothed, "true signals\n")
cat("  - No smoothing detected only", tp_raw, "true signals\n")
cat("  - Sensitivity improvement:", sprintf("%.1f%%", sensitivity_improvement * 100), "\n")
cat("  - FDR (Smoothed):", sprintf("%.1f%%", fdr_smoothed * 100), "\n\n")

if (sensitivity_smoothed > sensitivity_raw * 2 && fdr_smoothed < 0.2) {
  cat("✓✓✓ VALIDATION SUCCESSFUL ✓✓✓\n")
  cat("Phylogenetic smoothing dramatically improves sensitivity for weak signals!\n")
  cat("This approach is ready for full implementation in LinDA v2.\n")
} else if (sensitivity_smoothed > sensitivity_raw && fdr_smoothed < 0.3) {
  cat("✓ VALIDATION PROMISING ✓\n")
  cat("Phylogenetic smoothing shows improvement but may need parameter tuning.\n")
  cat("Consider optimizing lambda or using adaptive methods.\n")
} else {
  cat("⚠ WARNING: Results are inconclusive.\n")
  cat("Consider:\n")
  cat("  - Increasing signal clade size\n")
  cat("  - Adjusting effect size\n")
  cat("  - Tuning lambda parameter\n")
}
