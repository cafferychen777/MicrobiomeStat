#!/usr/bin/env Rscript
#' ==============================================================================
#' Pilot Experiment 1: Trajectory Spline Model Validation
#' ==============================================================================
#'
#' Objective: Prove that spline-based trajectory modeling can detect non-linear
#'            temporal patterns that linear models miss.
#'
#' Design:
#' - Simulate microbiome data with an "inverted-U" signal (spike at middle timepoints)
#' - Compare Method A (Linear LMM) vs Method B (Spline LMM)
#' - Measure power (true positive rate) and FDR control
#'
#' Author: Chen Yang (cafferychen777)
#' Date: 2025-12-24
#' ==============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(splines)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(42)

# ==============================================================================
# 1. Simulation Parameters
# ==============================================================================

n_subjects_per_group <- 50  # Subjects per group
n_timepoints <- 5           # Time points per subject
n_features <- 100           # Total number of features (taxa)
n_signal_features <- 10     # Features with true signals

# Signal parameters
signal_strength <- 2.0      # Effect size for the inverted-U pattern
noise_sd <- 1.0            # Residual standard deviation

cat("=== Pilot Experiment 1: Trajectory Spline Validation ===\n\n")
cat("Simulation Parameters:\n")
cat("  - Subjects per group:", n_subjects_per_group, "\n")
cat("  - Time points:", n_timepoints, "\n")
cat("  - Total features:", n_features, "\n")
cat("  - Signal features:", n_signal_features, "\n")
cat("  - Signal strength:", signal_strength, "\n\n")

# ==============================================================================
# 2. Data Simulation
# ==============================================================================

simulate_trajectory_data <- function() {

  # Create subject IDs
  n_total <- n_subjects_per_group * 2
  subject_ids <- paste0("S", sprintf("%03d", 1:n_total))

  # Assign groups
  groups <- rep(c("Control", "Treatment"), each = n_subjects_per_group)

  # Expand to all time points
  data_frame <- expand.grid(
    SubjectID = subject_ids,
    Time = 1:n_timepoints
  ) %>%
    mutate(
      Group = rep(groups, each = n_timepoints)
    ) %>%
    arrange(SubjectID, Time)

  # Initialize feature matrix
  n_obs <- nrow(data_frame)
  feature_matrix <- matrix(NA, nrow = n_obs, ncol = n_features)
  colnames(feature_matrix) <- paste0("Feature", sprintf("%03d", 1:n_features))

  # Simulate features
  for (i in 1:n_features) {

    # Random intercepts per subject
    subject_effects <- rnorm(n_total, mean = 0, sd = 0.5)
    subject_effect_vec <- rep(subject_effects, each = n_timepoints)

    # Base abundance (random)
    base_abundance <- rnorm(n_obs, mean = 0, sd = noise_sd)

    # Add signal to first n_signal_features
    if (i <= n_signal_features) {
      # Inverted-U pattern: peaks at Time = 3 (middle)
      # Formula: signal_strength * exp(-(Time - 3)^2 / 2)
      time_vec <- data_frame$Time
      group_vec <- data_frame$Group

      # Only Treatment group has the signal
      signal <- ifelse(
        group_vec == "Treatment",
        signal_strength * exp(-(time_vec - 3)^2 / 2),
        0
      )

      feature_matrix[, i] <- base_abundance + subject_effect_vec + signal

    } else {
      # Null features (no signal)
      feature_matrix[, i] <- base_abundance + subject_effect_vec
    }
  }

  # Combine into data frame
  result <- cbind(data_frame, as.data.frame(feature_matrix))

  return(result)
}

cat("Generating simulated data...\n")
sim_data <- simulate_trajectory_data()
cat("  Data dimensions:", nrow(sim_data), "observations x",
    (ncol(sim_data) - 3), "features\n\n")

# ==============================================================================
# 3. Visualize True Signal Pattern
# ==============================================================================

# Extract data for one signal feature
example_feature <- "Feature001"
plot_data <- sim_data %>%
  select(SubjectID, Time, Group, all_of(example_feature)) %>%
  group_by(Group, Time) %>%
  summarise(
    Mean = mean(!!sym(example_feature)),
    SE = sd(!!sym(example_feature)) / sqrt(n()),
    .groups = "drop"
  )

p_signal <- ggplot(plot_data, aes(x = Time, y = Mean, color = Group, fill = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = Mean - SE, ymax = Mean + SE), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Control" = "#3498db", "Treatment" = "#e74c3c")) +
  scale_fill_manual(values = c("Control" = "#3498db", "Treatment" = "#e74c3c")) +
  labs(
    title = "True Signal Pattern: Inverted-U Trajectory",
    subtitle = paste("Feature:", example_feature),
    x = "Time Point",
    y = "Abundance (Log Scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# Save plot
dir.create("pilot_experiments/results", recursive = TRUE, showWarnings = FALSE)
ggsave("pilot_experiments/results/01a_true_signal_pattern.pdf",
       p_signal, width = 8, height = 6)
cat("✓ Saved signal pattern plot\n\n")

# ==============================================================================
# 4. Method A: Linear LMM (Baseline)
# ==============================================================================

cat("Running Method A: Linear LMM...\n")

run_linear_lmm <- function(data, feature_name) {
  # Formula: Y ~ Group + Time + Group:Time + (1|SubjectID)
  formula_str <- paste0(feature_name, " ~ Group * Time + (1|SubjectID)")

  tryCatch({
    fit <- lmer(as.formula(formula_str), data = data)

    # Extract p-value for Group:Time interaction
    coef_summary <- summary(fit)$coefficients

    # Look for interaction term
    interaction_term <- grep("Group.*:Time|Time.*:Group", rownames(coef_summary))

    if (length(interaction_term) > 0) {
      p_value <- coef_summary[interaction_term, "Pr(>|t|)"]
    } else {
      p_value <- NA
    }

    return(p_value)

  }, error = function(e) {
    return(NA)
  })
}

# Test on all features
feature_names <- paste0("Feature", sprintf("%03d", 1:n_features))
results_linear <- sapply(feature_names, function(fname) {
  run_linear_lmm(sim_data, fname)
})

cat("  Completed:", sum(!is.na(results_linear)), "/", n_features, "features\n")
cat("  Significant features (P < 0.05):", sum(results_linear < 0.05, na.rm = TRUE), "\n")
cat("  True positives detected:",
    sum(results_linear[1:n_signal_features] < 0.05, na.rm = TRUE),
    "/", n_signal_features, "\n\n")

# ==============================================================================
# 5. Method B: Spline LMM (Proposed)
# ==============================================================================

cat("Running Method B: Spline LMM...\n")

run_spline_lmm <- function(data, feature_name) {
  # Formula: Y ~ Group * ns(Time, df=3) + (1|SubjectID)
  # Test the omnibus hypothesis of Group:Spline interaction

  tryCatch({
    # Fit full model with interaction
    formula_full <- paste0(feature_name, " ~ Group * ns(Time, df=3) + (1|SubjectID)")
    fit_full <- lmer(as.formula(formula_full), data = data)

    # Fit reduced model without interaction
    formula_reduced <- paste0(feature_name, " ~ Group + ns(Time, df=3) + (1|SubjectID)")
    fit_reduced <- lmer(as.formula(formula_reduced), data = data)

    # Likelihood ratio test
    lr_test <- anova(fit_reduced, fit_full)
    p_value <- lr_test$`Pr(>Chisq)`[2]

    return(p_value)

  }, error = function(e) {
    return(NA)
  })
}

# Test on all features
results_spline <- sapply(feature_names, function(fname) {
  run_spline_lmm(sim_data, fname)
})

cat("  Completed:", sum(!is.na(results_spline)), "/", n_features, "features\n")
cat("  Significant features (P < 0.05):", sum(results_spline < 0.05, na.rm = TRUE), "\n")
cat("  True positives detected:",
    sum(results_spline[1:n_signal_features] < 0.05, na.rm = TRUE),
    "/", n_signal_features, "\n\n")

# ==============================================================================
# 6. Performance Comparison
# ==============================================================================

cat("=== Performance Comparison ===\n\n")

# Create results table
results_df <- data.frame(
  Feature = feature_names,
  TrueSignal = c(rep(TRUE, n_signal_features), rep(FALSE, n_features - n_signal_features)),
  P_Linear = results_linear,
  P_Spline = results_spline
) %>%
  mutate(
    Sig_Linear = P_Linear < 0.05,
    Sig_Spline = P_Spline < 0.05
  )

# Calculate metrics
alpha <- 0.05

# True Positives
tp_linear <- sum(results_df$TrueSignal & results_df$Sig_Linear, na.rm = TRUE)
tp_spline <- sum(results_df$TrueSignal & results_df$Sig_Spline, na.rm = TRUE)

# False Positives
fp_linear <- sum(!results_df$TrueSignal & results_df$Sig_Linear, na.rm = TRUE)
fp_spline <- sum(!results_df$TrueSignal & results_df$Sig_Spline, na.rm = TRUE)

# Sensitivity (Power)
sensitivity_linear <- tp_linear / n_signal_features
sensitivity_spline <- tp_spline / n_signal_features

# False Discovery Rate
fdr_linear <- fp_linear / max(1, sum(results_df$Sig_Linear, na.rm = TRUE))
fdr_spline <- fp_spline / max(1, sum(results_df$Sig_Spline, na.rm = TRUE))

# Print comparison table
comparison_table <- data.frame(
  Metric = c("True Positives", "False Positives", "Sensitivity (Power)",
             "False Discovery Rate", "Total Significant"),
  Linear_LMM = c(
    tp_linear,
    fp_linear,
    sprintf("%.1f%%", sensitivity_linear * 100),
    sprintf("%.1f%%", fdr_linear * 100),
    sum(results_df$Sig_Linear, na.rm = TRUE)
  ),
  Spline_LMM = c(
    tp_spline,
    fp_spline,
    sprintf("%.1f%%", sensitivity_spline * 100),
    sprintf("%.1f%%", fdr_spline * 100),
    sum(results_df$Sig_Spline, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

print(comparison_table)
cat("\n")

# Power improvement
power_improvement <- (sensitivity_spline - sensitivity_linear) / max(0.01, sensitivity_linear)
cat("Power Improvement:", sprintf("%.1f%%", power_improvement * 100), "\n\n")

# ==============================================================================
# 7. Visualization: P-value Comparison
# ==============================================================================

# Scatter plot of p-values
plot_pvalues <- results_df %>%
  filter(!is.na(P_Linear) & !is.na(P_Spline)) %>%
  mutate(
    P_Linear_plot = pmax(P_Linear, 1e-10),
    P_Spline_plot = pmax(P_Spline, 1e-10)
  )

p_comparison <- ggplot(plot_pvalues, aes(x = P_Linear_plot, y = P_Spline_plot,
                                          color = TrueSignal)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 0.05, linetype = "dotted", color = "red") +
  scale_x_log10(
    breaks = c(1e-10, 1e-5, 0.01, 0.05, 0.1, 1),
    labels = c("1e-10", "1e-5", "0.01", "0.05", "0.1", "1")
  ) +
  scale_y_log10(
    breaks = c(1e-10, 1e-5, 0.01, 0.05, 0.1, 1),
    labels = c("1e-10", "1e-5", "0.01", "0.05", "0.1", "1")
  ) +
  scale_color_manual(
    values = c("TRUE" = "#e74c3c", "FALSE" = "#95a5a6"),
    labels = c("TRUE" = "True Signal", "FALSE" = "Null")
  ) +
  labs(
    title = "P-value Comparison: Linear vs Spline LMM",
    subtitle = paste0("Signal Features: ", n_signal_features, " | Null Features: ",
                     n_features - n_signal_features),
    x = "P-value (Linear LMM)",
    y = "P-value (Spline LMM)",
    color = "Feature Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(face = "bold")
  ) +
  annotate("text", x = 0.1, y = 1e-8, label = "Spline Wins",
           color = "#27ae60", fontface = "bold", size = 5) +
  annotate("text", x = 1e-8, y = 0.1, label = "Linear Wins",
           color = "#2980b9", fontface = "bold", size = 5)

ggsave("pilot_experiments/results/01b_pvalue_comparison.pdf",
       p_comparison, width = 10, height = 8)
cat("✓ Saved p-value comparison plot\n\n")

# ==============================================================================
# 8. Visualization: ROC-style Power Curve
# ==============================================================================

# Calculate power at different thresholds
thresholds <- c(0.001, 0.01, 0.05, 0.1, 0.2)

power_data <- lapply(thresholds, function(th) {
  data.frame(
    Threshold = th,
    Method = c("Linear LMM", "Spline LMM"),
    Power = c(
      sum(results_df$TrueSignal & results_df$P_Linear < th, na.rm = TRUE) / n_signal_features,
      sum(results_df$TrueSignal & results_df$P_Spline < th, na.rm = TRUE) / n_signal_features
    ),
    FDR = c(
      sum(!results_df$TrueSignal & results_df$P_Linear < th, na.rm = TRUE) /
        max(1, sum(results_df$P_Linear < th, na.rm = TRUE)),
      sum(!results_df$TrueSignal & results_df$P_Spline < th, na.rm = TRUE) /
        max(1, sum(results_df$P_Spline < th, na.rm = TRUE))
    )
  )
}) %>% bind_rows()

p_power <- ggplot(power_data, aes(x = Threshold, y = Power, color = Method, group = Method)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  scale_x_log10(breaks = thresholds) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("Linear LMM" = "#3498db", "Spline LMM" = "#e74c3c")) +
  labs(
    title = "Statistical Power Comparison",
    subtitle = "Ability to detect true non-linear signals",
    x = "Significance Threshold (α)",
    y = "Power (True Positive Rate)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

p_fdr <- ggplot(power_data, aes(x = Threshold, y = FDR, color = Method, group = Method)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_x_log10(breaks = thresholds) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("Linear LMM" = "#3498db", "Spline LMM" = "#e74c3c")) +
  labs(
    title = "False Discovery Rate Comparison",
    subtitle = "Type I error control",
    x = "Significance Threshold (α)",
    y = "False Discovery Rate"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

p_combined <- p_power / p_fdr
ggsave("pilot_experiments/results/01c_power_fdr_curves.pdf",
       p_combined, width = 10, height = 10)
cat("✓ Saved power/FDR curves\n\n")

# ==============================================================================
# 9. Save Results
# ==============================================================================

# Save detailed results
write.csv(results_df,
          "pilot_experiments/results/01_detailed_results.csv",
          row.names = FALSE)

# Save summary
summary_results <- list(
  simulation_params = list(
    n_subjects_per_group = n_subjects_per_group,
    n_timepoints = n_timepoints,
    n_features = n_features,
    n_signal_features = n_signal_features,
    signal_strength = signal_strength
  ),
  performance = comparison_table,
  power_improvement_pct = power_improvement * 100
)

saveRDS(summary_results, "pilot_experiments/results/01_summary.rds")

cat("=== Experiment Complete ===\n\n")
cat("Key Findings:\n")
cat("  - Spline LMM detected", tp_spline, "true signals\n")
cat("  - Linear LMM detected only", tp_linear, "true signals\n")
cat("  - Power improvement:", sprintf("%.1f%%", power_improvement * 100), "\n")
cat("  - FDR (Spline):", sprintf("%.1f%%", fdr_spline * 100), "\n\n")

if (sensitivity_spline > sensitivity_linear * 1.5 && fdr_spline < 0.2) {
  cat("✓✓✓ VALIDATION SUCCESSFUL ✓✓✓\n")
  cat("Spline-based trajectory modeling shows substantial improvement!\n")
  cat("This approach is ready for full implementation in LinDA v2.\n")
} else {
  cat("⚠ WARNING: Results are inconclusive.\n")
  cat("Consider adjusting simulation parameters or signal patterns.\n")
}
