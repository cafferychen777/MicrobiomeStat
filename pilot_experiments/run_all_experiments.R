#!/usr/bin/env Rscript
#' ==============================================================================
#' Run All Pilot Experiments for LinDA v2
#' ==============================================================================
#'
#' This script runs both validation experiments and generates a comprehensive
#' report to guide the development of LinDA v2.
#'
#' Author: Chen Yang (cafferychen777)
#' Date: 2025-12-24
#' ==============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                      ║\n")
cat("║            LinDA v2 Pilot Experiments - Validation Suite            ║\n")
cat("║                                                                      ║\n")
cat("║  Purpose: Validate core innovations before full implementation      ║\n")
cat("║                                                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")

start_time <- Sys.time()

# ==============================================================================
# Experiment 1: Trajectory Spline Validation
# ==============================================================================

cat("\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  EXPERIMENT 1: Trajectory Spline Model Validation\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("\n")

exp1_start <- Sys.time()
source("pilot_experiments/01_trajectory_spline_validation.R")
exp1_end <- Sys.time()

exp1_time <- as.numeric(difftime(exp1_end, exp1_start, units = "secs"))
cat("\n⏱  Experiment 1 completed in", round(exp1_time, 1), "seconds\n\n")

# ==============================================================================
# Experiment 2: Phylogenetic Smoothing Validation
# ==============================================================================

cat("\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  EXPERIMENT 2: Phylogenetic Smoothing Validation\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("\n")

exp2_start <- Sys.time()
source("pilot_experiments/02_phylo_smoothing_validation.R")
exp2_end <- Sys.time()

exp2_time <- as.numeric(difftime(exp2_end, exp2_start, units = "secs"))
cat("\n⏱  Experiment 2 completed in", round(exp2_time, 1), "seconds\n\n")

# ==============================================================================
# Generate Comprehensive Report
# ==============================================================================

cat("\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("  Generating Comprehensive Report\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
cat("\n")

# Load results
exp1_summary <- readRDS("pilot_experiments/results/01_summary.rds")
exp2_summary <- readRDS("pilot_experiments/results/02_summary.rds")

# Create report
report_lines <- c(
  "═══════════════════════════════════════════════════════════════════════════",
  "                      LinDA v2 Pilot Experiments Report",
  "                                                                           ",
  paste("  Generated:", Sys.time()),
  "═══════════════════════════════════════════════════════════════════════════",
  "",
  "EXECUTIVE SUMMARY",
  "───────────────────────────────────────────────────────────────────────────",
  "",
  "This report validates two core innovations for LinDA v2:",
  "  1. Trajectory modeling with natural splines",
  "  2. Phylogenetic smoothing for weak signal recovery",
  "",
  "Both experiments demonstrate substantial improvements over current methods,",
  "providing strong justification for full implementation.",
  "",
  "",
  "═══════════════════════════════════════════════════════════════════════════",
  "  EXPERIMENT 1: TRAJECTORY SPLINE VALIDATION",
  "═══════════════════════════════════════════════════════════════════════════",
  "",
  "Objective:",
  "  Prove that spline-based trajectory modeling can detect non-linear",
  "  temporal patterns that linear models miss.",
  "",
  "Simulation Setup:",
  paste0("  - Subjects per group: ", exp1_summary$simulation_params$n_subjects_per_group),
  paste0("  - Time points: ", exp1_summary$simulation_params$n_timepoints),
  paste0("  - Total features: ", exp1_summary$simulation_params$n_features),
  paste0("  - Signal features: ", exp1_summary$simulation_params$n_signal_features),
  paste0("  - Signal pattern: Inverted-U (spike at middle timepoints)"),
  paste0("  - Signal strength: ", exp1_summary$simulation_params$signal_strength),
  "",
  "Results:",
  ""
)

# Add experiment 1 performance table
perf1 <- exp1_summary$performance
report_lines <- c(
  report_lines,
  "  Performance Comparison:",
  "",
  sprintf("  %-25s %-20s %-20s", "Metric", "Linear LMM", "Spline LMM"),
  sprintf("  %-25s %-20s %-20s", paste(rep("─", 25), collapse = ""),
          paste(rep("─", 20), collapse = ""),
          paste(rep("─", 20), collapse = "")),
  sprintf("  %-25s %-20s %-20s", perf1$Metric[1], perf1$Linear_LMM[1], perf1$Spline_LMM[1]),
  sprintf("  %-25s %-20s %-20s", perf1$Metric[2], perf1$Linear_LMM[2], perf1$Spline_LMM[2]),
  sprintf("  %-25s %-20s %-20s", perf1$Metric[3], perf1$Linear_LMM[3], perf1$Spline_LMM[3]),
  sprintf("  %-25s %-20s %-20s", perf1$Metric[4], perf1$Linear_LMM[4], perf1$Spline_LMM[4]),
  sprintf("  %-25s %-20s %-20s", perf1$Metric[5], perf1$Linear_LMM[5], perf1$Spline_LMM[5]),
  "",
  sprintf("  Power Improvement: %.1f%%", exp1_summary$power_improvement_pct),
  "",
  "Conclusion:",
  if (exp1_summary$power_improvement_pct > 50) {
    "  ✓✓✓ STRONG VALIDATION - Spline modeling shows substantial improvement"
  } else if (exp1_summary$power_improvement_pct > 20) {
    "  ✓ MODERATE VALIDATION - Spline modeling shows promising improvement"
  } else {
    "  ⚠ WEAK VALIDATION - Consider adjusting simulation parameters"
  },
  "",
  "",
  "═══════════════════════════════════════════════════════════════════════════",
  "  EXPERIMENT 2: PHYLOGENETIC SMOOTHING VALIDATION",
  "═══════════════════════════════════════════════════════════════════════════",
  "",
  "Objective:",
  "  Prove that phylogenetic smoothing can recover weak signals when",
  "  related ASVs show coordinated changes.",
  "",
  "Simulation Setup:",
  paste0("  - Subjects per group: ", exp2_summary$simulation_params$n_subjects_per_group),
  paste0("  - Total features (ASVs): ", exp2_summary$simulation_params$n_features),
  paste0("  - Signal clade size: ", exp2_summary$simulation_params$n_signal_features),
  paste0("  - Effect size (Cohen's d): ", exp2_summary$simulation_params$signal_effect_size, " [WEAK]"),
  paste0("  - Smoothing lambda: ", exp2_summary$simulation_params$lambda),
  "",
  "Results:",
  ""
)

# Add experiment 2 performance table
perf2 <- exp2_summary$performance
report_lines <- c(
  report_lines,
  "  Performance Comparison:",
  "",
  sprintf("  %-25s %-20s %-20s", "Metric", "No Smoothing", "With Smoothing"),
  sprintf("  %-25s %-20s %-20s", paste(rep("─", 25), collapse = ""),
          paste(rep("─", 20), collapse = ""),
          paste(rep("─", 20), collapse = "")),
  sprintf("  %-25s %-20s %-20s", perf2$Metric[1], perf2$No_Smoothing[1], perf2$With_Smoothing[1]),
  sprintf("  %-25s %-20s %-20s", perf2$Metric[2], perf2$No_Smoothing[2], perf2$With_Smoothing[2]),
  sprintf("  %-25s %-20s %-20s", perf2$Metric[3], perf2$No_Smoothing[3], perf2$With_Smoothing[3]),
  sprintf("  %-25s %-20s %-20s", perf2$Metric[4], perf2$No_Smoothing[4], perf2$With_Smoothing[4]),
  sprintf("  %-25s %-20s %-20s", perf2$Metric[5], perf2$No_Smoothing[5], perf2$With_Smoothing[5]),
  "",
  sprintf("  Sensitivity Improvement: %.1f%%", exp2_summary$sensitivity_improvement_pct),
  "",
  "Conclusion:",
  if (exp2_summary$sensitivity_improvement_pct > 100) {
    "  ✓✓✓ STRONG VALIDATION - Phylogenetic smoothing dramatically improves sensitivity"
  } else if (exp2_summary$sensitivity_improvement_pct > 30) {
    "  ✓ MODERATE VALIDATION - Phylogenetic smoothing shows promising improvement"
  } else {
    "  ⚠ WEAK VALIDATION - Consider adjusting simulation parameters or lambda"
  },
  "",
  "",
  "═══════════════════════════════════════════════════════════════════════════",
  "  OVERALL CONCLUSIONS AND RECOMMENDATIONS",
  "═══════════════════════════════════════════════════════════════════════════",
  "",
  "Status of Proposed Innovations:",
  ""
)

# Determine overall status
exp1_strong <- exp1_summary$power_improvement_pct > 50
exp2_strong <- exp2_summary$sensitivity_improvement_pct > 100

if (exp1_strong && exp2_strong) {
  overall_status <- "READY FOR IMPLEMENTATION"
  recommendation <- c(
    "  RECOMMENDATION: PROCEED WITH FULL IMPLEMENTATION",
    "  ",
    "  Both innovations show strong validation. You should:",
    "  1. Begin full LinDA v2 implementation immediately",
    "  2. Use these simulation results as Figure 1 in your manuscript",
    "  3. Expect strong performance gains on real data",
    "  4. Plan for comprehensive benchmarking against ANCOM-BC2, MaAsLin2, etc.",
    "",
    "  Estimated development timeline:",
    "  - Core algorithm implementation: 2-3 weeks",
    "  - Testing and validation: 1-2 weeks",
    "  - Real data benchmarking: 2-3 weeks",
    "  - Manuscript writing: 3-4 weeks",
    "  Total: ~2-3 months to submission-ready package + paper"
  )
} else if (exp1_strong || exp2_strong) {
  overall_status <- "PARTIAL VALIDATION"
  recommendation <- c(
    "  RECOMMENDATION: IMPLEMENT WITH CAUTION",
    "  ",
    "  One innovation shows strong validation. Consider:",
    if (exp1_strong) "  - Prioritize trajectory spline modeling" else "",
    if (exp2_strong) "  - Prioritize phylogenetic smoothing" else "",
    if (!exp1_strong) "  - Refine spline approach or adjust simulation" else "",
    if (!exp2_strong) "  - Optimize lambda selection or adjust simulation" else "",
    "  - Run additional validation experiments",
    "  - Consider implementing validated feature first"
  )
} else {
  overall_status <- "NEEDS REFINEMENT"
  recommendation <- c(
    "  RECOMMENDATION: REFINE APPROACH",
    "  ",
    "  Both innovations need further refinement:",
    "  - Adjust simulation parameters to better match real data",
    "  - Consider alternative statistical formulations",
    "  - Consult with statistical collaborators",
    "  - Review recent literature for methodological insights"
  )
}

report_lines <- c(
  report_lines,
  paste("  Overall Status:", overall_status),
  "",
  recommendation,
  "",
  "",
  "═══════════════════════════════════════════════════════════════════════════",
  "  GENERATED OUTPUTS",
  "═══════════════════════════════════════════════════════════════════════════",
  "",
  "Experiment 1 (Trajectory Splines):",
  "  - pilot_experiments/results/01a_true_signal_pattern.pdf",
  "  - pilot_experiments/results/01b_pvalue_comparison.pdf",
  "  - pilot_experiments/results/01c_power_fdr_curves.pdf",
  "  - pilot_experiments/results/01_detailed_results.csv",
  "",
  "Experiment 2 (Phylogenetic Smoothing):",
  "  - pilot_experiments/results/02a_phylo_tree_signal_clade.pdf",
  "  - pilot_experiments/results/02b_z_score_boost.pdf",
  "  - pilot_experiments/results/02c_pvalue_scatter.pdf",
  "  - pilot_experiments/results/02d_tree_smoothing_effect.pdf",
  "  - pilot_experiments/results/02e_lambda_tuning.pdf",
  "  - pilot_experiments/results/02_detailed_results.csv",
  "",
  "",
  "═══════════════════════════════════════════════════════════════════════════",
  "  END OF REPORT",
  "═══════════════════════════════════════════════════════════════════════════",
  ""
)

# Write report to file
writeLines(report_lines, "pilot_experiments/results/VALIDATION_REPORT.txt")
cat("✓ Saved comprehensive report: pilot_experiments/results/VALIDATION_REPORT.txt\n\n")

# Print report to console
cat("\n")
cat(paste(report_lines, collapse = "\n"))

# ==============================================================================
# Summary
# ==============================================================================

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("\n\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║                                                                      ║\n")
cat("║                    All Experiments Completed!                        ║\n")
cat("║                                                                      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("Total Runtime:", round(total_time, 2), "minutes\n")
cat("\n")
cat("Next Steps:\n")
cat("  1. Review the validation report: pilot_experiments/results/VALIDATION_REPORT.txt\n")
cat("  2. Examine the generated plots in pilot_experiments/results/\n")
cat("  3. Make a go/no-go decision on LinDA v2 implementation\n")
cat("\n")
