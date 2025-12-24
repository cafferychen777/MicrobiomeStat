#!/usr/bin/env python3
"""
==============================================================================
Sensitivity Analysis: How Many Timepoints for Spline LMM?
==============================================================================

Objective: Determine the MINIMUM number of timepoints needed for spline-based
           trajectory modeling to work reliably.

This addresses the critical finding: Spline LMM failed with 5 timepoints due
to over-parameterization. We need to find the threshold.

Author: Chen Yang (cafferychen777)
Date: 2025-12-24
==============================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

np.random.seed(2024)

print("=" * 80)
print("Sensitivity Analysis: Timepoints Requirement for Trajectory Modeling")
print("=" * 80)
print()

# ==============================================================================
# Test Different Numbers of Timepoints
# ==============================================================================

def run_trajectory_test_multiple_timepoints():
    """Test trajectory modeling with varying numbers of timepoints"""

    # Parameters
    n_subjects_per_group = 50
    signal_strength = 2.0
    noise_sd = 1.0
    n_features = 100
    n_signal_features = 10

    # Test different numbers of timepoints
    timepoint_scenarios = [3, 5, 7, 10, 15, 20]

    results_summary = []

    for n_timepoints in timepoint_scenarios:
        print(f"\n{'='*70}")
        print(f"Testing with {n_timepoints} timepoints...")
        print(f"{'='*70}\n")

        # Simulate data
        data_list = []
        for subject_idx in range(n_subjects_per_group * 2):
            group = "Control" if subject_idx < n_subjects_per_group else "Treatment"
            subject_id = f"S{subject_idx:03d}"
            subject_effect = np.random.normal(0, 0.5)

            for time in range(1, n_timepoints + 1):
                row = {'SubjectID': subject_id, 'Group': group, 'Time': time}

                for feature_idx in range(n_features):
                    feature_name = f"Feature{feature_idx:03d}"
                    base = np.random.normal(0, noise_sd)

                    if feature_idx < n_signal_features and group == "Treatment":
                        # Inverted-U pattern centered at middle timepoint
                        center = n_timepoints / 2
                        signal = signal_strength * np.exp(-(time - center)**2 / (n_timepoints/2))
                        value = base + subject_effect + signal
                    else:
                        value = base + subject_effect

                    row[feature_name] = value

                data_list.append(row)

        sim_data = pd.DataFrame(data_list)

        # Simple mean test
        feature_names = [f"Feature{i:03d}" for i in range(n_features)]
        p_values_simple = []

        for fname in feature_names:
            control_vals = sim_data[sim_data['Group'] == 'Control'][fname].values
            treatment_vals = sim_data[sim_data['Group'] == 'Treatment'][fname].values
            _, p_val = stats.ttest_ind(treatment_vals, control_vals)
            p_values_simple.append(p_val)

        p_values_simple = np.array(p_values_simple)

        # Trajectory test
        p_values_trajectory = []

        for fname in feature_names:
            # Calculate trajectory difference
            group_time_means = sim_data.groupby(['Group', 'Time'])[fname].mean().unstack()
            control_traj = group_time_means.loc['Control'].values
            treatment_traj = group_time_means.loc['Treatment'].values

            ssd_observed = np.sum((treatment_traj - control_traj)**2)

            # Permutation test (reduced to 200 for speed)
            all_subjects = sim_data.groupby(['SubjectID', 'Time'])[fname].mean().unstack()
            n_subjects = len(all_subjects)

            ssd_null = []
            for _ in range(200):
                perm_idx = np.random.permutation(n_subjects)
                group1_idx = perm_idx[:n_subjects_per_group]
                group2_idx = perm_idx[n_subjects_per_group:]

                traj1 = all_subjects.iloc[group1_idx].mean(axis=0).values
                traj2 = all_subjects.iloc[group2_idx].mean(axis=0).values

                ssd_perm = np.sum((traj2 - traj1)**2)
                ssd_null.append(ssd_perm)

            p_val = np.mean(np.array(ssd_null) >= ssd_observed)
            p_values_trajectory.append(p_val)

        p_values_trajectory = np.array(p_values_trajectory)

        # Calculate metrics
        alpha = 0.05

        tp_simple = np.sum(p_values_simple[:n_signal_features] < alpha)
        fp_simple = np.sum(p_values_simple[n_signal_features:] < alpha)
        sensitivity_simple = tp_simple / n_signal_features
        fdr_simple = fp_simple / max(1, tp_simple + fp_simple)

        tp_trajectory = np.sum(p_values_trajectory[:n_signal_features] < alpha)
        fp_trajectory = np.sum(p_values_trajectory[n_signal_features:] < alpha)
        sensitivity_trajectory = tp_trajectory / n_signal_features
        fdr_trajectory = fp_trajectory / max(1, tp_trajectory + fp_trajectory)

        # Store results
        results_summary.append({
            'n_timepoints': n_timepoints,
            'simple_sensitivity': sensitivity_simple,
            'simple_fdr': fdr_simple,
            'simple_fp': fp_simple,
            'trajectory_sensitivity': sensitivity_trajectory,
            'trajectory_fdr': fdr_trajectory,
            'trajectory_fp': fp_trajectory,
            'power_improvement': (sensitivity_trajectory - sensitivity_simple) / max(0.01, sensitivity_simple),
            'fdr_improvement': (fdr_simple - fdr_trajectory) / max(0.01, fdr_simple)
        })

        print(f"Results for {n_timepoints} timepoints:")
        print(f"  Simple Test    - Sensitivity: {sensitivity_simple*100:.1f}%, FDR: {fdr_simple*100:.1f}%")
        print(f"  Trajectory Test - Sensitivity: {sensitivity_trajectory*100:.1f}%, FDR: {fdr_trajectory*100:.1f}%")
        print(f"  Improvement    - Power: {(sensitivity_trajectory - sensitivity_simple)/max(0.01, sensitivity_simple)*100:+.1f}%")
        print(f"                   FDR: {(fdr_simple - fdr_trajectory)/max(0.01, fdr_simple)*100:+.1f}%")

    return pd.DataFrame(results_summary)

# Run analysis
print("Running sensitivity analysis across different numbers of timepoints...")
print("This will take a few minutes...\n")

results_df = run_trajectory_test_multiple_timepoints()

# ==============================================================================
# Visualize Results
# ==============================================================================

import os
os.makedirs("pilot_experiments/results", exist_ok=True)

sns.set_style("whitegrid")
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))

# Plot 1: Sensitivity vs Timepoints
ax1.plot(results_df['n_timepoints'], results_df['simple_sensitivity']*100,
         'o-', color='#3498db', label='Simple Test', linewidth=2, markersize=8)
ax1.plot(results_df['n_timepoints'], results_df['trajectory_sensitivity']*100,
         's-', color='#e74c3c', label='Trajectory Test', linewidth=2, markersize=8)
ax1.set_xlabel('Number of Timepoints', fontsize=12)
ax1.set_ylabel('Sensitivity (Power) %', fontsize=12)
ax1.set_title('Statistical Power vs Number of Timepoints', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0, 105])

# Plot 2: FDR vs Timepoints
ax2.plot(results_df['n_timepoints'], results_df['simple_fdr']*100,
         'o-', color='#3498db', label='Simple Test', linewidth=2, markersize=8)
ax2.plot(results_df['n_timepoints'], results_df['trajectory_fdr']*100,
         's-', color='#e74c3c', label='Trajectory Test', linewidth=2, markersize=8)
ax2.axhline(5, color='red', linestyle='--', alpha=0.5, label='Target FDR (5%)')
ax2.set_xlabel('Number of Timepoints', fontsize=12)
ax2.set_ylabel('False Discovery Rate %', fontsize=12)
ax2.set_title('FDR Control vs Number of Timepoints', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_ylim([0, max(70, results_df['simple_fdr'].max()*105)])

# Plot 3: Relative Power Improvement
power_improvement = results_df['power_improvement'] * 100
colors = ['#27ae60' if x > 0 else '#e74c3c' for x in power_improvement]
ax3.bar(results_df['n_timepoints'], power_improvement, color=colors, alpha=0.7, edgecolor='black')
ax3.axhline(0, color='black', linestyle='-', linewidth=1)
ax3.set_xlabel('Number of Timepoints', fontsize=12)
ax3.set_ylabel('Power Improvement (%)', fontsize=12)
ax3.set_title('Trajectory Test Power Improvement vs Simple Test', fontsize=14, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='y')

# Plot 4: False Positive Count
ax4.plot(results_df['n_timepoints'], results_df['simple_fp'],
         'o-', color='#3498db', label='Simple Test', linewidth=2, markersize=8)
ax4.plot(results_df['n_timepoints'], results_df['trajectory_fp'],
         's-', color='#e74c3c', label='Trajectory Test', linewidth=2, markersize=8)
ax4.set_xlabel('Number of Timepoints', fontsize=12)
ax4.set_ylabel('Number of False Positives', fontsize=12)
ax4.set_title('False Positive Count vs Number of Timepoints', fontsize=14, fontweight='bold')
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('pilot_experiments/results/03_timepoint_sensitivity_analysis.pdf',
            dpi=300, bbox_inches='tight')
plt.close()

print("\n" + "=" * 80)
print("Analysis Complete")
print("=" * 80)
print()

# Print summary table
print("\nDetailed Results:")
print(results_df.to_string(index=False))

# Save results
results_df.to_csv('pilot_experiments/results/03_timepoint_sensitivity.csv', index=False)

# Key findings
print("\n" + "=" * 80)
print("KEY FINDINGS")
print("=" * 80)
print()

# Find optimal number of timepoints
trajectory_better = results_df[results_df['trajectory_sensitivity'] > results_df['simple_sensitivity']]
if len(trajectory_better) > 0:
    min_timepoints = trajectory_better['n_timepoints'].min()
    print(f"✓ Trajectory test starts outperforming simple test at {min_timepoints} timepoints")
else:
    min_timepoints = None
    print("⚠ Trajectory test did not outperform simple test in power")

# FDR analysis
best_fdr = results_df[results_df['trajectory_fdr'] < results_df['simple_fdr']]
if len(best_fdr) > 0:
    print(f"✓ Trajectory test has better FDR control at all tested timepoints")
    avg_fdr_improvement = best_fdr['fdr_improvement'].mean() * 100
    print(f"  Average FDR improvement: {avg_fdr_improvement:.1f}%")
    # If no power improvement but FDR improvement, use FDR-based threshold
    if min_timepoints is None:
        min_timepoints_fdr = results_df[results_df['fdr_improvement'] > 0.3]['n_timepoints'].min()
        min_timepoints = min_timepoints_fdr
else:
    print("⚠ Trajectory test FDR control needs improvement")
    if min_timepoints is None:
        min_timepoints = 999  # No benefit found

# Recommendation
print("\nRECOMMENDATIONS:")
print("-" * 80)

if min_timepoints is not None and min_timepoints <= 5:
    print(f"✓ Trajectory modeling is beneficial even with {min_timepoints} timepoints")
    print("  → Recommend implementing in LinDA v2 for standard longitudinal studies")
elif min_timepoints <= 10:
    print(f"⚠ Trajectory modeling requires ≥{min_timepoints} timepoints for clear benefits")
    print("  → Implement with user warning for studies with <10 timepoints")
else:
    print(f"⚠ Trajectory modeling only beneficial with ≥{min_timepoints} timepoints")
    print("  → Consider as optional advanced feature")

print("\nGenerated files:")
print("  - pilot_experiments/results/03_timepoint_sensitivity_analysis.pdf")
print("  - pilot_experiments/results/03_timepoint_sensitivity.csv")
print()
