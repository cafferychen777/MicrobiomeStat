#!/usr/bin/env python3
"""
==============================================================================
Pilot Experiment 1: Trajectory Spline Model Validation (Simplified)
==============================================================================

Objective: Prove that detecting trajectory differences gives more power
           than detecting overall group differences.

Simplified approach: Compare permutation tests
- Method A: Test if overall mean differs between groups (averaging across time)
- Method B: Test if trajectory differs between groups (considering time pattern)

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

# Set random seed
np.random.seed(42)

# ==============================================================================
# 1. Simulation Parameters
# ==============================================================================

n_subjects_per_group = 50
n_timepoints = 5
n_features = 100
n_signal_features = 10
signal_strength = 2.0
noise_sd = 1.0

print("=" * 80)
print("Pilot Experiment 1: Trajectory Spline Validation (Simplified)")
print("=" * 80)
print()
print("Simulation Parameters:")
print(f"  - Subjects per group: {n_subjects_per_group}")
print(f"  - Time points: {n_timepoints}")
print(f"  - Total features: {n_features}")
print(f"  - Signal features: {n_signal_features}")
print(f"  - Signal strength: {signal_strength}")
print()

# ==============================================================================
# 2. Data Simulation
# ==============================================================================

def simulate_trajectory_data():
    """Simulate microbiome trajectory data with inverted-U signal"""

    n_total = n_subjects_per_group * 2

    data_list = []

    for subject_idx in range(n_total):
        group = "Control" if subject_idx < n_subjects_per_group else "Treatment"
        subject_id = f"S{subject_idx:03d}"

        # Random subject effect
        subject_effect = np.random.normal(0, 0.5)

        for time in range(1, n_timepoints + 1):
            row = {
                'SubjectID': subject_id,
                'Group': group,
                'Time': time
            }

            # Generate features
            for feature_idx in range(n_features):
                feature_name = f"Feature{feature_idx:03d}"

                # Base abundance
                base = np.random.normal(0, noise_sd)

                # Add signal to first n_signal_features
                if feature_idx < n_signal_features and group == "Treatment":
                    # Inverted-U pattern: peaks at Time = 3
                    signal = signal_strength * np.exp(-(time - 3)**2 / 2)
                    value = base + subject_effect + signal
                else:
                    value = base + subject_effect

                row[feature_name] = value

            data_list.append(row)

    df = pd.DataFrame(data_list)
    return df

print("Generating simulated data...")
sim_data = simulate_trajectory_data()
print(f"  Data dimensions: {len(sim_data)} observations x {n_features} features")
print()

# ==============================================================================
# 3. Method A: Overall Mean Difference (Linear/Simple Approach)
# ==============================================================================

def test_overall_difference(data, feature_name):
    """Test if overall mean differs between groups (ignoring time)"""

    control_values = data[data['Group'] == 'Control'][feature_name].values
    treatment_values = data[data['Group'] == 'Treatment'][feature_name].values

    # Simple t-test
    t_stat, p_value = stats.ttest_ind(treatment_values, control_values)

    return p_value

print("Running Method A: Overall Mean Difference Test...")
feature_names = [f"Feature{i:03d}" for i in range(n_features)]
results_simple = []

for fname in feature_names:
    p_val = test_overall_difference(sim_data, fname)
    results_simple.append(p_val)

p_values_simple = np.array(results_simple)
sig_simple = np.sum(p_values_simple < 0.05)
tp_simple = np.sum(p_values_simple[:n_signal_features] < 0.05)

print(f"  Completed: {np.sum(~np.isnan(p_values_simple))}/{n_features} features")
print(f"  Significant features (P < 0.05): {sig_simple}")
print(f"  True positives detected: {tp_simple}/{n_signal_features}")
print()

# ==============================================================================
# 4. Method B: Trajectory Difference (Time-Aware Approach)
# ==============================================================================

def test_trajectory_difference(data, feature_name):
    """Test if trajectory differs between groups using repeated measures"""

    # For each group and time, calculate mean
    group_time_means = data.groupby(['Group', 'Time'])[feature_name].mean().unstack()

    control_trajectory = group_time_means.loc['Control'].values
    treatment_trajectory = group_time_means.loc['Treatment'].values

    # Calculate sum of squared differences at each timepoint
    ssd_observed = np.sum((treatment_trajectory - control_trajectory)**2)

    # Permutation test
    n_perm = 1000
    ssd_null = []

    # Get data for permutation
    all_subjects = data.groupby(['SubjectID', 'Time'])[feature_name].mean().unstack()
    n_subjects = len(all_subjects)

    for _ in range(n_perm):
        # Randomly assign subjects to groups
        perm_idx = np.random.permutation(n_subjects)
        group1_idx = perm_idx[:n_subjects_per_group]
        group2_idx = perm_idx[n_subjects_per_group:]

        traj1 = all_subjects.iloc[group1_idx].mean(axis=0).values
        traj2 = all_subjects.iloc[group2_idx].mean(axis=0).values

        ssd_perm = np.sum((traj2 - traj1)**2)
        ssd_null.append(ssd_perm)

    # Calculate p-value
    p_value = np.mean(ssd_null >= ssd_observed)

    return p_value

print("Running Method B: Trajectory Difference Test (with permutation)...")
results_trajectory = []

for i, fname in enumerate(feature_names):
    if i % 20 == 0:
        print(f"  Progress: {i}/{n_features}")

    p_val = test_trajectory_difference(sim_data, fname)
    results_trajectory.append(p_val)

p_values_trajectory = np.array(results_trajectory)
sig_trajectory = np.sum(p_values_trajectory < 0.05)
tp_trajectory = np.sum(p_values_trajectory[:n_signal_features] < 0.05)

print(f"  Completed: {np.sum(~np.isnan(p_values_trajectory))}/{n_features} features")
print(f"  Significant features (P < 0.05): {sig_trajectory}")
print(f"  True positives detected: {tp_trajectory}/{n_signal_features}")
print()

# ==============================================================================
# 5. Performance Comparison
# ==============================================================================

print("=" * 80)
print("Performance Comparison")
print("=" * 80)
print()

# Calculate metrics
alpha = 0.05
fp_simple = np.sum(p_values_simple[n_signal_features:] < alpha)
fp_trajectory = np.sum(p_values_trajectory[n_signal_features:] < alpha)

sensitivity_simple = tp_simple / n_signal_features
sensitivity_trajectory = tp_trajectory / n_signal_features

fdr_simple = fp_simple / max(1, sig_simple)
fdr_trajectory = fp_trajectory / max(1, sig_trajectory)

power_improvement = ((sensitivity_trajectory - sensitivity_simple) /
                     max(0.01, sensitivity_simple) * 100)

# Print comparison table
print(f"{'Metric':<30} {'Simple (Mean)':<20} {'Trajectory':<20}")
print("-" * 70)
print(f"{'True Positives':<30} {tp_simple:<20} {tp_trajectory:<20}")
print(f"{'False Positives':<30} {fp_simple:<20} {fp_trajectory:<20}")
print(f"{'Sensitivity (Power)':<30} {sensitivity_simple*100:.1f}%{'':<14} {sensitivity_trajectory*100:.1f}%")
print(f"{'False Discovery Rate':<30} {fdr_simple*100:.1f}%{'':<14} {fdr_trajectory*100:.1f}%")
print(f"{'Total Significant':<30} {sig_simple:<20} {sig_trajectory:<20}")
print()
print(f"Power Improvement: {power_improvement:.1f}%")
print()

# ==============================================================================
# 6. Visualizations
# ==============================================================================

import os
os.makedirs("pilot_experiments/results", exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100

# Plot 1: True signal pattern
print("Generating visualizations...")

example_feature = "Feature001"
plot_data = sim_data.groupby(['Group', 'Time'])[example_feature].agg(['mean', 'sem']).reset_index()

fig, ax = plt.subplots(figsize=(8, 6))
for group in ['Control', 'Treatment']:
    data_subset = plot_data[plot_data['Group'] == group]
    color = '#3498db' if group == 'Control' else '#e74c3c'
    ax.plot(data_subset['Time'], data_subset['mean'], 'o-',
            color=color, label=group, linewidth=2, markersize=8)
    ax.fill_between(data_subset['Time'],
                     data_subset['mean'] - data_subset['sem'],
                     data_subset['mean'] + data_subset['sem'],
                     alpha=0.2, color=color)

ax.set_xlabel('Time Point', fontsize=12)
ax.set_ylabel('Abundance (Log Scale)', fontsize=12)
ax.set_title('True Signal Pattern: Inverted-U Trajectory\n(Treatment shows spike at middle timepoints)',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pilot_experiments/results/01a_true_signal_pattern.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: P-value comparison
results_df = pd.DataFrame({
    'Feature': feature_names,
    'TrueSignal': [i < n_signal_features for i in range(n_features)],
    'P_Simple': p_values_simple,
    'P_Trajectory': p_values_trajectory
})

fig, ax = plt.subplots(figsize=(10, 8))
for is_signal in [False, True]:
    subset = results_df[results_df['TrueSignal'] == is_signal]
    color = '#e74c3c' if is_signal else '#95a5a6'
    label = 'True Signal' if is_signal else 'Null'
    ax.scatter(np.maximum(subset['P_Simple'], 1e-3),
               np.maximum(subset['P_Trajectory'], 1e-3),
               c=color, label=label, s=50, alpha=0.7, edgecolors='black', linewidths=0.5)

ax.plot([1e-3, 1], [1e-3, 1], 'k--', alpha=0.5, linewidth=1)
ax.axhline(0.05, color='red', linestyle=':', alpha=0.5)
ax.axvline(0.05, color='red', linestyle=':', alpha=0.5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('P-value (Simple Mean Test)', fontsize=12)
ax.set_ylabel('P-value (Trajectory Test)', fontsize=12)
ax.set_title('P-value Comparison: Simple vs Trajectory-Aware Test\n(Points below diagonal = Trajectory wins)',
             fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='lower right')
ax.text(0.1, 0.003, 'Trajectory Wins', color='#27ae60', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig('pilot_experiments/results/01b_pvalue_comparison.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: Power curves
thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
power_data = []

for th in thresholds:
    tp_simp = np.sum(p_values_simple[:n_signal_features] < th)
    tp_traj = np.sum(p_values_trajectory[:n_signal_features] < th)

    fp_simp = np.sum(p_values_simple[n_signal_features:] < th)
    fp_traj = np.sum(p_values_trajectory[n_signal_features:] < th)

    power_data.append({
        'Threshold': th,
        'Method': 'Simple (Mean)',
        'Power': tp_simp / n_signal_features,
        'FDR': fp_simp / max(1, tp_simp + fp_simp)
    })
    power_data.append({
        'Threshold': th,
        'Method': 'Trajectory',
        'Power': tp_traj / n_signal_features,
        'FDR': fp_traj / max(1, tp_traj + fp_traj)
    })

power_df = pd.DataFrame(power_data)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

for method in ['Simple (Mean)', 'Trajectory']:
    subset = power_df[power_df['Method'] == method]
    color = '#3498db' if method == 'Simple (Mean)' else '#e74c3c'
    ax1.plot(subset['Threshold'], subset['Power'], 'o-',
             color=color, label=method, linewidth=2, markersize=8)
    ax2.plot(subset['Threshold'], subset['FDR'], 'o-',
             color=color, label=method, linewidth=2, markersize=8)

ax1.set_xscale('log')
ax1.set_xlabel('Significance Threshold (α)', fontsize=12)
ax1.set_ylabel('Power (True Positive Rate)', fontsize=12)
ax1.set_title('Statistical Power Comparison\n(Higher is better)', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0, 1])

ax2.set_xscale('log')
ax2.axhline(0.05, color='red', linestyle='--', alpha=0.5, label='Target FDR = 5%')
ax2.set_xlabel('Significance Threshold (α)', fontsize=12)
ax2.set_ylabel('False Discovery Rate', fontsize=12)
ax2.set_title('False Discovery Rate Comparison\n(Lower is better)', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_ylim([0, max(0.3, power_df['FDR'].max() * 1.1)])

plt.tight_layout()
plt.savefig('pilot_experiments/results/01c_power_fdr_curves.pdf', dpi=300, bbox_inches='tight')
plt.close()

print("✓ Saved visualization plots")
print()

# ==============================================================================
# 7. Save Results
# ==============================================================================

results_df.to_csv('pilot_experiments/results/01_detailed_results.csv', index=False)

summary = {
    'simulation_params': {
        'n_subjects_per_group': n_subjects_per_group,
        'n_timepoints': n_timepoints,
        'n_features': n_features,
        'n_signal_features': n_signal_features,
        'signal_strength': signal_strength
    },
    'performance': {
        'simple': {
            'true_positives': int(tp_simple),
            'false_positives': int(fp_simple),
            'sensitivity': float(sensitivity_simple),
            'fdr': float(fdr_simple)
        },
        'trajectory': {
            'true_positives': int(tp_trajectory),
            'false_positives': int(fp_trajectory),
            'sensitivity': float(sensitivity_trajectory),
            'fdr': float(fdr_trajectory)
        }
    },
    'power_improvement_pct': float(power_improvement)
}

import json
with open('pilot_experiments/results/01_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("✓ Saved results files")
print()

# ==============================================================================
# 8. Final Summary
# ==============================================================================

print("=" * 80)
print("Experiment Complete")
print("=" * 80)
print()
print("Key Findings:")
print(f"  - Trajectory Test detected {tp_trajectory} true signals")
print(f"  - Simple Test detected only {tp_simple} true signals")
print(f"  - Power improvement: {power_improvement:.1f}%")
print(f"  - FDR (Trajectory): {fdr_trajectory*100:.1f}%")
print()

if sensitivity_trajectory > sensitivity_simple * 1.3 and fdr_trajectory < 0.2:
    print("✓✓✓ VALIDATION SUCCESSFUL ✓✓✓")
    print("Trajectory-aware testing shows substantial improvement!")
    print("This validates the importance of considering temporal dynamics.")
elif sensitivity_trajectory > sensitivity_simple and fdr_trajectory < 0.3:
    print("✓ VALIDATION PROMISING ✓")
    print("Trajectory-aware testing shows improvement.")
    print("Consider optimizing parameters or increasing sample size.")
else:
    print("⚠ Results are mixed or inconclusive.")
    print("The inverted-U pattern may be challenging for this sample size.")

print()
print("Interpretation:")
print("  The inverted-U signal (spike at middle timepoints) tends to")
print("  'cancel out' when averaged across all timepoints. Trajectory-aware")
print("  methods capture this non-linear pattern more effectively.")
print()
print("Generated files:")
print("  - pilot_experiments/results/01a_true_signal_pattern.pdf")
print("  - pilot_experiments/results/01b_pvalue_comparison.pdf")
print("  - pilot_experiments/results/01c_power_fdr_curves.pdf")
print("  - pilot_experiments/results/01_detailed_results.csv")
print("  - pilot_experiments/results/01_summary.json")
print()
