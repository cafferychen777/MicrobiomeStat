#!/usr/bin/env python3
"""
==============================================================================
Pilot Experiment 1: Trajectory Spline Model Validation (Python Version)
==============================================================================

Objective: Prove that spline-based trajectory modeling can detect non-linear
           temporal patterns that linear models miss.

Author: Chen Yang (cafferychen777)
Date: 2025-12-24
==============================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.interpolate import BSpline
from scipy.linalg import lstsq
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
print("Pilot Experiment 1: Trajectory Spline Validation")
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

    # Create data structure
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
# 3. Method A: Linear Model
# ==============================================================================

def run_linear_model(data, feature_name):
    """Run linear model: Y ~ Group * Time"""

    # Prepare design matrix
    group_indicator = (data['Group'] == 'Treatment').astype(int).values
    time = data['Time'].values

    # Design matrix: [intercept, group, time, group:time]
    X = np.column_stack([
        np.ones(len(data)),
        group_indicator,
        time,
        group_indicator * time
    ])

    y = data[feature_name].values

    # Fit model using least squares
    try:
        coeffs, residuals, rank, s = lstsq(X, y)

        # Calculate residual variance
        n = len(y)
        p = X.shape[1]
        sigma2 = residuals[0] / (n - p) if len(residuals) > 0 else np.var(y - X @ coeffs)

        # Standard errors
        XtX_inv = np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(XtX_inv) * sigma2)

        # T-statistic for interaction term (group:time)
        interaction_coeff = coeffs[3]
        interaction_se = se[3]
        t_stat = interaction_coeff / interaction_se

        # P-value (two-tailed)
        df = n - p
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df))

        return p_value

    except:
        return np.nan

print("Running Method A: Linear Model...")
feature_names = [f"Feature{i:03d}" for i in range(n_features)]
results_linear = {}

for fname in feature_names:
    results_linear[fname] = run_linear_model(sim_data, fname)

p_values_linear = np.array([results_linear[f] for f in feature_names])
sig_linear = np.sum(p_values_linear < 0.05)
tp_linear = np.sum(p_values_linear[:n_signal_features] < 0.05)

print(f"  Completed: {np.sum(~np.isnan(p_values_linear))}/{n_features} features")
print(f"  Significant features (P < 0.05): {sig_linear}")
print(f"  True positives detected: {tp_linear}/{n_signal_features}")
print()

# ==============================================================================
# 4. Method B: Spline Model (Approximation)
# ==============================================================================

def run_spline_model(data, feature_name):
    """Run spline model: Y ~ Group * spline(Time)"""

    # Create polynomial basis (approximation of splines)
    # Using quadratic terms to capture non-linearity
    group_indicator = (data['Group'] == 'Treatment').astype(int).values
    time = data['Time'].values
    time_centered = time - np.mean(time)

    # Design matrix with polynomial terms
    X = np.column_stack([
        np.ones(len(data)),
        group_indicator,
        time_centered,
        time_centered**2,
        group_indicator * time_centered,
        group_indicator * time_centered**2
    ])

    y = data[feature_name].values

    try:
        # Full model
        coeffs_full, residuals_full, rank, s = lstsq(X, y)
        rss_full = residuals_full[0] if len(residuals_full) > 0 else np.sum((y - X @ coeffs_full)**2)

        # Reduced model (no interaction)
        X_reduced = np.column_stack([
            np.ones(len(data)),
            group_indicator,
            time_centered,
            time_centered**2
        ])
        coeffs_reduced, residuals_reduced, rank, s = lstsq(X_reduced, y)
        rss_reduced = residuals_reduced[0] if len(residuals_reduced) > 0 else np.sum((y - X_reduced @ coeffs_reduced)**2)

        # Likelihood ratio test
        n = len(y)
        df_full = X.shape[1]
        df_reduced = X_reduced.shape[1]
        df_diff = df_full - df_reduced

        f_stat = ((rss_reduced - rss_full) / df_diff) / (rss_full / (n - df_full))
        p_value = 1 - stats.f.cdf(f_stat, df_diff, n - df_full)

        return p_value

    except:
        return np.nan

print("Running Method B: Spline Model...")
results_spline = {}

for fname in feature_names:
    results_spline[fname] = run_spline_model(sim_data, fname)

p_values_spline = np.array([results_spline[f] for f in feature_names])
sig_spline = np.sum(p_values_spline < 0.05)
tp_spline = np.sum(p_values_spline[:n_signal_features] < 0.05)

print(f"  Completed: {np.sum(~np.isnan(p_values_spline))}/{n_features} features")
print(f"  Significant features (P < 0.05): {sig_spline}")
print(f"  True positives detected: {tp_spline}/{n_signal_features}")
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
fp_linear = np.sum(p_values_linear[n_signal_features:] < alpha)
fp_spline = np.sum(p_values_spline[n_signal_features:] < alpha)

sensitivity_linear = tp_linear / n_signal_features
sensitivity_spline = tp_spline / n_signal_features

fdr_linear = fp_linear / max(1, sig_linear)
fdr_spline = fp_spline / max(1, sig_spline)

power_improvement = (sensitivity_spline - sensitivity_linear) / max(0.01, sensitivity_linear) * 100

# Print comparison table
print(f"{'Metric':<30} {'Linear Model':<20} {'Spline Model':<20}")
print("-" * 70)
print(f"{'True Positives':<30} {tp_linear:<20} {tp_spline:<20}")
print(f"{'False Positives':<30} {fp_linear:<20} {fp_spline:<20}")
print(f"{'Sensitivity (Power)':<30} {sensitivity_linear*100:.1f}%{'':<14} {sensitivity_spline*100:.1f}%")
print(f"{'False Discovery Rate':<30} {fdr_linear*100:.1f}%{'':<14} {fdr_spline*100:.1f}%")
print(f"{'Total Significant':<30} {sig_linear:<20} {sig_spline:<20}")
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
ax.set_title('True Signal Pattern: Inverted-U Trajectory', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pilot_experiments/results/01a_true_signal_pattern.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: P-value comparison
results_df = pd.DataFrame({
    'Feature': feature_names,
    'TrueSignal': [i < n_signal_features for i in range(n_features)],
    'P_Linear': p_values_linear,
    'P_Spline': p_values_spline
})

fig, ax = plt.subplots(figsize=(10, 8))
for is_signal in [False, True]:
    subset = results_df[results_df['TrueSignal'] == is_signal]
    color = '#e74c3c' if is_signal else '#95a5a6'
    label = 'True Signal' if is_signal else 'Null'
    ax.scatter(np.maximum(subset['P_Linear'], 1e-10),
               np.maximum(subset['P_Spline'], 1e-10),
               c=color, label=label, s=50, alpha=0.7)

ax.plot([1e-10, 1], [1e-10, 1], 'k--', alpha=0.5, linewidth=1)
ax.axhline(0.05, color='red', linestyle=':', alpha=0.5)
ax.axvline(0.05, color='red', linestyle=':', alpha=0.5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('P-value (Linear Model)', fontsize=12)
ax.set_ylabel('P-value (Spline Model)', fontsize=12)
ax.set_title('P-value Comparison: Linear vs Spline Model', fontsize=14, fontweight='bold')
ax.legend(fontsize=10, loc='lower right')
ax.text(0.1, 1e-8, 'Spline Wins', color='#27ae60', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('pilot_experiments/results/01b_pvalue_comparison.pdf', dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: Power curves
thresholds = [0.001, 0.01, 0.05, 0.1, 0.2]
power_data = []

for th in thresholds:
    tp_lin = np.sum(p_values_linear[:n_signal_features] < th)
    tp_spl = np.sum(p_values_spline[:n_signal_features] < th)

    fp_lin = np.sum(p_values_linear[n_signal_features:] < th)
    fp_spl = np.sum(p_values_spline[n_signal_features:] < th)

    power_data.append({
        'Threshold': th,
        'Method': 'Linear Model',
        'Power': tp_lin / n_signal_features,
        'FDR': fp_lin / max(1, tp_lin + fp_lin)
    })
    power_data.append({
        'Threshold': th,
        'Method': 'Spline Model',
        'Power': tp_spl / n_signal_features,
        'FDR': fp_spl / max(1, tp_spl + fp_spl)
    })

power_df = pd.DataFrame(power_data)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

for method in ['Linear Model', 'Spline Model']:
    subset = power_df[power_df['Method'] == method]
    color = '#3498db' if method == 'Linear Model' else '#e74c3c'
    ax1.plot(subset['Threshold'], subset['Power'], 'o-',
             color=color, label=method, linewidth=2, markersize=8)
    ax2.plot(subset['Threshold'], subset['FDR'], 'o-',
             color=color, label=method, linewidth=2, markersize=8)

ax1.set_xscale('log')
ax1.set_xlabel('Significance Threshold (α)', fontsize=12)
ax1.set_ylabel('Power (True Positive Rate)', fontsize=12)
ax1.set_title('Statistical Power Comparison', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0, 1])

ax2.set_xscale('log')
ax2.axhline(0.05, color='red', linestyle='--', alpha=0.5)
ax2.set_xlabel('Significance Threshold (α)', fontsize=12)
ax2.set_ylabel('False Discovery Rate', fontsize=12)
ax2.set_title('False Discovery Rate Comparison', fontsize=14, fontweight='bold')
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
        'linear': {
            'true_positives': int(tp_linear),
            'false_positives': int(fp_linear),
            'sensitivity': float(sensitivity_linear),
            'fdr': float(fdr_linear)
        },
        'spline': {
            'true_positives': int(tp_spline),
            'false_positives': int(fp_spline),
            'sensitivity': float(sensitivity_spline),
            'fdr': float(fdr_spline)
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
print(f"  - Spline Model detected {tp_spline} true signals")
print(f"  - Linear Model detected only {tp_linear} true signals")
print(f"  - Power improvement: {power_improvement:.1f}%")
print(f"  - FDR (Spline): {fdr_spline*100:.1f}%")
print()

if sensitivity_spline > sensitivity_linear * 1.5 and fdr_spline < 0.2:
    print("✓✓✓ VALIDATION SUCCESSFUL ✓✓✓")
    print("Spline-based trajectory modeling shows substantial improvement!")
    print("This approach is ready for full implementation in LinDA v2.")
else:
    print("⚠ WARNING: Results are inconclusive.")
    print("Consider adjusting simulation parameters or signal patterns.")

print()
print("Generated files:")
print("  - pilot_experiments/results/01a_true_signal_pattern.pdf")
print("  - pilot_experiments/results/01b_pvalue_comparison.pdf")
print("  - pilot_experiments/results/01c_power_fdr_curves.pdf")
print("  - pilot_experiments/results/01_detailed_results.csv")
print("  - pilot_experiments/results/01_summary.json")
print()
