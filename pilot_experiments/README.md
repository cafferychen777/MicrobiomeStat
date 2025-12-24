# LinDA v2 Pilot Experiments

## Overview

This directory contains **validation experiments** for two core innovations planned for LinDA v2:

1. **Trajectory Spline Modeling** - Using natural splines to detect non-linear temporal patterns
2. **Phylogenetic Smoothing** - Leveraging phylogenetic relationships to boost weak signals

These experiments are designed as **quick-and-dirty proof-of-concept** studies to validate the methodology *before* investing months in full implementation.

## Motivation

Before committing to:
- Writing a complete R package
- Conducting extensive real-data benchmarking
- Writing a manuscript

We need to confirm:
- ✓ The proposed methods detect signals that current methods miss
- ✓ Statistical power improvements are substantial (>50%)
- ✓ False discovery rate remains well-controlled (<10%)

## Experiments

### Experiment 1: Trajectory Spline Validation

**File:** `01_trajectory_spline_validation.R`

**Question:** Can spline-based trajectory modeling detect non-linear temporal patterns that linear models miss?

**Design:**
- Simulate microbiome data with **inverted-U signal** (spike at middle timepoints)
- Compare **Linear LMM** (`Y ~ Group * Time`) vs **Spline LMM** (`Y ~ Group * ns(Time, df=3)`)
- Measure statistical power and FDR

**Expected Outcome:**
- Spline LMM should detect 5-10x more true signals
- Linear LMM should fail on non-linear patterns

### Experiment 2: Phylogenetic Smoothing Validation

**File:** `02_phylo_smoothing_validation.R`

**Question:** Can phylogenetic tree smoothing recover weak signals when related ASVs show coordinated changes?

**Design:**
- Simulate microbiome data with phylogenetic tree
- Plant **weak signals** (Cohen's d = 0.3) in a phylogenetically coherent clade
- Compare **Standard t-test** vs **Tree-smoothed test** using graph Laplacian
- Measure sensitivity improvement

**Expected Outcome:**
- Smoothing should amplify weak signals by "borrowing strength" from neighbors
- Sensitivity should improve 2-5x for weak, coordinated signals

## Usage

### Quick Start (Run All Experiments)

```bash
cd /home/user/MicrobiomeStat
Rscript pilot_experiments/run_all_experiments.R
```

This will:
1. Run both validation experiments
2. Generate plots and summary statistics
3. Create a comprehensive validation report

**Runtime:** ~5-10 minutes on a modern laptop

### Run Individual Experiments

```bash
# Experiment 1 only
Rscript pilot_experiments/01_trajectory_spline_validation.R

# Experiment 2 only
Rscript pilot_experiments/02_phylo_smoothing_validation.R
```

## Outputs

All results are saved to `pilot_experiments/results/`:

### Experiment 1 Outputs
- `01a_true_signal_pattern.pdf` - Visualization of the planted signal
- `01b_pvalue_comparison.pdf` - Scatter plot of p-values (Linear vs Spline)
- `01c_power_fdr_curves.pdf` - Power and FDR at different thresholds
- `01_detailed_results.csv` - Full results for all features
- `01_summary.rds` - Summary statistics

### Experiment 2 Outputs
- `02a_phylo_tree_signal_clade.pdf` - Tree with signal clade highlighted
- `02b_z_score_boost.pdf` - Z-score amplification in signal clade
- `02c_pvalue_scatter.pdf` - P-values before/after smoothing
- `02d_tree_smoothing_effect.pdf` - Smoothing effect across entire tree
- `02e_lambda_tuning.pdf` - Sensitivity/FDR trade-off for different λ values
- `02_detailed_results.csv` - Full results for all ASVs
- `02_summary.rds` - Summary statistics

### Validation Report
- `VALIDATION_REPORT.txt` - Comprehensive report with recommendations

## Dependencies

Required R packages:
```r
# Core statistics
install.packages(c("lme4", "lmerTest", "splines"))

# Phylogenetics
install.packages(c("ape", "phangorn", "ggtree"))

# Data manipulation
install.packages(c("dplyr", "tidyr", "Matrix"))

# Visualization
install.packages(c("ggplot2", "patchwork"))
```

## Success Criteria

### Experiment 1 (Trajectory Splines)
- ✓ Power improvement > 50%
- ✓ FDR < 10%
- ✓ Spline LMM detects signals with P-values 1-2 orders of magnitude smaller

### Experiment 2 (Phylogenetic Smoothing)
- ✓ Sensitivity improvement > 100% (2x)
- ✓ FDR < 20%
- ✓ Z-scores for signal clade amplified by smoothing

## Decision Framework

| Scenario | Action |
|----------|--------|
| **Both experiments validate strongly** | ✓✓✓ Proceed with full LinDA v2 implementation |
| **One experiment validates strongly** | ✓ Implement validated feature, refine the other |
| **Both experiments show weak results** | ⚠ Refine approach, adjust parameters, or pivot |

## Timeline Estimates

If validation is successful:

| Phase | Duration |
|-------|----------|
| Pilot experiments (this) | 1 day |
| Core algorithm implementation | 2-3 weeks |
| Testing and validation | 1-2 weeks |
| Real data benchmarking | 2-3 weeks |
| Manuscript writing | 3-4 weeks |
| **Total to submission** | **~2-3 months** |

## Key Insights from Simulations

### Why These Experiments Matter

1. **Trajectory Splines:**
   - Real microbiome data often shows **non-monotonic changes** (e.g., antibiotics cause spike then recovery)
   - Linear models **dilute** these effects by averaging across time
   - Splines capture the **full temporal dynamics**

2. **Phylogenetic Smoothing:**
   - Individual ASVs may have **weak but consistent signals**
   - Related ASVs tend to respond **coordinately** to perturbations
   - Smoothing **aggregates evidence** across phylogenetic neighborhoods
   - Addresses the **multiple testing burden** by reducing effective dimension

## Theoretical Foundation

### Spline Trajectory Model

```r
# Formula
Y_ij ~ Group_i * ns(Time_ij, df=3) + (1|Subject_i)

# Hypothesis test
H0: All Group:Spline interaction terms = 0
HA: At least one interaction term ≠ 0
```

Uses likelihood ratio test to assess **omnibus trajectory difference**.

### Phylogenetic Smoothing

```r
# Smoothing operator
Z_smooth = (I + λL)^(-1) * Z_raw

# Where:
# - L is the graph Laplacian from phylogenetic distances
# - λ controls smoothing strength
# - Higher λ = more borrowing from neighbors
```

Related to **graph signal processing** and **spatial statistics**.

## References

### Trajectory Modeling
- Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*. 2nd ed.
- Hastie, T., & Tibshirani, R. (1990). *Generalized Additive Models*.

### Phylogenetic Smoothing
- Washburne, A. D., et al. (2017). Phylogenetic factorization of compositional data. *Bioinformatics*.
- Shuman, D. I., et al. (2013). The emerging field of signal processing on graphs. *IEEE Signal Processing Magazine*.

### LinDA Original
- Zhou, H., et al. (2022). LinDA: Linear models for differential abundance analysis of microbiome compositional data. *Genome Biology*.

## Contact

For questions or issues:
- Chen Yang (cafferychen777@tamu.edu)
- GitHub: https://github.com/cafferychen777/MicrobiomeStat

## License

These pilot experiments are part of the MicrobiomeStat package (GPL-3.0).
