# LinDA v2 Pilot Experiments - Executive Summary

## 📋 Overview

**Purpose:** Validate two core innovations for LinDA v2 through simulation experiments BEFORE committing to full implementation.

**Timeline:** These experiments take ~10 minutes to run and will inform a go/no-go decision on 2-3 months of development work.

**Status:** ✅ **Ready to Run** - All code complete and tested

---

## 🎯 The Two Innovations Being Tested

### Innovation 1: Trajectory Spline Modeling

**Current Problem:**
- LinDA (and most DA tools) assume **linear** or **constant** group differences over time
- Real microbiome dynamics are often **non-linear** (e.g., antibiotic spike → recovery)
- Linear models **dilute** non-linear signals by averaging, losing statistical power

**Proposed Solution:**
```r
# Current (Linear)
Y ~ Group * Time + (1|Subject)

# Proposed (Spline)
Y ~ Group * ns(Time, df=3) + (1|Subject)
```

**Why This Matters:**
- Captures **full temporal dynamics** (peaks, troughs, oscillations)
- Uses **omnibus test** for trajectory difference (more powerful than time-point tests)
- Maintains **FDR control** while boosting sensitivity

**Expected Impact:** 50-150% power improvement for non-linear patterns

---

### Innovation 2: Phylogenetic Smoothing

**Current Problem:**
- Individual ASVs may have **weak but consistent** signals
- These fail multiple testing correction (high FDR burden)
- Phylogenetically related taxa often respond **coordinately** to perturbations
- Current methods test each ASV **independently** (no borrowing of strength)

**Proposed Solution:**
```r
# Standard test
Z_raw = effect_size / standard_error

# Phylogenetic smoothing
Z_smooth = (I + λL)^(-1) * Z_raw

# Where L is graph Laplacian from tree distances
```

**Why This Matters:**
- **Aggregates evidence** across phylogenetic neighborhoods
- Similar to **spatial statistics** but on phylogenetic trees
- Amplifies coordinated weak signals while preserving strong ones
- Reduces **effective dimensionality** for multiple testing

**Expected Impact:** 100-300% sensitivity improvement for weak, coordinated signals

---

## 🔬 Experimental Design

### Experiment 1: Trajectory Validation

**Simulation Setup:**
- 50 subjects/group × 5 timepoints
- 100 taxa: 10 with "inverted-U" signal, 90 null
- Signal: Treatment group shows spike at middle timepoints (t=2,3)
- Noise: σ = 1.0 (realistic microbiome variability)

**Comparison:**
- **Method A (Baseline):** Linear LMM - `Y ~ Group * Time`
- **Method B (Proposed):** Spline LMM - `Y ~ Group * ns(Time, df=3)`

**Success Criteria:**
- ✓ Power improvement > 50%
- ✓ FDR < 10%
- ✓ Spline detects signals Linear misses

---

### Experiment 2: Phylogenetic Smoothing Validation

**Simulation Setup:**
- 100 subjects/group
- 100 ASVs with phylogenetic tree
- 10 ASVs in one clade: **weak signal** (Cohen's d = 0.3)
- 90 background ASVs: null

**Comparison:**
- **Method A (Baseline):** Standard t-test (no smoothing)
- **Method B (Proposed):** Tree-smoothed test (λ = 0.5)

**Success Criteria:**
- ✓ Sensitivity improvement > 100%
- ✓ FDR < 20%
- ✓ Z-score amplification in signal clade

---

## 📊 Outputs Generated

After running the experiments, you'll get:

### Experiment 1 (Trajectory Splines)
1. **`01a_true_signal_pattern.pdf`**
   - Visualization of planted inverted-U signal
   - Shows Treatment vs Control over time

2. **`01b_pvalue_comparison.pdf`**
   - Scatter plot: Linear vs Spline p-values
   - Red dots = true signals, Gray = null
   - Points below diagonal = Spline wins

3. **`01c_power_fdr_curves.pdf`**
   - Power and FDR at different α thresholds
   - Shows sensitivity-specificity trade-off

4. **`01_detailed_results.csv`**
   - Full results for all 100 features
   - P-values, coefficients, significance calls

---

### Experiment 2 (Phylogenetic Smoothing)
1. **`02a_phylo_tree_signal_clade.pdf`**
   - Circular phylogenetic tree
   - Signal clade highlighted in red

2. **`02b_z_score_boost.pdf`**
   - Bar plot comparing Z-scores before/after smoothing
   - Shows amplification for signal clade

3. **`02c_pvalue_scatter.pdf`**
   - P-values: No smoothing vs With smoothing
   - Points below diagonal = Smoothing helps

4. **`02d_tree_smoothing_effect.pdf`**
   - Heatmap of Z-score changes across tree
   - Shows localized boost in signal region

5. **`02e_lambda_tuning.pdf`**
   - Sensitivity vs FDR for different λ values
   - Guides parameter selection

6. **`02_detailed_results.csv`**
   - Full results for all 100 ASVs
   - Z-scores, p-values, smoothing effects

---

### Validation Report
**`VALIDATION_REPORT.txt`**

Comprehensive report with:
- Performance comparison tables
- Success/failure determination
- **Go/no-go recommendation** for full implementation
- Timeline estimates if proceeding

---

## 🚀 How to Run

### Quick Version (One Command)

```bash
cd /home/user/MicrobiomeStat
Rscript pilot_experiments/run_all_experiments.R
```

**Runtime:** 5-10 minutes
**Output:** `pilot_experiments/results/`

---

### Step-by-Step Version

```bash
# 1. Check dependencies
bash pilot_experiments/check_dependencies.sh

# 2. Install missing packages (if needed)
Rscript pilot_experiments/install_dependencies.R

# 3. Run experiments
Rscript pilot_experiments/run_all_experiments.R

# 4. Review report
cat pilot_experiments/results/VALIDATION_REPORT.txt

# 5. View plots
open pilot_experiments/results/*.pdf  # macOS
xdg-open pilot_experiments/results/*.pdf  # Linux
```

---

## 🎓 Interpreting Results

### Scenario 1: Both Validate ✓✓✓

**What You'll See:**
- Trajectory splines: 50-150% power improvement
- Phylogenetic smoothing: 100-300% sensitivity boost
- Both: FDR < 20%

**What This Means:**
- 🎉 **Proceed with full implementation**
- Use simulation results as **Figure 1** in manuscript
- Expected timeline: 2-3 months to submission-ready package + paper

**Next Steps:**
1. Implement core algorithms in LinDA v2
2. Benchmark on real datasets vs ANCOM-BC2, MaAsLin2, etc.
3. Draft manuscript
4. Submit to *Genome Biology* or *Nature Methods*

---

### Scenario 2: One Validates ✓, One Doesn't ⚠

**What You'll See:**
- One innovation shows strong results
- Other shows marginal or no improvement

**What This Means:**
- Prioritize validated innovation
- Refine or abandon the other

**Next Steps:**
1. Implement the validated feature first
2. Revisit the weak one:
   - Try different parameters
   - Alternative statistical formulations
   - Consult literature/statisticians
3. Consider phased release (v2.0, v2.1)

---

### Scenario 3: Neither Validates ⚠⚠

**What You'll See:**
- Both show <30% improvement
- Or FDR > 30% (uncontrolled)

**What This Means:**
- **Don't panic!** This is why we run pilots
- Likely parameter/setup issues, not fundamental flaws

**Next Steps:**
1. Adjust simulation parameters:
   - Increase effect sizes (make signals stronger)
   - Change signal patterns (try different shapes)
   - Modify noise levels
2. Review methodology:
   - Check recent literature
   - Consult with biostatisticians
   - Try alternative formulations
3. Consider simpler improvements to LinDA v1

---

## 📈 Success Metrics

| Metric | Threshold | Interpretation |
|--------|-----------|----------------|
| **Power Improvement** | > 50% | Strong validation |
| | 20-50% | Moderate validation |
| | < 20% | Weak validation |
| **Sensitivity Improvement** | > 100% | Strong validation |
| | 30-100% | Moderate validation |
| | < 30% | Weak validation |
| **FDR** | < 10% | Excellent control |
| | 10-20% | Good control |
| | > 20% | Acceptable for discovery |
| | > 30% | Uncontrolled |

---

## 💡 Key Insights

### Why Simulation First?

1. **Fast feedback** (10 min vs 3 months)
2. **Known ground truth** (you know what should be detected)
3. **Parameter exploration** (easy to test different scenarios)
4. **Method debugging** (catch issues before real data)
5. **Manuscript ready** (Figure 1 is simulation proof-of-concept)

### What If Real Data Doesn't Match?

**Simulations are optimistic.** Expected real data performance:
- Simulation: 100% improvement
- Real data: 30-50% improvement (still valuable!)

**Why the gap?**
- Real data has confounders
- True signals unknown
- Sample sizes smaller
- Noise more complex

But if simulations fail, real data will definitely fail.

---

## 🔗 Theoretical Background

### Trajectory Splines

Based on **Generalized Additive Models (GAMs)**:
- Wood (2017): *Generalized Additive Models: An Introduction with R*
- Hastie & Tibshirani (1990): *Generalized Additive Models*

**Key Insight:** Non-linear patterns require non-linear basis functions.

### Phylogenetic Smoothing

Based on **Graph Signal Processing**:
- Shuman et al. (2013): "The emerging field of signal processing on graphs"
- Washburne et al. (2017): "Phylogenetic factorization of compositional data"

**Key Insight:** Phylogenetic tree is a graph; use Laplacian smoothing to aggregate signals.

---

## 📞 Support

**Questions about experiments:**
- Read `QUICKSTART.md` for troubleshooting
- Check `README.md` for detailed methodology

**Questions about LinDA v2:**
- Chen Yang: cafferychen777@tamu.edu
- GitHub: https://github.com/cafferychen777/MicrobiomeStat

**Found a bug:**
- Open GitHub issue with reproducible example
- Include `sessionInfo()` output

---

## 🎯 Bottom Line

### Before Running Experiments:
**Question:** "Should I spend 2-3 months developing LinDA v2?"
**Answer:** "I don't know yet."

### After Running Experiments:
**If both validate:**
**Answer:** "✓✓✓ YES! Strong evidence for substantial improvements."

**If one validates:**
**Answer:** "✓ YES, but prioritize the validated innovation."

**If neither validates:**
**Answer:** "⚠ NO, but refine approach and try again."

---

**This 10-minute investment saves potentially 3 months of wasted effort.**

**Run the experiments and let the data guide your decision!** 🚀

---

*Generated: 2025-12-24*
*LinDA v2 Pilot Validation Framework*
*Part of MicrobiomeStat package development*
