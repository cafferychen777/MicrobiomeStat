# Quick Start Guide - LinDA v2 Pilot Experiments

## 🎯 Goal

Run **2 validation experiments** (~10 minutes) to determine if LinDA v2's proposed innovations are worth implementing.

## 📋 Prerequisites

1. **R version ≥ 4.0** installed
2. Internet connection (for package installation)
3. ~500 MB disk space

## 🚀 Three-Step Setup

### Step 1: Check Dependencies

```bash
cd /home/user/MicrobiomeStat
bash pilot_experiments/check_dependencies.sh
```

**If you see errors**, proceed to Step 2. Otherwise, skip to Step 3.

### Step 2: Install Missing Packages

```bash
Rscript pilot_experiments/install_dependencies.R
```

This will automatically install all required R packages (~5 minutes).

### Step 3: Run Experiments

```bash
Rscript pilot_experiments/run_all_experiments.R
```

**Runtime:** ~5-10 minutes

## 📊 What You Get

After completion, check `pilot_experiments/results/`:

### Key Files

1. **`VALIDATION_REPORT.txt`** ← **READ THIS FIRST**
   - Executive summary
   - Performance metrics
   - Go/no-go recommendation

2. **PDF Plots** (10 files)
   - Signal patterns
   - P-value comparisons
   - Power curves
   - Phylogenetic trees

3. **CSV Data** (2 files)
   - Full statistical results
   - Can be imported for further analysis

## 🎓 Understanding the Results

### Experiment 1: Trajectory Splines

**Question:** Can we detect non-linear temporal patterns?

**What to Look For:**
- **Power improvement > 50%** → Strong validation ✓
- **FDR < 10%** → Good error control ✓

**Example Good Result:**
```
Power Improvement: 150%
Linear LMM:  20% sensitivity
Spline LMM:  50% sensitivity
```

### Experiment 2: Phylogenetic Smoothing

**Question:** Can we recover weak signals in related taxa?

**What to Look For:**
- **Sensitivity improvement > 100%** → Strong validation ✓
- **FDR < 20%** → Acceptable error control ✓

**Example Good Result:**
```
Sensitivity Improvement: 200%
No Smoothing:    10% sensitivity
With Smoothing:  30% sensitivity
```

## 🔍 Troubleshooting

### Problem: `Rscript: command not found`

**Solution:**
```bash
# Find R installation
which R

# If found, use full path
/usr/bin/Rscript pilot_experiments/run_all_experiments.R

# If not found, install R
# Ubuntu/Debian:
sudo apt-get install r-base r-base-dev

# macOS:
brew install r
```

### Problem: Package installation fails

**Common Issues:**

1. **ggtree installation fails**
   ```r
   # In R console:
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install("ggtree")
   ```

2. **Compilation errors**
   ```bash
   # Ubuntu/Debian: Install development tools
   sudo apt-get install build-essential gfortran

   # macOS: Install Xcode command line tools
   xcode-select --install
   ```

3. **Permission errors**
   ```r
   # In R console: Use personal library
   install.packages("package_name", lib = "~/R/packages")
   ```

### Problem: Experiments run but produce poor results

**Possible Reasons:**

1. **Random seed variation**
   - Rerun with different seeds
   - Results should be qualitatively similar

2. **Parameter mismatch**
   - Check if signal strength is too weak/strong
   - Adjust in the R scripts

## 📖 What's Next?

### If Both Experiments Validate ✓✓✓

**You should:**
1. **Celebrate!** 🎉 Your innovations are solid.
2. **Proceed with implementation** immediately.
3. **Use simulation results as Figure 1** in your manuscript.
4. **Plan comprehensive benchmarking** on real datasets.

**Timeline:**
- Full implementation: 2-3 weeks
- Testing: 1-2 weeks
- Benchmarking: 2-3 weeks
- Manuscript: 3-4 weeks
- **Total: 2-3 months to submission**

### If One Experiment Validates ✓

**You should:**
1. **Prioritize the validated innovation** for implementation.
2. **Refine the other approach:**
   - Adjust simulation parameters
   - Try alternative formulations
   - Consult statistical literature
3. **Consider a phased rollout:**
   - LinDA v2.0: Trajectory splines only
   - LinDA v2.1: Add phylogenetic smoothing later

### If Neither Experiment Validates ⚠

**Don't panic!** This is why we run pilots. Options:

1. **Adjust simulation parameters:**
   - Increase effect sizes
   - Change signal patterns
   - Modify noise levels

2. **Revisit methodology:**
   - Review recent literature
   - Consult with statisticians
   - Consider alternative approaches

3. **Focus on incremental improvements:**
   - Optimize existing LinDA features
   - Better diagnostics/visualization
   - Integration with other tools

## 💡 Pro Tips

### Running Individual Experiments

If you only want to test one innovation:

```bash
# Trajectory splines only
Rscript pilot_experiments/01_trajectory_spline_validation.R

# Phylogenetic smoothing only
Rscript pilot_experiments/02_phylo_smoothing_validation.R
```

### Customizing Parameters

Edit the R scripts directly to adjust:

**Experiment 1 (Lines 20-28):**
```r
n_subjects_per_group <- 50   # More subjects = more power
signal_strength <- 2.0        # Stronger signal = easier detection
n_timepoints <- 5             # More timepoints = better trajectory
```

**Experiment 2 (Lines 20-27):**
```r
n_subjects_per_group <- 100  # More subjects = detect weaker signals
signal_effect_size <- 0.3    # Smaller = more realistic weak signal
lambda <- 0.5                # Smoothing strength (0-2)
```

### Viewing Results

**On a remote server?** Download the PDFs:

```bash
# Compress results
tar -czf linda2_pilot_results.tar.gz pilot_experiments/results/

# Download via scp
scp user@server:/path/to/linda2_pilot_results.tar.gz .

# Extract
tar -xzf linda2_pilot_results.tar.gz
```

**Want interactive plots?** Modify scripts to use `plotly` instead of `ggplot2`.

## 📚 Further Reading

- **LinDA Paper:** Zhou et al. (2022), Genome Biology
- **Splines:** Wood (2017), Generalized Additive Models
- **Graph Signal Processing:** Shuman et al. (2013), IEEE Signal Processing

## 🆘 Getting Help

**Questions?** Contact:
- Chen Yang: cafferychen777@tamu.edu
- GitHub Issues: https://github.com/cafferychen777/MicrobiomeStat/issues

**Found a bug?** Please report with:
- R version (`R.version.string`)
- Package versions (`sessionInfo()`)
- Error messages (full output)

---

**Good luck with your validation experiments!** 🚀
