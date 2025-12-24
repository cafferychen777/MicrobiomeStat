#!/usr/bin/env Rscript
#' ==============================================================================
#' Install Dependencies for LinDA v2 Pilot Experiments
#' ==============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat("  Installing R Packages for LinDA v2 Pilot Experiments\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat("\n")

# List of required packages
packages_needed <- c(
  # Core statistics
  "lme4",         # Linear mixed-effects models
  "lmerTest",     # Tests for lmer models
  "splines",      # Spline basis functions (usually included with R)
  "modeest",      # Mode estimation (for LinDA bias correction)

  # Phylogenetics
  "ape",          # Phylogenetic tree manipulation
  "phangorn",     # Phylogenetic analysis
  "ggtree",       # Tree visualization

  # Data manipulation
  "dplyr",        # Data wrangling
  "tidyr",        # Data tidying
  "Matrix",       # Sparse matrices

  # Visualization
  "ggplot2",      # Grammar of graphics
  "patchwork"     # Combine plots
)

cat("Packages to install:\n")
cat(paste0("  - ", packages_needed, collapse = "\n"))
cat("\n\n")

# Check what's already installed
already_installed <- sapply(packages_needed, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

packages_to_install <- packages_needed[!already_installed]

if (length(packages_to_install) == 0) {
  cat("✓ All packages are already installed!\n\n")
} else {
  cat("Installing", length(packages_to_install), "packages...\n\n")

  # Set CRAN mirror
  options(repos = c(CRAN = "https://cloud.r-project.org/"))

  # Install packages
  for (pkg in packages_to_install) {
    cat("Installing", pkg, "...\n")

    tryCatch({
      # Special handling for ggtree (Bioconductor package)
      if (pkg == "ggtree") {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", quiet = TRUE)
        }
        BiocManager::install("ggtree", update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, quiet = TRUE)
      }
      cat("  ✓ Installed", pkg, "\n")
    }, error = function(e) {
      cat("  ✗ Failed to install", pkg, "\n")
      cat("    Error:", e$message, "\n")
    })
  }

  cat("\n")
}

# Verify installation
cat("Verifying installation...\n\n")

all_installed <- sapply(packages_needed, function(pkg) {
  is_installed <- requireNamespace(pkg, quietly = TRUE)

  if (is_installed) {
    cat(sprintf("  ✓ %-20s installed\n", pkg))
  } else {
    cat(sprintf("  ✗ %-20s MISSING\n", pkg))
  }

  return(is_installed)
})

cat("\n")

if (all(all_installed)) {
  cat("════════════════════════════════════════════════════════════════════\n")
  cat("  ✓✓✓ Installation Complete!\n")
  cat("════════════════════════════════════════════════════════════════════\n")
  cat("\n")
  cat("You can now run the pilot experiments:\n")
  cat("  Rscript pilot_experiments/run_all_experiments.R\n")
  cat("\n")
} else {
  cat("════════════════════════════════════════════════════════════════════\n")
  cat("  ⚠ Some packages failed to install\n")
  cat("════════════════════════════════════════════════════════════════════\n")
  cat("\n")
  cat("Please install them manually:\n")
  missing <- packages_needed[!all_installed]
  cat(sprintf('  install.packages(c("%s"))\n', paste(missing, collapse = '", "')))
  cat("\n")
}
