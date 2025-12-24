#!/bin/bash
#' ==============================================================================
#' Check Dependencies for LinDA v2 Pilot Experiments
#' ==============================================================================

echo "════════════════════════════════════════════════════════════════════"
echo "  LinDA v2 Pilot Experiments - Dependency Checker"
echo "════════════════════════════════════════════════════════════════════"
echo ""

# Check R installation
echo "Checking R installation..."
if command -v R &> /dev/null; then
    R_VERSION=$(R --version | head -1)
    echo "  ✓ R is installed: $R_VERSION"
else
    echo "  ✗ R is NOT installed"
    echo ""
    echo "Please install R:"
    echo "  - Ubuntu/Debian: sudo apt-get install r-base r-base-dev"
    echo "  - MacOS: brew install r"
    echo "  - Windows: Download from https://cran.r-project.org/"
    echo ""
    exit 1
fi

echo ""

# Check required R packages
echo "Checking R packages..."
R --quiet --vanilla <<EOF
packages_needed <- c(
  "lme4", "lmerTest", "splines",
  "ape", "phangorn", "ggtree",
  "dplyr", "tidyr", "Matrix",
  "ggplot2", "patchwork", "modeest"
)

installed <- sapply(packages_needed, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

cat("\n")
for (i in seq_along(packages_needed)) {
  pkg <- packages_needed[i]
  if (installed[i]) {
    cat(sprintf("  ✓ %s is installed\n", pkg))
  } else {
    cat(sprintf("  ✗ %s is MISSING\n", pkg))
  }
}

missing <- packages_needed[!installed]
if (length(missing) > 0) {
  cat("\nTo install missing packages, run:\n")
  cat(sprintf('  install.packages(c("%s"))\n', paste(missing, collapse = '", "')))
  quit(status = 1)
} else {
  cat("\n✓ All required packages are installed!\n")
  quit(status = 0)
}
EOF

if [ $? -eq 0 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════════════════"
    echo "  ✓✓✓ All dependencies satisfied!"
    echo "════════════════════════════════════════════════════════════════════"
    echo ""
    echo "You can now run the experiments:"
    echo "  Rscript pilot_experiments/run_all_experiments.R"
    echo ""
else
    echo ""
    echo "════════════════════════════════════════════════════════════════════"
    echo "  ⚠ Please install missing packages first"
    echo "════════════════════════════════════════════════════════════════════"
    echo ""
    exit 1
fi
