# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MicrobiomeStat is a comprehensive R package for advanced longitudinal microbiome and multi-omics data analysis. It specializes in in-depth longitudinal microbiome analysis across time and provides tools for cross-sectional studies, paired sample analysis, and automated report generation.

## Development Commands

### Build and Check
```bash
# Build package
R CMD build .

# Check package (basic)
R CMD check MicrobiomeStat_*.tar.gz

# Check without vignettes and examples (faster)
R CMD check --no-build-vignettes --no-manual --no-examples MicrobiomeStat_*.tar.gz
```

### Testing
```bash
# Run tests using testthat
Rscript -e 'devtools::test()'

# Run a specific test file
Rscript -e 'testthat::test_file("tests/testthat/test-linda.R")'
```

### Documentation
```bash
# Generate documentation
Rscript -e 'devtools::document()'

# Build pkgdown site
Rscript -e 'pkgdown::build_site()'
```

### Installation
```bash
# Install from source (development)
Rscript -e 'devtools::install()'

# Install dependencies
Rscript -e 'remotes::install_deps(dependencies = c("Depends", "Imports", "LinkingTo", "Suggests"))'
```

## Code Architecture

### Core Data Structure
The package uses a `data.obj` structure - a list containing:
- **feature.tab**: Matrix of feature abundances (rows: features, columns: samples)
- **feature.ann**: Matrix/data.frame of feature annotations (taxonomy)
- **meta.dat**: Data.frame of sample metadata
- **tree**: Phylogenetic tree (optional)
- **feature.agg.list**: List of aggregated features at different taxonomic levels

### Key Function Categories

1. **Data Import/Conversion** (`mStat_convert_*`, `mStat_import_*`):
   - Convert from various formats (phyloseq, QIIME2, DADA2, etc.)
   - Located in R/mStat_convert_*.R and R/mStat_import_*.R

2. **Data Processing** (`mStat_*`):
   - Normalization, filtering, rarefaction, aggregation
   - Core utilities in R/mStat_*.R files

3. **Diversity Analysis**:
   - Alpha diversity: `generate_alpha_*` functions
   - Beta diversity: `generate_beta_*` functions
   - Calculate diversity metrics in R/mStat_calculate_*.R

4. **Statistical Testing**:
   - LinDA (Linear models for Differential Abundance): R/linda.R
   - Various test functions: `generate_*_test_*` pattern

5. **Visualization**:
   - Extensive plotting functions: `generate_*_boxplot`, `generate_*_heatmap`, etc.
   - Report generation: `mStat_generate_report_*` for different study designs

6. **Study Design Support**:
   - Cross-sectional: `*_single` suffix functions
   - Paired samples: `*_pair` suffix functions  
   - Longitudinal: `*_long` suffix functions

### Important Validation
All data objects should be validated using `mStat_validate_data()` which ensures:
- Proper list structure
- meta.dat is a data.frame
- Row/column name alignment between components

### Key Dependencies
The package heavily relies on:
- ggplot2 for visualization
- lmerTest for mixed models
- vegan for ecological analyses
- dplyr/tidyr for data manipulation
- ggtree/ggtreeExtra for phylogenetic visualization

## Testing Approach
Tests are located in tests/testthat/ and use the testthat framework. Test files follow the pattern test-*.R and include unit tests for core functions, especially the LinDA statistical method.