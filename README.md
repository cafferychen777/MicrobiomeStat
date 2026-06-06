# MicrobiomeStat <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![R-CMD-check](https://github.com/cafferychen777/MicrobiomeStat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cafferychen777/MicrobiomeStat/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

MicrobiomeStat is an R package for statistical analysis and visualization of microbiome and multi-omics data. It supports cross-sectional, paired-sample, and longitudinal study designs through a unified interface, with built-in automated report generation.

## Key Features

- **Unified interface** -- `test_alpha()`, `test_beta()`, `test_taxa()` and their `plot_*()` counterparts automatically dispatch to the correct method based on study design.
- **LinDA** -- Linear models for differential abundance analysis of compositional data ([Zhou et al. 2022](https://doi.org/10.1186/s13059-022-02655-5)), with an enhanced `linda2()` variant supporting detection-depth weighting and tree-guided smoothing.
- **Longitudinal analysis** -- Trend tests, volatility assessment, spaghetti plots, and per-time-point comparisons for time-series microbiome data.
- **Automated reports** -- `mStat_generate_report_single()`, `_pair()`, and `_long()` produce publication-ready PDF/HTML reports from a single function call.
- **Broad data import** -- Native support for phyloseq, QIIME2, DADA2, mothur, BIOM, DESeq2, edgeR, SummarizedExperiment, and MultiAssayExperiment objects.
- **120+ functions** covering alpha diversity, beta diversity, taxonomic composition, and their change/trend/volatility variants across all three study designs.

## Installation

Install from CRAN:

```r
install.packages("MicrobiomeStat")
```

Or install the development version from GitHub for the full feature set:

```r
# install.packages("devtools")
devtools::install_github("cafferychen777/MicrobiomeStat")
```

## Quick Start

```r
library(MicrobiomeStat)

# Load example data
data(peerj32.obj)

# Alpha diversity
alpha.obj <- mStat_calculate_alpha_diversity(
  peerj32.obj$feature.tab,
  alpha.name = c("shannon", "observed_species")
)

generate_alpha_boxplot_single(
  data.obj = peerj32.obj,
  alpha.obj = alpha.obj,
  alpha.name = "shannon",
  subject.var = "subject",
  time.var = "time",
  t.level = "2",
  group.var = "group",
  strata.var = "sex"
)

# Beta diversity ordination
generate_beta_ordination_single(
  data.obj = peerj32.obj,
  dist.obj = NULL,
  pc.obj = NULL,
  subject.var = "subject",
  time.var = "time",
  t.level = "2",
  group.var = "group",
  strata.var = "sex",
  dist.name = "BC"
)

# Differential abundance with LinDA
generate_taxa_test_single(
  data.obj = peerj32.obj,
  subject.var = "subject",
  time.var = "time",
  t.level = "2",
  group.var = "group",
  adj.vars = "sex",
  feature.level = "Genus",
  feature.dat.type = "count"
)
```

## Data Import

```r
# From phyloseq
data.obj <- mStat_convert_phyloseq_to_data_obj(phyloseq_obj)

# From QIIME2
data.obj <- mStat_import_qiime2_as_data_obj(
  otu.qza = "table.qza",
  taxa.qza = "taxonomy.qza",
  sam.tab = "metadata.txt"
)

# From DADA2
data.obj <- mStat_import_dada2_as_data_obj(
  seq.tab = seqtab, tax.tab = taxtab, sam.tab = metadata
)
```

See `?mStat_import_biom_as_data_obj`, `?mStat_import_mothur_as_data_obj`, and the Bioconductor converters for other formats.

## Documentation

- **Tutorials**: [www.microbiomestat.wiki](https://www.microbiomestat.wiki)
- **Function reference**: [cafferychen777.github.io/MicrobiomeStat](https://cafferychen777.github.io/MicrobiomeStat/)
- **Vignettes**: `vignette("MicrobiomeStat")` and `vignette("longitudinal-analysis")`

## Citation

If you use MicrobiomeStat in your research, please cite:

> Zhou, H., He, K., Chen, J., & Zhang, X. (2022). LinDA: linear models for differential abundance analysis of microbiome compositional data. *Genome Biology*, 23(1), 95. <https://doi.org/10.1186/s13059-022-02655-5>

```r
citation("MicrobiomeStat")
```

## Community

- [GitHub Issues](https://github.com/cafferychen777/MicrobiomeStat/issues) -- bug reports and feature requests
- [Discord](https://discord.gg/U6vNamyk) -- questions and discussion
