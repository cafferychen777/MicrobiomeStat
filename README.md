MicrobiomeStat: Your Best Choice for Microbiome Analysis Ever! 🧬🔬
================

<!-- badges: start -->

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MicrobiomeStat)](https://cran.r-project.org/package=MicrobiomeStat)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`MicrobiomeStat` is a powerful and versatile R package, designed
specifically for analyzing longitudinal microbiome data. It extends to
handle multi-omics data, making it an asset for broader biological
research.

Crucially, `MicrobiomeStat` shines not only in longitudinal studies, but
is also adept at handling cross-sectional studies. For more detailed
tutorials and information, please visit our main page:
[www.microbiomestat.wiki](http://www.microbiomestat.wiki).

## Why Choose MicrobiomeStat? 👍

Choosing the right tools for your microbiome data analysis can make a
significant difference in the results. We have compared `MicrobiomeStat`
with other similar packages, and we believe our package stands out in
many aspects. Below are two comparative tables, showing the distinct
advantages of `MicrobiomeStat` over other packages.

## Comparative Table 1: Longitudinal Packages

|                    Features                    |                                      MicrobiomeStat                                      |              q2-longitudinal               |                               SplinectomeR                               |                             coda4microbiome                              |
|:----------------------------------------------:|:----------------------------------------------------------------------------------------:|:------------------------------------------:|:------------------------------------------------------------------------:|:------------------------------------------------------------------------:|
|               Data Input Formats               | QIIME2, DADA2, BIOM, Mothur, Phyloseq, DGEList, DESeqDataSet, SummarizedExperiment, etc. |                   QIIME2                   | Custom format (No other formats import or conversion functions provided) | Custom format (No other formats import or conversion functions provided) |
|       Supports Alpha Diversity Analysis        |                                            ✔️                                            |  Metrics need to be calculated externally  |                 Metrics need to be calculated externally                 |                                    ❌                                    |
|        Supports Beta Diversity Analysis        |                                            ✔️                                            | Distances need to be calculated externally |                                    ❌                                    |                                    ❌                                    |
|             Paired Sample Analysis             |                                            ✔️                                            |                     ✔️                     |                                    ❌                                    |                                    ❌                                    |
|         Multiple Visualization Styles          |   Boxplot, Violinplot, Spaghettiplot, Heatmap, Stacked Barplot, Areaplot, Volcano Plot   |           Spaghettiplot, Boxplot           |                                 Lineplot                                 |                 Boxplot, Densityplot, Barplot, Lineplot                  |
|          Key Features Identification           |                                            ✔️                                            |                     ✔️                     |                                    ✔️                                    |                                    ✔️                                    |
|                Diff. Prevalence                |                                            ✔️                                            |                     ❌                     |                                    ❌                                    |                                    ❌                                    |
|    Paired Samples Support in Taxa Analysis     |                                            ✔️                                            |                     ✔️                     |                                    ❌                                    |                                    ❌                                    |
| Time Points Subgroup Analysis in Taxa Analysis |                                            ✔️                                            |                     ❌                     |                                    ❌                                    |                                    ❌                                    |
|           Similar Abundance Pattern            |                                            ✔️                                            |                     ❌                     |                                    ❌                                    |                                    ❌                                    |
|             Similar Change Pattern             |                                            ✔️                                            |                     ❌                     |                                    ❌                                    |                                    ❌                                    |
|    Similar Changes in Composition Over Time    |                                            ✔️                                            |                     ❌                     |                                    ❌                                    |                                    ❌                                    |
|               Publication Ready                |                                            ✔️                                            |                     ✔️                     |                                    ❌                                    |                                    ❌                                    |
|              Ongoing Development               |                                            ✔️                                            |                     ✔️                     |                                    ✔️                                    |                                    ✔️                                    |

As shown in Comparative Table 1, `MicrobiomeStat` provides a variety of
data input formats, supports both alpha and beta diversity analysis,
paired sample analysis, multiple visualization styles, and many other
features. This makes `MicrobiomeStat` more comprehensive and flexible
than many other packages dedicated to longitudinal data analysis.

## Comparative Table 2: Comprehensive Packages

|             Features              |                                                    MicrobiomeStat                                                     |                              microbiomeutilities                               |                         phyloseq                         |               microbiomemarker                |                                                               MicrobiomeAnalyst                                                                |                                             microbiomeeco                                             |                                                     EasyAmplicon                                                     |        STAMP         |                qiime2                 |                                            MicrobiotaProcess                                            |
|:---------------------------------:|:---------------------------------------------------------------------------------------------------------------------:|:------------------------------------------------------------------------------:|:--------------------------------------------------------:|:---------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------:|:--------------------:|:-------------------------------------:|:-------------------------------------------------------------------------------------------------------:|
|        Data Input Formats         |                  QIIME2, DADA2, BIOM, Mothur, Phyloseq, DGEList, DESeqDataSet, SummarizedExperiment                   |                                    phyloseq                                    |               QIIME2, DADA2, BIOM, Mothur                |     QIIME2, DADA2, BIOM, Mothur, Phyloseq     |                                                            txt, .csv, biom, mothur                                                             |                               QIIME, QIIME2, HUMAnN, Kraken2, phyloseq                                |                                                Paired fastq/fq files                                                 |         .tsv         |              fastq, biom              | qiime2, dada2, MetaPhlAn, Returns MPSE object, Phyloseq, TreeSummarizedExperiment, SummarizedExperiment |
| Supports Alpha Diversity Analysis |                                                          ✔️                                                           |                                       ✔️                                       |                            ✔️                            |                      ✔️                       |                                                                       ✔️                                                                       |                                                  ✔️                                                   |                                                          ✔️                                                          |          ❌          |                  ✔️                   |                                                   ✔️                                                    |
| Supports Beta Diversity Analysis  |                                                          ✔️                                                           |                                       ✔️                                       |                            ✔️                            |                      ✔️                       |                                                                       ✔️                                                                       |                                                  ✔️                                                   |                                                          ✔️                                                          |          ❌          |                  ✔️                   |                                                   ✔️                                                    |
|      Paired Sample Analysis       |                                                          ✔️                                                           |                                       ✔️                                       |                            ❌                            |                      ❌                       |                                                                       ❌                                                                       |                                                  ❌                                                   |                                                          ❌                                                          |          ❌          |                  ✔️                   |                                                   ❌                                                    |
|    Longitudinal Data Analysis     |                                                          ✔️                                                           |                                       ✔️                                       |                            ❌                            |                      ❌                       |                                                                       ❌                                                                       |                                                  ❌                                                   |                                                          ❌                                                          |          ❌          |                  ✔️                   |                                                   ❌                                                    |
|   Multiple Visualization Styles   | Boxplot, Violinplot, Scatterplot, Barplot, Heatmap, Dotplot, Ordination plot, Area plot, Volcano plot, Spaghetti plot | Density plot, Barplot, Stripchart, Boxplot, Heatmap, Area plot, Spaghetti plot | Boxplot, Barplot, Heatmap, Network plot, Ordination plot | Boxplot, Heatmap, Dotplot, Barplot, Cladogram | Stacked bar/area plot, Interactive pie chart, Rarefaction curve, Phylogenetic tree, Heat tree, Boxplot, Heatmap, Network plot, Ordination plot | Barplot, Boxplot, Heatmap, Pie chart, Venn diagram, Ordination plot, Circular heatmap, Sankey diagram | Barplot, Cladogram, Phylogenetic tree, Heatmap, Sankey diagram, Venn diagram, Ordination plot, Volcano plot, Boxplot | Barplot, Scatterplot | Scatterplot, Ordination plot, Boxplot |                   Barplot, Boxplot, Heatmap, Ordination plot, Tree plot, Volcano plot                   |
|         Publication Ready         |                                                          ✔️                                                           |                                       ✔️                                       |                            ❌                            |                      ❌                       |                                                                       ✔️                                                                       |                                                  ❌                                                   |                                                          ❌                                                          |          ✔️          |                  ❌                   |                                                   ❌                                                    |
|        Ongoing Development        |                                                          ✔️                                                           |                                       ❌                                       |                            ✔️                            |                      ❌                       |                                                                       ✔️                                                                       |                                                  ✔️                                                   |                                                          ✔️                                                          |          ❌          |                  ✔️                   |                                                   ✔️                                                    |

Comparative Table 2, comparing comprehensive packages, further
illustrates the strengths of `MicrobiomeStat`. Not only does
`MicrobiomeStat` support a wide range of data input formats and
diversity analyses, but it also excels in paired sample analysis and
longitudinal data analysis, among others.

In conclusion, `MicrobiomeStat` offers a robust and comprehensive
solution for your microbiome data analysis needs, outperforming many
other packages in both versatility and user-friendliness. Why not give
it a try and see the difference for yourself?

## Features at a glance 🌟

| Feature                        | Description                                                                                      |
|--------------------------------|--------------------------------------------------------------------------------------------------|
| Data Import and Conversion     | Supports numerous input formats from popular tools like QIIME2, Mothur, DADA2, Phyloseq and more |
| Cross-sectional Study Analysis | Performs comprehensive analysis of cross-sectional studies                                       |
| Paired Sample Analysis         | Excellent tool for analyzing paired samples                                                      |
| Longitudinal Study Analysis    | Allows for exploring the temporal dynamics of the microbiome                                     |
| Visualization Capabilities     | Offers a wide variety of visualization styles                                                    |
| Data Export                    | Supports export of analysis results in diverse formats                                           |
| Ongoing Development            | Continual feature refinement and new functionality addition                                      |

## Output interpretation 📝

To learn how to interpret the output objects and visualizations, visit
our wiki: [www.microbiomestat.wiki](http://www.microbiomestat.wiki).

## Known issues and limitations ❗️

We strive to make `MicrobiomeStat` as comprehensive as possible, but
there may still be limitations. All known issues can be found in the
issues section of our GitHub repository.

## Support & Contact 📬

| Name         | Email                       |
|--------------|-----------------------------|
| Jun Chen     | <Chen.Jun2@mayo.edu>        |
| Caffery Yang | <cafferychen7850@gmail.com> |

Stay updated with the latest `MicrobiomeStat` features by visiting our
wiki: [www.microbiomestat.wiki](http://www.microbiomestat.wiki)
