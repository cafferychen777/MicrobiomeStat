% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_change_heatmap_pair.R
\name{generate_taxa_change_heatmap_pair}
\alias{generate_taxa_change_heatmap_pair}
\title{Generate Taxa Change Heatmap Pair}
\usage{
generate_taxa_change_heatmap_pair(
  data.obj,
  subject.var,
  time.var,
  group.var = NULL,
  strata.var = NULL,
  change.base = NULL,
  feature.change.func = "relative change",
  feature.level = NULL,
  feature.dat.type = c("count", "proportion", "other"),
  features.plot = NULL,
  top.k.plot = NULL,
  top.k.func = NULL,
  prev.filter = 0.01,
  abund.filter = 0.01,
  base.size = 10,
  palette = NULL,
  cluster.rows = NULL,
  cluster.cols = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5,
  ...
)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.}

\item{subject.var}{A character string specifying the subject variable in the metadata.}

\item{time.var}{A character string specifying the time variable in the metadata.}

\item{group.var}{A character string specifying the grouping variable in the metadata. Default is NULL.}

\item{strata.var}{A character string specifying the stratification variable in the metadata. Default is NULL.}

\item{change.base}{A numeric value specifying the baseline time point for computing change. This should match one of the time points in the time variable. Default is 1, which assumes the first time point is the baseline.}

\item{feature.change.func}{Specifies the method or function to compute the change between two time points.
The following options are available:

- A custom function: If you provide a user-defined function, it should take two numeric arguments corresponding to the values at the two time points (`value_time_1` and `value_time_2`) and return the computed change. This custom function will be applied directly.

- "log fold change": Computes the log2 fold change between the two time points. For zero values, imputation is performed using half of the minimum nonzero value for each feature level at the respective time point before taking the logarithm.

- "relative change": Computes the relative change as `(value_time_2 - value_time_1) / (value_time_2 + value_time_1)`. If both time points have a value of 0, the change is defined as 0.

- "absolute change": Computes the difference between the values at the two time points.

- Any other value (or if the parameter is omitted): By default, the function computes the absolute change as described above.}

\item{feature.level}{The column name in the feature annotation matrix (feature.ann) of data.obj
to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
is also NULL.}

\item{feature.dat.type}{The type of the feature data, which determines how the data is handled in downstream analyses.
Should be one of:
- "count": Raw count data, will be normalized by the function.
- "proportion": Data that has already been normalized to proportions/percentages.
- "other": Custom abundance data that has unknown scaling. No normalization applied.
The choice affects preprocessing steps as well as plot axis labels.
Default is "count", which assumes raw OTU table input.}

\item{features.plot}{A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.}

\item{top.k.plot}{Integer specifying number of top k features to plot, when `features.plot` is NULL.
Default is NULL, which case all features passing filters will be plotted.}

\item{top.k.func}{Function to use for selecting top k features, when `features.plot` is NULL.
Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
case features will be selected by abundance.}

\item{prev.filter}{Numeric value specifying the minimum prevalence threshold for filtering
taxa before analysis. Taxa with prevalence below this value will be removed.
Prevalence is calculated as the proportion of samples where the taxon is present.
Default 0 removes no taxa by prevalence filtering.}

\item{abund.filter}{Numeric value specifying the minimum abundance threshold for filtering
taxa before analysis. Taxa with mean abundance below this value will be removed.
Abundance refers to counts or proportions depending on \code{feature.dat.type}.
Default 0 removes no taxa by abundance filtering.}

\item{base.size}{Base font size for the generated plots.}

\item{palette}{Specifies the color palette to be used for annotating groups and strata in the heatmap.
The parameter can be provided in several ways:
- As a character string denoting a predefined palette name.
  Available predefined palettes include 'npg', 'aaas', 'nejm',
  'lancet', 'jama', 'jco', and 'ucscgb', sourced from the `mStat_get_palette` function.
- As a vector of color codes in a format accepted by ggplot2
  (e.g., hexadecimal color codes).
If `palette` is NULL or an unrecognized string, a default color palette will be used.
The function assigns colors from this palette to the unique levels of
`group.var` and, if provided, `strata.var`. When both `group.var` and
`strata.var` are present, `group.var` levels are colored using the
beginning of the palette, while `strata.var` levels are colored using
the reversed palette, ensuring a distinct color representation for each.
If only `group.var` is provided, its levels are assigned colors from the
palette sequentially. If neither `group.var` nor `strata.var` is provided,
no annotation colors are applied.}

\item{cluster.rows}{A logical variable indicating if rows should be clustered. Default is TRUE.}

\item{cluster.cols}{A logical variable indicating if columns should be clustered. Default is NULL.}

\item{pdf}{If TRUE, save the plot as a PDF file (default: TRUE)}

\item{file.ann}{(Optional) A character string specifying a file annotation to include in the generated PDF file's name.}

\item{pdf.wid}{Width of the PDF plots.}

\item{pdf.hei}{Height of the PDF plots.}

\item{...}{Additional parameters to be passed to pheatmap function}
}
\value{
If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the pheatmap::pheatmap plot. If `pdf` is set to FALSE, the function will return the pheatmap plot without creating a PDF file.
}
\description{
This function generates a heatmap showing the pairwise changes in relative abundances of taxa between different time points.
The data used in this visualization will be first filtered based on prevalence and abundance thresholds.
The plot can either be displayed interactively or saved as a PDF file.
}
\examples{
\dontrun{
# Load required libraries and example data
library(pheatmap)
data(peerj32.obj)
generate_taxa_change_heatmap_pair(
  data.obj = peerj32.obj,
  subject.var = "subject",
  time.var = "time",
  group.var = "group",
  strata.var = "sex",
  change.base = "1",
  feature.change.func = "relative change",
  feature.level = c("Genus"),
  feature.dat.type = "count",
  features.plot = NULL,
  top.k.plot = 10,
  top.k.func = "sd",
  prev.filter = 0.1,
  abund.filter = 0.001,
  base.size = 10,
  palette = NULL,
  cluster.rows = NULL,
  cluster.cols = FALSE,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)

generate_taxa_change_heatmap_pair(
  data.obj = peerj32.obj,
  subject.var = "subject",
  time.var = "time",
  group.var = "group",
  strata.var = "sex",
  change.base = "1",
  feature.change.func = "relative change",
  feature.level = c("Genus"),
  feature.dat.type = "count",
  features.plot = NULL,
  top.k.plot = 10,
  top.k.func = "sd",
  prev.filter = 0.1,
  abund.filter = 0.001,
  base.size = 10,
  palette = NULL,
  cluster.rows = FALSE,
  cluster.cols = FALSE,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)

data(subset_pairs.obj)
generate_taxa_change_heatmap_pair(
  data.obj = subset_pairs.obj,
  subject.var = "MouseID",
  time.var = "Antibiotic",
  group.var = "Sex",
  strata.var = NULL,
  change.base = "Baseline",
  feature.change.func = "relative change",
  feature.level = c("Genus"),
  feature.dat.type = "count",
  features.plot = NULL,
  top.k.plot = 10,
  top.k.func = "sd",
  prev.filter = 0.1,
  abund.filter = 0.001,
  base.size = 10,
  palette = NULL,
  cluster.rows = NULL,
  cluster.cols = FALSE,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)

generate_taxa_change_heatmap_pair(
  data.obj = subset_pairs.obj,
  subject.var = "MouseID",
  time.var = "Antibiotic",
  group.var = "Sex",
  strata.var = NULL,
  change.base = "Baseline",
  feature.change.func = "relative change",
  feature.level = c("Genus"),
  feature.dat.type = "count",
  features.plot = NULL,
  top.k.plot = 10,
  top.k.func = "sd",
  prev.filter = 0.1,
  abund.filter = 0.001,
  base.size = 10,
  palette = NULL,
  cluster.rows = FALSE,
  cluster.cols = FALSE,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
}
}
