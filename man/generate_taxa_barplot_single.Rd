% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_barplot_single.R
\name{generate_taxa_barplot_single}
\alias{generate_taxa_barplot_single}
\title{Generate Stacked Taxa Barplots for Longitudinal Data}
\usage{
generate_taxa_barplot_single(
  data.obj,
  subject.var,
  time.var = NULL,
  t.level = NULL,
  group.var = NULL,
  strata.var = NULL,
  feature.level = "original",
  feature.dat.type = c("count", "proportion", "other"),
  feature.number = 20,
  base.size = 10,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5,
  ...
)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.}

\item{subject.var}{A string indicating the variable for subject identifiers.}

\item{time.var}{A string indicating the variable for time points.}

\item{t.level}{Character string specifying the time level/value to subset data to,
if a time variable is provided. Default NULL does not subset data.}

\item{group.var}{A string indicating the variable for group identifiers. Default is NULL.}

\item{strata.var}{A string indicating the variable for strata identifiers. Default is NULL.}

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

\item{feature.number}{A numeric value indicating the number of top abundant features to retain in the plot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 20.}

\item{base.size}{A numeric value indicating the base font size for the plot.}

\item{theme.choice}{Plot theme choice. Can be one of:
  - "prism": ggprism::theme_prism()
  - "classic": theme_classic()
  - "gray": theme_gray()
  - "bw": theme_bw()
Default is "bw".}

\item{custom.theme}{A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. To use a custom theme, first create a theme object with ggplot2::theme(), then pass it to this argument. For example:

```r
my_theme <- ggplot2::theme(
  axis.title = ggplot2::element_text(size=16, color="red"),
  legend.position = "none"
)
```

Then pass `my_theme` to `custom.theme`. Default is NULL, which will use the default theme based on `theme.choice`.}

\item{palette}{A character vector specifying the color palette. Default is NULL.}

\item{pdf}{A logical value indicating whether to save the plot as a PDF. Default is TRUE.}

\item{file.ann}{A string for additional annotation to the file name. Default is NULL.}

\item{pdf.wid}{A numeric value specifying the width of the PDF file. Default is 11.}

\item{pdf.hei}{A numeric value specifying the height of the PDF file. Default is 8.5.}

\item{...}{Additional arguments to be passed to the function.}
}
\value{
A list of ggplot objects of the taxa barplots.
}
\description{
This function generates stacked barplots to visualize the taxonomic composition of samples for longitudinal data.
It also provides options for grouping and stratifying data.
}
\details{
This function generates a stacked barplot of taxa proportions for longitudinal data.
The barplot can be stratified by a group variable and/or other variables.
It also allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
}
\examples{
\dontrun{
data(peerj32.obj)
generate_taxa_barplot_single(
  data.obj = peerj32.obj,
  subject.var = "subject",
  time.var = NULL,
  t.level = NULL,
  group.var = "group",
  strata.var = "sex",
  feature.level = "Family",
  feature.dat.type = "count",
  feature.number = 30,
  base.size = 10,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
generate_taxa_barplot_single(
  data.obj = peerj32.obj,
  subject.var = "subject",
  time.var = NULL,
  t.level = NULL,
  group.var = "group",
  strata.var = NULL,
  feature.level = "Family",
  feature.dat.type = "count",
  feature.number = 30,
  base.size = 10,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
data("subset_T2D.obj")
generate_taxa_barplot_single(
  data.obj = subset_T2D.obj,
  subject.var = "subject_id",
  time.var = "visit_number_num",
  t.level = 1,
  group.var = "subject_race",
  strata.var = "subject_gender",
  feature.level = c("Phylum", "Family", "Genus"),
  feature.dat.type = "count",
  feature.number = 30,
  base.size = 10,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
generate_taxa_barplot_single(
  data.obj = subset_T2D.obj,
  subject.var = "subject_id",
  time.var = "visit_number_num",
  t.level = 1,
  group.var = "subject_race",
  strata.var = NULL,
  feature.level = c("Phylum", "Family", "Genus"),
  feature.dat.type = "count",
  feature.number = 30,
  base.size = 10,
  theme.choice = "bw",
  custom.theme = NULL,
  palette = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
}
}
