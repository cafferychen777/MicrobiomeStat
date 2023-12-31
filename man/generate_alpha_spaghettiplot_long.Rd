% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_alpha_spaghettiplot_long.R
\name{generate_alpha_spaghettiplot_long}
\alias{generate_alpha_spaghettiplot_long}
\title{Generate an alpha diversity line plot for longitudinal data}
\usage{
generate_alpha_spaghettiplot_long(
  data.obj,
  alpha.obj = NULL,
  alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
  depth = NULL,
  subject.var,
  time.var,
  t0.level,
  ts.levels,
  group.var = NULL,
  strata.var = NULL,
  adj.vars = NULL,
  base.size = 16,
  palette = NULL,
  theme.choice = "bw",
  custom.theme = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5,
  ...
)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.}

\item{alpha.obj}{An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.}

\item{alpha.name}{The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou". Previously named as `alpha.index`.}

\item{depth}{An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.}

\item{subject.var}{The name of the subject variable.}

\item{time.var}{The name of the time variable.}

\item{t0.level}{Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.}

\item{ts.levels}{Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.}

\item{group.var}{The name of the group variable.}

\item{strata.var}{The name of the strata variable (default is NULL).}

\item{adj.vars}{A character vector of variable names to be used for adjustment.}

\item{base.size}{The base font size for the plot.}

\item{palette}{An optional parameter specifying the color palette to be used for the plot.
It can be either a character string specifying the name of a predefined
palette or a vector of color codes in a format accepted by ggplot2
(e.g., hexadecimal color codes). Available predefined palettes include
'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
by various scientific publications and the `ggsci` package. If `palette`
is not provided or an unrecognized palette name is given, a default color
palette will be used. Ensure the number of colors in the palette is at
least as large as the number of groups being plotted.}

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

\item{pdf}{Logical, whether to save the plot as a PDF (default is TRUE).}

\item{file.ann}{The annotation to be added to the PDF file name (default is NULL).}

\item{pdf.wid}{The width of the output PDF file. Default is 11.}

\item{pdf.hei}{The height of the output PDF file. Default is 8.5.}

\item{...}{Additional arguments passed to ggplot().}
}
\value{
A ggplot object of the alpha diversity line plot.
}
\description{
This function creates a ggplot object of alpha diversity (e.g., Shannon index) line plot for longitudinal data,
showing individual subject trajectories and the mean trajectory for each group.
}
\examples{
\dontrun{
data("subset_T2D.obj")
T2D.alpha.obj <- mStat_calculate_alpha_diversity(subset_T2D.obj$feature.tab,"shannon")
generate_alpha_spaghettiplot_long(
  data.obj = subset_T2D.obj,
  alpha.obj = T2D.alpha.obj,
  alpha.name = c("shannon"),
  subject.var = "subject_id",
  time.var = "visit_number",
  t0.level = sort(unique(subset_T2D.obj$meta.dat$visit_number))[1],
  ts.levels = sort(unique(subset_T2D.obj$meta.dat$visit_number))[-1],
  group.var = "subject_gender",
  strata.var = "subject_race",
  adj.vars = "sample_body_site",
  theme.choice = "bw",
  palette = "lancet",
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)

data("ecam.obj")
generate_alpha_spaghettiplot_long(
  data.obj = ecam.obj,
  alpha.obj = NULL,
  alpha.name = c("shannon","simpson", "observed_species"),
  subject.var = "subject.id",
  time.var = "month_num",
  t0.level = NULL,
  ts.levels = NULL,
  group.var = "delivery",
  strata.var = "diet",
  adj.vars = NULL,
  theme.choice = "bw",
  palette = "lancet",
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
}
}
