% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_alpha_change_boxplot_pair.R
\name{generate_alpha_change_boxplot_pair}
\alias{generate_alpha_change_boxplot_pair}
\title{Generate boxplot comparing change in specified alpha diversity index}
\usage{
generate_alpha_change_boxplot_pair(
  data.obj,
  alpha.obj = NULL,
  alpha.name = c("shannon", "observed_species"),
  depth = NULL,
  subject.var,
  time.var,
  group.var = NULL,
  strata.var = NULL,
  adj.vars = NULL,
  change.base = NULL,
  alpha.change.func = c("log fold change"),
  base.size = 16,
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

\item{alpha.obj}{An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.}

\item{alpha.name}{The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou". Previously named as `alpha.index`.}

\item{depth}{An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.}

\item{subject.var}{The variable in the metadata table that represents the subject.}

\item{time.var}{The variable in the metadata table that represents the time.}

\item{group.var}{(Optional) The variable in the metadata table that represents the grouping factor.}

\item{strata.var}{(Optional) The variable in the metadata table that represents the stratification factor.}

\item{adj.vars}{A character vector of variable names to be used for adjustment.}

\item{change.base}{The base time for calculating the change in alpha diversity.}

\item{alpha.change.func}{Function or method for calculating change in alpha diversity
  between two timepoints. This allows flexible options to quantify change:

  - If a function is provided: The function will be applied to compare alpha diversity
    at timepoint t vs baseline t0. The function should take two arguments
    representing the alpha diversity values at t and t0. For instance, a custom function to
    calculate the percentage change might look like:
    \preformatted{
      percentage_change <- function(t, t0) {
        return ((t - t0) / t0) * 100
      }
    }
    You can then pass this function as the value for `alpha.change.func`.

  - If a string is provided, the following options are supported:
    - 'log fold change': Calculates the log2 fold change of alpha diversity at t compared to t0.
    - 'absolute change': Calculates the absolute difference in alpha diversity at t compared to t0.
    - Any other value: A warning will be given that the provided method is not recognized,
      and the default method ('absolute change') will be used.

  - Default behavior (if no recognized string or function is provided) is to compute the absolute difference between t and t0.}

\item{base.size}{The base font size for the plot.}

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

\item{palette}{An optional parameter specifying the color palette to be used for the plot.
It can be either a character string specifying the name of a predefined
palette or a vector of color codes in a format accepted by ggplot2
(e.g., hexadecimal color codes). Available predefined palettes include
'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
by various scientific publications and the `ggsci` package. If `palette`
is not provided or an unrecognized palette name is given, a default color
palette will be used. Ensure the number of colors in the palette is at
least as large as the number of groups being plotted.}

\item{pdf}{(Optional) A boolean indicating whether to save the output as a PDF file.}

\item{file.ann}{(Optional) A string for annotating the output file name.}

\item{pdf.wid}{The width of the output PDF file. Default is 11.}

\item{pdf.hei}{The height of the output PDF file. Default is 8.5.}

\item{...}{(Optional) Additional arguments to pass to the plotting function.}
}
\value{
A boxplot displaying the change in the specified alpha diversity index between two time points, stratified by the specified grouping and/or strata variables (if provided). The boxplot will be saved as a PDF if `pdf` is set to `TRUE`.
}
\description{
This function generates a boxplot comparing the change in a specified alpha diversity index between two time points. The change can be calculated as the log or the absolute value. Several optional arguments are available for customizing the output, such as strata or group variables.
}
\details{
This parameter allows flexible quantification of how alpha diversity
  changes from baseline. Log-ratio is commonly used to compare relative
  difference. Absolute difference indicates the magnitude of change.
  Custom functions can also be supplied to calculate change as needed.

This function extracts metadata, calculates alpha diversity for specified
  indices, compares values between two timepoints, applies log2 fold change by default,
  visualizes the change using boxplots and optional faceting and saving as PDF.
}
\examples{
\dontrun{
library(vegan)
data(peerj32.obj)
generate_alpha_change_boxplot_pair(
  data.obj = peerj32.obj,
  alpha.obj = NULL,
  alpha.name = c("simpson"),
  subject.var = "subject",
  time.var = "time",
  group.var = "group",
  strata.var = "sex",
  adj.vars = "sex",
  change.base = "1",
  alpha.change.func = "absolute change",
  base.size = 16,
  theme.choice = "bw",
  palette = "lancet",
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
data("subset_pairs.obj")
generate_alpha_change_boxplot_pair(
  data.obj = subset_pairs.obj,
  alpha.obj = NULL,
  alpha.name = c("simpson"),
  subject.var = "MouseID",
  time.var = "Antibiotic",
  group.var = "Sex",
  strata.var = NULL,
  adj.vars = NULL,
  change.base = "Baseline",
  alpha.change.func = "log fold change",
  base.size = 16,
  theme.choice = "bw",
  palette = "lancet",
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 11,
  pdf.hei = 8.5
)
}
}
