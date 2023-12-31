% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_change_dotplot_pair.R
\name{generate_taxa_change_dotplot_pair}
\alias{generate_taxa_change_dotplot_pair}
\title{Generate Taxa Stack Dotplot Pair}
\usage{
generate_taxa_change_dotplot_pair(
  data.obj,
  subject.var,
  time.var,
  group.var = NULL,
  strata.var = NULL,
  change.base = "1",
  feature.change.func = "log fold change",
  feature.level = NULL,
  feature.dat.type = c("count", "proportion", "other"),
  features.plot = NULL,
  top.k.plot = NULL,
  top.k.func = NULL,
  prev.filter = 0.001,
  abund.filter = 0.001,
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
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).}

\item{subject.var}{A character string defining subject variable in meta_tab}

\item{time.var}{A character string defining time variable in meta_tab}

\item{group.var}{A character string defining group variable in meta_tab used for sorting and facetting}

\item{strata.var}{(Optional) A character string defining strata variable in meta_tab used for sorting and facetting}

\item{change.base}{A numeric value setting base for the change (usually 1)}

\item{feature.change.func}{A method or function specifying how to compute the change in feature abundance or prevalence between two time points.
The following options are available:

- A custom function: If you provide a user-defined function, it should take two numeric arguments corresponding to the values at the two time points and return the computed change. This function can be applied to compute changes both in abundance (`time1_mean_abundance` and `time2_mean_abundance`) and prevalence (`time1_prevalence` and `time2_prevalence`).

- "log fold change": Computes the log2 fold change between the two time points. To handle zeros, a small offset (0.00001) is added before taking the logarithm. This method can be applied for both abundance and prevalence changes.

- "relative change": Computes the relative change as `(time2_value - time1_value) / (time2_value + time1_value)`. If both time points have a value of 0, the change is defined as 0. This method can be applied for both abundance and prevalence changes.

- "absolute change": Computes the difference between the values at the two time points. This method can be applied for both abundance and prevalence changes.

- Any other value (or if the parameter is omitted): By default, the function will compute the absolute change as described above, regardless of whether it is abundance or prevalence data.}

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

\item{top.k.plot}{Integer specifying number of top abundant features to plot, when `features.plot` is NULL.
Default is NULL, in which case all features passing filters will be plotted.}

\item{top.k.func}{Function to use for selecting top abundant features, when `features.plot` is NULL.
Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
case features will be selected by mean abundance.}

\item{prev.filter}{Numeric value specifying the minimum prevalence threshold for filtering
taxa before analysis. Taxa with prevalence below this value will be removed.
Prevalence is calculated as the proportion of samples where the taxon is present.
Default 0 removes no taxa by prevalence filtering.}

\item{abund.filter}{Numeric value specifying the minimum abundance threshold for filtering
taxa before analysis. Taxa with mean abundance below this value will be removed.
Abundance refers to counts or proportions depending on \code{feature.dat.type}.
Default 0 removes no taxa by abundance filtering.}

\item{base.size}{Base font size for the generated plots.}

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

\item{palette}{Color palette used for the plots.}

\item{pdf}{If TRUE, save the plot as a PDF file (default: TRUE)}

\item{file.ann}{(Optional) A character string specifying a file annotation to include in the generated PDF file's name}

\item{pdf.wid}{Width of the PDF plots.}

\item{pdf.hei}{Height of the PDF plots.}

\item{...}{Additional parameters to be passed}
}
\value{
If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
}
\description{
This function generates a stacked dotplot of specified taxa level with paired samples. The data used in this
visualization will be first filtered based on prevalence and abundance thresholds. The plot can either be displayed
interactively or saved as a PDF file.
}
\examples{
\dontrun{

# Note: In the RStudio viewer, the plot might appear cluttered if there are many taxa.
# It's recommended to view the generated PDF for better clarity. If it still feels
# overcrowded in the PDF, consider increasing the 'pdf.wid' value to adjust the width of the plot.

data(peerj32.obj)
generate_taxa_change_dotplot_pair(
  data.obj = peerj32.obj,
  subject.var = "subject",
  time.var = "time",
  group.var = "group",
  strata.var = "sex",
  change.base = "1",
  feature.change.func = "log fold change",
  feature.level = "Family",
  feature.dat.type = "count",
  features.plot = NULL,
  top.k.plot = 20,
  top.k.func = "mean",
  prev.filter = 0.01,
  abund.filter = 1e-4,
  base.size = 16,
  theme.choice = "bw",
  custom.theme = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 30,
  pdf.hei = 10
)

data("subset_pairs.obj")
generate_taxa_change_dotplot_pair(
  data.obj = subset_pairs.obj,
  subject.var = "MouseID",
  time.var = "Antibiotic",
  group.var = "Sex",
  strata.var = NULL,
  change.base = "Baseline",
  feature.change.func = "log fold change",
  feature.level = "Family",
  feature.dat.type = "count",
  features.plot = NULL,
  top.k.plot = 20,
  top.k.func = "mean",
  prev.filter = 0.01,
  abund.filter = 1e-4,
  base.size = 16,
  theme.choice = "bw",
  custom.theme = NULL,
  pdf = TRUE,
  file.ann = NULL,
  pdf.wid = 30,
  pdf.hei = 10
)
}
# View the result
}
