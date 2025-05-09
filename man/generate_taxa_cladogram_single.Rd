% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_cladogram_single.R
\name{generate_taxa_cladogram_single}
\alias{generate_taxa_cladogram_single}
\title{Generate a Circular Cladogram with Heatmap for Taxa}
\usage{
generate_taxa_cladogram_single(
  data.obj,
  test.list = NULL,
  group.var = NULL,
  feature.level,
  feature.mt.method = "none",
  cutoff = 1,
  color.group.level = NULL,
  palette = NULL,
  pdf = FALSE,
  pdf.width = 10,
  pdf.height = 10,
  time.var = NULL,
  t.level = NULL,
  adj.vars = NULL,
  prev.filter = 0.1,
  abund.filter = 1e-04,
  feature.dat.type = "count"
)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which includes components like feature.tab, feature.ann, meta.dat, etc.}

\item{test.list}{A list of test results. If NULL, it will be generated using generate_taxa_test_single.}

\item{group.var}{The name of the grouping variable in meta.dat.}

\item{feature.level}{A character vector specifying taxonomic levels to be analyzed.}

\item{feature.mt.method}{Multiple testing method for features, "none" (default), "fdr", or other methods supported by p.adjust.}

\item{cutoff}{The p-value cutoff for significance.}

\item{color.group.level}{The taxonomic level used to color-code the branches of the cladogram.}

\item{palette}{An optional vector of colors to be used for the plot. If NULL, a default color palette will be used.}

\item{pdf}{Boolean indicating whether to save the plot as a PDF.}

\item{pdf.width}{The width of the PDF file if saved.}

\item{pdf.height}{The height of the PDF file if saved.}

\item{time.var}{Character string specifying the column name in metadata containing time variable. Used when test.list is NULL.}

\item{t.level}{Character string specifying the time level/value to subset data to. Used when test.list is NULL.}

\item{adj.vars}{Character vector specifying column names in metadata containing covariates. Used when test.list is NULL.}

\item{prev.filter}{Numeric value specifying the minimum prevalence threshold for filtering taxa. Used when test.list is NULL.}

\item{abund.filter}{Numeric value specifying the minimum abundance threshold for filtering taxa. Used when test.list is NULL.}

\item{feature.dat.type}{The type of the feature data, which determines data handling. Should be one of "count", "proportion", or "other". For CLR-transformed data, use "other". Used when test.list is NULL.}
}
\value{
A ggplot object representing the circular heatmap with phylogenetic tree.

A phylogenetic tree object
}
\description{
This function generates a circular cladogram with an integrated heatmap for taxonomic data.
It visualizes the phylogenetic relationships between different taxa and their abundances or other
coefficients across different taxonomic levels using a tree-like structure (cladogram).
}
\examples{
\dontrun{
data(subset_T2D.obj)

test.list <- generate_taxa_test_single(
    data.obj = subset_T2D.obj,
    time.var = "visit_number",
    t.level = NULL,
    group.var = "subject_race",
    adj.vars = "subject_gender",
    feature.level = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    feature.dat.type = "count",
    prev.filter = 0.1,
    abund.filter = 0.0001,
)

plot.list <- generate_taxa_cladogram_single(
  data.obj = subset_T2D.obj,
  test.list = test.list,
  group.var = "subject_gender",
  feature.level = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  feature.mt.method = "none",
  cutoff = 0.9,
  color.group.level = "Order"
)

test.list <- generate_taxa_test_single(
    data.obj = subset_T2D.obj,
    time.var = "visit_number",
    t.level = NULL,
    group.var = "subject_race",
    adj.vars = "subject_gender",
    feature.level = c("Order"),
    feature.dat.type = "count",
    prev.filter = 0.1,
    abund.filter = 0.0001,
)

plot.list <- generate_taxa_cladogram_single(
  data.obj = subset_T2D.obj,
  test.list = test.list,
  group.var = "subject_gender",
  feature.level = c("Order"),
  feature.mt.method = "none",
  cutoff = 0.9,
  color.group.level = "Order"
)

data(peerj32.obj)

test.list <- generate_taxa_test_single(
    data.obj = peerj32.obj,
    time.var = "time",
    t.level = NULL,
    group.var = "group",
    adj.vars = "sex",
    feature.level = c("Phylum","Family","Genus"),
    feature.dat.type = "count",
    prev.filter = 0.1,
    abund.filter = 0.0001,
)

plot.list <- generate_taxa_cladogram_single(
  data.obj = peerj32.obj,
  test.list = test.list,
  group.var = "group",
  feature.level = c("Phylum", "Family", "Genus"),
  cutoff = 0.3,
  color.group.level = "Family"
)
}

