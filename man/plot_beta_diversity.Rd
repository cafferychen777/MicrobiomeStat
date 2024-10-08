% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_beta_diversity.R
\name{plot_beta_diversity}
\alias{plot_beta_diversity}
\title{Plot Beta Diversity}
\usage{
plot_beta_diversity(
  dist.obj,
  pc.obj = NULL,
  meta.dat,
  measure,
  group.var = NULL,
  strata.var = NULL,
  subject.var = NULL,
  time.var = NULL,
  adj.vars = NULL,
  time.point.plot,
  is.plot.change = TRUE,
  plot.type = c("PcoA", "boxplot", "spaghettiplot"),
  ...
)
}
\arguments{
\item{dist.obj}{A distance matrix or object containing beta diversity measures.}

\item{pc.obj}{A matrix or data frame containing principal coordinates.}

\item{meta.dat}{A data frame containing metadata for the samples.}

\item{measure}{Character string specifying the beta diversity measure to plot.}

\item{group.var}{Character string specifying the grouping variable in meta.dat.}

\item{strata.var}{Character string specifying the stratification variable in meta.dat.}

\item{subject.var}{Character string specifying the subject variable in meta.dat.}

\item{time.var}{Character string specifying the time variable in meta.dat.}

\item{adj.vars}{Character vector specifying variables to adjust for. If specified, adjusted PCoA will be performed.}

\item{time.point.plot}{Character vector specifying time points to plot. The first point will be the reference.}

\item{is.plot.change}{Logical, whether to plot change from baseline. Only effective when there are multiple time points.}

\item{plot.type}{Character string specifying the plot type, one of c("PCoA", "boxplot", "spaghettiplot").}

\item{...}{Additional arguments passed to the underlying plotting functions.}
}
\value{
A ggplot object representing the beta diversity plot.
}
\description{
This function, part of the MicrobiomeStat package, generates various types of plots
to visualize beta diversity measures in microbiome data. It supports different
visualization schemes for single time point, paired time points, and longitudinal data.
}
\details{
This function provides a flexible framework for visualizing beta diversity in microbiome data.
It can handle single time point, paired time points, and longitudinal data. The function offers
different plot types including PCoA plots, box plots, and spaghetti plots, and can visualize
changes over time.

The behavior of the function depends on the number of time points specified:
- For a single time point or when time variable is not specified, only PCoA plot is available.
- For two time points, both PCoA and change boxplot are available.
- For more than two time points, PCoA, change boxplot, and change spaghettiplot are available.

When `is.plot.change` is TRUE and multiple time points are specified, the function will calculate
changes from the baseline (first time point).
}
\examples{
data(peerj32.obj)
dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, c("BC", "Jaccard"))
plot_beta_diversity(dist.obj = dist.obj,
                    meta.dat = peerj32.obj$meta.dat,
                    measure = c("BC", "Jaccard"),
                    group.var = "group",
                    time.var = "time",
                    time.point.plot = c("1"),
                    plot.type = "PCoA")

plot_beta_diversity(dist.obj = dist.obj,
                    meta.dat = peerj32.obj$meta.dat,
                    measure = c("BC", "Jaccard"),
                    group.var = "group",
                    subject.var = "subject",
                    time.var = "time",
                    time.point.plot = c("1","2"),
                    plot.type = "PCoA",
                    is.plot.change = FALSE)

plot_beta_diversity(dist.obj = dist.obj,
                    meta.dat = peerj32.obj$meta.dat,
                    measure = c("BC", "Jaccard"),
                    group.var = "group",
                    subject.var = "subject",
                    time.var = "time",
                    time.point.plot = c("1","2"),
                    plot.type = "boxplot",
                    is.plot.change = TRUE)

data(ecam.obj)
dist.obj <- mStat_calculate_beta_diversity(ecam.obj, c("BC", "Jaccard"))
plot_beta_diversity(dist.obj = dist.obj,
                    meta.dat = ecam.obj$meta.dat,
                    measure = c("BC", "Jaccard"),
                    group.var = "delivery",
                    subject.var = "subject.id",
                    time.var = "month_num",
                    time.point.plot = unique(ecam.obj$meta.dat$month_num),
                    plot.type = "boxplot",
                    is.plot.change = TRUE)

plot_beta_diversity(dist.obj = dist.obj,
                    meta.dat = ecam.obj$meta.dat,
                    measure = c("BC", "Jaccard"),
                    group.var = "delivery",
                    subject.var = "subject.id",
                    time.var = "month_num",
                    time.point.plot = unique(ecam.obj$meta.dat$month_num),
                    plot.type = "spaghettiplot",
                    is.plot.change = TRUE)

plot_beta_diversity(dist.obj = dist.obj,
                    meta.dat = ecam.obj$meta.dat,
                    measure = c("BC", "Jaccard"),
                    group.var = "delivery",
                    subject.var = "subject.id",
                    time.var = "month_num",
                    time.point.plot = unique(ecam.obj$meta.dat$month_num),
                    plot.type = "PCoA",
                    is.plot.change = FALSE)

}
