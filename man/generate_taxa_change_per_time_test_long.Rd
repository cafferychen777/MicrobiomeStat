% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_change_per_time_test_long.R
\name{generate_taxa_change_per_time_test_long}
\alias{generate_taxa_change_per_time_test_long}
\title{Longitudinal Taxa Test in Microbiome Data}
\usage{
generate_taxa_change_per_time_test_long(
  data.obj,
  subject.var,
  time.var = NULL,
  t0.level = NULL,
  ts.levels = NULL,
  group.var = NULL,
  adj.vars = NULL,
  feature.level,
  feature.change.func = "relative change",
  feature.dat.type = c("count", "proportion", "other"),
  prev.filter = 0.001,
  abund.filter = 0.001,
  ...
)
}
\arguments{
\item{data.obj}{A MicrobiomeStat data object containing microbiome data and metadata.}

\item{subject.var}{A string specifying the column name in meta.dat that uniquely identifies each subject.}

\item{time.var}{Optional; a string representing the time variable in the meta.dat. If provided, enables longitudinal analysis.}

\item{t0.level}{Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.}

\item{ts.levels}{Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.}

\item{group.var}{Optional; a string specifying the group variable in meta.dat for between-group comparisons.}

\item{adj.vars}{Optional; a vector of strings representing covariates in meta.dat for adjustment in the analysis.}

\item{feature.level}{A string or vector of strings indicating the taxonomic level(s) for analysis (e.g., "Phylum", "Class").}

\item{feature.change.func}{A function or character string specifying how to calculate
the change from baseline value. This allows flexible options:
- If a function is provided, it will be applied to each row to calculate change.
  The function should take 2 arguments: value at timepoint t and value at baseline t0.
- If a character string is provided, following options are supported:
  - 'relative change': (value_t - value_t0) / (value_t + value_t0)
  - 'absolute change': value_t - value_t0
  - 'log fold change': log2(value_t + 1e-5) - log2(value_t0 + 1e-5)
- Default is 'relative change'.

If none of the above options are matched, an error will be thrown indicating
the acceptable options or prompting the user to provide a custom function.}

\item{feature.dat.type}{Character; "count" or "proportion", indicating the type of feature data.}

\item{prev.filter}{Numeric; a minimum prevalence threshold for taxa inclusion in the analysis.}

\item{abund.filter}{Numeric; a minimum abundance threshold for taxa inclusion in the analysis.}

\item{...}{Additional arguments passed to other methods.}
}
\value{
A nested list structure. The top level of the list corresponds to different time points, and each element contains a list of dataframes for each taxonomic level. Each dataframe provides statistical analysis results for taxa at that level and time point.

This function is especially useful for longitudinal microbiome studies, facilitating the exploration of temporal patterns in microbial communities. By analyzing different time points against a baseline, it helps to uncover significant temporal shifts in the abundance of various taxa.

The function is tailored for investigations that aim to monitor changes in microbial communities over time, such as in response to treatments or environmental changes. The structured output assists in interpreting temporal trends and identifying key taxa that contribute to these changes.
}
\description{
This function performs a comprehensive analysis of microbiome data, focusing on the longitudinal trends of various taxa. It is specifically designed to work with data where the primary interest is to understand how the abundance of microbial taxa changes over time, across different groups, or under different conditions.
}
\details{
Utilizing a standard linear model approach, the function is adept at identifying significant temporal variations in taxa abundance. It provides a robust framework for comparing microbial communities at different time points, thereby offering valuable insights into the dynamics of these communities over extended periods.

The function's ability to handle both count and proportion data types, along with its features for adjusting covariates, makes it a versatile tool for microbiome research. It is particularly beneficial for studies investigating the effects of treatments, environmental changes, or other interventions on the microbiome's composition and behavior over time.


The function integrates various data manipulations, normalization procedures, and statistical tests to assess the significance of taxa changes over time or between groups. It allows for the adjustment of covariates and is capable of handling both count and proportion data types.

The function uses a standard linear model (lm) to analyze the data. It handles fixed effects to account for the influence of different variables on the taxa. Filtering is performed based on prevalence and abundance thresholds, and normalization and aggregation procedures are applied as necessary.

A key feature of the function is its ability to conduct differential abundance analysis separately for each time point in the longitudinal data. This method is particularly effective for identifying significant changes in taxa at specific time points, offering insights into the temporal dynamics of the microbiome.
}
\examples{
\dontrun{
# Example1: Analyzing the Type 2 Diabetes dataset
data("subset_T2D.obj")
# Longitudinal analysis of microbial changes in different racial groups
result <- generate_taxa_change_per_time_test_long(
  data.obj = subset_T2D.obj,
  subject.var = "subject_id",
  time.var = "visit_number",
  t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
  ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
  group.var = "subject_race",
  adj.vars = "subject_gender",
  prev.filter = 0.1,
  abund.filter = 0.001,
  feature.level = c("Genus", "Family"),
  feature.dat.type = "count"
)
# Visualizing the results for the Type 2 Diabetes dataset
dotplot_T2D <- generate_taxa_per_time_dotplot_long(
  data.obj = subset_T2D.obj,
  test.list = result,
  t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
  ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
  group.var = "subject_race",
  time.var = "visit_number",
  feature.level = c("Genus", "Family")
)
result2 <- generate_taxa_change_per_time_test_long(
  data.obj = subset_T2D.obj,
  subject.var = "subject_id",
  time.var = "visit_number",
  t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
  ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
  group.var = "subject_race",
  adj.vars = NULL,
  prev.filter = 0.1,
  abund.filter = 0.001,
  feature.level = c("Genus", "Family"),
  feature.dat.type = "count"
)
}
}
