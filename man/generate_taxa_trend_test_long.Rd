% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_trend_test_long.R
\name{generate_taxa_trend_test_long}
\alias{generate_taxa_trend_test_long}
\title{Longitudinal Taxa Trend Test Generation}
\usage{
generate_taxa_trend_test_long(
  data.obj,
  subject.var,
  time.var = NULL,
  group.var = NULL,
  adj.vars = NULL,
  feature.level,
  prev.filter = 0,
  abund.filter = 0,
  feature.dat.type = c("count", "proportion", "other"),
  ...
)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.}

\item{subject.var}{A character string that indicates the column name in the metadata which uniquely identifies each subject or sample.}

\item{time.var}{A character string representing the time variable column in the metadata. Time points should be numeric. If not, the function will convert it to numeric. Default is NULL.}

\item{group.var}{A character string specifying the grouping variable column in the metadata. This variable differentiates between different experimental or observational groups.}

\item{adj.vars}{A vector of character strings. Each string should denote a column name in the metadata that will serve as a covariate in the analysis. These variables might account for potential confounding influences. Default is NULL.}

\item{feature.level}{A character string indicating the taxonomic resolution for analysis (e.g., "Phylum", "Class"). This choice will determine the granularity of the analysis.}

\item{prev.filter}{Numeric value specifying the minimum prevalence threshold for filtering
taxa before analysis. Taxa with prevalence below this value will be removed.
Prevalence is calculated as the proportion of samples where the taxon is present.}

\item{abund.filter}{Numeric value specifying the minimum abundance threshold for filtering
taxa before analysis. Taxa with mean abundance below this value will be removed.
Abundance refers to counts or proportions depending on \code{feature.dat.type}.}

\item{feature.dat.type}{A character string, either "count", "proportion", or "other", indicating the nature of the data in the `data.obj`. 
This helps the function to determine if normalization is required. 
- "count": Raw count data that will be automatically normalized using TSS.
- "proportion": Pre-normalized data (e.g., relative abundance).
- "other": Custom pre-processed data. Use with caution; requires appropriate pre-processing for compositional data analysis.
Default is "count".}

\item{...}{Additional arguments to cater to any specialized requirements. For now, these are placeholder and not used.}
}
\value{
A list of dataframes, with each dataframe representing a specific taxonomic level (as specified in `feature.level`). These dataframes contain essential statistics, including taxa changes, p-values, and other metrics derived from the linear model.
}
\description{
This function is designed to conduct a longitudinal trend test on microbiome data. The primary aim is to discern how the abundance of various microbial taxa changes over time and/or in response to different experimental or observational groups. The function delivers robust statistical insights that enable researchers to draw meaningful conclusions about the dynamics of microbial populations.
}
\details{
Based on whether group.var, adj.vars, and time.var are NULL, the formula tests:

- When time.var is NULL:
  - When group.var is NULL and adj.vars is NOT NULL:
    - Tests adj.vars main effects only.
    - Adjusted for adj.vars but not for group.var or time.var.
  - When group.var is NOT NULL and adj.vars is NOT NULL:
    - Tests adj.vars and group.var main effects.
    - Adjusted for adj.vars but not for time.var.
  - When group.var is NOT NULL and adj.vars is NULL:
    - Tests group.var main effects only.
    - Unadjusted for adj.vars but adjusted for group.var.
  - When both group.var and adj.vars are NULL:
    - Tests the intercept only.
    - Unadjusted analysis.

- When time.var is NOT NULL:
  - When group.var is NULL and adj.vars is NOT NULL:
    - Tests adj.vars and time.var main effects.
    - Adjusted for adj.vars but not for group.var.
  - When group.var is NOT NULL and adj.vars is NOT NULL:
    - Tests adj.vars main effects.
    - Tests group.var and time.var main effects.
    - Tests group.var x time.var interaction.
    - Adjusted for adj.vars.
  - When group.var is NOT NULL and adj.vars is NULL:
    - Tests group.var and time.var main effects.
    - Tests group.var x time.var interaction.
    - Unadjusted analysis.
  - When both group.var and adj.vars are NULL:
    - Tests time.var main effect only.
    - Unadjusted analysis.

The formula combines the appropriate terms based on which variables are NULL.
Subject variability is accounted for through random effects.

When group.var = NULL and adj.vars = NULL and time.var is NOT NULL,
the slope of time.var is tested without adjusting for any additional covariates.
}
\examples{
\dontrun{
# Example 1
data("ecam.obj")
generate_taxa_trend_test_long(
  data.obj = ecam.obj,
  subject.var = "studyid",
  time.var = "month_num",
  group.var = "delivery",
  adj.vars = "diet",
  feature.level = c("Phylum","Class"),
  feature.dat.type = c("proportion")
)
generate_taxa_trend_test_long(
  data.obj = ecam.obj,
  subject.var = "studyid",
  time.var = "month_num",
  group.var = "delivery",
  feature.level = c("Phylum","Class"),
  feature.dat.type = c("proportion")
)
generate_taxa_trend_test_long(
  data.obj = ecam.obj,
  subject.var = "studyid",
  time.var = "month_num",
  group.var = NULL,
  feature.level = c("Phylum","Class"),
  feature.dat.type = c("proportion")
)

# Example 2
data("subset_T2D.obj")
test.list <- generate_taxa_trend_test_long(
  data.obj = subset_T2D.obj,
  subject.var = "subject_id",
  time.var = "visit_number",
  group.var = "subject_race",
  adj.vars = "sample_body_site",
  prev.filter = 0.1,
  abund.filter = 0.001,
  feature.level = c("Genus","Family"),
  feature.dat.type = c("count")
)

plot.list <- generate_taxa_trend_volcano_long(
  data.obj = subset_T2D.obj,
  group.var = "subject_race",
  time.var = "visit_number_num",
  test.list = test.list,
  feature.mt.method = "none")

}
}
