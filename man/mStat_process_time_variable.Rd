% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mStat_process_time_variable.R
\name{mStat_process_time_variable}
\alias{mStat_process_time_variable}
\title{Handle Time Variable in Meta Data Table}
\usage{
mStat_process_time_variable(
  data.obj,
  time.var,
  t0.level = NULL,
  ts.levels = NULL
)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
table (`feature.tab`), feature names, a full feature name list, and a feature aggregation list
(`feature.agg.list`).}

\item{time.var}{A character string that specifies the name of the time variable in the meta data table.}

\item{t0.level}{Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.}

\item{ts.levels}{Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.}
}
\value{
A data frame that is the modified version of the input meta data table.
}
\description{
This function handles a time variable in a meta data table. It sets the levels of the time variable based on provided parameters and subsets the data if necessary.
}
\details{
The function first gets the unique values of the time variable. Then, it sets the levels of the time variable based on the provided t0.level and ts.levels. If these parameters are not provided, the levels are set to the unique values of the time variable. The function then modifies the time variable in the meta data table based on its type (factor, numeric, or character). If the time variable is a character, it is converted to numeric. If there are any NA values generated during this conversion, a message is displayed. Finally, if the levels set do not include all the unique values of the time variable, the function subsets the data to exclude the missing levels and displays a message.
}
\examples{
# Assuming 'meta_data' is your meta data table and 'time' is your time variable
# meta_data_handled <- mStat_process_time_variable(meta_data, 'time')
data(ecam.obj)
mStat_process_time_variable(data.obj = ecam.obj,time.var = "month",
t0.level = "0",ts.levels = c("1","2"))

}
\seealso{
\code{\link[dplyr]{mutate}}, \code{\link[rlang]{sym}}
}
\author{
Your Name
}
