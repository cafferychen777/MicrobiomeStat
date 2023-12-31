% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_alpha_test_pair.R
\name{is_rarefied}
\alias{is_rarefied}
\title{Check if the data has been rarefied}
\usage{
is_rarefied(data.obj)
}
\arguments{
\item{data.obj}{A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.}
}
\value{
A boolean value indicating whether the data is rarefied. It returns
  TRUE if the data is rarefied, FALSE otherwise.
}
\description{
This function checks if the data has been rarefied by inspecting if the sum
of each column (which represents each sample in the OTU table) is equal.
}
\examples{
# Assuming peerj32.obj is a data object with OTU and metadata tables
data(peerj32.obj)
is_rarefied(peerj32.obj)

}
