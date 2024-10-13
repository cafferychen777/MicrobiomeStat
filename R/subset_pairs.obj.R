#' Paired Sample 16S Sequencing Data from the MiPair Study
#'
#' The subset_pairs.obj dataset contains paired sample 16S sequencing data from the MiPair study. This data set is designed for design-based comparative analysis with paired microbiome data, focusing on the effects of antibiotic treatment on mouse gut microbiome.
#'
#' @format
#' A MicrobiomeStat Data Object containing the following data matrices as described in detail by Jang H, Koh H, Gu W, Kang B. (2022):
#' \describe{
#'   \item{feature.tab}{A matrix of microbial abundances, where rows represent taxa and columns represent samples.}
#'   \item{meta.dat}{A data frame with 80 rows and 3 variables:
#'     \describe{
#'       \item{Antibiotic}{Factor with 2 levels: "Baseline" and "Week 2", indicating the time point of sample collection relative to antibiotic treatment.}
#'       \item{MouseID}{Factor with unique identifiers for each mouse, ranging from 1 to 224.}
#'       \item{Sex}{Factor with 2 levels: "M" (male) and "F" (female), indicating the sex of each mouse.}
#'     }
#'   }
#'   \item{feature.ann}{A data frame containing taxonomic annotations for the microbial features.}
#'   \item{tree}{A phylogenetic tree object representing the evolutionary relationships among the microbial taxa.}
#' }
#'
#' The dataset includes samples from 80 mice, with paired measurements taken at baseline and 2 weeks after antibiotic treatment. There are 40 unique mice, each sampled at two time points.
#'
#' @usage
#' data(subset_pairs.obj)
#'
#' @author
#' Jang H, Koh H, Gu W, Kang B
#'
#' @references
#' Jang H, Koh H, Gu W, Kang B. (2022) Integrative web cloud computing and analytics using MiPair for design-based comparative analysis with paired microbiome data. Scientific Reports 12(20465).
#'
#' @source
#' Data source: https://github.com/YJ7599/MiPairGit
#'
#' @examples
#' data(subset_pairs.obj)
#' # View the first few rows of the metadata
#' head(subset_pairs.obj$meta.dat)
#' 
#' # Check the dimensions of the feature table
#' dim(subset_pairs.obj$feature.tab)
#' 
#' # Summarize the sex distribution
#' table(subset_pairs.obj$meta.dat$Sex)
#' 
#' # Count the number of samples at each time point
#' table(subset_pairs.obj$meta.dat$Antibiotic)
"subset_pairs.obj"
