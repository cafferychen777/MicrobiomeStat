#' Probiotics Intervention Data from the PeerJ32 Study
#'
#' The peerj32.obj dataset contains high-throughput profiling data from 389 human blood serum lipids 
#' and 130 intestinal genus-level bacteria from 44 samples (22 subjects from 2 time points; 
#' before and after probiotic/placebo intervention). The data set is designed to investigate 
#' associations between intestinal bacteria and host lipid metabolism.
#'
#' @format A MicrobiomeStat Data Object containing the following components:
#' \describe{
#'   \item{feature.tab}{A matrix of microbial abundances (130 intestinal genus-level bacteria)}
#'   \item{meta.dat}{A data frame with 44 observations and 4 variables:
#'     \describe{
#'       \item{time}{Character. Time point of sample collection ("1" for before, "2" for after intervention)}
#'       \item{sex}{Character. Sex of the subject}
#'       \item{subject}{Character. Unique identifier for each subject (e.g., "S1", "S2", etc.)}
#'       \item{group}{Character. Intervention group ("Placebo" or "Probiotic")}
#'     }
#'   }
#'   \item{feature.ann}{Annotation data for the 389 human blood serum lipids}
#'   \item{feature.tab.lipids}{A matrix of lipid abundances (389 human blood serum lipids)}
#' }
#'
#' @usage
#' data(peerj32.obj)
#'
#' @author
#' Leo Lahti <microbiome-admin@googlegroups.com>
#'
#' @references
#' Lahti L, Salonen A, Kekkonen RA, Salojärvi J, Jalanka-Tuovinen J, Palva A, Orešič M, de Vos WM. (2013) 
#' Associations between the human intestinal microbiota, Lactobacillus rhamnosus GG and serum lipids 
#' indicated by integrated analysis of high-throughput profiling data. 
#' PeerJ 1:e32 https://doi.org/10.7717/peerj.32
#'
#' @source
#' Data source: https://doi.org/10.7717/peerj.32
#'
#' @examples
#' data(peerj32.obj)
#' str(peerj32.obj$meta.dat)
#'
"peerj32.obj"
