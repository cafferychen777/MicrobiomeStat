#' Subset T2D Microbiome Data from the HMP2 Study
#'
#' The subset_T2D.obj dataset contains 16S rRNA data for 79 patients with type 2 diabetes (T2D) 
#' over multiple visits, gathered as part of the Integrative Human Microbiome Project (iHMP) or HMP2 study. 
#' The dataset includes a mStat object with 12,062 taxa and 2,208 samples.
#'
#' @format A MicrobiomeStat Data Object with the following components:
#' \describe{
#'   \item{feature.tab}{A matrix of microbial abundances}
#'   \item{meta.dat}{A data frame with 575 observations and 14 variables:
#'     \describe{
#'       \item{file_id}{Character. Unique identifier for each file}
#'       \item{md5}{Character. MD5 hash of the file}
#'       \item{size}{Character. File size in bytes}
#'       \item{urls}{Character. URL for file download}
#'       \item{sample_id}{Character. Unique identifier for each sample}
#'       \item{file_name}{Character. Name of the file}
#'       \item{subject_id}{Character. Unique identifier for each subject}
#'       \item{sample_body_site}{Character. Body site of sample collection (e.g., "feces")}
#'       \item{visit_number}{Character. Visit number as a string}
#'       \item{subject_gender}{Character. Gender of the subject}
#'       \item{subject_race}{Character. Race/ethnicity of the subject}
#'       \item{study_full_name}{Character. Full name of the study (e.g., "prediabetes")}
#'       \item{project_name}{Character. Name of the project}
#'       \item{visit_number_num}{Numeric. Visit number as a numeric value}
#'     }
#'   }
#'   \item{tree}{Phylogenetic tree of the taxa}
#'   \item{tax.tab}{Taxonomic classification of the taxa}
#' }
#'
#' @source Data source: Integrative Human Microbiome Project (iHMP)
#'   \url{https://www.hmpdacc.org/ihmp/}
#'
#' @references
#' Zhou, W., Sailani, M.R., Contrepois, K. et al. Longitudinal multi-omics of 
#' host–microbe dynamics in prediabetes. Nature 569, 663–671 (2019). 
#' https://doi.org/10.1038/s41586-019-1236-x
#'
#' @usage
#' data(subset_T2D.obj)
#'
#' @examples
#' data(subset_T2D.obj)
#' str(subset_T2D.obj$meta.dat)
#'
"subset_T2D.obj"
