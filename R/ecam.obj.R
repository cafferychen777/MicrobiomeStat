#' Data from the Early Childhood Antibiotics and the Microbiome (ECAM) study
#'
#' The ecam.obj dataset is derived from a subset of the Early Childhood
#' Antibiotics and the Microbiome (ECAM) study, which tracked the microbiome
#' composition and development of 43 infants in the United States from birth
#' to 2 years of age, identifying associations between the microbiome and
#' antibiotic exposure, delivery mode, and diet.
#'
#' ecam.obj is a MicrobiomeStat Data Object.
#'
#' The meta.dat component of ecam.obj is a data frame with 875 observations
#' and 14 variables:
#' \itemize{
#'   \item X.SampleID: Character. Unique identifier for each sample.
#'   \item antiexposedall: Character. Indicates antibiotic exposure ('y' or 'n').
#'   \item day_of_life: Character. Day of life when the sample was collected.
#'   \item delivery: Character. Mode of delivery (e.g., 'Vaginal').
#'   \item diet: Character. Diet type (e.g., 'bd' for breast-fed).
#'   \item diet_3: Character. Categorized diet type (e.g., 'eb' for exclusively breast-fed).
#'   \item mom_child: Character. Indicates whether the sample is from mother or child ('C' for child).
#'   \item month: Character. Month when the sample was collected.
#'   \item month_of_life: Character. Precise month of life when the sample was collected.
#'   \item sample_summary: Character. Summary of sample characteristics.
#'   \item sex: Character. Sex of the infant.
#'   \item studyid: Character. Study identifier.
#'   \item subject.id: Character. Unique identifier for each subject.
#'   \item month_num: Numeric. Month number.
#' }
#'
#' @references
#' Bokulich, Nicholas A., et al. "Antibiotics, birth mode, and diet shape
#' microbiome maturation during early life." Science translational medicine
#' 8.343 (2016): 343ra82-343ra82.
#'
#' @source
#' Data source: https://github.com/FrederickHuangLin/ANCOM/tree/master/data
#'
"ecam.obj"
