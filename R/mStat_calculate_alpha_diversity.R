#' Calculate alpha diversity indices
#'
#' This function calculates several alpha diversity indices (Shannon, Simpson, Observed Species, Chao1, ACE, Pielou, and Faith's PD) using the vegan and picante packages. The function takes a feature table (x) as input and returns a list containing the requested alpha diversity indices.
#' @name mStat_calculate_alpha_diversity
#' @param x Feature table (e.g., OTU/ASV table) with taxa in rows and samples in columns.
#' @param alpha.name character vector containing the names of alpha diversity indices to calculate. Possible values are: "shannon", "simpson", "observed_species", "chao1", "ace", "pielou", and "faith_pd".
#' @param tree Phylogenetic tree object of class "phylo". Required for Faith's phylogenetic diversity ("faith_pd") calculation.
#' @return A list containing the requested alpha diversity indices.
#' @examples
#' \dontrun{
#' # Create example feature table
#' otu.tab <- matrix(data = rpois(100, 5), nrow = 10, ncol = 10)
#' rownames(otu.tab) <- paste0("Taxon_", 1:10)
#' colnames(otu.tab) <- paste0("Sample_", 1:10)
#'
#' # Calculate alpha diversity indices
#' alpha.obj <- mStat_calculate_alpha_diversity(x = otu.tab, alpha.name =
#' c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou", "faith_pd"), tree = tree)
#' }
#' @export
mStat_calculate_alpha_diversity <- function(x, alpha.name, tree = NULL) {

  # Check if alpha.name is NULL and return early if so
  if (is.null(alpha.name)){
    return()
  }
  
  # Check if Faith's PD is requested but tree is missing
  if ("faith_pd" %in% alpha.name && is.null(tree)) {
    stop("Phylogenetic tree is required for Faith's phylogenetic diversity calculation. Please provide a tree object.")
  }

  # Check if the data has been rarefied by comparing the column sums
  # Rarefaction ensures equal sampling depth across all samples
  if (length(unique(colSums(x))) != 1) {
    warning("It appears the data may not have been rarefied. Please verify.")
  }

  # Transpose the input matrix to have samples as rows and taxa as columns
  # This is the required format for vegan package functions
  x_transpose <- t(x)

  # Calculate alpha diversity indices using lapply for efficient iteration
  alpha.obj <- lapply(alpha.name, function(index) {
    # Inform the user about the current diversity index being calculated
    message(paste("Calculating", index, "diversity..."))

    # Use switch to select the appropriate diversity calculation based on the index
    result <- switch(index,
                     # Shannon diversity: measure of entropy, sensitive to rare species
                     shannon = vegan::diversity(x_transpose, index = "shannon"),
                     # Simpson diversity: measure of dominance, less sensitive to rare species
                     simpson = vegan::diversity(x_transpose, index = "simpson"),
                     # Observed species: simple count of non-zero species
                     observed_species = vegan::specnumber(x_transpose),
                     # Chao1: non-parametric estimator of species richness
                     chao1 = vegan::estimateR(x_transpose)[2, ],
                     # ACE: Abundance-based Coverage Estimator of species richness
                     ace = vegan::estimateR(x_transpose)[4, ],
                     # Pielou's evenness: Shannon diversity divided by log of species richness
                     # Handle edge cases: when specnumber <= 1, log is 0 or undefined
                     pielou = {
                       shannon <- vegan::diversity(x_transpose, index = "shannon")
                       spec_num <- vegan::specnumber(x_transpose)
                       log_spec <- log(spec_num, exp(1))
                       # Set pielou to NA when species count <= 1 (undefined evenness)
                       pielou_val <- shannon / log_spec
                       pielou_val[spec_num <= 1] <- NA_real_
                       pielou_val
                     },
                     # Faith's phylogenetic diversity: sum of phylogenetic branch lengths
                     faith_pd = {
                       if (!requireNamespace("picante", quietly = TRUE)) {
                         stop("Package 'picante' is required for Faith's phylogenetic diversity calculation. Please install it using: install.packages('picante')")
                       }
                       # Calculate Faith's PD using picante package
                       pd_result <- picante::pd(x_transpose, tree, include.root = TRUE)
                       pd_result$PD
                     }
    )

    # Create a tibble with the diversity index and sample names
    # Convert to a named vector (row names as sample names)
    tibble(!!index := result, sample = rownames(x_transpose)) %>% column_to_rownames("sample")
  })

  # Assign names to the list elements based on the input alpha.name
  names(alpha.obj) <- alpha.name

  # Inform the user that all diversity calculations are complete
  message("Diversity calculations complete.")
  
  # Return the list of calculated alpha diversity indices
  return(alpha.obj)
}
