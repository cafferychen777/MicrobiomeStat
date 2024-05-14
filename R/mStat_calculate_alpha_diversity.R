#' Calculate alpha diversity indices
#'
#' This function calculates several alpha diversity indices (Shannon, Simpson, Observed Species, Chao1, ACE, and Pielou) using the vegan package. The function takes an OTU table (x) as input and returns a list containing the requested alpha diversity indices.
#' @name mStat_calculate_alpha_diversity
#' @param x OTU table with taxa in rows and samples in columns.
#' @param alpha.name character vector containing the names of alpha diversity indices to calculate. Possible values are: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @return A list containing the requested alpha diversity indices.
#' @examples
#' \dontrun{
#' # Create example OTU table
#' otu.tab <- matrix(data = rpois(100, 5), nrow = 10, ncol = 10)
#' rownames(otu.tab) <- paste0("Taxon_", 1:10)
#' colnames(otu.tab) <- paste0("Sample_", 1:10)
#'
#' # Calculate alpha diversity indices
#' alpha.obj <- mStat_calculate_alpha_diversity(x = otu.tab, alpha.name =
#' c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"))
#' }
#' @export
mStat_calculate_alpha_diversity <- function(x, alpha.name) {

  if (is.null(alpha.name)){
    return()
  }

  # Check if the data has been rarefied
  if (length(unique(colSums(x))) != 1) {
    warning("It appears the data may not have been rarefied. Please verify.")
  }

  x_transpose <- t(x)

    alpha.obj <- lapply(alpha.name, function(index) {
      message(paste("Calculating", index, "diversity..."))

      result <- switch(index,
                       shannon = vegan::diversity(x_transpose, index = "shannon"),
                       simpson = vegan::diversity(x_transpose, index = "simpson"),
                       observed_species = vegan::specnumber(x_transpose),
                       chao1 = vegan::estimateR(x_transpose)[2, ],
                       ace = vegan::estimateR(x_transpose)[4, ],
                       pielou = vegan::diversity(x_transpose, index = "shannon") / log(vegan::specnumber(x_transpose), exp(1))
      )

      tibble(!!index := result, sample = rownames(x_transpose)) %>% column_to_rownames("sample")
    })

    names(alpha.obj) <- alpha.name

  message("Diversity calculations complete.")
  return(alpha.obj)
}
