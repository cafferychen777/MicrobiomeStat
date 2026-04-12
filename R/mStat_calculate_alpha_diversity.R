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
#' # Calculate non-phylogenetic alpha diversity indices
#' alpha.obj <- mStat_calculate_alpha_diversity(
#'   x = otu.tab,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou")
#' )
#'
#' # Calculate Faith's PD when picante is available
#' if (requireNamespace("picante", quietly = TRUE)) {
#'   tree <- ape::rtree(n = nrow(otu.tab), tip.label = rownames(otu.tab))
#'   faith_pd <- mStat_calculate_alpha_diversity(
#'     x = otu.tab,
#'     alpha.name = "faith_pd",
#'     tree = tree
#'   )
#' }
#' }
#' @export
mStat_calculate_alpha_diversity <- function(x, alpha.name, tree = NULL) {
  if (is.null(alpha.name)){
    return()
  }

  valid_alpha <- c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou", "faith_pd")
  invalid_alpha <- setdiff(alpha.name, valid_alpha)
  if (length(invalid_alpha) != 0) {
    stop(
      "Unsupported alpha diversity indices: ",
      paste(invalid_alpha, collapse = ", "),
      ". Supported options are: ",
      paste(valid_alpha, collapse = ", ")
    )
  }

  if ("faith_pd" %in% alpha.name && is.null(tree)) {
    stop("Phylogenetic tree is required for Faith's phylogenetic diversity calculation. Please provide a tree object.")
  }

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
                     pielou = {
                       shannon <- vegan::diversity(x_transpose, index = "shannon")
                       spec_num <- vegan::specnumber(x_transpose)
                       log_spec <- log(spec_num, exp(1))
                       pielou_val <- shannon / log_spec
                       pielou_val[spec_num <= 1] <- NA_real_
                       pielou_val
                     },
                     faith_pd = {
                       if (!requireNamespace("picante", quietly = TRUE)) {
                         stop("Package 'picante' is required for Faith's phylogenetic diversity calculation. Please install it using: install.packages('picante')")
                       }
                       pd_result <- picante::pd(x_transpose, tree, include.root = TRUE)
                       pd_result$PD
                     }
    )

    tibble(!!index := result, sample = rownames(x_transpose)) %>% tibble::column_to_rownames("sample")
  })

  names(alpha.obj) <- alpha.name

  message("Diversity calculations complete.")
  return(alpha.obj)
}
