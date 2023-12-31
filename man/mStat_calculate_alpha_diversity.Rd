% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mStat_calculate_alpha_diversity.R
\name{mStat_calculate_alpha_diversity}
\alias{mStat_calculate_alpha_diversity}
\title{Calculate alpha diversity indices}
\usage{
mStat_calculate_alpha_diversity(x, alpha.name)
}
\arguments{
\item{x}{OTU table with taxa in rows and samples in columns.}

\item{alpha.name}{character vector containing the names of alpha diversity indices to calculate. Possible values are: "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".}
}
\value{
A list containing the requested alpha diversity indices.
}
\description{
This function calculates several alpha diversity indices (Shannon, Simpson, Observed Species, Chao1, ACE, and Pielou) using the vegan package. The function takes an OTU table (x) as input and returns a list containing the requested alpha diversity indices.
}
\examples{
\dontrun{
# Create example OTU table
otu.tab <- matrix(data = rpois(100, 5), nrow = 10, ncol = 10)
rownames(otu.tab) <- paste0("Taxon_", 1:10)
colnames(otu.tab) <- paste0("Sample_", 1:10)

# Calculate alpha diversity indices
alpha.obj <- mStat_calculate_alpha_diversity(x = otu.tab, alpha.name =
c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"))
}
}
