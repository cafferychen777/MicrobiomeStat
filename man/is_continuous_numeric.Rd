% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_taxa_indiv_change_scatterplot_pair.R
\name{is_continuous_numeric}
\alias{is_continuous_numeric}
\title{Check if a Variable is Continuous Numeric}
\usage{
is_continuous_numeric(x)
}
\arguments{
\item{x}{A variable that you want to check.}
}
\value{
Logical. Returns TRUE if the variable is continuous numeric, FALSE otherwise.
}
\description{
This function checks if a given variable is continuous numeric by checking if it is numeric and has at least 10 unique values.
}
\examples{
is_continuous_numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) # TRUE
is_continuous_numeric(c(1, 2, 3, 4, 5)) # FALSE
is_continuous_numeric(c('a', 'b', 'c')) # FALSE

}
