#' @title Parameter Documentation Templates
#'
#' @description
#' This file provides centralized parameter documentation for MicrobiomeStat.
#' Other functions should use `@inheritParams` to inherit these definitions,
#' ensuring consistency across the package.
#'
#' @details
#' ## Usage
#'
#' To inherit parameter documentation in another function:
#' ```
#' #' @inheritParams mStat_data_obj_doc
#' ```
#'
#' This will inherit all matching parameter names from the template function.
#'
#' ## Available Templates
#'
#' - `mStat_data_obj_doc`: Core data object and common analysis parameters
#' - `mStat_plot_params_doc`: Visualization-related parameters
#' - `mStat_test_params_doc`: Statistical testing parameters
#'
#' @name mStat-param-docs
#' @keywords internal
NULL


#' Core Data Object Parameter Documentation
#'
#' Template function defining the standard documentation for data.obj and
#' common analysis parameters used throughout MicrobiomeStat.
#'
#' @param data.obj A MicrobiomeStat data object, which is a list containing
#'   at minimum the following components:
#'   \itemize{
#'     \item \code{feature.tab}: A matrix of feature abundances (taxa/genes as
#'       rows, samples as columns)
#'     \item \code{meta.dat}: A data frame of sample metadata (samples as rows)
#'   }
#'   Optional components include:
#'   \itemize{
#'     \item \code{feature.ann}: A matrix/data frame of feature annotations
#'       (e.g., taxonomy)
#'     \item \code{tree}: A phylogenetic tree object (class "phylo")
#'     \item \code{feature.agg.list}: Pre-aggregated feature tables by taxonomy
#'   }
#'   Data objects can be created using converters like
#'   \code{\link{mStat_convert_phyloseq_to_data_obj}} or importers like
#'   \code{\link{mStat_import_qiime2_as_data_obj}}.
#'
#' @param subject.var Character string specifying the column name in meta.dat
#'   that uniquely identifies each subject or sample unit. Required for
#'   longitudinal and paired designs to track repeated measurements.
#'
#' @param time.var Character string specifying the column name in meta.dat
#'   containing the time variable. Required for longitudinal and paired analyses.
#'   Should be a factor or character with meaningful time point labels.
#'
#' @param group.var Character string specifying the column name in meta.dat
#'   containing the grouping variable (e.g., treatment, condition, phenotype).
#'   Used for between-group comparisons.
#'
#' @param strata.var Character string specifying the column name in meta.dat
#'   for stratification. When provided, analyses and visualizations will be
#'   performed separately within each stratum (e.g., by site, batch, or sex).
#'
#' @param adj.vars Character vector specifying column names in meta.dat to be
#'   used as covariates for adjustment in statistical models. These variables
#'   will be included as fixed effects.
#'
#' @param feature.level Character vector specifying the taxonomic or annotation
#'   level(s) for analysis. Should match column names in feature.ann, such as
#'   "Phylum", "Family", "Genus", etc. Use "original" to analyze at the
#'   original feature level without aggregation.
#'
#' @param feature.dat.type Character string specifying the data type of
#'   feature.tab. One of:
#'   \itemize{
#'     \item "count": Raw count data (will be normalized)
#'     \item "proportion": Relative abundance data (should sum to 1 per sample)
#'     \item "other": Pre-transformed data (no transformation applied)
#'   }
#'
#' @param prev.filter Numeric value between 0 and 1. Features with prevalence
#'   (proportion of non-zero samples) below this threshold will be excluded
#'   from analysis. Default is usually 0 (no filtering).
#'
#' @param abund.filter Numeric value. Features with mean abundance below this
#'   threshold will be excluded from analysis. Default is usually 0 (no filtering).
#'
#' @param t0.level Character or numeric value specifying the baseline time
#'   point for longitudinal or paired analyses. Should match a value in the
#'   time.var column.
#'
#' @param ts.levels Character vector specifying the follow-up time points for
#'   longitudinal or paired analyses. Should match values in the time.var column.
#'
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return NULL (documentation-only function)
#'
#' @keywords internal
mStat_data_obj_doc <- function(data.obj,
                                subject.var = NULL,
                                time.var = NULL,
                                group.var = NULL,
                                strata.var = NULL,
                                adj.vars = NULL,
                                feature.level = NULL,
                                feature.dat.type = c("count", "proportion", "other"),
                                prev.filter = 0,
                                abund.filter = 0,
                                t0.level = NULL,
                                ts.levels = NULL,
                                ...) {
 NULL
}


#' Plot Parameter Documentation
#'
#' Template function defining standard documentation for visualization parameters.
#'
#' @param base.size Numeric value specifying the base font size for plot text
#'   elements. Default is typically 16.
#'
#' @param theme.choice Character string specifying the ggplot2 theme to use.
#'   Options include:
#'   \itemize{
#'     \item "bw": Black and white theme (theme_bw)
#'     \item "classic": Classic theme (theme_classic)
#'     \item "minimal": Minimal theme (theme_minimal)
#'     \item "prism": GraphPad Prism-like theme
#'     \item "nature": Nature journal style
#'     \item "light": Light theme (theme_light)
#'   }
#'   Can also use a custom ggplot2 theme object via custom.theme.
#'
#' @param custom.theme A custom ggplot2 theme object to override theme.choice.
#'   Should be created using ggplot2::theme() or a complete theme function.
#'
#' @param palette Character vector of colors or a named palette for the plot.
#'   If NULL, uses default MicrobiomeStat color scheme. Can be:
#'   \itemize{
#'     \item A vector of color codes (e.g., c("#E41A1C", "#377EB8"))
#'     \item A palette name recognized by the plotting function
#'   }
#'
#' @param pdf Logical. If TRUE, saves the plot(s) to PDF file(s) in the
#'   current working directory. Default is TRUE.
#'
#' @param file.ann Character string for additional annotation to append to
#'   output filenames. Useful for distinguishing multiple outputs.
#'
#' @param pdf.wid Numeric value specifying the width of PDF output in inches.
#'   Default is typically 11.
#'
#' @param pdf.hei Numeric value specifying the height of PDF output in inches.
#'   Default is typically 8.5.
#'
#' @return NULL (documentation-only function)
#'
#' @keywords internal
mStat_plot_params_doc <- function(base.size = 16,
                                   theme.choice = "bw",
                                   custom.theme = NULL,
                                   palette = NULL,
                                   pdf = TRUE,
                                   file.ann = NULL,
                                   pdf.wid = 11,
                                   pdf.hei = 8.5) {
 NULL
}


#' Statistical Test Parameter Documentation
#'
#' Template function defining standard documentation for statistical testing parameters.
#'
#' @param alpha.obj A list containing pre-calculated alpha diversity indices.
#'   If NULL and alpha diversity is needed, it will be calculated automatically.
#'   Names should match the alpha.name parameter (e.g., "shannon", "simpson").
#'   See \code{\link{mStat_calculate_alpha_diversity}}.
#'
#' @param alpha.name Character vector specifying which alpha diversity indices
#'   to analyze. Options include:
#'   \itemize{
#'     \item "shannon": Shannon diversity index
#'     \item "simpson": Simpson diversity index
#'     \item "observed_species": Observed species richness
#'     \item "chao1": Chao1 richness estimator
#'     \item "ace": ACE richness estimator
#'     \item "pielou": Pielou's evenness
#'   }
#'
#' @param dist.obj A list of pre-calculated distance matrices. If NULL and
#'   distances are needed, they will be calculated automatically.
#'   List names should match dist.name (e.g., "BC" for Bray-Curtis).
#'   See \code{\link{mStat_calculate_beta_diversity}}.
#'
#' @param dist.name Character vector specifying which distance metrics to use.
#'   Options depend on available methods:
#'   \itemize{
#'     \item "BC": Bray-Curtis dissimilarity
#'     \item "Jaccard": Jaccard distance
#'     \item "UniFrac": Unweighted UniFrac (requires tree)
#'     \item "GUniFrac": Generalized UniFrac (requires tree)
#'     \item "WUniFrac": Weighted UniFrac (requires tree)
#'     \item "JS": Jensen-Shannon divergence
#'   }
#'
#' @param depth Numeric value or NULL. Rarefaction depth for diversity
#'   calculations. If NULL, uses minimum sample depth or no rarefaction.
#'
#' @return NULL (documentation-only function)
#'
#' @keywords internal
mStat_test_params_doc <- function(alpha.obj = NULL,
                                   alpha.name = c("shannon", "simpson", "observed_species"),
                                   dist.obj = NULL,
                                   dist.name = c("BC", "Jaccard"),
                                   depth = NULL) {
 NULL
}
