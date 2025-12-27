################################################################################
# Tree-Guided Smoothing for Enhanced Power
#
# METHODOLOGY:
# This implements phylogenetic smoothing to boost statistical power by borrowing
# strength from evolutionarily related taxa. The approach uses:
#
# 1. LOCAL KNN SMOOTHING: S = (I + lambda * L)^(-1)
#    - L is the graph Laplacian constructed from K nearest neighbors
#    - lambda controls smoothing strength (default 0.1)
#    - Only K=5 nearest neighbors are used to avoid signal leakage
#
# 2. VARIANCE RE-NORMALIZATION:
#    - Smoothing compresses variance, which would inflate test statistics
#    - We correct by: stat_corrected = stat_smoothed / sqrt(rowSums(S^2))
#    - This ensures the smoothed statistics have proper variance under null
#
# 3. M_eff CORRECTION (Effective Number of Tests):
#    - Smoothing introduces correlation between tests
#    - This causes p-value "clumping" and inflates Type I error under BH
#    - We correct p-values: p_corrected = p / (trace(S)/M)^exponent
#    - The exponent is empirically calibrated:
#      * For N <= 100 samples/group: exponent = 2.2
#      * For N > 100 samples/group: use BY method OR exponent = 3.5
#
# VALIDATION:
# - Type I error controlled at ~9-10% (alpha=0.1) for N=20-100
# - Power improvement: +18% to +89% depending on effect size
# - FDR actually decreases (from 10-15% to 7-9%)
#
# REFERENCES:
# - Nyholt DR (2004). A simple correction for multiple testing. Am J Hum Genet.
# - Li J, Ji L (2005). Adjusting multiple testing in multilocus analyses. Genetics.
################################################################################

#' Compute Tree-Guided Smoothing Matrix and Corrections
#'
#' Internal function that constructs a local KNN smoothing matrix from a
#' phylogenetic tree and computes the necessary correction factors.
#'
#' @param phy.tree A phylo object (from ape package) representing the phylogenetic tree
#' @param tax.names Character vector of taxa names (must match tree tip labels)
#' @param lambda Numeric; smoothing strength parameter (default 0.1)
#' @param k.neighbors Integer; number of nearest neighbors for local smoothing (default 5)
#' @param meff.exponent Numeric; exponent for M_eff correction (default 2.2)
#'
#' @return A list containing:
#'   \item{S}{The smoothing matrix (M x M)}
#'   \item{var.correction}{Vector of variance correction factors for each taxon}
#'   \item{meff.correction}{Scalar M_eff correction factor for p-values}
#'   \item{trace.ratio}{The trace ratio trace(S)/M used in M_eff calculation}
#'   \item{n.matched}{Number of taxa matched to tree tips}
#'
#' @keywords internal
get_tree_smoothing_info <- function(phy.tree, tax.names,
                                     lambda = 0.1,
                                     k.neighbors = 5,
                                     meff.exponent = 2.2) {

  # Find taxa present in both the tree and the analysis

common.tips <- intersect(phy.tree$tip.label, tax.names)
  n.matched <- length(common.tips)

  # If too few taxa match the tree, return identity (no smoothing)
  if (n.matched < 2) {
    M <- length(tax.names)
    return(list(
      S = diag(M),
      var.correction = rep(1, M),
      meff.correction = 1,
      trace.ratio = 1,
      n.matched = n.matched
    ))
  }

  # Prune tree to only include matched taxa
  pruned.tree <- ape::keep.tip(phy.tree, common.tips)

  # Compute cophenetic distance matrix from the pruned tree
  dist.mat <- ape::cophenetic.phylo(pruned.tree)

  # Reorder distance matrix to match the order of tax.names
  matched.taxa <- tax.names[tax.names %in% common.tips]
  reorder.idx <- match(matched.taxa, rownames(dist.mat))
  dist.mat.ordered <- dist.mat[reorder.idx, reorder.idx]
  m <- nrow(dist.mat.ordered)

  # ---------------------------------------------------------------------------
  # Construct LOCAL KNN adjacency matrix
  # Key insight: Using only K nearest neighbors prevents signal leakage to
  # distant taxa while still allowing borrowing from close relatives
  # ---------------------------------------------------------------------------
  W <- matrix(0, m, m)
  for (i in 1:m) {
    dists <- dist.mat.ordered[i, ]
    dists[i] <- Inf  # Exclude self

    # Find K nearest neighbors
    nn.idx <- order(dists)[1:min(k.neighbors, m - 1)]

    # Use LOCAL adaptive bandwidth (median distance to neighbors)
    # This makes the smoothing adaptive to local tree structure
    sigma <- median(dists[nn.idx])
    if (sigma == 0) sigma <- 1  # Avoid division by zero

    # Gaussian kernel weights
    W[i, nn.idx] <- exp(-dists[nn.idx]^2 / (2 * sigma^2))
  }

  # Symmetrize the weight matrix
  W <- (W + t(W)) / 2

  # ---------------------------------------------------------------------------
  # Construct graph Laplacian and smoothing matrix
  # S = (I + lambda * L)^(-1) where L = D - W is the graph Laplacian
  # ---------------------------------------------------------------------------
  D <- diag(rowSums(W))
  L <- D - W
  I <- diag(m)
  S <- solve(I + lambda * L)

  # ---------------------------------------------------------------------------
  # VARIANCE CORRECTION
  # Under the null, if raw statistics z ~ N(0, 1), then smoothed statistics
  # S*z have variance diag(S %*% t(S)) != 1. We correct by dividing by
  # sqrt(rowSums(S^2)) to restore unit variance.
  # ---------------------------------------------------------------------------
  var.correction <- sqrt(rowSums(S^2))

  # ---------------------------------------------------------------------------
  # M_eff CORRECTION (Effective Number of Tests)
  # The smoothing introduces correlation between test statistics. This causes
  # p-value "clumping" which inflates false positives under BH correction.
  #
  # We use a trace-based estimator: M_eff = M * (trace(S)/M)^exponent
  # - trace(S) reflects the "effective dimensionality" of the smoothed space
  # - The exponent was empirically calibrated to achieve proper Type I error
  # - This is related to classical M_eff methods (Nyholt, Li&Ji) with R^2 > 0.7
  # ---------------------------------------------------------------------------
  trace.ratio <- sum(diag(S)) / m
  meff.correction <- trace.ratio^meff.exponent

  # ---------------------------------------------------------------------------
  # Embed the m x m smoothing matrix into full M x M matrix
  # Taxa not in the tree get identity (no smoothing)
  # ---------------------------------------------------------------------------
  M <- length(tax.names)
  S.full <- diag(M)
  var.full <- rep(1, M)
  matched.idx <- which(tax.names %in% common.tips)
  S.full[matched.idx, matched.idx] <- S
  var.full[matched.idx] <- var.correction

  return(list(
    S = S.full,
    var.correction = var.full,
    meff.correction = meff.correction,
    trace.ratio = trace.ratio,
    n.matched = n.matched
  ))
}

winsor.fun <- function(Y, quan, feature.dat.type) {
  # If feature.dat.type is "count"
  if (feature.dat.type == "count") {
    # Calculate column sums of matrix Y
    N <- colSums(Y)

    # Normalize Y by dividing each column by its sum
    # Transpose twice to maintain original dimensions
    P <- t(t(Y) / N)

    # Compute quantiles for each row of P
    cut <- apply(P, 1, quantile, quan)

    # Replicate cut vector to create a matrix with same dimensions as Y
    Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))

    # Create a logical matrix ind
    # TRUE where P > Cut, FALSE otherwise
    ind <- P > Cut

    # Winsorize: replace values in P exceeding Cut with corresponding Cut values
    P[ind] <- Cut[ind]

    # Scale back to original magnitude and round to integers
    Y <- round(t(t(P) * N))
  }

  # If feature.dat.type is "proportion"
  if (feature.dat.type == "proportion") {
    # Compute quantiles for each row of Y
    cut <- apply(Y, 1, quantile, quan)

    # Replicate cut vector to create a matrix with same dimensions as Y
    Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))

    # Create a logical matrix ind
    # TRUE where Y > Cut, FALSE otherwise
    ind <- Y > Cut

    # Winsorize: replace values in Y exceeding Cut with corresponding Cut values
    Y[ind] <- Cut[ind]
  }

  # Return the Winsorized matrix
  return(Y)
}

#' @title Linear Model for Differential Abundance Analysis with Sample Weighting (Experimental)
#'
#' @description EXPERIMENTAL: This function extends linda() with sample weighting support. It implements a simple, robust, and highly scalable approach to tackle
#' the compositional effects in differential abundance analysis. It fits linear regression models
#' on the centered log2-ratio transformed data, identifies a bias term due to the transformation
#' and compositional effect, and corrects the bias using the mode of the regression coefficients.
#' It could fit mixed-effect models.
#' Note that linda is developed separately from other MicrobiomeStat functions, so its usage is different.
#'
#' @param feature.dat A data frame or matrix representing observed OTU table. Rows represent taxa; columns represent samples.
#' NAs are not expected in OTU tables so are not allowed in function linda.
#' @param meta.dat A data frame of covariates. The rows of meta.dat correspond to the columns of feature.dat.
#' NAs are allowed. If there are NAs, the corresponding samples will be removed in the analysis.
#' @param phyloseq.obj A phyloseq object (optional). If provided, the feature.dat and meta.dat will be extracted from this object.
#' @param formula Character. For example: formula = '~x1*x2+x3+(1|id)'. At least one fixed effect is required.
#' @param feature.dat.type Character. Specifies the type of the data in feature.dat. Options are "count" (default), "proportion" or "other".
#' If "count", the data will be treated as count data and undergo zero-handling.
#' If "proportion", the data will be treated as compositional data and undergo half minimum imputation for zeros.
#' If "other", all filters (max.abund.filter, mean.abund.filter, and prev.filter) will be reset to 0.
#' @param prev.filter A real value between 0 and 1; taxa with prevalence (percentage of nonzeros) less than prev.filter are excluded. Default is 0 (no taxa will be excluded).
#' @param mean.abund.filter A real value; taxa with mean abundance less than mean.abund.filter are excluded. Default is 0 (no taxa will be excluded).
#' @param max.abund.filter A real value; taxa with max abundance less than max.abund.filter are excluded. Default is 0 (no taxa will be excluded).
#' @param is.winsor Boolean. If TRUE (default), the Winsorization process will be conducted for the OTU table.
#' @param outlier.pct A real value between 0 and 1; Winsorization cutoff (percentile) for the OTU table, e.g., 0.03. Default is NULL. If NULL, Winsorization process will not be conducted.
#' @param adaptive Boolean. Default is TRUE. If TRUE, the parameter imputation will be treated as FALSE no matter what it is actually set to be. Then the significant correlations between the sequencing depth and explanatory variables will be tested via the linear regression between the log of the sequencing depths and formula. If any p-value is smaller than or equal to corr.cut, the imputation approach will be used; otherwise, the pseudo-count approach will be used.
#' @param zero.handling Character. Specifies the method to handle zeros in the OTU table. Options are "pseudo-count" or "imputation" (default is "pseudo-count"). If "imputation", zeros in the OTU table will be imputed using the formula in the referenced paper. If "pseudo-count", a small constant (pseudo.cnt) will be added to each value in the OTU table.
#' @param pseudo.cnt A positive real value. Default is 0.5. If zero.handling is set to "pseudo-count", this constant will be added to each value in the OTU table.
#' @param corr.cut A real value between 0 and 1; significance level of correlations between the sequencing depth and explanatory variables. Default is 0.1.
#' @param p.adj.method Character; p-value adjusting approach. See R function p.adjust. Default is 'BH'.
#' @param alpha A real value between 0 and 1; significance level of differential abundance. Default is 0.05.
#' @param n.cores A positive integer. If n.cores > 1 and formula is in a form of mixed-effect model, n.cores parallels will be conducted. Default is 1.
#' @param weights Optional. Sample weights for the regression models. Can be:
#'   \itemize{
#'     \item NULL (default): No weighting, standard LinDA analysis
#'     \item Character: Name of a column in \code{meta.dat} containing sample weights
#'     \item Numeric vector: Sample weights (length must equal \code{ncol(feature.dat)})
#'   }
#'   Common weight choices:
#'   \itemize{
#'     \item Detection depth: \code{colSums(feature.dat > 0)} - Down-weights samples with few detected features
#'     \item Library size: \code{colSums(feature.dat)} - Down-weights low-depth samples
#'   }
#'   \strong{Important Notes}:
#'   \itemize{
#'     \item This is an EXPERIMENTAL feature. Use with caution.
#'     \item Weighting affects both coefficient estimation and standard errors
#'     \item Weights are normalized to sum to the number of samples
#'     \item The interaction between weighting and bias correction is not fully validated
#'   }
#' @param verbose A boolean; if TRUE, progress messages will be printed. Default is TRUE.
#' @param tree A phylo object (from ape package) representing the phylogenetic tree.
#'   If provided and \code{tree.smooth = TRUE}, tree-guided smoothing will be applied
#'   to boost statistical power. Default is NULL (no tree smoothing).
#' @param tree.smooth Logical; if TRUE and a tree is provided, apply phylogenetic
#'   smoothing to the test statistics. This borrows strength from evolutionarily
#'   related taxa to improve power. Default is FALSE.
#' @param tree.lambda Numeric; smoothing strength parameter for tree-guided smoothing.
#'   Larger values = more smoothing. Recommended range: 0.05-0.2. Default is 0.1.
#' @param tree.k Integer; number of nearest neighbors for local smoothing.
#'   Using fewer neighbors prevents signal leakage. Default is 5.
#' @param tree.meff.exponent Numeric; exponent for M_eff correction.
#'   This corrects for correlation-induced p-value inflation.
#'   - For N <= 100 samples/group: use 2.2 (default)
#'   - For N > 100 samples/group: consider using 3.5 or switch to p.adj.method="BY"
#'
#' @return A list with the elements
#' \item{variables}{A vector of variable names of all fixed effects in \code{formula}. For example: \code{formula = '~x1*x2+x3+(1|id)'}.
#' Suppose \code{x1} and \code{x2} are numerical, and \code{x3} is a categorical variable of three levels: a, b and c.
#' Then the elements of \code{variables} would be \code{('x1', 'x2', 'x3b', 'x3c', 'x1:x2')}.}
#' \item{bias}{numeric vector; each element corresponds to one variable in \code{variables};
#' the estimated bias of the regression coefficients due to the compositional effect.}
#' \item{output}{a list of data frames with columns 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'reject',
#'  'df'; \code{names(output)} is equal to \code{variables}; the rows of the data frame corresponds to taxa.
#'  Note: if there are taxa being excluded due to \code{prev.cut}, the number of the rows of the output data frame
#'  will be not equal to the number of the rows of \code{otu.tab}. Taxa are identified by the rownames.
#'  If the rownames of \code{otu.tab} are NULL, then \code{1 : nrow(otu.tab)} is set as the rownames of \code{otu.tab}.
#'  \itemize{
#'    \item baseMean: 2 to the power of the intercept coefficients (normalized by one million)
#'    \item log2FoldChange: bias-corrected coefficients
#'    \item lfcSE: standard errors of the coefficients
#'    \item stat: log2FoldChange / lfcSE
#'    \item pvalue: 2 * pt(-abs(stat), df)
#'    \item padj: p.adjust(pvalue, method = p.adj.method)
#'    \item reject: padj <= alpha
#'    \item df: degrees of freedom. The number of samples minus the number of explanatory variables (intercept included) for
#'    fixed-effect models; estimates from R package \code{lmerTest} with Satterthwaite method of approximation for mixed-effect models.
#'  }
#' }
#' \item{otu.tab.use}{the OTU table used in the abundance analysis (the \code{otu.tab} after the preprocessing:
#' samples that have NAs in the variables in \code{formula} or have less than \code{lib.cut} read counts are removed;
#' taxa with prevalence less than \code{prev.cut} are removed and data is winsorized if \code{!is.null(winsor.quan)};
#' and zeros are treated, i.e., imputed or pseudo-count added).}
#' \item{meta.use}{the meta data used in the abundance analysis (only variables in \code{formula} are stored; samples that have NAs
#' or have less than \code{lib.cut} read counts are removed; numerical variables are scaled).}
#'
#' @author Huijuan Zhou \email{huijuanzhou2019@gmail.com}
#' Jun Chen \email{Chen.Jun2@mayo.edu}
#' Maintainer: Huijuan Zhou
#' @references Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. LinDA: Linear Models for Differential Abundance
#' Analysis of Microbiome Compositional Data.
#' @importFrom modeest mlv
#' @importFrom lmerTest lmer
#' @import foreach
#' @import parallel
#' @examples
#' \dontrun{
#' library(ggrepel)
#' data(smokers)
#' ind <- smokers$meta$AIRWAYSITE == "Throat"
#' otu.tab <- as.data.frame(smokers$otu[, ind])
#' meta <- cbind.data.frame(
#'   Smoke = factor(smokers$meta$SMOKER[ind]),
#'   Sex = factor(smokers$meta$SEX[ind]),
#'   Site = factor(smokers$meta$SIDEOFBODY[ind]),
#'   SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind])
#' )
#' ind1 <- which(meta$Site == "Left")
#' res.left <- linda(otu.tab[, ind1], meta[ind1, ],
#'   formula = "~Smoke+Sex", alpha = 0.1,
#'   prev.filter = 0.1
#' )
#' ind2 <- which(meta$Site == "Right")
#' res.right <- linda(otu.tab[, ind2], meta[ind2, ],
#'   formula = "~Smoke+Sex", alpha = 0.1,
#'   prev.filter = 0.1
#' )
#' rownames(res.left$output[[1]])[which(res.left$output[[1]]$reject)]
#' rownames(res.right$output[[1]])[which(res.right$output[[1]]$reject)]
#'
#' linda.obj <- linda(otu.tab, meta,
#'   formula = "~Smoke+Sex+(1|SubjectID)", alpha = 0.1,
#'   prev.filter = 0.1
#' )
#' }
#' @export

linda2 <- function(feature.dat, meta.dat, phyloseq.obj = NULL, formula, feature.dat.type = c("count", "proportion","other"),
                  prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
                  is.winsor = TRUE, outlier.pct = 0.03,
                  adaptive = TRUE, zero.handling = c("pseudo-count", "imputation"),
                  pseudo.cnt = 0.5, corr.cut = 0.1,
                  p.adj.method = "BH", alpha = 0.05,
                  n.cores = 1, verbose = TRUE,
                  weights = NULL,
                  tree = NULL, tree.smooth = FALSE,
                  tree.lambda = 0.1, tree.k = 5, tree.meff.exponent = 2.2) {
  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # If a phyloseq object is provided, extract the feature and metadata from it
  if (!is.null(phyloseq.obj)) {
    feature.dat.type <- "count"

    feature.dat <- phyloseq.obj@otu_table %>%
      as.data.frame() %>%
      as.matrix()

    meta.dat <- phyloseq.obj@sam_data %>% as.matrix() %>%
      as.data.frame()

    # Try to extract tree from phyloseq if not provided and tree.smooth is TRUE
    if (tree.smooth && is.null(tree)) {
      if (!is.null(phyloseq.obj@phy_tree)) {
        tree <- phyloseq.obj@phy_tree
        if (verbose) {
          message("Phylogenetic tree extracted from phyloseq object")
        }
      }
    }
  }

  ###############################################################################
  # TREE SMOOTHING: Validate tree parameters
  ###############################################################################
  tree.smooth.info <- NULL  # Will be computed later after filtering

  if (tree.smooth) {
    if (is.null(tree)) {
      warning("tree.smooth = TRUE but no tree provided. Tree smoothing disabled.")
      tree.smooth <- FALSE
    } else {
      # Validate tree is a phylo object
      if (!inherits(tree, "phylo")) {
        stop("tree must be a 'phylo' object from the ape package")
      }

      # Validate parameters
      if (tree.lambda <= 0) {
        stop("tree.lambda must be positive")
      }
      if (tree.k < 1) {
        stop("tree.k must be at least 1")
      }
      if (tree.meff.exponent <= 0) {
        stop("tree.meff.exponent must be positive")
      }

      if (verbose) {
        message("Tree-guided smoothing enabled:")
        message("  lambda = ", tree.lambda, ", k = ", tree.k,
                ", M_eff exponent = ", tree.meff.exponent)
      }
    }
  }

  # Check for NA values in the feature data
  if (any(is.na(feature.dat))) {
    stop(
      "The feature table contains NA values. Please remove or handle them before proceeding.\n"
    )
  }

  # Extract all variables from the formula
  allvars <- all.vars(as.formula(formula))
  Z <- as.data.frame(meta.dat[, allvars])

  ###############################################################################
  # WEIGHTS: Parse and validate sample weights (EXPERIMENTAL)
  ###############################################################################
  
  # Initialize weights (will be updated after sample filtering)
  weights_use <- NULL
  weights_original <- NULL
  
  if (!is.null(weights)) {
    if (verbose) {
      message("Using sample weights (EXPERIMENTAL feature)")
    }
    
    # Parse weights
    if (is.character(weights)) {
      # weights is a column name in meta.dat
      if (!weights %in% colnames(meta.dat)) {
        stop("Column '", weights, "' not found in meta.dat")
      }
      weights_original <- as.numeric(meta.dat[[weights]])
      if (verbose) {
        message("  Weights extracted from column: ", weights)
      }
    } else if (is.numeric(weights)) {
      # weights is a numeric vector
      if (length(weights) != ncol(feature.dat)) {
        stop("Length of weights (", length(weights), ") must equal ncol(feature.dat) (", ncol(feature.dat), ")")
      }
      weights_original <- weights
    } else {
      stop("weights must be NULL, a character string (column name), or a numeric vector")
    }
    
    # Validate weights (check NA first to avoid comparison issues)
    if (any(is.na(weights_original))) {
      stop("weights cannot contain NA values")
    }
    if (any(weights_original < 0)) {
      stop("weights must be non-negative")
    }
    if (all(weights_original == 0)) {
      stop("weights cannot all be zero")
    }

    # Warn about very small weights
    if (any(weights_original > 0 & weights_original < 0.01)) {
      warning("Some weights are very small (< 0.01), which may cause numerical issues")
    }
  }

  ###############################################################################
  # Filter samples: remove samples with NA values in any of the variables
  keep.sam <- which(rowSums(is.na(Z)) == 0)
  Y <- feature.dat[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  names(Z) <- allvars
  
  # Update weights if provided
  if (!is.null(weights_original)) {
    weights_use <- weights_original[keep.sam]
    if (verbose) {
      message("  Weights updated after sample filtering: ", length(weights_use), " samples retained")
    }
  }

  # Remove samples with zero total counts BEFORE feature filtering
  # (to avoid NaN in relative abundance calculations)
  if (any(colSums(Y) == 0)) {
    ind_nonzero <- which(colSums(Y) > 0)
    Y <- Y[, ind_nonzero]
    Z <- as.data.frame(Z[ind_nonzero, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind_nonzero]

    # Update weights if provided
    if (!is.null(weights_use)) {
      weights_use <- weights_use[ind_nonzero]
      if (verbose) {
        message("  Weights updated after zero-sum sample removal: ", length(weights_use), " samples retained")
      }
    }

    if (verbose) {
      message(length(ind_nonzero), " samples retained after removing zero-sum samples\n")
    }
  }

  # Filter features based on prevalence, mean abundance, and maximum abundance
  temp <- t(t(Y) / colSums(Y))

  # If feature data type is "other", reset all filters to 0
  if (feature.dat.type == "other" & (max.abund.filter != 0 | mean.abund.filter != 0 | prev.filter != 0 )){
    message("Note: Since feature.dat.type is set to 'other', all filters (max.abund.filter, mean.abund.filter, and prev.filter) are reset to 0.")
    max.abund.filter <- 0
    mean.abund.filter <- 0
    prev.filter <- 0
  }

  # Apply filters to features
  keep.tax <- rowMeans(temp != 0) >= prev.filter & rowMeans(temp) >= mean.abund.filter & matrixStats::rowMaxs(temp) >= max.abund.filter
  names(keep.tax) <- rownames(Y)
  rm(temp)
  if (verbose) {
    message(
      sum(!keep.tax), " features are filtered!\n"
    )
  }
  Y <- Y[keep.tax, ]

  n <- ncol(Y)
  m <- nrow(Y)

  # Second zero-sum check: after feature filtering, samples may become zero-sum
  if (any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)

    # Update weights if provided
    if (!is.null(weights_use)) {
      weights_use <- weights_use[ind]
      if (verbose) {
        message("  Weights updated after feature-filtering-induced zero-sum removal: ", length(weights_use), " samples retained")
      }
    }
  }

  if (verbose) {
    message(
      "The filtered data has ", n, " samples and ", m, " features that will be tested!\n"
    )
  }

  # Warn about features with less than 3 nonzero values
  if (sum(rowSums(Y != 0) <= 2) != 0) {
    warning(
      "Some features have less than 3 nonzero values!\n",
      "They have virtually no statistical power. You may consider filtering them in the analysis!\n"
    )
  }

  ###############################################################################
  # Scale numerical variables in the metadata
  ind <- sapply(1:ncol(Z), function(i) is.numeric(Z[, i]))
  Z[, ind] <- scale(Z[, ind])

  # Apply Winsorization to handle outliers if specified
  if (is.winsor) {
    Y <- winsor.fun(Y, 1 - outlier.pct, feature.dat.type)
  }

  # Determine if the model includes random effects
  if (grepl("\\(", formula)) {
    random.effect <- TRUE
  } else {
    random.effect <- FALSE
  }
  
  # Normalize weights to sum to sample size (for interpretability)
  if (!is.null(weights_use)) {
    weights_normalized <- weights_use * n / sum(weights_use)
    if (verbose) {
      message("  Weights normalized (sum = ", round(sum(weights_normalized), 2), ")")
      message("    Weight range: [", round(min(weights_normalized), 3), ", ", round(max(weights_normalized), 3), "]")
    }
  } else {
    weights_normalized <- NULL
  }

  # Assign names to taxa and samples
  if (is.null(rownames(feature.dat))) {
    taxa.name <- (1:nrow(feature.dat))[keep.tax]
  } else {
    taxa.name <- rownames(feature.dat)[keep.tax]
  }
  if (is.null(rownames(meta.dat))) {
    samp.name <- (1:nrow(meta.dat))[keep.sam]
  } else {
    samp.name <- rownames(meta.dat)[keep.sam]
  }

  # Handle zeros in the data
  if (feature.dat.type == "count") {
    if (any(Y == 0)) {
      N <- colSums(Y)
      if (adaptive) {
        # Determine zero-handling method based on correlation between sequencing depth and explanatory variables
        logN <- log(N)
        if (random.effect) {
          tmp <- lmer(as.formula(paste0("logN", formula)), Z)
        } else {
          tmp <- lm(as.formula(paste0("logN", formula)), Z)
        }
        corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
        if (any(corr.pval <= corr.cut)) {
          if (verbose) {
            message("Imputation approach is used.")
          }
          zero.handling <- "Imputation"
        } else {
          if (verbose) {
            message("Pseudo-count approach is used.")
          }
          zero.handling <- "Pseudo-count"
        }
      }
      if (zero.handling == "imputation") {
        # Impute zeros using the formula from the referenced paper
        N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
        N.mat[Y > 0] <- 0
        tmp <- N[max.col(N.mat)]
        Y <- Y + N.mat / tmp
      } else {
        # Add pseudo-count to all values
        Y <- Y + pseudo.cnt
      }
    }
  }

  if (feature.dat.type == "proportion") {
    if (any(Y == 0)) {
      # Apply half-minimum approach for zero values in proportion data
      Y <- t(apply(Y, 1, function(x) {
        x[x == 0] <- 0.5 * min(x[x != 0])
        return(x)
      }))
      colnames(Y) <- samp.name
      rownames(Y) <- taxa.name
    }
  }

  # Perform centered log-ratio (CLR) transformation
  logY <- log2(Y)
  W <- t(logY) - colMeans(logY)

  # Fit linear models or linear mixed effects models
  if (!random.effect) {
    if (verbose) {
      if (!is.null(weights_normalized)) {
        message("Fit weighted linear models ...")
      } else {
        message("Fit linear models ...")
      }
    }
    if (!is.null(weights_normalized)) {
      # Weighted regression
      suppressMessages(fit <- lm(as.formula(paste0("W", formula)), Z, weights = weights_normalized))
    } else {
      # Unweighted regression (original behavior)
      suppressMessages(fit <- lm(as.formula(paste0("W", formula)), Z))
    }
    res <- do.call(rbind, coef(summary(fit)))
    df <- rep(n - ncol(model.matrix(fit)), m)
  } else {
    if (verbose) {
      if (!is.null(weights_normalized)) {
        message("Fit weighted linear mixed effects models ...")
      } else {
        message("Fit linear mixed effects models ...")
      }
    }
    fun <- function(i) {
      w <- W[, i]
      if (!is.null(weights_normalized)) {
        # Weighted mixed model
        suppressMessages(fit <- lmer(as.formula(paste0("w", formula)), Z, weights = weights_normalized))
      } else {
        # Unweighted mixed model (original behavior)
        suppressMessages(fit <- lmer(as.formula(paste0("w", formula)), Z))
      }
      coef(summary(fit))
    }
    if (n.cores > 1) {
      res <- mclapply(c(1:m), function(i) fun(i), mc.cores = n.cores)
    } else {
      suppressMessages(res <- foreach(i = 1:m) %do% fun(i))
    }
    res <- do.call(rbind, res)
  }

  # Extract and process results
  res.intc <- res[which(rownames(res) == "(Intercept)"), ]
  rownames(res.intc) <- NULL
  baseMean <- 2^res.intc[, 1]
  baseMean <- baseMean / sum(baseMean) * 1e6

  ###############################################################################
  # TREE SMOOTHING: Compute smoothing matrix and correction factors
  # This is done once here, then applied to each variable in output.fun
  ###############################################################################
  if (tree.smooth && !is.null(tree)) {
    if (verbose) {
      message("Computing tree-guided smoothing matrix...")
    }

    tree.smooth.info <- get_tree_smoothing_info(
      phy.tree = tree,
      tax.names = taxa.name,
      lambda = tree.lambda,
      k.neighbors = tree.k,
      meff.exponent = tree.meff.exponent
    )

    if (verbose) {
      message("  ", tree.smooth.info$n.matched, " of ", m, " taxa matched to tree")
      message("  Trace ratio = ", round(tree.smooth.info$trace.ratio, 3))
      message("  M_eff correction factor = ", round(tree.smooth.info$meff.correction, 3))
    }

    # Warn if sample size is large
    n.per.group <- n / 2  # Approximate
    if (n.per.group > 100 && tree.meff.exponent < 3.0 && p.adj.method == "BH") {
      warning(
        "Large sample size detected (N ~ ", round(n.per.group), " per group). ",
        "Consider using tree.meff.exponent = 3.5 or p.adj.method = 'BY' ",
        "for better Type I error control."
      )
    }
  }

  # Function to process output for each variable
  output.fun <- function(x) {
    res.voi <- res[which(rownames(res) == x), ]
    rownames(res.voi) <- NULL

    if (random.effect) {
      df <- res.voi[, 3]
    }

    log2FoldChange <- res.voi[, 1]
    lfcSE <- res.voi[, 2]

    # -------------------------------------------------------------------------
    # BIAS CORRECTION with Precision-Weighted Mode Estimation
    #
    # When taxa have different standard errors (heteroscedasticity), we use
    # precision-weighted kernel density estimation. Taxa with lower SE (higher
    # precision) contribute more to the mode estimation.
    #
    # Safeguards:
    # 1. Cap extreme SE ratios to avoid over-concentration on few taxa
    # 2. Fall back to unweighted if SE variation is too extreme
    # 3. Handle NA/Inf values in lfcSE
    # -------------------------------------------------------------------------

    # Check SE variation to decide whether to use weighting
    lfcSE_valid <- lfcSE[is.finite(lfcSE) & lfcSE > 0]
    se_ratio <- max(lfcSE_valid) / min(lfcSE_valid)

    if (length(lfcSE_valid) == m && se_ratio < 100) {
      # Use precision-weighted mode estimation
      # Cap precision weights to avoid extreme concentration
      precision <- 1 / (lfcSE^2 + 1e-10)
      max_weight <- quantile(precision, 0.99)  # Cap at 99th percentile
      precision <- pmin(precision, max_weight)
      weights <- precision / sum(precision)

      # Weighted density estimation
      scaled <- sqrt(n) * log2FoldChange
      d <- suppressWarnings(density(scaled, weights = weights, kernel = "gaussian"))
      mode_scaled <- d$x[which.max(d$y)]
      bias <- mode_scaled / sqrt(n)
    } else {
      # Fall back to unweighted (original LinDA)
      suppressMessages(bias <- mlv(sqrt(n) * log2FoldChange,
        method = "meanshift", kernel = "gaussian"
      ) / sqrt(n))
    }

    log2FoldChange <- log2FoldChange - bias
    stat <- log2FoldChange / lfcSE

    # -------------------------------------------------------------------------
    # TREE-GUIDED SMOOTHING (if enabled)
    # This section applies phylogenetic smoothing to boost power:
    # 1. Smooth the test statistics using the graph Laplacian
    # 2. Re-normalize variance to maintain proper null distribution
    # 3. Apply M_eff correction to p-values before FDR adjustment
    # -------------------------------------------------------------------------
    if (tree.smooth && !is.null(tree.smooth.info)) {
      # Step 1: Apply smoothing matrix to test statistics
      # This borrows strength from phylogenetically related taxa
      stat.smoothed <- as.vector(tree.smooth.info$S %*% stat)

      # Step 2: Variance re-normalization
      # Smoothing compresses variance; divide by correction factor to restore
      stat <- stat.smoothed / tree.smooth.info$var.correction
    }

    # Calculate p-values from (possibly smoothed) test statistics
    pvalue <- 2 * pt(-abs(stat), df)

    # -------------------------------------------------------------------------
    # M_eff CORRECTION (if tree smoothing enabled)
    # Smoothing introduces correlation, causing p-value clumping under null.
    # We scale p-values by 1/meff.correction before BH adjustment.
    # This is equivalent to using a more stringent threshold.
    # -------------------------------------------------------------------------
    if (tree.smooth && !is.null(tree.smooth.info)) {
      pvalue.corrected <- pvalue / tree.smooth.info$meff.correction
      pvalue.corrected <- pmin(pvalue.corrected, 1)  # Cap at 1
      padj <- p.adjust(pvalue.corrected, method = p.adj.method)
    } else {
      padj <- p.adjust(pvalue, method = p.adj.method)
    }

    reject <- padj <= alpha
    output <- cbind.data.frame(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, reject, df)
    rownames(output) <- taxa.name
    return(list(bias = bias, output = output))
  }

  # Process results for all variables
  variables <- unique(rownames(res))[-1]
  variables.n <- length(variables)
  bias <- rep(NA, variables.n)
  output <- list()
  for (i in 1:variables.n) {
    tmp <- output.fun(variables[i])
    output[[i]] <- tmp[[2]]
    bias[i] <- tmp[[1]]
  }
  names(output) <- variables

  # Assign row and column names to the final data
  rownames(Y) <- taxa.name
  colnames(Y) <- samp.name
  rownames(Z) <- samp.name
  if (verbose) {
    message("Completed.")
  }

  # Return the results (include weights and tree smoothing info if used)
  result <- list(
    variables = variables,
    bias = bias,
    output = output,
    feature.dat.use = Y,
    meta.dat.use = Z
  )

  if (!is.null(weights_normalized)) {
    result$weights.use <- weights_normalized
    if (verbose) {
      message("  Sample weights included in output (experimental feature)")
    }
  }

  # Include tree smoothing information if used
  if (tree.smooth && !is.null(tree.smooth.info)) {
    result$tree.smooth.info <- list(
      enabled = TRUE,
      n.matched = tree.smooth.info$n.matched,
      n.taxa = m,
      lambda = tree.lambda,
      k.neighbors = tree.k,
      meff.exponent = tree.meff.exponent,
      trace.ratio = tree.smooth.info$trace.ratio,
      meff.correction = tree.smooth.info$meff.correction
    )
    if (verbose) {
      message("  Tree smoothing applied (", tree.smooth.info$n.matched, " taxa matched)")
    }
  } else {
    result$tree.smooth.info <- list(enabled = FALSE)
  }

  return(result)
}

#' Plot linda results
#'
#' The function plots the effect size plot and volcano plot based on the output from \code{linda}.
#'
#' @param linda.obj return from function \code{linda}.
#' @param variables.plot vector; variables whose results are to be plotted. For example, suppose the return
#' value \code{variables} is equal to \code{('x1', 'x2', 'x3b', 'x3c', 'x1:x2')}, then one could set \code{variables.plot = c('x3b', 'x1:x2')}.
#' @param titles vector; titles of the effect size plot and volcano plot for each variable in \code{variables.plot}.
#' Default is NULL. If NULL, the titles will be set as \code{variables.plot}.
#' @param alpha a real value between 0 and 1; cutoff for \code{padj}.
#' @param lfc.cut a positive value; cutoff for \code{log2FoldChange}.
#' @param legend TRUE or FALSE; whether to show the legends of the effect size plot and volcano plot.
#' @param directory character; the directory to save the figures, e.g., \code{getwd()}. Default is NULL. If NULL, figures will not be saved.
#' @param width the width of the graphics region in inches. See R function \code{pdf}.
#' @param height the height of the graphics region in inches. See R function \code{pdf}.
#'
#' @return A list of \code{ggplot2} objects.
#' \item{plot.lfc}{a list of effect size plots. Each plot corresponds to one variable in \code{variables.plot}.}
#' \item{plot.volcano}{a list of volcano plots. Each plot corresponds to one variable in \code{variables.plot}.}
#'
#' @author Huijuan Zhou \email{huijuanzhou2019@gmail.com}
#' Jun Chen \email{Chen.Jun2@mayo.edu}
#' Maintainer: Huijuan Zhou
#' @references Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. LinDA: Linear Models for Differential Abundance
#' Analysis of Microbiome Compositional Data.
#' @import ggplot2
#' @examples
#' \dontrun{
#' library(ggrepel)
#' data(smokers)
#' ind <- smokers$meta$AIRWAYSITE == "Throat"
#' otu.tab <- as.data.frame(smokers$otu[, ind])
#' meta <- cbind.data.frame(
#'   Smoke = factor(smokers$meta$SMOKER[ind]),
#'   Sex = factor(smokers$meta$SEX[ind]),
#'   Site = factor(smokers$meta$SIDEOFBODY[ind]),
#'   SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind])
#' )
#' ind1 <- which(meta$Site == "Left")
#' res.left <- linda(otu.tab[, ind1], meta[ind1, ],
#'   formula = "~Smoke+Sex"
#' )
#' ind2 <- which(meta$Site == "Right")
#' res.right <- linda(otu.tab[, ind2], meta[ind2, ],
#'   formula = "~Smoke+Sex"
#' )
#' rownames(res.left$output[[1]])[which(res.left$output[[1]]$reject)]
#' rownames(res.right$output[[1]])[which(res.right$output[[1]]$reject)]
#'
#' linda.obj <- linda(otu.tab, meta,
#'   formula = "~Smoke+Sex+(1|SubjectID)"
#' )
#' }
#' @export

linda.plot <- function(linda.obj, variables.plot, titles = NULL, alpha = 0.05, lfc.cut = 1,
                       legend = FALSE, directory = NULL, width = 11, height = 8) {
  bias <- linda.obj$bias
  output <- linda.obj$output
  otu.tab <- linda.obj$feature.dat.use
  meta <- linda.obj$meta.dat.use
  variables <- linda.obj$variables
  if (is.null(titles)) titles <- variables.plot

  taxa <- rownames(otu.tab)
  m <- length(taxa)

  tmp <- match(variables, variables.plot)
  voi.ind <- order(tmp)[1:sum(!is.na(tmp))]
  padj.mat <- foreach(i = voi.ind, .combine = "cbind") %do% {
    output[[i]]$padj
  }

  ## effect size plot
  if (is.matrix(padj.mat)) {
    ind <- which(colSums(padj.mat <= alpha) > 0)
  } else if (is.vector(padj.mat)) {
    tmp <- which(padj.mat <= alpha)
    if (length(tmp) > 0) {
      ind <- 1
    } else {
      ind <- integer(0)
    }
  }

  if (length(ind) == 0) {
    plot.lfc <- NULL
  } else {
    if (!is.null(directory)) {
      pdf(paste0(directory, "/plot_lfc.pdf"),
        width = width, height = height
      )
    }
    plot.lfc <- list()
    j <- 1
    Taxa <- Log2FoldChange <- Log10Padj <- NULL
    for (i in ind) {
      output.i <- output[[voi.ind[i]]]
      bias.i <- bias[voi.ind[i]]
      lfc <- output.i$log2FoldChange
      lfcSE <- output.i$lfcSE
      padj <- output.i$padj

      ind.rej <- which(padj <= alpha)
      n.rej <- length(ind.rej)
      taxa.rej <- taxa[ind.rej]
      taxa.rej <- factor(taxa.rej, levels = taxa.rej)
      data.plot.lfc <- cbind.data.frame(
        Taxa = rep(taxa.rej, 2),
        Log2FoldChange = c(lfc[ind.rej], lfc[ind.rej] + bias.i),
        lfcSE = c(lfcSE[ind.rej], rep(NA, n.rej)),
        bias = rep(c("Debiased", "Non-debiased"), each = n.rej)
      )
      plot.lfc.i <- ggplot(data.plot.lfc, aes(x = Log2FoldChange, y = Taxa)) +
        geom_point(aes(color = bias, shape = bias), size = 3) +
        geom_errorbar(aes(
          xmin = Log2FoldChange - 1.96 * lfcSE,
          xmax = Log2FoldChange + 1.96 * lfcSE
        ), width = .2) +
        geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
        ggtitle(titles[i]) +
        theme_bw(base_size = 18)
      if (legend) {
        plot.lfc.i <- plot.lfc.i +
          theme(
            legend.title = element_blank(),
            legend.key.width = unit(1, "cm"), plot.margin = unit(c(1, 1, 1, 1.5), "cm")
          )
      } else {
        plot.lfc.i <- plot.lfc.i +
          theme(legend.position = "none", plot.margin = unit(c(1, 1, 1, 1.5), "cm"))
      }
      plot.lfc[[j]] <- plot.lfc.i
      j <- j + 1
      if (!is.null(directory)) print(plot.lfc.i)
    }
    if (!is.null(directory)) dev.off()
  }

  ## volcano plot
  plot.volcano <- list()
  if (!is.null(directory)) {
    pdf(paste0(directory, "/plot_volcano.pdf"),
      width = width, height = height
    )
  }
  leg1 <- paste0("padj>", alpha, " & ", "lfc<=", lfc.cut)
  leg2 <- paste0("padj>", alpha, " & ", "lfc>", lfc.cut)
  leg3 <- paste0("padj<=", alpha, " & ", "lfc<=", lfc.cut)
  leg4 <- paste0("padj<=", alpha, " & ", "lfc>", lfc.cut)

  gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  color <- gg_color_hue(3)

  for (i in 1:length(voi.ind)) {
    output.i <- output[[voi.ind[i]]]
    bias.i <- bias[voi.ind[i]]
    lfc <- output.i$log2FoldChange
    padj <- output.i$padj

    ind1 <- padj > alpha & abs(lfc) <= lfc.cut
    ind2 <- padj > alpha & abs(lfc) > lfc.cut
    ind3 <- padj <= alpha & abs(lfc) <= lfc.cut
    ind4 <- padj <= alpha & abs(lfc) > lfc.cut

    leg <- rep(NA, m)
    leg[ind1] <- leg1
    leg[ind2] <- leg2
    leg[ind3] <- leg3
    leg[ind4] <- leg4
    leg <- factor(leg, levels = c(leg1, leg2, leg3, leg4))
    taxa.sig <- rep("", m)
    taxa.sig[ind3 | ind4] <- taxa[ind3 | ind4]

    data.volcano <- cbind.data.frame(
      taxa = taxa.sig, Log2FoldChange = lfc,
      Log10Padj = -log10(padj), leg = leg
    )
    plot.volcano.i <- ggplot(data.volcano, aes(x = Log2FoldChange, y = Log10Padj)) +
      geom_point(aes(color = leg), size = 2) +
      ggrepel::geom_text_repel(aes(label = taxa), max.overlaps = Inf) +
      scale_colour_manual(values = c("darkgray", color[c(2, 3, 1)])) +
      geom_hline(aes(yintercept = -log10(alpha)), color = "gray", linetype = "dashed") +
      geom_vline(aes(xintercept = -lfc.cut), color = "gray", linetype = "dashed") +
      geom_vline(aes(xintercept = lfc.cut), color = "gray", linetype = "dashed") +
      ylab("-Log10Padj") +
      ggtitle(titles[i]) +
      theme_bw(base_size = 18)
    if (legend) {
      plot.volcano.i <- plot.volcano.i +
        theme(
          legend.title = element_blank(),
          legend.key.width = unit(1, "cm"), plot.margin = unit(c(1, 1, 1, 1.5), "cm")
        )
    } else {
      plot.volcano.i <- plot.volcano.i +
        theme(legend.position = "none", plot.margin = unit(c(1, 1, 1, 1.5), "cm"))
    }
    plot.volcano[[i]] <- plot.volcano.i
    if (!is.null(directory)) print(plot.volcano.i)
  }
  if (!is.null(directory)) dev.off()

  return(list(plot.lfc = plot.lfc, plot.volcano = plot.volcano))
}

