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

#' @title Linear (Lin) model for differential abundance (DA) analysis
#'
#' @description This function implements a simple, robust, and highly scalable approach to tackle
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
#' @param verbose A boolean; if TRUE, progress messages will be printed. Default is TRUE.
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

linda <- function(feature.dat, meta.dat, phyloseq.obj = NULL, formula, feature.dat.type = c("count", "proportion","other"),
                  prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
                  is.winsor = TRUE, outlier.pct = 0.03,
                  adaptive = TRUE, zero.handling = c("pseudo-count", "imputation"),
                  pseudo.cnt = 0.5, corr.cut = 0.1,
                  p.adj.method = "BH", alpha = 0.05,
                  n.cores = 1, verbose = TRUE) {
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
  # Filter samples: remove samples with NA values in any of the variables
  keep.sam <- which(rowSums(is.na(Z)) == 0)
  Y <- feature.dat[, keep.sam]
  Z <- as.data.frame(Z[keep.sam, ])
  names(Z) <- allvars

  # Remove samples with zero total counts BEFORE computing proportions
  # This prevents division by zero in the next step
  col_sums <- colSums(Y)
  if (any(col_sums == 0)) {
    zero_samples <- which(col_sums == 0)
    if (verbose) {
      message(
        length(zero_samples), " samples with zero total counts removed before filtering!\n"
      )
    }
    Y <- Y[, col_sums > 0, drop = FALSE]
    Z <- as.data.frame(Z[col_sums > 0, , drop = FALSE])
    names(Z) <- allvars
    keep.sam <- keep.sam[col_sums > 0]
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

  # Remove samples with zero total counts after feature filtering
  if (any(colSums(Y) == 0)) {
    ind <- which(colSums(Y) > 0)
    Y <- Y[, ind]
    Z <- as.data.frame(Z[ind, ])
    names(Z) <- allvars
    keep.sam <- keep.sam[ind]
    n <- ncol(Y)
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
      message("Fit linear models ...")
    }
    suppressMessages(fit <- lm(as.formula(paste0("W", formula)), Z))
    res <- do.call(rbind, coef(summary(fit)))
    df <- rep(n - ncol(model.matrix(fit)), m)
  } else {
    if (verbose) {
      message("Fit linear mixed effects models ...")
    }
    fun <- function(i) {
      w <- W[, i]
      suppressMessages(fit <- lmer(as.formula(paste0("w", formula)), Z))
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

  # Function to process output for each variable
  output.fun <- function(x) {
    res.voi <- res[which(rownames(res) == x), ]
    rownames(res.voi) <- NULL

    if (random.effect) {
      df <- res.voi[, 3]
    }

    log2FoldChange <- res.voi[, 1]
    lfcSE <- res.voi[, 2]
    
    # Estimate and correct for bias using the mode of the regression coefficients
    suppressMessages(bias <- mlv(sqrt(n) * log2FoldChange,
      method = "meanshift", kernel = "gaussian"
    ) / sqrt(n))
    log2FoldChange <- log2FoldChange - bias
    stat <- log2FoldChange / lfcSE

    # Calculate p-values and adjust for multiple testing
    pvalue <- 2 * pt(-abs(stat), df)
    padj <- p.adjust(pvalue, method = p.adj.method)
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
  
  # Return the results
  return(list(variables = variables, bias = bias, output = output, feature.dat.use = Y, meta.dat.use = Z))
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