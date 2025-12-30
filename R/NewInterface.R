#' @title MicrobiomeStat Unified Interface
#' @description Simplified, unified API for microbiome data analysis.
#'
#' This module provides 6 core functions that automatically route to the
#' appropriate底层 functions based on study design and parameters:
#'
#' - `plot_taxa()`: Visualize taxonomic composition and abundance
#' - `plot_alpha()`: Visualize alpha diversity metrics
#' - `plot_beta()`: Visualize beta diversity and ordination
#' - `test_taxa()`: Statistical tests for taxonomic features
#' - `test_beta()`: Statistical tests for beta diversity
#' - `test_alpha()`: Statistical tests for alpha diversity
#'
#' @name MicrobiomeStat-unified
NULL


# =============================================================================
# Taxa Visualization
# =============================================================================

#' Unified Taxa Visualization
#'
#' One function for all taxa/feature visualizations. Automatically detects study
#' design (cross-sectional, longitudinal, or paired) and routes to the appropriate
#' underlying function.
#'
#' @param data.obj A MicrobiomeStat data object containing feature.tab, meta.dat,
#'   and optionally feature.ann, tree, and feature.agg.list.
#' @param plot.type Type of visualization:
#'   - "barplot": Stacked bar plot of relative abundances
#'   - "boxplot": Box plots comparing groups/time points
#'   - "heatmap": Heatmap of feature abundances
#'   - "dotplot": Dot plot with effect sizes
#'   - "areaplot": Stacked area plot (longitudinal only)
#'   - "spaghettiplot": Individual trajectories (longitudinal only)
#'   - "cladogram": Phylogenetic cladogram (single time point only)
#' @param subject.var Character. Name of the subject/sample ID variable in meta.dat.
#'   Required for longitudinal and paired designs.
#' @param time.var Character. Name of the time variable in meta.dat.
#'   NULL for cross-sectional studies.
#' @param group.var Character. Name of the grouping variable (e.g., treatment, condition).
#' @param strata.var Character. Name of the stratification variable for faceting.
#' @param time.points Time point specification. Can be:
#'   - NULL: use all available time points (auto-detect design)
#'   - Single value: cross-sectional at that time point
#'   - Vector of 2: paired design (baseline, followup)
#'   - Vector of >2: longitudinal with specific time points
#'   - List: list(baseline = "T0", followup = c("T1", "T2"))
#' @param feature.level Character vector. Taxonomic level(s) for aggregation.
#'   Options include column names from feature.ann or "original" for no aggregation.
#' @param feature.select Feature selection. Can be:
#'   - Integer: top N features by mean abundance (default: 20)
#'   - Character vector: specific feature names to display
#'   - NULL: use function defaults
#' @param feature.dat.type Data type: "count", "proportion", or "other".
#' @param change.type For paired/longitudinal comparisons:
#'   - "none": show raw values (default)
#'   - "relative": relative change (x1-x0)/(x1+x0)
#'   - "log_fold": log2 fold change
#'   - "absolute": absolute difference
#' @param prev.filter Numeric. Minimum prevalence threshold (0-1) for filtering.
#' @param abund.filter Numeric. Minimum abundance threshold for filtering.
#' @param transform Transformation for visualization: "identity", "sqrt", or "log".
#' @param theme Theme specification. Can be:
#'   - Character: preset name ("bw", "classic", "minimal", "prism")
#'   - List: list(base.size = 12, choice = "bw", palette = NULL)
#' @param ... Additional arguments passed to the underlying function.
#'
#' @return A list of ggplot objects, one per feature.level.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Cross-sectional barplot
#' plot_taxa(peerj32.obj, "barplot",
#'           group.var = "group",
#'           feature.level = "Phylum")
#'
#' # Longitudinal boxplot with change
#' plot_taxa(peerj32.obj, "boxplot",
#'           subject.var = "subject",
#'           time.var = "time",
#'           group.var = "group",
#'           change.type = "relative")
#'
#' # Paired heatmap
#' plot_taxa(peerj32.obj, "heatmap",
#'           subject.var = "subject",
#'           time.var = "time",
#'           time.points = c("1", "2"),
#'           feature.level = "Genus")
#' }
#'
#' @export
plot_taxa <- function(data.obj,
                      plot.type = c("barplot", "boxplot", "heatmap", "dotplot",
                                    "areaplot", "spaghettiplot", "cladogram"),
                      subject.var = NULL,
                      time.var = NULL,
                      group.var = NULL,
                      strata.var = NULL,
                      time.points = NULL,
                      feature.level = "original",
                      feature.select = 20,
                      feature.dat.type = c("count", "proportion", "other"),
                      change.type = c("none", "relative", "log_fold", "absolute"),
                      prev.filter = 0,
                      abund.filter = 0,
                      transform = c("identity", "sqrt", "log"),
                      theme = "bw",
                      ...) {

  # Match arguments
  plot.type <- match.arg(plot.type)
  feature.dat.type <- match.arg(feature.dat.type)
  change.type <- match.arg(change.type)
  transform <- match.arg(transform)

  # Detect study design
  design_info <- detect_study_design(data.obj, subject.var, time.var, time.points)

  # Determine if this is a change plot
  is_change <- change.type != "none"

  # Validate inputs
  validation <- validate_inputs(
    data.obj = data.obj,
    category = "taxa",
    plot_type = if (is_change) paste0("change_", plot.type) else plot.type,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    design_info = design_info
  )

  # If change plot not available, try regular plot

if (!validation$is_valid && is_change) {
    validation <- validate_inputs(
      data.obj = data.obj,
      category = "taxa",
      plot_type = plot.type,
      subject.var = subject.var,
      time.var = time.var,
      group.var = group.var,
      design_info = design_info
    )
    if (validation$is_valid) {
      message("Note: change plot not available, using regular ", plot.type)
      is_change <- FALSE
    }
  }

  if (!validation$is_valid) {
    stop(format_validation_errors(validation, "plot_taxa"))
  }

  # Get target function
  target_func <- validation$target_func

  # Build argument list
  args <- list(
    data.obj = data.obj,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    strata.var = strata.var,
    feature.level = feature.level,
    feature.dat.type = feature.dat.type,
    prev.filter = prev.filter,
    abund.filter = abund.filter,
    feature.select = feature.select,
    theme = theme,
    change.type = change.type
  )

  # Add transform for boxplot
  if (plot.type == "boxplot") {
    args$transform <- transform
  }

  # Add extra arguments
  extra_args <- list(...)
  args <- c(args, extra_args)

  # Resolve parameters for target function
  resolved <- resolve_params(args, target_func, design_info)

  # Clean parameters (remove NULLs except key params, remove subject.var for single design)
  resolved <- clean_resolved_params(resolved, design_info$design)

  # Call the target function
  result <- tryCatch({
    do.call(target_func, resolved)
  }, error = function(e) {
    stop("Error in ", target_func, ": ", e$message, call. = FALSE)
  })

  # Add metadata
  attr(result, "source") <- "plot_taxa"
  attr(result, "plot.type") <- plot.type
  attr(result, "design") <- design_info$design

  return(result)
}


# =============================================================================
# Alpha Diversity Visualization
# =============================================================================

#' Unified Alpha Diversity Visualization
#'
#' One function for all alpha diversity visualizations.
#'
#' @inheritParams plot_taxa
#' @param plot.type Type of visualization:
#'   - "boxplot": Box plots comparing groups/time points
#'   - "spaghettiplot": Individual trajectories (longitudinal)
#'   - "dotplot": Dot plot of per-time statistics (longitudinal)
#' @param alpha.name Character vector. Alpha diversity indices to calculate.
#'   Options: "shannon", "simpson", "observed_species", "chao1", "ace", "pielou".
#' @param alpha.obj Optional pre-calculated alpha diversity object.
#' @param depth Rarefaction depth. NULL for no rarefaction.
#'
#' @return A list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Cross-sectional alpha boxplot
#' plot_alpha(peerj32.obj, "boxplot",
#'            alpha.name = c("shannon", "observed_species"),
#'            group.var = "group")
#'
#' # Longitudinal spaghettiplot
#' plot_alpha(peerj32.obj, "spaghettiplot",
#'            subject.var = "subject",
#'            time.var = "time",
#'            alpha.name = "shannon",
#'            group.var = "group")
#' }
#'
#' @export
plot_alpha <- function(data.obj,
                       plot.type = c("boxplot", "spaghettiplot", "dotplot"),
                       alpha.name = c("shannon", "observed_species"),
                       alpha.obj = NULL,
                       depth = NULL,
                       subject.var = NULL,
                       time.var = NULL,
                       group.var = NULL,
                       strata.var = NULL,
                       time.points = NULL,
                       change.type = c("none", "relative", "log_fold", "absolute"),
                       theme = "bw",
                       ...) {

  # Match arguments
  plot.type <- match.arg(plot.type)
  change.type <- match.arg(change.type)

  # Detect study design
  design_info <- detect_study_design(data.obj, subject.var, time.var, time.points)

  # Determine if this is a change plot
  is_change <- change.type != "none"

  # Validate inputs
  validation <- validate_inputs(
    data.obj = data.obj,
    category = "alpha",
    plot_type = if (is_change) paste0("change_", plot.type) else plot.type,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    design_info = design_info
  )

  # Fallback for change plots
  if (!validation$is_valid && is_change) {
    validation <- validate_inputs(
      data.obj = data.obj,
      category = "alpha",
      plot_type = plot.type,
      subject.var = subject.var,
      time.var = time.var,
      group.var = group.var,
      design_info = design_info
    )
    if (validation$is_valid) {
      message("Note: change plot not available, using regular ", plot.type)
      is_change <- FALSE
    }
  }

  if (!validation$is_valid) {
    stop(format_validation_errors(validation, "plot_alpha"))
  }

  # Get target function
  target_func <- validation$target_func

  # Build argument list
  args <- list(
    data.obj = data.obj,
    alpha.obj = alpha.obj,
    alpha.name = alpha.name,
    depth = depth,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    strata.var = strata.var,
    theme = theme,
    change.type = if (is_change) change.type else NULL
  )

  # Add extra arguments
  extra_args <- list(...)
  args <- c(args, extra_args)

  # Resolve parameters
  resolved <- resolve_params(args, target_func, design_info)
  resolved <- clean_resolved_params(resolved, design_info$design)

  # Call the target function
  result <- tryCatch({
    do.call(target_func, resolved)
  }, error = function(e) {
    stop("Error in ", target_func, ": ", e$message, call. = FALSE)
  })

  attr(result, "source") <- "plot_alpha"
  attr(result, "plot.type") <- plot.type
  attr(result, "design") <- design_info$design

  return(result)
}


# =============================================================================
# Beta Diversity Visualization
# =============================================================================

#' Unified Beta Diversity Visualization
#'
#' One function for all beta diversity visualizations.
#'
#' @inheritParams plot_taxa
#' @param plot.type Type of visualization:
#'   - "ordination": PCoA/MDS ordination plot
#'   - "boxplot": Box plots of beta diversity distances
#'   - "spaghettiplot": Individual distance trajectories
#'   - "dotplot": Dot plot of per-time statistics
#' @param dist.name Character vector. Distance metrics to use.
#'   Options: "BC" (Bray-Curtis), "Jaccard", "UniFrac", "GUniFrac", "WUniFrac", "JS".
#' @param dist.obj Optional pre-calculated distance object.
#' @param pc.obj Optional pre-calculated PC/MDS object (for ordination).
#'
#' @return A list of ggplot objects.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Cross-sectional ordination
#' plot_beta(peerj32.obj, "ordination",
#'           dist.name = c("BC", "Jaccard"),
#'           group.var = "group")
#'
#' # Longitudinal ordination
#' plot_beta(peerj32.obj, "ordination",
#'           subject.var = "subject",
#'           time.var = "time",
#'           dist.name = "BC",
#'           group.var = "group")
#' }
#'
#' @export
plot_beta <- function(data.obj,
                      plot.type = c("ordination", "boxplot", "spaghettiplot", "dotplot"),
                      dist.name = c("BC", "Jaccard"),
                      dist.obj = NULL,
                      pc.obj = NULL,
                      subject.var = NULL,
                      time.var = NULL,
                      group.var = NULL,
                      strata.var = NULL,
                      time.points = NULL,
                      change.type = c("none", "relative"),
                      theme = "bw",
                      ...) {

  # Match arguments
  plot.type <- match.arg(plot.type)
  change.type <- match.arg(change.type)

  # Detect study design
  design_info <- detect_study_design(data.obj, subject.var, time.var, time.points)

  # Validate inputs
  validation <- validate_inputs(
    data.obj = data.obj,
    category = "beta",
    plot_type = plot.type,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    design_info = design_info
  )

  if (!validation$is_valid) {
    stop(format_validation_errors(validation, "plot_beta"))
  }

  # Get target function
  target_func <- validation$target_func

  # Build argument list
  args <- list(
    data.obj = data.obj,
    dist.obj = dist.obj,
    pc.obj = pc.obj,
    dist.name = dist.name,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    strata.var = strata.var,
    theme = theme
  )

  # Add extra arguments
  extra_args <- list(...)
  args <- c(args, extra_args)

  # Resolve parameters
  resolved <- resolve_params(args, target_func, design_info)
  resolved <- clean_resolved_params(resolved, design_info$design)

  # Call the target function
  result <- tryCatch({
    do.call(target_func, resolved)
  }, error = function(e) {
    stop("Error in ", target_func, ": ", e$message, call. = FALSE)
  })

  attr(result, "source") <- "plot_beta"
  attr(result, "plot.type") <- plot.type
  attr(result, "design") <- design_info$design

  return(result)
}


# =============================================================================
# Taxa Statistical Testing
# =============================================================================

#' Unified Taxa Statistical Testing
#'
#' One function for all taxonomic feature statistical tests.
#'
#' @inheritParams plot_taxa
#' @param test.type Type of test:
#'   - "difference": Test for group differences (cross-sectional or paired)
#'   - "trend": Test for temporal trends (longitudinal)
#'   - "volatility": Test for temporal instability (longitudinal)
#'   - "association": Test for associations with continuous variables
#' @param adj.vars Character vector. Names of adjustment/covariate variables.
#' @param feature.mt.method Multiple testing correction: "fdr", "bonferroni", or "none".
#' @param feature.sig.level Significance threshold (default 0.05).
#'
#' @return A list containing statistical results (p-values, coefficients, etc.).
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Cross-sectional group comparison
#' test_taxa(peerj32.obj, "difference",
#'           group.var = "group",
#'           feature.level = "Genus")
#'
#' # Longitudinal trend test
#' test_taxa(peerj32.obj, "trend",
#'           subject.var = "subject",
#'           time.var = "time",
#'           group.var = "group",
#'           feature.level = "Genus")
#' }
#'
#' @export
test_taxa <- function(data.obj,
                      test.type = c("difference", "trend", "volatility", "association"),
                      subject.var = NULL,
                      time.var = NULL,
                      group.var = NULL,
                      adj.vars = NULL,
                      time.points = NULL,
                      feature.level = "original",
                      feature.dat.type = c("count", "proportion", "other"),
                      prev.filter = 0.1,
                      abund.filter = 0.0001,
                      feature.mt.method = c("fdr", "bonferroni", "none"),
                      feature.sig.level = 0.05,
                      ...) {

  # Match arguments
  test.type <- match.arg(test.type)
  feature.dat.type <- match.arg(feature.dat.type)
  feature.mt.method <- match.arg(feature.mt.method)

  # Detect study design
  design_info <- detect_study_design(data.obj, subject.var, time.var, time.points)

  # Validate inputs
  validation <- validate_inputs(
    data.obj = data.obj,
    category = "taxa_test",
    plot_type = test.type,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    design_info = design_info
  )

  if (!validation$is_valid) {
    stop(format_validation_errors(validation, "test_taxa"))
  }

  # Get target function
  target_func <- validation$target_func

  # Build argument list
  args <- list(
    data.obj = data.obj,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    adj.vars = adj.vars,
    feature.level = feature.level,
    feature.dat.type = feature.dat.type,
    prev.filter = prev.filter,
    abund.filter = abund.filter
  )

  # Add extra arguments
  extra_args <- list(...)
  args <- c(args, extra_args)

  # Resolve parameters
  resolved <- resolve_params(args, target_func, design_info)
  resolved <- clean_resolved_params(resolved, design_info$design)

  # Call the target function
  result <- tryCatch({
    do.call(target_func, resolved)
  }, error = function(e) {
    stop("Error in ", target_func, ": ", e$message, call. = FALSE)
  })

  attr(result, "source") <- "test_taxa"
  attr(result, "test.type") <- test.type
  attr(result, "design") <- design_info$design

  return(result)
}


# =============================================================================
# Beta Diversity Statistical Testing
# =============================================================================

#' Unified Beta Diversity Statistical Testing
#'
#' One function for all beta diversity statistical tests (PERMANOVA, etc.).
#'
#' @inheritParams plot_beta
#' @param test.type Type of test:
#'   - "difference": PERMANOVA for group differences
#'   - "trend": Test for temporal trends in community structure
#'   - "volatility": Test for temporal instability
#' @param adj.vars Character vector. Names of adjustment/covariate variables.
#'
#' @return A list containing statistical results.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Cross-sectional PERMANOVA
#' test_beta(peerj32.obj, "difference",
#'           dist.name = c("BC", "Jaccard"),
#'           group.var = "group")
#' }
#'
#' @export
test_beta <- function(data.obj,
                      test.type = c("difference", "trend", "volatility"),
                      dist.name = c("BC", "Jaccard"),
                      dist.obj = NULL,
                      subject.var = NULL,
                      time.var = NULL,
                      group.var = NULL,
                      adj.vars = NULL,
                      time.points = NULL,
                      ...) {

  # Match arguments
  test.type <- match.arg(test.type)

  # Detect study design
  design_info <- detect_study_design(data.obj, subject.var, time.var, time.points)

  # Validate inputs
  validation <- validate_inputs(
    data.obj = data.obj,
    category = "beta_test",
    plot_type = test.type,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    design_info = design_info
  )

  if (!validation$is_valid) {
    stop(format_validation_errors(validation, "test_beta"))
  }

  # Get target function
  target_func <- validation$target_func

  # Build argument list
  args <- list(
    data.obj = data.obj,
    dist.obj = dist.obj,
    dist.name = dist.name,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    adj.vars = adj.vars
  )

  # Add extra arguments
  extra_args <- list(...)
  args <- c(args, extra_args)

  # Resolve parameters
  resolved <- resolve_params(args, target_func, design_info)
  resolved <- clean_resolved_params(resolved, design_info$design)

  # Call the target function
  result <- tryCatch({
    do.call(target_func, resolved)
  }, error = function(e) {
    stop("Error in ", target_func, ": ", e$message, call. = FALSE)
  })

  attr(result, "source") <- "test_beta"
  attr(result, "test.type") <- test.type
  attr(result, "design") <- design_info$design

  return(result)
}


# =============================================================================
# Alpha Diversity Statistical Testing
# =============================================================================

#' Unified Alpha Diversity Statistical Testing
#'
#' One function for all alpha diversity statistical tests.
#'
#' @inheritParams plot_alpha
#' @param test.type Type of test:
#'   - "difference": Test for group differences in alpha diversity
#'   - "trend": Test for temporal trends in alpha diversity
#'   - "volatility": Test for temporal instability in alpha diversity
#'   - "per_time": Test at each time point separately
#' @param adj.vars Character vector. Names of adjustment/covariate variables.
#'
#' @return A list containing statistical results for each alpha diversity metric.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#'
#' # Cross-sectional test
#' test_alpha(peerj32.obj, "difference",
#'            alpha.name = c("shannon", "simpson"),
#'            group.var = "group")
#'
#' # Longitudinal trend test
#' test_alpha(peerj32.obj, "trend",
#'            alpha.name = "shannon",
#'            subject.var = "subject",
#'            time.var = "time",
#'            group.var = "group")
#' }
#'
#' @export
test_alpha <- function(data.obj,
                       test.type = c("difference", "trend", "volatility", "per_time"),
                       alpha.name = c("shannon", "simpson", "observed_species"),
                       alpha.obj = NULL,
                       depth = NULL,
                       subject.var = NULL,
                       time.var = NULL,
                       group.var = NULL,
                       adj.vars = NULL,
                       time.points = NULL,
                       ...) {

  # Match arguments
  test.type <- match.arg(test.type)

  # Detect study design
  design_info <- detect_study_design(data.obj, subject.var, time.var, time.points)

  # Validate inputs
  validation <- validate_inputs(
    data.obj = data.obj,
    category = "alpha_test",
    plot_type = test.type,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    design_info = design_info
  )

  if (!validation$is_valid) {
    stop(format_validation_errors(validation, "test_alpha"))
  }

  # Get target function
  target_func <- validation$target_func

  # Build argument list
  args <- list(
    data.obj = data.obj,
    alpha.obj = alpha.obj,
    alpha.name = alpha.name,
    depth = depth,
    subject.var = subject.var,
    time.var = time.var,
    group.var = group.var,
    adj.vars = adj.vars
  )

  # Add extra arguments
  extra_args <- list(...)
  args <- c(args, extra_args)

  # Resolve parameters
  resolved <- resolve_params(args, target_func, design_info)
  resolved <- clean_resolved_params(resolved, design_info$design)

  # Call the target function
  result <- tryCatch({
    do.call(target_func, resolved)
  }, error = function(e) {
    stop("Error in ", target_func, ": ", e$message, call. = FALSE)
  })

  attr(result, "source") <- "test_alpha"
  attr(result, "test.type") <- test.type
  attr(result, "design") <- design_info$design

  return(result)
}
