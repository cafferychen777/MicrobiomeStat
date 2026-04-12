#' @title Interface Utilities for MicrobiomeStat Unified API
#' @description Internal utility functions for the unified interface layer.
#' These functions handle study design detection, parameter resolution,
#' function routing, and input validation.
#'
#' @name interface_utils
#' @keywords internal
NULL

# =============================================================================
# Study Design Detection
# =============================================================================

#' Detect Study Design Type
#'
#' Automatically determines whether the data represents a cross-sectional (single),
#' longitudinal (long), or paired (pair) study design based on the provided parameters.
#'
#' @param data.obj MicrobiomeStat data object
#' @param subject.var Subject identifier variable name
#' @param time.var Time variable name (NULL for cross-sectional)
#' @param time.points Time point specification. Can be:
#'   - NULL: use all available time points
#'   - Single value: cross-sectional at that time point
#'   - Vector of length 2: paired design (baseline, followup)
#'   - Vector of length > 2: longitudinal design
#'   - List with `baseline` and `followup` elements: paired design
#'
#' @return A list containing:
#'   - `design`: "single", "long", or "pair"
#'   - `t0.level`: baseline time point (NULL for single)
#'   - `ts.levels`: follow-up time points (NULL for single)
#'   - `t.level`: single time point for subsetting (for single design)
#'
#' @keywords internal
detect_study_design <- function(data.obj,
                                 subject.var = NULL,
                                 time.var = NULL,
                                 time.points = NULL) {



  result <- list(
    design = "single",
    t0.level = NULL,
    ts.levels = NULL,
    t.level = NULL
  )

  # No time variable = cross-sectional
  if (is.null(time.var)) {
    return(result)
  }

  # Get available time points from data
  available_times <- NULL
  if (!is.null(data.obj$meta.dat) && time.var %in% colnames(data.obj$meta.dat)) {
    available_times <- mStat_order_time_values(data.obj$meta.dat[[time.var]])
  }

  # Parse time.points specification
  if (is.null(time.points)) {
    # Use all available time points
    if (is.null(available_times) || length(available_times) < 2) {
      # Only one time point or unknown - treat as single
      result$t.level <- if (!is.null(available_times)) available_times[1] else NULL
      return(result)
    }

    if (length(available_times) == 2) {
      result$design <- "pair"
      result$t0.level <- available_times[1]
      result$ts.levels <- available_times[2]
    } else {
      result$design <- "long"
      result$t0.level <- available_times[1]
      result$ts.levels <- available_times[-1]
    }

  } else if (is.list(time.points)) {
    # List format: list(baseline = "T0", followup = c("T1", "T2"))
    if (!is.null(time.points$baseline)) {
      result$t0.level <- time.points$baseline
    }
    if (!is.null(time.points$followup)) {
      result$ts.levels <- time.points$followup
    }

    if (!is.null(result$t0.level) && !is.null(result$ts.levels)) {
      result$design <- if (length(result$ts.levels) == 1) "pair" else "long"
    } else if (!is.null(result$t0.level)) {
      result$design <- "single"
      result$t.level <- result$t0.level
    }

  } else if (length(time.points) == 1) {
    # Single time point - cross-sectional
    result$design <- "single"
    result$t.level <- time.points

  } else if (length(time.points) == 2) {
    # Two time points - paired
    result$design <- "pair"
    result$t0.level <- time.points[1]
    result$ts.levels <- time.points[2]

  } else {
    # Multiple time points - longitudinal
    result$design <- "long"
    result$t0.level <- time.points[1]
    result$ts.levels <- time.points[-1]
  }

  return(result)
}


#' @keywords internal
mStat_design_time_values <- function(design_info) {
  if (is.null(design_info) || is.null(design_info$design)) {
    return(character())
  }

  if (identical(design_info$design, "single") && !is.null(design_info$t.level)) {
    return(as.character(design_info$t.level))
  }

  if (design_info$design %in% c("pair", "long")) {
    return(as.character(c(design_info$t0.level, design_info$ts.levels)))
  }

  character()
}


#' @keywords internal
mStat_scope_unified_context <- function(data.obj,
                                        design_info,
                                        time.var = NULL,
                                        alpha.obj = NULL,
                                        dist.obj = NULL,
                                        pc.obj = NULL) {
  scoped <- list(
    data.obj = data.obj,
    alpha.obj = alpha.obj,
    dist.obj = dist.obj,
    pc.obj = pc.obj
  )

  if (is.null(data.obj$meta.dat) || is.null(time.var) || !time.var %in% colnames(data.obj$meta.dat)) {
    return(scoped)
  }

  candidate_samIDs <- rownames(data.obj$meta.dat)

  if (!is.null(alpha.obj) && length(alpha.obj) > 0 && !is.null(alpha.obj[[1]])) {
    alpha_samIDs <- rownames(alpha.obj[[1]])
    if (is.null(alpha_samIDs)) {
      stop("alpha.obj must have sample row names.", call. = FALSE)
    }
    candidate_samIDs <- intersect(candidate_samIDs, alpha_samIDs)
  }

  if (!is.null(dist.obj) && length(dist.obj) > 0 && !is.null(dist.obj[[1]])) {
    dist_samIDs <- mStat_get_dist_labels(dist.obj[[1]])
    candidate_samIDs <- intersect(candidate_samIDs, dist_samIDs)
  }

  if (!is.null(pc.obj) && length(pc.obj) > 0 && !is.null(pc.obj[[1]]$points)) {
    pc_samIDs <- rownames(pc.obj[[1]]$points)
    candidate_samIDs <- intersect(candidate_samIDs, pc_samIDs)
  }

  requested_time_values <- unique(mStat_design_time_values(design_info))

  if (length(requested_time_values) > 0) {
    scoped_meta <- data.obj$meta.dat[candidate_samIDs, , drop = FALSE]
    resolved_time <- mStat_resolve_time_levels(
      values = scoped_meta[[time.var]],
      time.var = time.var,
      requested_levels = requested_time_values,
      context = "unified interface time scoping"
    )
    candidate_samIDs <- rownames(
      scoped_meta[
        mStat_match_metadata_values(scoped_meta[[time.var]], resolved_time$kept_levels),
        ,
        drop = FALSE
      ]
    )
  }

  if (length(candidate_samIDs) == 0) {
    stop("No samples remain after applying the requested time scope.", call. = FALSE)
  }

  if (!identical(candidate_samIDs, rownames(data.obj$meta.dat))) {
    scoped$data.obj <- mStat_subset_data(data.obj, samIDs = candidate_samIDs)
  }

  if (!is.null(alpha.obj)) {
    scoped$alpha.obj <- mStat_subset_alpha(alpha.obj, candidate_samIDs)
  }

  if (!is.null(dist.obj) || !is.null(pc.obj)) {
    aligned <- mStat_subset_beta_objects(dist.obj, pc.obj, candidate_samIDs)
    scoped$dist.obj <- aligned$dist.obj
    scoped$pc.obj <- aligned$pc.obj
  }

  scoped
}


# =============================================================================
# Parameter Resolution
# =============================================================================

#' Resolve Theme Parameters
#'
#' Converts theme specification to individual parameters expected by underlying functions.
#'
#' @param theme Theme specification. Can be:
#'   - Character: preset theme name ("bw", "classic", "minimal", etc.)
#'   - List: with elements base.size, choice, custom, palette
#'
#' @return List with base.size, theme.choice, custom.theme, palette
#'
#' @keywords internal
resolve_theme <- function(theme = NULL) {
  defaults <- list(
    base.size = 12,
    theme.choice = "bw",
    custom.theme = NULL,
    palette = NULL
  )

  if (is.null(theme)) {
    return(defaults)
  }

  if (is.character(theme) && length(theme) == 1) {
    defaults$theme.choice <- theme
    return(defaults)
  }

  if (is.list(theme)) {
    if (!is.null(theme$base.size)) defaults$base.size <- theme$base.size
    if (!is.null(theme$choice)) defaults$theme.choice <- theme$choice
    if (!is.null(theme$custom)) defaults$custom.theme <- theme$custom
    if (!is.null(theme$palette)) defaults$palette <- theme$palette
  }

  return(defaults)
}


#' Resolve Feature Selection Parameters
#'
#' Converts unified feature.select parameter to function-specific parameters.
#'
#' @param feature.select Feature selection specification. Can be:
#'   - Integer: top N features by abundance
#'   - Character vector: specific feature names
#'   - NULL: use default or all features
#' @param target_type Target function type: "barplot", "boxplot", "heatmap", etc.
#'
#' @return List with appropriate parameters (feature.number or top.k.plot, features.plot)
#'
#' @keywords internal
resolve_feature_select <- function(feature.select = NULL, target_type = "default") {
  result <- list(
    feature.number = NULL,
    top.k.plot = NULL,
    features.plot = NULL
  )

  if (is.null(feature.select)) {
    # Use defaults based on function type
    if (target_type %in% c("barplot", "areaplot")) {
      result$feature.number <- 20
    }
    return(result)
  }

  if (is.numeric(feature.select) && length(feature.select) == 1) {
    # Integer: top N features
    if (target_type %in% c("barplot", "areaplot")) {
      result$feature.number <- as.integer(feature.select)
    } else {
      result$top.k.plot <- as.integer(feature.select)
    }
  } else if (is.character(feature.select)) {
    # Character vector: specific features
    result$features.plot <- feature.select
  }

  return(result)
}


#' Resolve Change Type Parameters
#'
#' Converts unified change.type to function-specific parameters.
#'
#' @param change.type Change type: "none", "relative", "log_fold", "absolute"
#' @param context Change context: "taxa" or "alpha"
#'
#' @return List with the resolved change function and related parameters
#'
#' @keywords internal
resolve_change_type <- function(change.type = "none", context = c("taxa", "alpha")) {
  context <- match.arg(context)
  mapping <- list(
    none = NULL,
    relative = .CHANGE_RELATIVE,
    log_fold = .CHANGE_LOG_FOLD,
    absolute = .CHANGE_ABSOLUTE
  )

  func <- mapping[[change.type]]
  if (is.null(func) && change.type != "none") {
    warning("Unknown change.type '", change.type, "', using 'relative change'")
    func <- .CHANGE_RELATIVE
  }

  if (context == "alpha" && identical(change.type, "relative")) {
    warning("Alpha diversity change plots do not support relative change; using 'absolute change' instead.")
    func <- .CHANGE_ABSOLUTE
  }

  result <- list(is_change_plot = change.type != "none")
  result[[if (context == "alpha") "alpha.change.func" else "feature.change.func"]] <- func
  return(result)
}


#' Resolve All Parameters for Target Function
#'
#' Main parameter resolution function that prepares all parameters for the target function.
#'
#' @param args List of user-provided arguments
#' @param target_func Name of target function
#' @param design_info Study design information from detect_study_design()
#'
#' @return List of resolved parameters ready for do.call()
#'
#' @keywords internal
resolve_params <- function(args, target_func, design_info) {
  # Start with passed arguments
 resolved <- args

  # Determine function category
  is_barplot <- grepl("barplot", target_func)
  is_areaplot <- grepl("areaplot", target_func)

  # Resolve theme
  if (!is.null(args$theme)) {
    theme_params <- resolve_theme(args$theme)
    resolved$base.size <- theme_params$base.size
    resolved$theme.choice <- theme_params$theme.choice
    resolved$custom.theme <- theme_params$custom.theme
    resolved$palette <- theme_params$palette
    resolved$theme <- NULL  # Remove unified param
  }

  # Resolve feature selection
  if (!is.null(args$feature.select) || is_barplot || is_areaplot) {
    target_type <- if (is_barplot) "barplot" else if (is_areaplot) "areaplot" else "default"
    feature_params <- resolve_feature_select(args$feature.select, target_type)

    if (!is.null(feature_params$feature.number)) {
      resolved$feature.number <- feature_params$feature.number
    }
    if (!is.null(feature_params$top.k.plot)) {
      resolved$top.k.plot <- feature_params$top.k.plot
    }
    if (!is.null(feature_params$features.plot)) {
      resolved$features.plot <- feature_params$features.plot
    }
    resolved$feature.select <- NULL  # Remove unified param
  }

  # Resolve change type
  if (!is.null(args$change.type)) {
    change_context <- if (grepl("^generate_alpha_", target_func)) "alpha" else "taxa"
    change_params <- resolve_change_type(args$change.type, context = change_context)
    if (!is.null(change_params$feature.change.func)) {
      resolved$feature.change.func <- change_params$feature.change.func
    }
    if (!is.null(change_params$alpha.change.func)) {
      resolved$alpha.change.func <- change_params$alpha.change.func
    }
    resolved$change.type <- NULL  # Remove unified param
  }

  # Add time point parameters from design info
  if (design_info$design == "single" && !is.null(design_info$t.level)) {
    resolved$t.level <- design_info$t.level
  } else if (design_info$design %in% c("pair", "long")) {
    # Pair-specific change functions use change.base instead of t0.level/ts.levels
    uses_change_base <- grepl("_pair$", target_func) &&
                        (grepl("change", target_func) || grepl("alpha_test_pair", target_func))

    # Functions that don't need t0.level/ts.levels at all (trend/volatility tests)
    no_time_params <- grepl("_trend_test_", target_func) ||
                      grepl("_volatility_test_", target_func)

    if (uses_change_base) {
      resolved$change.base <- design_info$t0.level
      # Don't add t0.level/ts.levels for change functions
    } else if (!no_time_params) {
      resolved$t0.level <- design_info$t0.level
      resolved$ts.levels <- design_info$ts.levels
    }
  }

  # Remove time.points (unified param)
  resolved$time.points <- NULL

  # Only add pdf = FALSE for plot functions (not test functions)
  if (is.null(resolved$pdf) &&
      !grepl("_test_", target_func) && !grepl("_trend_", target_func) &&
      !grepl("_volatility_", target_func) && !grepl("_association_", target_func)) {
    resolved$pdf <- FALSE
  }

  return(resolved)
}


# =============================================================================
# Function Routing
# =============================================================================

#' Get Available Function Variants
#'
#' Returns information about available function variants for a given category and type.
#'
#' @return Named list of function patterns
#'
#' @keywords internal
get_function_registry <- function() {
  list(
    taxa = list(
      barplot = c(single = "generate_taxa_barplot_single",
                  long = "generate_taxa_barplot_long",
                  pair = "generate_taxa_barplot_pair"),
      areaplot = c(long = "generate_taxa_areaplot_long"),
      boxplot = c(single = "generate_taxa_boxplot_single",
                  long = "generate_taxa_boxplot_long"),
      heatmap = c(single = "generate_taxa_heatmap_single",
                  long = "generate_taxa_heatmap_long",
                  pair = "generate_taxa_heatmap_pair"),
      dotplot = c(single = "generate_taxa_dotplot_single",
                  pair = "generate_taxa_dotplot_pair"),
      spaghettiplot = c(long = "generate_taxa_spaghettiplot_long"),
      cladogram = c(single = "generate_taxa_cladogram_single"),
      # Change variants
      change_boxplot = c(pair = "generate_taxa_change_boxplot_pair"),
      change_heatmap = c(long = "generate_taxa_change_heatmap_long",
                         pair = "generate_taxa_change_heatmap_pair"),
      change_dotplot = c(pair = "generate_taxa_change_dotplot_pair"),
      change_scatterplot = c(pair = "generate_taxa_change_scatterplot_pair")
    ),
    alpha = list(
      boxplot = c(single = "generate_alpha_boxplot_single",
                  long = "generate_alpha_boxplot_long"),
      spaghettiplot = c(long = "generate_alpha_spaghettiplot_long"),
      dotplot = c(long = "generate_alpha_per_time_dotplot_long"),
      change_boxplot = c(pair = "generate_alpha_change_boxplot_pair")
    ),
    beta = list(
      ordination = c(single = "generate_beta_ordination_single",
                     long = "generate_beta_ordination_long",
                     pair = "generate_beta_ordination_pair"),
      boxplot = c(long = "generate_beta_change_boxplot_long",
                  pair = "generate_beta_change_boxplot_pair"),
      spaghettiplot = c(long = "generate_beta_change_spaghettiplot_long"),
      dotplot = c(long = "generate_beta_per_time_dotplot_long"),
      pc_boxplot = c(long = "generate_beta_pc_boxplot_long",
                     pair = "generate_beta_pc_change_boxplot_pair"),
      pc_spaghettiplot = c(long = "generate_beta_pc_spaghettiplot_long")
    ),
    # Test functions
    taxa_test = list(
      difference = c(single = "generate_taxa_test_single",
                     pair = "generate_taxa_test_pair"),
      trend = c(long = "generate_taxa_trend_test_long"),
      volatility = c(long = "generate_taxa_volatility_test_long"),
      association = c(long = "generate_taxa_association_test_long"),
      per_time = c(long = "generate_taxa_per_time_test_long")
    ),
    beta_test = list(
      difference = c(single = "generate_beta_test_single",
                     long = "generate_beta_change_per_time_test_long",
                     pair = "generate_beta_change_test_pair"),
      trend = c(long = "generate_beta_trend_test_long"),
      volatility = c(long = "generate_beta_volatility_test_long")
    ),
    alpha_test = list(
      difference = c(single = "generate_alpha_test_single",
                     long = "generate_alpha_change_per_time_test_long",
                     pair = "generate_alpha_change_test_pair"),
      trend = c(long = "generate_alpha_trend_test_long"),
      volatility = c(long = "generate_alpha_volatility_test_long"),
      per_time = c(long = "generate_alpha_per_time_test_long")
    )
  )
}


#' Route to Appropriate Function
#'
#' Determines the correct underlying function based on category, plot type, and study design.
#'
#' @param category Analysis category: "taxa", "alpha", "beta", "taxa_test", "beta_test", "alpha_test"
#' @param plot_type Type of visualization or test
#' @param design Study design: "single", "long", "pair"
#' @param is_change Whether this is a change/delta analysis
#'
#' @return Function name as character, or NULL if not found
#'
#' @keywords internal
route_function <- function(category, plot_type, design, is_change = FALSE) {
  registry <- get_function_registry()

  # Handle change plots
  if (is_change && !grepl("^change_", plot_type)) {
    change_type <- paste0("change_", plot_type)
    if (!is.null(registry[[category]][[change_type]])) {
      plot_type <- change_type
    }
  }

  # Look up in registry
  if (is.null(registry[[category]])) {
    return(NULL)
  }

  type_variants <- registry[[category]][[plot_type]]
  if (is.null(type_variants)) {
    return(NULL)
  }

  if (!design %in% names(type_variants)) {
    return(NULL)
  }

  type_variants[[design]]
}


# =============================================================================
# Input Validation
# =============================================================================

#' Validate Input Parameters
#'
#' Performs comprehensive validation of input parameters before routing to underlying functions.
#'
#' @param data.obj MicrobiomeStat data object
#' @param category Analysis category
#' @param plot_type Plot or test type
#' @param subject.var Subject variable
#' @param time.var Time variable
#' @param group.var Group variable
#' @param strata.var Stratification variable
#' @param design_info Study design info
#'
#' @return List with is_valid (logical), errors (character vector), warnings (character vector)
#'
#' @keywords internal
validate_inputs <- function(data.obj,
                            category,
                            plot_type,
                            subject.var = NULL,
                            time.var = NULL,
                            group.var = NULL,
                            strata.var = NULL,
                            design_info = NULL) {
  errors <- character()
  warnings <- character()

  # Check data.obj structure
  if (is.null(data.obj)) {
    errors <- c(errors, "data.obj is required")
  } else {
    if (!is.list(data.obj)) {
      errors <- c(errors, "data.obj must be a list")
    } else {
      if (is.null(data.obj$feature.tab)) {
        errors <- c(errors, "data.obj must contain feature.tab")
      }
      if (is.null(data.obj$meta.dat)) {
        errors <- c(errors, "data.obj must contain meta.dat")
      }
    }
  }

  # If basic validation failed, return early
  if (length(errors) > 0) {
    return(list(is_valid = FALSE, errors = errors, warnings = warnings))
  }

  # Check metadata variables exist
  meta_cols <- colnames(data.obj$meta.dat)

  if (!is.null(subject.var) && !subject.var %in% meta_cols) {
    errors <- c(errors, paste0("subject.var '", subject.var, "' not found in meta.dat"))
  }

  if (!is.null(time.var) && !time.var %in% meta_cols) {
    errors <- c(errors, paste0("time.var '", time.var, "' not found in meta.dat"))
  }

  if (!is.null(group.var) && !group.var %in% meta_cols) {
    errors <- c(errors, paste0("group.var '", group.var, "' not found in meta.dat"))
  }

  if (!is.null(strata.var) && !strata.var %in% meta_cols) {
    errors <- c(errors, paste0("strata.var '", strata.var, "' not found in meta.dat"))
  }

  requires_group_var <-
    category %in% c("taxa_test", "beta_test", "alpha_test") &&
    plot_type %in% c("difference", "per_time", "association", "volatility")
  if (requires_group_var && is.null(group.var)) {
    errors <- c(errors, paste0("group.var is required for ", category, "/", plot_type))
  }

  # Check design compatibility
  if (!is.null(design_info)) {
    if (design_info$design %in% c("long", "pair") && is.null(subject.var)) {
      errors <- c(errors, "subject.var is required for longitudinal/paired designs")
    }

    if (design_info$design %in% c("long", "pair") && is.null(time.var)) {
      errors <- c(errors, "time.var is required for longitudinal/paired designs")
    }
  }

  # Check if target function exists
  requested_design <- if (!is.null(design_info)) design_info$design else "single"
  func_name <- route_function(category, plot_type, requested_design)
  if (is.null(func_name)) {
    registry <- get_function_registry()
    supported_designs <- names(registry[[category]][[plot_type]])
    errors <- c(
      errors,
      paste0(
        "No function available for ", category, "/", plot_type,
        " with ", requested_design, " design. Supported designs: ",
        paste(supported_designs, collapse = ", ")
      )
    )
  } else if (!exists(func_name, mode = "function")) {
    errors <- c(errors, paste0("Function '", func_name, "' not found"))
  }

  return(list(
    is_valid = length(errors) == 0,
    errors = errors,
    warnings = warnings,
    target_func = func_name
  ))
}


#' Clean Resolved Parameters for Function Call
#'
#' Removes NULL values from resolved parameters while preserving
#' key parameters that underlying functions may require even when NULL.
#' For single design functions, subject.var is removed since those functions
#' no longer accept this parameter.
#'
#' @param resolved Named list of resolved parameters
#' @param design Study design: "single", "pair", or "long"
#'
#' @return Cleaned parameter list
#'
#' @keywords internal
clean_resolved_params <- function(resolved, design = "single") {
  # Parameters that should be passed even when NULL
  # (some underlying functions have these as required params without defaults)
  keep_if_null <- c("time.var", "group.var", "strata.var",
                    "adj.vars", "t.level", "t0.level", "ts.levels")

  # For pair and long designs, subject.var is still required
  if (design %in% c("pair", "long")) {
    keep_if_null <- c("subject.var", keep_if_null)
  } else {
    # For single design, remove subject.var entirely (functions don't accept it)
    resolved$subject.var <- NULL
  }

  resolved[!sapply(names(resolved), function(n) {
    is.null(resolved[[n]]) && !n %in% keep_if_null
  })]
}


#' Drop Unsupported Parameters for Target Function
#'
#' Removes arguments that the target function cannot accept unless it
#' explicitly declares `...`.
#'
#' @param resolved Named list of resolved parameters
#' @param target_func Function name as character
#'
#' @return Filtered parameter list
#'
#' @keywords internal
filter_target_params <- function(resolved, target_func = NULL) {
  if (is.null(target_func) || !exists(target_func, mode = "function")) {
    return(resolved)
  }

  target_formals <- names(formals(get(target_func, mode = "function")))
  if ("..." %in% target_formals) {
    return(resolved)
  }

  resolved[names(resolved) %in% target_formals]
}


#' Format Validation Errors for Display
#'
#' Creates user-friendly error messages from validation results.
#'
#' @param validation Validation result from validate_inputs()
#' @param func_name Name of the calling function
#'
#' @return Formatted error message string
#'
#' @keywords internal
format_validation_errors <- function(validation, func_name = "plot function") {
  if (validation$is_valid) {
    return(NULL)
  }

  msg <- paste0("\n", func_name, " validation failed:\n")

  for (err in validation$errors) {
    msg <- paste0(msg, "  - ", err, "\n")
  }

  if (length(validation$warnings) > 0) {
    msg <- paste0(msg, "\nWarnings:\n")
    for (warn in validation$warnings) {
      msg <- paste0(msg, "  - ", warn, "\n")
    }
  }

  return(msg)
}
