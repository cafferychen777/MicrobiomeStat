#' @title Change Calculation Utilities
#' @description Constants and helper functions for computing abundance/diversity
#'   changes between time points. Centralizes the logic that was previously
#'   duplicated across 10+ visualization and testing functions.
#' @name change_utils
#' @keywords internal
NULL

# =============================================================================
# Change Method Constants
# =============================================================================

#' @noRd
.CHANGE_LOG_FOLD  <- "log fold change"

#' @noRd
.CHANGE_RELATIVE  <- "relative change"

#' @noRd
.CHANGE_ABSOLUTE  <- "absolute change"

#' @noRd
.CHANGE_METHODS   <- c(.CHANGE_LOG_FOLD, .CHANGE_RELATIVE, .CHANGE_ABSOLUTE)


# =============================================================================
# Taxa Change Computation
# =============================================================================

#' Compute taxa abundance change between two time points
#'
#' Applies the specified change method element-wise. For log fold change, uses
#' a per-feature half-minimum pseudocount calculated across both time points to
#' avoid bias from asymmetric zero handling.
#'
#' @param value_after Numeric vector of abundances at the later time point.
#' @param value_before Numeric vector of abundances at the earlier time point.
#' @param method A character string (\code{"log fold change"},
#'   \code{"relative change"}, \code{"absolute change"}) or a function
#'   accepting two arguments (after, before).
#' @param feature_id Optional character vector of feature identifiers (same
#'   length as \code{value_after}). Required for per-feature pseudocount
#'   calculation when \code{method = "log fold change"}.
#' @param verbose Logical; if \code{TRUE}, emit a message describing the
#'   zero-handling strategy. Default \code{TRUE}.
#'
#' @return Numeric vector of computed changes, same length as the inputs.
#' @keywords internal
compute_taxa_change <- function(value_after,
                                value_before,
                                method,
                                feature_id = NULL,
                                verbose = TRUE) {

  # Custom function: pass through directly

  if (is.function(method)) {
    return(method(value_after, value_before))
  }

  switch(
    method,
    "log fold change" = {
      va <- value_after
      vb <- value_before

      if (!is.null(feature_id)) {
        # Per-feature half-minimum pseudocount across BOTH time points
        uid <- unique(feature_id)
        for (feat in uid) {
          idx <- feature_id == feat
          all_nonzero <- c(vb[idx & vb > 0], va[idx & va > 0])
          pseudo <- if (length(all_nonzero) > 0) min(all_nonzero) / 2 else 1e-10
          vb[idx & vb == 0] <- pseudo
          va[idx & va == 0] <- pseudo
        }
      } else {
        # Scalar half-minimum across all values (single-feature context)
        all_nonzero <- c(vb[vb > 0], va[va > 0])
        pseudo <- if (length(all_nonzero) > 0) min(all_nonzero) / 2 else 1e-10
        vb[vb == 0] <- pseudo
        va[va == 0] <- pseudo
      }

      if (verbose) {
        message(
          "Zero-handling: Per-taxon half-minimum pseudocount across BOTH time points.\n",
          "  This ensures unbiased log fold change calculations."
        )
      }

      log2(va) - log2(vb)
    },
    "relative change" = {
      dplyr::case_when(
        value_after == 0 & value_before == 0 ~ 0,
        TRUE ~ (value_after - value_before) / (value_after + value_before)
      )
    },
    "absolute change" = {
      value_after - value_before
    },
    # Default fallback
    {
      value_after - value_before
    }
  )
}


# =============================================================================
# Alpha Diversity Change Computation
# =============================================================================

#' Compute alpha diversity change between two time points
#'
#' Alpha diversity indices (Shannon, Simpson, etc.) are always positive, so no
#' zero-value handling is needed. Only supports \code{"log fold change"} and
#' \code{"absolute change"}.
#'
#' @param value_after Numeric vector of diversity values at the later time point.
#' @param value_before Numeric vector of diversity values at the earlier time point.
#' @param method A character string (\code{"log fold change"},
#'   \code{"absolute change"}) or a function accepting two arguments (after, before).
#'
#' @return Numeric vector of computed changes, same length as the inputs.
#' @keywords internal
compute_alpha_change <- function(value_after, value_before, method) {

  if (is.function(method)) {
    return(method(value_after, value_before))
  }

  switch(
    method,
    "log fold change" = {
      log2(value_after / value_before)
    },
    "absolute change" = {
      value_after - value_before
    },
    # Default fallback
    {
      message("No valid change method provided. Defaulting to 'absolute change'.")
      value_after - value_before
    }
  )
}


# =============================================================================
# Change Method Description
# =============================================================================

#' Generate a human-readable description of a change method
#'
#' Used by report-generation functions to describe how changes were calculated.
#'
#' @param method A character string or function specifying the change method.
#' @param context Character string: \code{"taxa"} (default) for taxa abundance
#'   context, \code{"alpha"} for alpha diversity context.
#'
#' @return A single character string describing the method.
#' @keywords internal
describe_change_method <- function(method, context = "taxa") {
  if (is.function(method)) {
    return("The changes were computed using a custom function provided by the user.")
  }

  if (context == "alpha") {
    prefix <- "The changes from t0.level were"
  } else {
    prefix <- "The changes were"
  }

  switch(
    method,
    "relative change" = paste0(
      prefix,
      " relative changes, which were computed as (after.abund - before.abund) / (after.abund + before.abund) so the values lie between [-1, 1]."
    ),
    "absolute change" = paste0(
      prefix,
      " absolute changes, computed as the difference between after.abund and before.abund."
    ),
    "log fold change" = paste0(
      prefix,
      " log2 fold changes, computed as the logarithm of the ratio of after.abund to before.abund, with a per-taxon half-minimum pseudocount to avoid taking the log of zero."
    ),
    paste0(prefix, " computed using method: ", method, ".")
  )
}
