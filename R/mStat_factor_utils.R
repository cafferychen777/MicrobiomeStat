# Internal separator used by interaction()/separate() pairing for
# group.var and strata.var combination.  Using a multi-character token
# avoids false splits when level names contain dots or underscores.
.STRATA_SEP <- ":::"

#' Capture and Enforce Factor Levels in Metadata
#'
#' Reads factor levels from \code{data.obj$meta.dat} for the requested
#' variables, converts character columns to factors (preserving the
#' first-occurrence order), and returns both the updated data object
#' and the captured levels for later restoration.
#'
#' @param data.obj A MicrobiomeStat data object whose \code{meta.dat}
#'   component will be modified in place.
#' @param group.var,strata.var Optional variable names (strings or NULL).
#' @return A list with two elements:
#'   \describe{
#'     \item{data.obj}{The (possibly modified) data object.}
#'     \item{levels}{A named list with elements \code{group} and
#'       \code{strata}, each containing the captured level vector
#'       or NULL.}
#'   }
#' @noRd
mStat_capture_factor_levels <- function(data.obj,
                                        group.var = NULL,
                                        strata.var = NULL) {
  captured <- list(group = NULL, strata = NULL)
  vars     <- list(list(var = group.var,  key = "group"),
                   list(var = strata.var, key = "strata"))

  for (v in vars) {
    if (!is.null(v$var)) {
      col <- data.obj$meta.dat[[v$var]]
      lvls <- if (is.factor(col)) levels(col) else unique(col)
      data.obj$meta.dat[[v$var]] <- factor(col, levels = lvls)
      captured[[v$key]] <- lvls
    }
  }

  list(data.obj = data.obj, levels = captured)
}

#' Restore Factor Levels After \code{tidyr::separate()}
#'
#' \code{tidyr::separate()} converts its output columns to character.
#' This function re-applies the original factor levels captured by
#' \code{mStat_capture_factor_levels()}.
#'
#' @param df A data frame whose columns need factor restoration.
#' @param levels_info The \code{levels} element returned by
#'   \code{mStat_capture_factor_levels()}.
#' @param group.col Column name to receive group levels.
#'   Defaults to the original variable name, but callers may pass
#'   e.g. \code{paste0(group.var, "2")} when the separated column
#'   has been renamed.
#' @param strata.col Column name to receive strata levels.
#' @return The data frame with factor columns restored.
#' @noRd
mStat_restore_factor_levels <- function(df,
                                        levels_info,
                                        group.col = NULL,
                                        strata.col = NULL) {
  pairs <- list(list(col = group.col,  lvls = levels_info$group),
                list(col = strata.col, lvls = levels_info$strata))

  for (p in pairs) {
    if (!is.null(p$col) && !is.null(p$lvls) && p$col %in% names(df)) {
      df[[p$col]] <- factor(df[[p$col]], levels = p$lvls)
    }
  }

  df
}
