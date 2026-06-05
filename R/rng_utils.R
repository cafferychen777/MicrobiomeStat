#' Run an expression under a temporary RNG seed
#'
#' Sets the random seed to \code{seed}, evaluates \code{expr}, then restores
#' the previous global RNG state (or removes \code{.Random.seed} if none
#' existed). This makes stochastic operations reproducible without affecting
#' the caller's random stream.
#'
#' @param seed Integer seed passed to \code{\link{set.seed}}.
#' @param expr Expression to evaluate.
#' @return The value of \code{expr}.
#' @noRd
mStat_with_local_seed <- function(seed, expr) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expr)
}
