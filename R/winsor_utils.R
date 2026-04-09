#' @title Winsorization Utilities
#' @description Internal helper for LinDA-style winsorization. This is shared
#'   by `linda()` and `linda2()` to keep preprocessing behavior identical.
#' @name winsor_utils
#' @keywords internal
NULL


#' @keywords internal
winsor_feature_table <- function(Y, quan, feature.dat.type) {
  switch(feature.dat.type,
    count = {
      N <- colSums(Y)
      P <- t(t(Y) / N)
      cut <- apply(P, 1, quantile, quan)
      Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
      ind <- P > Cut
      P[ind] <- Cut[ind]
      Y <- round(t(t(P) * N))
    },
    proportion = {
      cut <- apply(Y, 1, quantile, quan)
      Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
      ind <- Y > Cut
      Y[ind] <- Cut[ind]
    }
  )

  Y
}


#' @keywords internal
winsor.fun <- winsor_feature_table


#' @keywords internal
mStat_has_random_effect_term <- function(formula) {
  formula_text <- paste(as.character(formula), collapse = " ")
  grepl("|", formula_text, fixed = TRUE)
}
