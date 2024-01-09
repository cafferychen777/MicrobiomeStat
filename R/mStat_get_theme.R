#' Get the ggplot2 Theme for Plots
#'
#' This function returns a ggplot2 theme based on the user's choice or a custom theme if provided.
#'
#' @param theme.choice A character string specifying the pre-defined theme choice.
#' Supported themes are "prism", "classic", "gray", and "bw". If an unsupported theme is provided,
#' the "bw" theme will be used as the default.
#' @param custom.theme An optional ggplot2 theme object to be used instead of the pre-defined themes.
#'
#' @return A ggplot2 theme object to be used in plotting functions.
#'
#' @details
#' The function simplifies the process of theme selection in plotting functions
#' by providing a standardized way to choose between several common themes or
#' to use a custom theme.
#'
#' @examples
#' \dontrun{
#' # Using a pre-defined theme
#' theme_to_use <- mStat_get_theme(theme.choice = "bw")
#'
#' # Using a custom theme
#' my_custom_theme <- theme_minimal() + theme(axis.text.x = element_text(angle = 45))
#' theme_to_use <- mStat_get_theme(custom.theme = my_custom_theme)
#' }
#'
#' @export
mStat_get_theme <- function(theme.choice = "bw", custom.theme = NULL) {

  if (!is.null(custom.theme)) {
    return(custom.theme)
  }

  theme_function <- switch(
    theme.choice,
    prism = ggprism::theme_prism(),   # Prism theme from ggprism package
    classic = theme_classic(),        # Classic theme
    gray = theme_gray(),              # Gray theme
    bw = theme_bw(),                  # Black and white theme
    light = theme_light(),            # Light theme
    dark = theme_dark(),              # Dark theme
    minimal = theme_minimal(),        # Minimal theme
    void = theme_void(),              # Void theme (no axes, labels, etc.)
    theme_bw()                        # Default case to black and white theme
  )

  return(theme_function)
}
