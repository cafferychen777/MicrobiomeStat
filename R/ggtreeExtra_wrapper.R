#' Wrapper for ggtreeExtra::geom_fruit to handle NSE correctly
#'
#' @param data Data to use for the layer
#' @param geom Name of the geom function to use (e.g., geom_tile)
#' @param mapping Aesthetic mapping
#' @param ... Additional arguments to pass to ggtreeExtra::geom_fruit
#'
#' @return A layer to be added to a ggtree plot
#' @export
#'
#' @keywords internal
mStat_geom_fruit <- function(data, geom = NULL, mapping = NULL, ...) {
  # Make sure required packages are available
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    stop("Package 'ggtree' is required but not installed.")
  }
  if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
    stop("Package 'ggtreeExtra' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  # Load ggplot2 locally to ensure geom_* functions are available
  # This is a standard practice for wrappers dealing with NSE
  if (is.null(geom)) {
    # If no geom is provided, use geom_tile by default
    geom <- ggplot2::geom_tile
  }
  
  # Call the actual function with proper NSE handling
  ggtreeExtra::geom_fruit(
    data = data,
    geom = geom,
    mapping = mapping,
    ...
  )
}
