#' Retrieve Custom Color Palettes for Data Visualization
#'
#' This function provides a selection of color palettes inspired by various scientific publications,
#' particularly the `ggsci` package, which offers a wide range of color palettes used in scientific journals and publications.
#' It allows users to choose a predefined palette from `ggsci` or specify their own set of colors.
#'
#' @param palette A character string specifying the name of a predefined palette from `ggsci` or
#'                a vector of color codes. If `NULL`, a default palette is used.
#'                Available predefined palettes include 'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb',
#'                all derived from the `ggsci` package.
#'
#' @return A vector of color codes. If a predefined palette name is provided, the corresponding
#'         palette from the `ggsci` package is returned. If a vector of color codes is provided, it is returned as-is.
#'         If the input is NULL or an unrecognized string, a default palette is returned.
#'
#' @examples
#' # Default palette
#' mStat_get_palette()
#'
#' # Predefined palette from `ggsci`: 'nejm'
#' mStat_get_palette("nejm")
#'
#' # Other predefined palettes
#' # 'npg' palette
#' mStat_get_palette("npg")
#'
#' # 'aaas' palette
#' mStat_get_palette("aaas")
#'
#' # 'lancet' palette
#' mStat_get_palette("lancet")
#'
#' # 'jama' palette
#' mStat_get_palette("jama")
#'
#' # 'jco' palette
#' mStat_get_palette("jco")
#'
#' # 'ucscgb' palette
#' mStat_get_palette("ucscgb")
#'
#' # Custom color vector
#' mStat_get_palette(c("#123456", "#654321"))
#'
#' # Unrecognized input returns default palette
#' mStat_get_palette("unknown")
#'
#' # Using a longer custom color vector
#' mStat_get_palette(c("#FF5733", "#33FF57", "#3357FF", "#FF33FF"))
#'
#' # Using an empty string (should return default palette)
#' mStat_get_palette("")
#'
#' # Using a single color code (treated as custom vector)
#' mStat_get_palette("#FF5733")
#' @export
mStat_get_palette <- function(palette = NULL) {
  # Define the custom color palettes
  ggsci_db <- list(
    npg = list(default = c(
      "Cinnabar" = "#E64B35", "Shakespeare" = "#4DBBD5",
      "PersianGreen" = "#00A087", "Chambray" = "#3C5488",
      "Apricot" = "#F39B7F", "WildBlueYonder" = "#8491B4",
      "MonteCarlo" = "#91D1C2", "Monza" = "#DC0000",
      "RomanCoffee" = "#7E6148", "Sandrift" = "#B09C85"
    )),
    aaas = list(default = c(
      "Chambray" = "#3B4992", "Red" = "#EE0000",
      "FunGreen" = "#008B45", "HoneyFlower" = "#631879",
      "Teal" = "#008280", "Monza" = "#BB0021",
      "ButterflyBush" = "#5F559B", "FreshEggplant" = "#A20056",
      "Stack" = "#808180", "CodGray" = "#1B1919"
    )),
    nejm = list(default = c(
      "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
      "Zest" = "#E18727", "Eucalyptus" = "#20854E",
      "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
      "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"
    )),
    lancet = list(default= c(
      "CongressBlue" = "#00468B", "Red" = "#ED0000",
      "Apple" = "#42B540", "BondiBlue" = "#0099B4",
      "TrendyPink" = "#925E9F", "MonaLisa" = "#FDAF91",
      "Carmine" = "#AD002A", "Edward" = "#ADB6B6",
      "CodGray" = "#1B1919"
    )),
    jama = list(default = c(
      "Limed Spruce" = "#374E55", "Anzac" = "#DF8F44",
      "Cerulean" = "#00A1D5", "Apple Blossom" = "#B24745",
      "Acapulco" = "#79AF97", "Kimberly" = "#6A6599",
      "Makara" = "#80796B"
    )),
    jco = list(default = c(
      "Lochmara" = "#0073C2", "Corn" = "#EFC000",
      "Gray" = "#868686", "ChestnutRose" = "#CD534C",
      "Danube" = "#7AA6DC", "RegalBlue" = "#003C67",
      "Olive" = "#8F7700", "MineShaft" = "#3B3B3B",
      "WellRead" = "#A73030", "KashmirBlue" = "#4A6990"
    )),
    ucscgb = list(default = c(
      "chr5" = "#FF0000", "chr8" = "#FF9900", "chr9" = "#FFCC00",
      "chr12" = "#00FF00", "chr15" = "#6699FF", "chr20" = "#CC33FF",
      "chr3" = "#99991E", "chrX" = "#999999", "chr6" = "#FF00CC",
      "chr4" = "#CC0000", "chr7" = "#FFCCCC", "chr10" = "#FFFF00",
      "chr11" = "#CCFF00", "chr13" = "#358000", "chr14" = "#0000CC",
      "chr16" = "#99CCFF", "chr17" = "#00FFFF", "chr18" = "#CCFFFF",
      "chr19" = "#9900CC", "chr21" = "#CC99FF", "chr1" = "#996600",
      "chr2" = "#666600", "chr22" = "#666666", "chrY" = "#CCCCCC",
      "chrUn" = "#79CC3D", "chrM" = "#CCCC99"
    ))
  )

  # Default palette
  default_palette <-
    c("#E31A1C",
    "#1F78B4",
    "#FB9A99",
    "#33A02C",
    "#FDBF6F",
    "#B2DF8A",
    "#A6CEE3",
    "#BA7A70",
    "#9D4E3F",
    "#829BAB"
  )

  # Determine which palette to return
  if (is.null(palette)) {
    return(default_palette)
  } else if (is.character(palette) && length(palette) == 1) {
    if (palette %in% names(ggsci_db)) {
      return(as.vector(ggsci_db[[palette]]$default))
    } else {
      return(default_palette)
    }
  } else if (is.vector(palette)) {
    return(palette)
  } else {
    return(default_palette)
  }

}
