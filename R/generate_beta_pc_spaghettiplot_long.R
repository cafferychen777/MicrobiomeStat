#' Generate spaghetti plots of beta diversity over time
#'
#' This function generates spaghetti plots of beta diversity metrics over time for longitudinal microbiome studies.
#' It allows plotting principal coordinates from different dimension reduction methods such as PCoA, NMDS, tSNE or UMAP.
#' By default it uses PCoA for dimension reduction.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param pc.obj A list containing the results of dimension reduction/Principal Component Analysis.
#' This should be the output from functions like \code{\link[MicrobiomeStat]{mStat_calculate_PC}},
#' containing the PC coordinates and other metadata. If NULL (default), dimension reduction
#' will be automatically performed using metric multidimensional scaling (MDS) via
#' \code{\link[MicrobiomeStat]{mStat_calculate_PC}}. The pc.obj list structure should contain:
#' \describe{
#'   \item{points}{A matrix with samples as rows and PCs as columns containing the coordinates.}
#'   \item{eig}{Eigenvalues for each PC dimension.}
#'   \item{vectors}{Loadings vectors for features onto each PC.}
#'   \item{Other metadata}{like method, dist.name, etc.}
#' }
#' See \code{\link[MicrobiomeStat]{mStat_calculate_PC}} function for details on output format.
#' @param pc.ind Numeric vector indicating which Principal Coordinates to plot.
#' @param subject.var Character string indicating the variable for subject identifiers.
#' @param time.var Character string indicating the variable for time points.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var Character string indicating the variable for group identifiers.
#' @param strata.var Character string indicating the variable for stratum identifiers.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param dist.name Vector of character strings indicating the distance measures to include in the plot.
#' @param base.size Numeric value indicating the base font size for the plot.
#' @param theme.choice
#' Plot theme choice. Specifies the visual style of the plot. Can be one of the following pre-defined themes:
#'   - "prism": Utilizes the ggprism::theme_prism() function from the ggprism package, offering a polished and visually appealing style.
#'   - "classic": Applies theme_classic() from ggplot2, providing a clean and traditional look with minimal styling.
#'   - "gray": Uses theme_gray() from ggplot2, which offers a simple and modern look with a light gray background.
#'   - "bw": Employs theme_bw() from ggplot2, creating a classic black and white plot, ideal for formal publications and situations where color is best minimized.
#'   - "light": Implements theme_light() from ggplot2, featuring a light theme with subtle grey lines and axes, suitable for a fresh, modern look.
#'   - "dark": Uses theme_dark() from ggplot2, offering a dark background, ideal for presentations or situations where a high-contrast theme is desired.
#'   - "minimal": Applies theme_minimal() from ggplot2, providing a minimalist theme with the least amount of background annotations and colors.
#'   - "void": Employs theme_void() from ggplot2, creating a blank canvas with no axes, gridlines, or background, ideal for custom, creative plots.
#' Each theme option adjusts various elements like background color, grid lines, and font styles to match the specified aesthetic.
#' Default is "bw", offering a universally compatible black and white theme suitable for a wide range of applications.
#' @param custom.theme
#' A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. Custom themes are useful for creating publication-ready figures with specific formatting requirements.
#'
#' To use a custom theme, create a theme object with ggplot2::theme(), including any desired customizations. Common customizations for publication-ready figures might include adjusting text size for readability, altering line sizes for clarity, and repositioning or formatting the legend. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=14, face="bold"),        # Bold axis titles with larger font
#'   axis.text = ggplot2::element_text(size=12),                      # Slightly larger axis text
#'   legend.position = "top",                                         # Move legend to the top
#'   legend.background = ggplot2::element_rect(fill="lightgray"),     # Light gray background for legend
#'   panel.background = ggplot2::element_rect(fill="white", colour="black"), # White panel background with black border
#'   panel.grid.major = ggplot2::element_line(colour = "grey90"),     # Lighter color for major grid lines
#'   panel.grid.minor = ggplot2::element_blank(),                     # Remove minor grid lines
#'   plot.title = ggplot2::element_text(size=16, hjust=0.5)           # Centered plot title with larger font
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. If `custom.theme` is NULL (the default), the theme is determined by `theme.choice`. This flexibility allows for both easy theme selection for general use and detailed customization for specific presentation or publication needs.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf Logical indicating whether to save the plot as a PDF file.
#' @param file.ann Optional character string to append to the file name.
#' @param pdf.wid Numeric value indicating the width of the PDF file.
#' @param pdf.hei Numeric value indicating the height of the PDF file.
#' @param ... Additional arguments to be passed to the plotting function.
#'
#' @return A list of ggplot objects for each distance measure and Principal Coordinate.
#' @examples
#' \dontrun{
#'
#' # Load required libraries and data
#' library(vegan)
#'
#' # Example with ecam.obj dataset
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, "BC")
#' pc.obj <- mStat_calculate_PC(dist.obj)
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:10],
#'   group.var = "diet",
#'   strata.var = "delivery",
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:10],
#'   group.var = "diet",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' # Example with peerj32.obj dataset
#' data(peerj32.obj)
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_beta_pc_spaghettiplot_long <- function(data.obj = NULL,
                                                dist.obj = NULL,
                                                pc.obj = NULL,
                                                pc.ind = c(1, 2),
                                                subject.var,
                                                time.var,
                                                t0.level = NULL,
                                                ts.levels = NULL,
                                                group.var = NULL,
                                                strata.var = NULL,
                                                adj.vars = NULL,
                                                dist.name = c("BC"),
                                                base.size = 16,
                                                theme.choice = "bw",
                                                custom.theme = NULL,
                                                palette = NULL,
                                                pdf = TRUE,
                                                file.ann = NULL,
                                                pdf.wid = 11,
                                                pdf.hei = 8.5,
                                                ...) {

  # Check if distance metrics are provided
  if (is.null(dist.name)){
    return()
  }

  # If distance object is not provided, calculate it from the data object
  if (is.null(dist.obj)) {
    # Process time variable and extract relevant metadata
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      )))
    
    # Calculate beta diversity
    dist.obj <-
      mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
    
    # If adjustment variables are provided, calculate adjusted distances
    if (!is.null(adj.vars)) {
      dist.obj <-
        mStat_calculate_adjusted_distance(
          data.obj = data.obj,
          dist.obj = dist.obj,
          adj.vars = adj.vars,
          dist.name = dist.name
        )
    }
  } else {
    # If data object is provided with metadata, process time variable
    if (!is.null(data.obj) & !is.null(data.obj$meta.dat)) {
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <-
        data.obj$meta.dat %>% select(all_of(c(
          subject.var, time.var, group.var, strata.var
        )))
    } else {
      # If no data object, extract metadata from distance object
      meta_tab <-
        attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(
          subject.var, time.var, group.var, strata.var
        )))
      data.obj <- list(meta.dat = meta_tab)
      data.obj <-
        mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat
      dist.obj <- mStat_subset_dist(dist.obj, colnames(meta_tab))
    }
  }

  # Get color palette
  col <- mStat_get_palette(palette)

  # Calculate new sizes based on base.size for consistent plot aesthetics
  title.size = base.size * 1.25
  axis.title.size = base.size * 0.75
  axis.text.size = base.size * 0.5
  legend.title.size = base.size * 1
  legend.text.size = base.size * 0.75

  # Get the appropriate theme
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Create a list to store plots for each distance metric
  plot_list <- lapply(dist.name, function(dist.name) {
    # If principal component object is not provided, calculate it using MDS
    if (is.null(pc.obj)) {
      message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
      message(
        "If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method."
      )
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj[dist.name],
          method = "mds",
          k = max(pc.ind),
          dist.name = dist.name
        )
    }

    # Extract principal component coordinates
    pc.mat <- pc.obj[[dist.name]]$points

    colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

    # Ensure PC matrix rows match metadata
    pc.mat <- pc.mat[rownames(meta_tab[, c(subject.var, time.var, group.var, strata.var)]),]

    pc.mat <- pc.mat %>% as_tibble()

    # Combine PC coordinates with metadata
    df <-
      cbind(pc.mat[, paste0("PC", pc.ind)], meta_tab[, c(subject.var, time.var, group.var, strata.var)])

    # Reshape data from wide to long format
    df <-
      df %>%
      as_tibble() %>%
      tidyr::gather(
        key = "PC",
        value = "value",
        -one_of(subject.var, group.var, time.var, strata.var)
      )

    # If no group variable is provided, create a dummy group
    if (is.null(group.var)){
      df <- df %>% dplyr::mutate("ALL" = "ALL")
      group.var = "ALL"
    }

    # Calculate mean values for each group and time point
    if (is.null(strata.var)) {
      df.mean <- df %>%
        dplyr::group_by(!!sym(time.var),!!sym(group.var), PC) %>%
        dplyr::summarize(mean_value = mean(value, na.rm = TRUE))
      df <-
        dplyr::left_join(df, df.mean, by = c(time.var, group.var, "PC"))
    } else {
      df.mean <- df %>%
        dplyr::group_by(!!sym(time.var),!!sym(group.var),!!sym(strata.var), PC) %>%
        dplyr::summarize(mean_value = mean(value, na.rm = TRUE))
      df <-
        dplyr::left_join(df, df.mean, by = c(time.var, group.var, strata.var, "PC"))
    }

    # Create a list to store plots for each principal component
    sub_plot_list <- lapply(unique(df$PC), function(pc.index) {
      sub_df <- df %>% filter(PC == pc.index)

      # Create the spaghetti plot
      p <- ggplot() +
        geom_point(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "value",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.3,
          size = 3
        ) +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_value",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        geom_point(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_value",
            group = group.var,
            color = group.var
          ),
          size = 5
        ) +
        labs(x = time.var, y = paste("Distance:",
                                     dist.name,
                                     " - Axis",
                                     gsub("PC", "", pc.index)), color = group.var) +
        scale_color_manual(
          values = col
        ) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = axis.title.size*2),
          axis.title.y = element_text(size = axis.title.size*2),
          axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size * 2),
          axis.text.y = element_text(size = base.size*2),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16 * 2),
          legend.title = ggplot2::element_text(size = 16 * 2),
          legend.key.size = unit(10, "mm"),
          legend.key.spacing = unit(2, "mm")
        )

      # Remove legend if there's only one group
      if (group.var == "ALL"){
        p <- p + theme(legend.position = "none")
      }

      # Add faceting if strata variable is provided
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(
          cols = vars(!!sym(strata.var)),
          scale = "free",
          space = "free"
        )
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "beta_pc_spaghettiplot_long_",
          dist.name,
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var
        )
        if (!is.null(group.var)) {
          pdf_name <- paste0(pdf_name, "_", "group_", group.var)
        }
        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
        }
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }
        pdf_name <- paste0(pdf_name, ".pdf")
        ggsave(
          filename = pdf_name,
          plot = p,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      }
      return(p)
    })

    names(sub_plot_list) <- unique(df$PC)
    return(sub_plot_list)
  })

  names(plot_list) <- dist.name

  return(plot_list)
}