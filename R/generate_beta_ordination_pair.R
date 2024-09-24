#' @title Generate beta ordination plots for paired samples
#'
#' @description This function, part of the MicrobiomeStat package, generates Principle Coordinate Analysis (PCoA) plots using the provided distance object. It is specifically designed for the analysis of microbiome data. This function is tailored for paired samples, such as those from longitudinal studies or experiments with multiple time points.
#' @name generate_beta_ordination_pair
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
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param group.var (Optional) The variable in the metadata table that represents the grouping factor.
#' @param strata.var (Optional) The variable in the metadata table that represents the stratification factor.
#'                   When provided, the plot will be faceted based on this variable. Each facet will display
#'                   all data points, with points relevant to the current facet shown in color and points
#'                   from other facets shown in gray. This allows for easy comparison across strata while
#'                   maintaining context of the overall data distribution.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param base.size (Optional) Base font size for the plot (default is 16).
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
#' @param pdf (Optional) A boolean indicating whether to save the output as a PDF file (default is TRUE).
#' @param file.ann (Optional) A string for annotating the output file name.
#' @param pdf.wid (Optional) The width of the PDF file if `pdf` is set to `TRUE` (default is 11).
#' @param pdf.hei (Optional) The height of the PDF file if `pdf` is set to `TRUE` (default is 8.5).
#' @param ... (Optional) Additional arguments to pass to the plotting function.
#'
#' @details The function is flexible and allows for various modifications, including the choice of distance measure and stratification factor, providing a comprehensive tool for microbiome beta diversity exploration. It integrates well with other MicrobiomeStat functions and takes their output as input.
#'
#' @return A PCoA plot displaying the beta diversity ordination, stratified by the specified grouping and/or strata variables (if provided). The plot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object, \code{\link[MicrobiomeStat]{mStat_calculate_PC}} for computing the principal coordinates, and \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{geom_boxplot}} for the underlying plot functions used, and \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_DESeqDataSet_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_phyloseq_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_SummarizedExperiment_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_qiime2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_mothur_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_dada2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_biom_as_data_obj}} for data conversion.
#'
#'
#' @examples
#' \dontrun{
#'
#' # Load necessary libraries and data
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Perform beta ordination pair analysis using `generate_beta_ordination_pair`
#' generate_beta_ordination_pair(
#'   data.obj      = peerj32.obj,
#'   dist.obj      = NULL,
#'   pc.obj        = NULL,
#'   subject.var   = "subject",
#'   time.var      = "time",
#'   group.var     = "group",
#'   strata.var    = "sex",
#'   adj.vars      = "sex",
#'   dist.name     = c("BC"),
#'   base.size     = 16,
#'   theme.choice  = "bw",
#'   custom.theme  = NULL,
#'   palette       = NULL,
#'   pdf           = TRUE,
#'   file.ann      = NULL,
#'   pdf.wid       = 11,
#'   pdf.hei       = 8.5
#' )
#'
#' data(subset_pairs.obj)
#'
#' # Perform beta ordination pair analysis using `generate_beta_ordination_pair`
#' generate_beta_ordination_pair(
#'   data.obj      = subset_pairs.obj,
#'   dist.obj      = NULL,
#'   pc.obj        = NULL,
#'   subject.var   = "MouseID",
#'   time.var      = "Antibiotic",
#'   group.var     = "Sex",
#'   strata.var    = NULL,
#'   adj.vars      = NULL,
#'   dist.name     = c("BC"),
#'   base.size     = 16,
#'   theme.choice  = "bw",
#'   custom.theme  = NULL,
#'   palette       = NULL,
#'   pdf           = TRUE,
#'   file.ann      = NULL,
#'   pdf.wid       = 11,
#'   pdf.hei       = 8.5
#' )
#' }
#' @export
generate_beta_ordination_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           dist.name = c('BC', 'Jaccard'),
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    if (is.null(dist.name)){
      return()
    }

    if (is.null(dist.obj)) {
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
      print(dist.obj)
    } else {
      if (is.null(data.obj)) {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      } else {
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      }
    }

    if (is.null(pc.obj)) {
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = 2,
          dist.name = dist.name
        )
    }

    col <- mStat_get_palette(palette)

    aes_function <- if (!is.null(group.var)) {
      aes(color = !!sym(group.var),
          shape = !!sym(time.var))
    } else {
      aes(color = !!sym(time.var))
    }

    # Assuming mStat_get_theme function is already defined
    # Replace the existing theme selection code with this:
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    plot_list <- lapply(dist.name, function(dist.name) {
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]
      df <- as.data.frame(pc.mat) %>%
        setNames(c("PC1", "PC2")) %>%
        dplyr::bind_cols(meta_tab[, c(subject.var, time.var, group.var, strata.var)]) %>%
        dplyr::mutate(
          x_start = PC1,
          y_start = PC2,
          x_end = NA,
          y_end = NA
        )

      Time_choices <-
        df %>% dplyr::select(all_of(time.var)) %>% dplyr::pull() %>% unique()

      df <- df %>%
        dplyr::arrange(!!sym(subject.var),!!sym(time.var)) %>% # Sort to ensure the correct order of time.
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::mutate(x_end = dplyr::lead(PC1),
                      y_end = dplyr::lead(PC2)) %>%
        dplyr::ungroup()

      # Create a dataset with all points for each facet
      if (!is.null(strata.var)) {
        all_strata <- unique(df[[strata.var]])
        df_all <- do.call(rbind, lapply(all_strata, function(s) {
          df_temp <- df
          df_temp$facet <- s
          df_temp$is_facet <- df_temp[[strata.var]] == s
          return(df_temp)
        }))
      } else {
        df_all <- df
        df_all$facet <- "All"
        df_all$is_facet <- TRUE
      }

      # Create base plot
      p <- ggplot2::ggplot(df_all, ggplot2::aes(PC1, PC2))

      # Add gray points and lines for all data
      p <- p +
        ggplot2::geom_point(data = subset(df_all, !is_facet),
                            aes_function,
                            size = 10, alpha = 0.5, color = "gray80") +
        ggplot2::geom_segment(
          data = subset(df_all, !is_facet),
          aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
          color = "gray80", size = 1, alpha = 0.5,
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open")
        )

      # Define aesthetic mapping
      aes_function <- if (!is.null(group.var)) {
        aes(shape = !!sym(time.var), color = !!sym(group.var))
      } else {
        aes(shape = !!sym(time.var), color = !!sym(time.var))
      }

      p <- p +
        ggplot2::geom_point(data = subset(df_all, is_facet), aes_function, size = 10, show.legend = TRUE) +
        ggplot2::geom_segment(
          data = subset(df_all, is_facet),
          aes(
            x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = if (!is.null(group.var)) !!sym(group.var) else !!sym(time.var)
          ),
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
          size = 1
        )

      p <- p +
        ggplot2::labs(
          x = ifelse(
            !is.null(pc.obj[[dist.name]]$eig),
            paste0("Axis 1 (", round(
              pc.obj[[dist.name]]$eig[1] / sum(pc.obj[[dist.name]]$eig) * 100, 2
            ), "%)"),
            "Axis 1"
          ),
          y = ifelse(
            !is.null(pc.obj[[dist.name]]$eig),
            paste0("Axis 2 (", round(
              pc.obj[[dist.name]]$eig[2] / sum(pc.obj[[dist.name]]$eig) * 100, 2
            ), "%)"),
            "Axis 2"
          )
        ) +
        scale_color_manual(values = col) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        theme_to_use +
        ggplot2::theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          axis.line.x = ggplot2::element_line(size = 1, colour = "black"),
          axis.line.y = ggplot2::element_line(size = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title = ggplot2::element_text(color = "black"),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color = "black", size = base.size),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(as.formula(paste("~ facet")))
      }


      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "beta_ordination_pair_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "dist.name_",
          dist.name
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
    names(plot_list) <- dist.name
    return(plot_list)
  }
