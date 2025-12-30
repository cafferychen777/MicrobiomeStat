#' Generate Beta Ordination Plots
#'
#' This function generates beta ordination plots using the specified distance measure.
#' The plots can be stratified by a given variable and the color and shape aesthetics
#' can be mapped to groups and time variables, respectively. The function also supports
#' saving the plots as PDF files.
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
#' @param subject.var Character string specifying the column name in metadata
#'                    containing the subject IDs. This should uniquely identify
#'                    each subject in the study. Required for connecting samples
#'                    from the same subject.
#' @param time.var Character string specifying the column name in metadata containing
#'                the time variable. This should contain time points for each
#'                sample. Required for ordering and connecting samples for the
#'                same subject over time.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var Character string specifying the column name in metadata containing
#'                 the grouping variable. This will be mapped to color aesthetic in the
#'                 plot. Optional, can be NULL.
#' @param strata.var Character string specifying the column name in metadata containing
#'                  the stratification variable. This will be used for nested faceting
#'                  in the plots. Optional, can be NULL.
#' @param adj.vars Character vector specifying column names in metadata to include as
#'                covariates in multivariate adjustment of the distance matrix before
#'                ordination. Can be empty or NULL if no adjustment needed.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#' @param base.size a numeric value specifying the base size for the plot text.
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
#' @param pdf a logical value indicating whether to save the plots as PDF files. Defaults to TRUE.
#' @param file.ann a character string specifying an annotation to add to the file names of the saved plots.
#' @param pdf.wid a numeric value specifying the width of the saved PDF files.
#' @param pdf.hei a numeric value specifying the height of the saved PDF files.
#' @param ... further arguments to be passed to the underlying functions.
#'
#' @return A list of ggplot2 objects representing the beta ordination plots.
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object, \code{\link[MicrobiomeStat]{mStat_calculate_PC}} for computing the principal coordinates, and \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{geom_boxplot}} for the underlying plot functions used, and \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_DESeqDataSet_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_phyloseq_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_SummarizedExperiment_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_qiime2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_mothur_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_dada2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_biom_as_data_obj}} for data conversion.
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = NULL,
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, "BC")
#' pc.obj <- mStat_calculate_PC(dist.obj)
#' generate_beta_ordination_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:10],
#'   group.var = "diet",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_ordination_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = NULL,
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 16,
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
generate_beta_ordination_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
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

    # Check if distance metrics are provided
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    if (is.null(dist.obj)) {
      # Process time variable and extract relevant metadata
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      
      # Adjust distances if adjustment variables are provided
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # If distance object is provided, process metadata accordingly
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      } else {
        # Extract metadata from distance object if data object is not provided
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
        data.obj <- list(meta.dat = meta_tab)
        data.obj <- mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat
      }
    }

    # Calculate principal coordinates if not provided
    if (is.null(pc.obj)) {
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = 2,
          dist.name = dist.name
        )
    }

    # Get color palette
    col <- mStat_get_palette(palette)

    # Define aesthetic mapping based on presence of group variable
    aes_function <- if (!is.null(group.var)) {
      aes(color = !!sym(group.var),
          alpha = !!sym(time.var))
    } else {
      aes(color = !!sym(time.var))
    }

    # Get appropriate theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each distance metric
    plot_list <- lapply(dist.name, function(dist.name) {
      # Extract principal coordinates
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]

      # Prepare data frame for plotting
      df <- as.data.frame(pc.mat) %>%
        setNames(c("PC1", "PC2")) %>%
        rownames_to_column("sample") %>%
        dplyr::inner_join(meta_tab %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var))) %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::mutate(x_start = PC1,
               y_start = PC2,
               x_end = NA,
               y_end = NA)

      # Calculate end points for arrows
      df <- df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::arrange(!!sym(time.var)) %>%
        dplyr::mutate(
          end_condition = if (is.factor(!!sym(time.var))) {
            levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
          } else {
            max(!!sym(time.var))
          },
          x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(PC1)),
          y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(PC2))
        ) %>%
        dplyr::ungroup()

      # Calculate mean positions for different grouping scenarios
      if (!is.null(strata.var) & !is.null(group.var)) {
        # Case: Both strata and group variables are present
        df_mean <- df %>%
          dplyr::group_by(!!sym(time.var), !!sym(group.var), !!sym(strata.var)) %>%
          dplyr::summarise(mean_PC1 = mean(PC1, na.rm = TRUE),
                           mean_PC2 = mean(PC2, na.rm = TRUE)) %>%
          dplyr::ungroup()

        df_mean <- df_mean %>%
          dplyr::group_by(!!sym(group.var), !!sym(strata.var)) %>%
          dplyr::arrange(!!sym(time.var)) %>%
          dplyr::mutate(
            end_condition = if (is.factor(!!sym(time.var))) {
              levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
            } else {
              max(!!sym(time.var))
            },
            x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC1)),
            y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC2))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(x_start = mean_PC1,
                        y_start = mean_PC2)

      } else if (!is.null(group.var)){
        # Case: Only group variable is present
        df_mean <- df %>%
          dplyr::group_by(!!sym(time.var), !!sym(group.var)) %>%
          dplyr::summarise(mean_PC1 = mean(PC1, na.rm = TRUE),
                           mean_PC2 = mean(PC2, na.rm = TRUE)) %>%
          dplyr::ungroup()

        df_mean <- df_mean %>%
          dplyr::group_by(!!sym(group.var)) %>%
          dplyr::arrange(!!sym(time.var)) %>%
          dplyr::mutate(
            end_condition = if (is.factor(!!sym(time.var))) {
              levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
            } else {
              max(!!sym(time.var))
            },
            x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC1)),
            y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC2))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(x_start = mean_PC1,
                        y_start = mean_PC2)

      } else {
        # Case: No grouping variables
        df_mean <- df %>%
          dplyr::group_by(!!sym(time.var)) %>%
          dplyr::summarise(mean_PC1 = mean(PC1, na.rm = TRUE),
                           mean_PC2 = mean(PC2, na.rm = TRUE)) %>%
          dplyr::ungroup()

        df_mean <- df_mean %>%
          dplyr::arrange(!!sym(time.var)) %>%
          dplyr::mutate(
            end_condition = if (is.factor(!!sym(time.var))) {
              levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
            } else {
              max(!!sym(time.var))
            },
            x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC1)),
            y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC2))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(x_start = mean_PC1,
                        y_start = mean_PC2)

      }

      # Create the plot
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(
          size = 5,
          aes_function,
          show.legend = T
        ) +
        ggplot2::geom_segment(
          aes(
            x = x_start,
            y = y_start,
            xend = x_end,
            yend = y_end
          ),
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
          size = 0.7,
          color = "gray70",
          alpha = 0.3
        ) +
        {
          # Add group-specific or time-specific arrows
          if (!is.null(group.var)){
            ggplot2::geom_segment(data = df_mean,
                                  aes(x = x_start,
                                      y = y_start,
                                      xend = x_end,
                                      yend = y_end,
                                      color = !!sym(group.var)),
                                  arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
                                  size = 1.5)
          } else {
            ggplot2::geom_segment(data = df_mean,
                                  aes(x = x_start,
                                      y = y_start,
                                      xend = x_end,
                                      yend = y_end,
                                      color = !!sym(time.var)),
                                  arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
                                  size = 1.5)
          }
        } +
        # Add axis labels with explained variance
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
        # Set color scale based on grouping variables
        {
          if (!is.null(group.var) | !is.null(strata.var)){
            scale_color_manual(values = col)
          } else {
            ggplot2::scale_color_gradientn(colors = c("#92c5de", "#0571b0", "#f4a582", "#ca0020"))
          }
        } +
        # Add reference lines
        ggplot2::geom_vline(xintercept = 0,
                            linetype = "dashed",
                            color = "black") +
        ggplot2::geom_hline(yintercept = 0,
                            linetype = "dashed",
                            color = "black") +
        theme_to_use  +
        # Customize theme elements
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
          axis.text = ggplot2::element_text(color = "black",
                                            size = base.size),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      # Add faceting if strata variable is present
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested_wrap(as.formula(paste(".~", strata.var)), ncol = 3)
      }

      # Save the plot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "beta_ordination_long_",
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

    # Assign names to the plot list
    names(plot_list) <- dist.name

    return(plot_list)
  }