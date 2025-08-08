#' Generate boxplot for specified alpha diversity index
#'
#' This function generates a boxplot of a specified alpha diversity index dplyr::across different groupings and time points, with optional stratification. The output can be saved as a PDF.
#' @name generate_alpha_boxplot_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", "pielou", and "faith_pd".
#' @param depth An integer specifying the sequencing depth for the "Rarefy" and "Rarefy-TSS" methods.
#' If NULL, no rarefaction is performed.
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param group.var An optional variable in the metadata table that represents the grouping factor.
#' @param strata.var An optional variable in the metadata table that represents the stratification factor.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param base.size The base font size for the plot.
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
#' @param pdf A boolean indicating whether to save the output as a PDF file.
#' @param file.ann A string for annotating the output file name.
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @return A boxplot displaying the specified alpha diversity index dplyr::across different groupings and time points, stratified by the specified stratification variable (if provided). The boxplot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' generate_alpha_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = c("sample_body_site"),
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#'
#' data("ecam.obj")
#' generate_alpha_boxplot_long(
#'   data.obj = ecam.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   adj.vars = NULL,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#'
#' data(peerj32.obj)
#' generate_alpha_boxplot_long(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5)
#'
#' data(subset_pairs.obj)
#' alpha.obj <- mStat_calculate_alpha_diversity(subset_pairs.obj$feature.tab,
#' c("shannon", "observed_species"))
#' generate_alpha_boxplot_long(
#'   data.obj = subset_pairs.obj,
#'   alpha.obj = alpha.obj,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   t0.level = "Baseline",
#'   ts.levels = NULL,
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#' }
#' @export
generate_alpha_boxplot_long <- function (data.obj,
                                         alpha.obj = NULL,
                                         alpha.name = c("shannon", "observed_species"),
                                         depth = NULL,
                                         subject.var,
                                         time.var,
                                         t0.level = NULL,
                                         ts.levels = NULL,
                                         group.var = NULL,
                                         strata.var = NULL,
                                         adj.vars = NULL,
                                         base.size = 12,
                                         theme.choice = "bw",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         ...) {
  # Check if alpha.name is provided
  if (is.null(alpha.name)){
    return()
  }

  # Check if alpha.obj is a list
  if (!is.null(alpha.obj) &&
      !is(alpha.obj, "list"))
    stop("`alpha.obj` should be a list or NULL.")
  if (!is.character(subject.var))
    stop("`subject.var` should be a character string.")
  if (!is.character(time.var))
    stop("`time.var` should be a character string.")
  if (!is.null(group.var) &&
      !is.character(group.var))
    stop("`group.var` should be a character string or NULL.")
  if (!is.null(strata.var) &&
      !is.character(strata.var))
    stop("`strata.var` should be a character string or NULL.")

  # Calculate alpha diversity if not provided
  if (is.null(alpha.obj)) {
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    if (!is.null(depth)) {
      message(
        "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <-
        mStat_rarefy_data(data.obj = data.obj, depth = depth)
    }
    otu_tab <- data.obj$feature.tab
    
    # Extract tree if faith_pd is requested
    tree <- NULL
    if ("faith_pd" %in% alpha.name) {
      tree <- data.obj$tree
    }
    
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name, tree = tree)
  } else {
    # Verify that all alpha.name are present in alpha.obj
    if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x)
      colnames(x))))) {
      missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
      stop(
        "The following alpha diversity indices are not available in alpha.obj: ",
        paste(missing_alphas, collapse = ", "),
        call. = FALSE
      )
    }
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
  }

  # Prepare metadata and alpha diversity data
  meta_tab <- data.obj$meta.dat %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(c(subject.var, group.var, time.var, strata.var, adj.vars)))

  time.levels <- meta_tab %>%
    dplyr::select(all_of(c(time.var))) %>%
    pull() %>%
    as.factor() %>%
    levels() %>%
    length()

  # Convert the alpha.obj list to a data frame
  alpha_df <-
    dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    dplyr::inner_join(meta_tab %>% rownames_to_column(var = "sample"),
                      by = c("sample"))

  # Set up theme and color palette
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Use mStat_get_palette to set the color palette
  col <- mStat_get_palette(palette)

  # Handle case when no grouping variable is provided
  if (is.null(group.var)) {
    alpha_df <- alpha_df %>% dplyr::mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  # Create a plot for each alpha diversity index
  plot_list <- lapply(alpha.name, function(index) {
    line_aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        group = !!sym(subject.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        group = !!sym(subject.var)
      )
    }

    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        fill = !!sym(time.var)
      )
    }

    # Adjust for covariates if specified
    if (!is.null(adj.vars)) {
      data_subset <- alpha_df %>%
        dplyr::select(all_of(adj.vars)) %>%
        dplyr::mutate(dplyr::across(where(is.character) &
                                      !is.factor, factor))

      M <-
        model.matrix(
          ~ 0 + .,
          data = data_subset,
          contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
        )

      # Center the covariates
      M_centered <- scale(M, scale = FALSE)

      # Fit regression model
      fit <- lm(alpha_df[[index]] ~ M_centered)

      # Compute the adjusted value
      adjusted_value <- fit$coefficients[1] + residuals(fit)

      # Update the alpha_df
      alpha_df[[index]] <- adjusted_value

      message(
        "Alpha diversity has been adjusted for the following covariates: ",
        paste(adj.vars, collapse = ", "),
        "."
      )
    }

    # Calculate average alpha diversity for large datasets
    average_alpha_df <- NULL
    if (length(unique(alpha_df[[time.var]])) > 10 ||
        length(unique(alpha_df[[subject.var]])) > 25) {
      if (!is.null(strata.var) & !is.null(group.var)) {
        average_alpha_df <- alpha_df %>%
          dplyr::group_by(!!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
          dplyr::summarise(dplyr::across(!!sym(index), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL") %>%
          dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
      } else if (!is.null(group.var)) {
        average_alpha_df <- alpha_df %>%
          dplyr::group_by(!!sym(group.var), !!sym(time.var)) %>%
          dplyr::summarise(dplyr::across(!!sym(index), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL") %>%
          dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
      } else {
        average_alpha_df <- alpha_df %>%
          dplyr::group_by(!!sym(time.var)) %>%
          dplyr::summarise(dplyr::across(!!sym(index), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL") %>%
          dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
      }
    }

    # Set up y-axis label
    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      y_label <- paste0(index, " index (adjusted by: ", covariates, ")")
    } else {
      y_label <- paste0(index, " index")
    }

    # Ensure time variable is a factor
    alpha_df <- alpha_df %>% dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))

    boxplot <- ggplot(alpha_df,
                      aes_function) +
      #geom_violin(trim = FALSE, alpha = 0.8) +
      stat_boxplot(geom = "errorbar",
                   position = position_dodge(width = 0.2),
                   width = 0.3) +
      geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.3,
        #fill = "white"
      ) +
      geom_line(
        line_aes_function,
        alpha = 0.8,
        linewidth = 0.6,
        color = "black",
        linetype = "dashed",
        data = if (!is.null(average_alpha_df))
          average_alpha_df
        else
          alpha_df
      ) +
      scale_fill_manual(values = col) +
      {
        if (time.levels > 2) {
          if (!is.null(strata.var) & !is.null(group.var)) {
            ggh4x::facet_nested(
              cols = vars(!!sym(group.var)),
              rows = vars(!!sym(strata.var)),
              scale = "free",
              space = "free"
            )
          } else {
            if (group.var != "ALL") {
              ggh4x::facet_nested(
                cols = vars(!!sym(group.var)),
                scale = "free",
                space = "free"
              )
            }
          }
        } else {
          if (!is.null(strata.var) & !is.null(group.var)) {
            ggh4x::facet_nested(
              cols = vars(!!sym(strata.var), !!sym(group.var)),
              scale = "free",
              space = "free"
            )
          } else {
            if (group.var != "ALL") {
              ggh4x::facet_nested(
                cols = vars(!!sym(group.var)),
                scale = "free",
                space = "free"
              )
            }
          }
        }
      } +
      labs(x = time.var,
           y = y_label)  +
      theme_to_use +
      theme(
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0, "cm"),
        strip.text.x = element_text(size = base.size, color = "black"),
        strip.text.y = element_text(size = base.size, color = "black"),
        axis.text.x = element_text(
          angle = 90,
          color = "black",
          vjust = 0.5,
          size = base.size * 0.75
        ),
        axis.text.y = element_text(color = "black", size = base.size),
        axis.title.x = element_text(size = base.size),
        axis.title.y = element_text(size = base.size),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        legend.text = ggplot2::element_text(size = 16 * 2),
        legend.title = ggplot2::element_text(size = 16 * 2),
        legend.key.size = unit(10, "mm"),
        legend.key.spacing = unit(2, "mm")
      ) + {
        if (group.var == "ALL") {
          guides(fill = "none")
        }
      }

    # Add geom_jitter() if the number of unique time points or subjects is greater than 10
    if (length(unique(alpha_df[[time.var]])) > 10 ||
        length(unique(alpha_df[[subject.var]])) > 10) {
      boxplot <- boxplot + geom_jitter(width = 0.1,
                                       alpha = 0.1,
                                       size = 1)
    }

    return(boxplot)
  })

  # Save plots as PDF if requested
  if (pdf) {
    for (plot_index in seq_along(plot_list)) {
      plot <- plot_list[[plot_index]]
      current_alpha_name <- alpha.name[plot_index]

      pdf_name <- paste0(
        "alpha_boxplot_long_",
        current_alpha_name,
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

      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      print(plot)
      dev.off()
    }
  }

  # Name the plots in the list
  names(plot_list) <- alpha.name

  return(plot_list)
}
