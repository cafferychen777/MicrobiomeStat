#' @title Generate boxplot for alpha diversity index at a single time point
#'
#' @description This function generates a boxplot of specified alpha diversity indices at a single time point dplyr::across different groupings, with optional stratification. The output can be saved as a PDF.
#' @name generate_alpha_boxplot_single
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param depth An integer specifying the sequencing depth for the "Rarefy" and "Rarefy-TSS" methods.
#' If NULL, no rarefaction is performed.
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var An optional variable in the metadata table that represents the grouping factor.
#' @param strata.var An optional variable in the metadata table that represents the stratification factor.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param base.size The base font size for the plot. Default is 16.
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
#' @param pdf A boolean indicating whether to save the output as a PDF file. Default is TRUE.
#' @param file.ann A string for annotating the output file name.
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @return A list of boxplots displaying the specified alpha diversity indices at the specified time point dplyr::across different groupings, stratified by the specified stratification variable (if provided). Each boxplot in the list corresponds to one of the alpha diversity indices specified in `alpha.name`. The boxplots will be saved as PDF files if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(ggh4x)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # First example with peerj32.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("simpson"),
#'   subject.var  = "subject",
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' alpha.obj <- mStat_calculate_alpha_diversity(peerj32.obj$feature.tab, "simpson")
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = alpha.obj,
#'   alpha.name   = c("simpson"),
#'   subject.var  = "subject",
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' # Load another dataset
#' data("subset_T2D.obj")
#'
#' # Second example with subset_T2D.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = subset_T2D.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("shannon"),
#'   subject.var  = "subject_id",
#'   time.var     = "visit_number",
#'   t.level      = "   3",
#'   group.var    = "subject_race",
#'   strata.var   = "subject_gender",
#'   adj.vars     = "sample_body_site",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 20,
#'   pdf.hei      = 8.5
#' )
#'
#' }
#' library(vegan)
#' library(ggh4x)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # First example with peerj32.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("simpson"),
#'   subject.var  = "subject",
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = "lancet",
#'   pdf          = FALSE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' @export
generate_alpha_boxplot_single <- function (data.obj,
                                           alpha.obj = NULL,
                                           alpha.name = c("shannon",
                                                          "observed_species"),
                                           depth = NULL,
                                           subject.var = NULL,
                                           time.var = NULL,
                                           t.level = NULL,
                                           group.var = NULL,
                                           strata.var = NULL,
                                           adj.vars = NULL,
                                           base.size = 16,
                                           theme.choice = "bw",
                                           custom.theme = NULL,
                                           palette = NULL,
                                           pdf = TRUE,
                                           file.ann = NULL,
                                           pdf.wid = 11,
                                           pdf.hei = 8.5,
                                           ...) {

  if (is.null(alpha.name)){
    return()
  }

  if (is.null(alpha.obj)) {
    if (!is.null(depth)) {
      message(
        "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj, depth = depth)
    }

    if (!is.null(time.var) & !is.null(t.level)){
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    otu_tab <- data.obj$feature.tab
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
  }

  meta_tab <- data.obj$meta.dat

  # Convert the alpha.obj list to a data frame
  alpha_df <-
    dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    dplyr::inner_join(meta_tab %>% rownames_to_column(var = "sample"),
                      by = c("sample"))

  if (is.null(group.var)) {
    alpha_df <- alpha_df %>% dplyr::mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  # Replace the existing theme selection code with this:
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  col <- mStat_get_palette(palette)

  # Create a plot for each alpha diversity index
  plot_list <- lapply(alpha.name, function(index) {
    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(group.var),
        y = !!sym(index),
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(group.var),
        y = !!sym(index),
        fill = !!sym(time.var)
      )
    }

    if (!is.null(adj.vars)) {
      # Factor transformation for non-numeric covariates
      data_subset <- alpha_df %>%
        dplyr::select(all_of(adj.vars)) %>%
        dplyr::mutate(dplyr::across(where(is.character) & !is.factor, factor))

      # Create a model matrix and set contrasts for non-numeric covariates.
      M <- model.matrix(
        ~ 0 + .,
        data = data_subset,
        contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
      )

      # Center the covariates (no scaling)
      M_centered <- scale(M, scale = FALSE)

      # Fit the regression model
      fit <- lm(alpha_df[[index]] ~ M_centered)

      # Calculate the adjusted alpha diversity value.
      adjusted_value <- fit$coefficients[1] + residuals(fit)

      # Update alpha diversity values in alpha_df.
      alpha_df[[index]] <- adjusted_value

      # Display message indicating that alpha diversity has been adjusted for specific covariates.
      message(
        "Alpha diversity has been adjusted for the following covariates: ",
        paste(adj.vars, collapse = ", "),
        "."
      )
    }

    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      y_label <-
        paste0(index, " index (adjusted by: ", covariates, ")")
    } else {
      y_label <- paste0(index, " index")
    }

    boxplot <- ggplot(alpha_df,
                      aes_function) +
      #geom_violin(trim = FALSE, alpha = 0.8) +
      geom_jitter(width = 0.3,
                  alpha = 0.5,
                  size = 2) +
      stat_boxplot(geom = "errorbar",
                   position = position_dodge(width = 0.2),
                   width = 0.3) +
      geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.3,
        #fill = "white"
      ) +
      scale_fill_manual(values = col) +
      {
        if (!is.null(strata.var) & !is.null(group.var)) {
          ggh4x::facet_nested(
            as.formula(paste(". ~", strata.var, "+", group.var)),
            drop = T,
            scale = "free",
            space = "free"
          )
        } else {
          if (group.var != "ALL") {
            ggh4x::facet_nested(
              as.formula(paste(". ~", group.var)),
              drop = T,
              scale = "free",
              space = "free"
            )
          }
        }
      } +
      labs(y = y_label,
           title = dplyr::if_else(
             !is.null(time.var) &
               !is.null(t.level),
             paste0(time.var, " = ", t.level),
             ""
           ))  +
      theme_to_use +
      theme(
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0, "cm"),
        strip.text = element_text(size = base.size),
        axis.text.y = element_text(color = "black", size = base.size),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = base.size),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 20)
      ) + {
        if (group.var == "ALL") {
          guides(fill = "none")
        }
      }

    return(boxplot)
  })

  # Save the plots as a PDF file
  if (pdf) {
    plot_list <- lapply(seq_along(plot_list), function(plot_index) {
      plot <- plot_list[[plot_index]]
      current_alpha_name <- alpha.name[plot_index]

      pdf_name <- paste0(
        "alpha_boxplot_single_",
        current_alpha_name,
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

      return(plot)
    })

  }

  names(plot_list) <- alpha.name

  return(plot_list)
}
