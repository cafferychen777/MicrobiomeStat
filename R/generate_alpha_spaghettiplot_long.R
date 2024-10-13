#' Generate an alpha diversity line plot for longitudinal data
#'
#' This function creates a ggplot object of alpha diversity (e.g., Shannon index) line plot for longitudinal data,
#' showing individual subject trajectories and the mean trajectory for each group.
#' @name generate_alpha_spaghettiplot_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou". Previously named as `alpha.index`.
#' @param depth An integer specifying the sequencing depth for the "Rarefy" and "Rarefy-TSS" methods.
#' If NULL, no rarefaction is performed.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param base.size The base font size for the plot.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
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
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param subject.var The name of the subject variable.
#' @param time.var The name of the time variable.
#' @param group.var The name of the group variable.
#' @param strata.var The name of the strata variable (default is NULL).
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param pdf Logical, whether to save the plot as a PDF (default is TRUE).
#' @param file.ann The annotation to be added to the PDF file name (default is NULL).
#' @param ... Additional arguments passed to ggplot().
#'
#' @return A ggplot object of the alpha diversity line plot.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' T2D.alpha.obj <- mStat_calculate_alpha_diversity(subset_T2D.obj$feature.tab,"shannon")
#' generate_alpha_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = T2D.alpha.obj,
#'   alpha.name = c("shannon"),
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = sort(unique(subset_T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = sort(unique(subset_T2D.obj$meta.dat$visit_number))[-1],
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("ecam.obj")
#' generate_alpha_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","simpson", "observed_species"),
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = NULL,
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_alpha_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","simpson", "observed_species"),
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "delivery",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_alpha_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","simpson", "observed_species"),
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   adj.vars = NULL,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_alpha_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","simpson", "observed_species"),
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   adj.vars = "antiexposedall",
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_alpha_spaghettiplot_long <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = c("shannon",
                          "simpson",
                          "observed_species",
                          "chao1",
                          "ace",
                          "pielou"),
           depth = NULL,
           subject.var,
           time.var,
           t0.level,
           ts.levels,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           base.size = 16,
           palette = NULL,
           theme.choice = "bw",
           custom.theme = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    
    # Exit the function if no alpha diversity indices are specified
    if (is.null(alpha.name)){
      return()
    }

    # Validate input data types to ensure proper function execution
    if (!is(data.obj, "list"))
      stop("`data.obj` should be a list.")
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
    # This ensures we have the necessary diversity metrics for visualization
    if (is.null(alpha.obj)) {
      # Process time variable to ensure proper ordering in longitudinal analysis
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

      # Perform rarefaction if depth is specified
      # Rarefaction standardizes sampling effort across all samples
      if (!is.null(depth)) {
        message(
          "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj, depth = depth)
      }

      otu_tab <- data.obj$feature.tab

      alpha.obj <-
        mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
    } else {
      # Verify that all requested alpha diversity indices are available
      if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x) colnames(x))))) {
        missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
        stop("The following alpha diversity indices are not available in alpha.obj: ",
             paste(missing_alphas, collapse = ", "), call. = FALSE)
      }
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    }

    # Extract relevant metadata for the analysis
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
        subject.var, group.var, time.var, strata.var, adj.vars
      )))

    # Combine alpha diversity data with metadata
    # This step creates a comprehensive dataset for visualization
    alpha.df <-
      dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
      dplyr::inner_join(
        meta_tab %>% dplyr::select(all_of(
          c(subject.var, time.var, group.var, strata.var, adj.vars)
        )) %>% rownames_to_column(var = "sample"),
        by = c("sample")
      )

    # If no group variable is specified, create a single group for all samples
    if (is.null(group.var)){
      alpha.df <- alpha.df %>% dplyr::mutate("ALL" = "ALL")
      group.var <- "ALL"
    }

    # Select the appropriate theme based on user input or default to "bw"
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Get the color palette for the plot
    col <- mStat_get_palette(palette)

    # Create a plot for each alpha diversity index
    plot_list <- lapply(alpha.name, function(index) {
      # Subset data for the current alpha diversity index
      sub_alpha.df <- alpha.df %>% dplyr::select(all_of(c(index,subject.var, time.var, group.var, strata.var, adj.vars)))

      # Adjust alpha diversity for covariates if specified
      # This step helps to control for potential confounding factors
      if (!is.null(adj.vars)) {

        # Convert non-numeric covariates to factors
        data_subset <- sub_alpha.df %>%
          dplyr::select(all_of(adj.vars)) %>%
          dplyr::mutate(dplyr::across(where(is.character) & !is.factor, factor))

        # Create a model matrix for non-numeric covariates
        M <- model.matrix(
          ~ 0 + .,
          data = data_subset,
          contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
        )

        # Center the covariates (no scaling)
        M_centered <- scale(M, scale = FALSE)

        # Fit a linear regression model to adjust for covariates
        fit <- lm(sub_alpha.df[[index]] ~ M_centered)

        # Compute the adjusted alpha diversity value
        adjusted_value <- fit$coefficients[1] + residuals(fit)

        # Update the alpha diversity values with the adjusted values
        sub_alpha.df[[index]] <- adjusted_value

        # Inform the user about the adjustment
        message("Alpha diversity has been adjusted for the following covariates: ",
                paste(adj.vars, collapse = ", "), ".")
      }

      # Calculate mean alpha diversity for each time point and group
      if (is.null(strata.var)) {
        sub_alpha.df.mean <- sub_alpha.df %>%
          dplyr::group_by(!!sym(time.var), !!sym(group.var)) %>%
          dplyr::summarize(mean_alpha = mean(!!sym(index), na.rm = TRUE))
        sub_alpha.df <-
          dplyr::left_join(sub_alpha.df, sub_alpha.df.mean, by = c(time.var, group.var))
      } else {
        sub_alpha.df.mean <- sub_alpha.df %>%
          dplyr::group_by(!!sym(time.var), !!sym(group.var), !!sym(strata.var)) %>%
          dplyr::summarize(mean_alpha = mean(!!sym(index), na.rm = TRUE))
        sub_alpha.df <-
          dplyr::left_join(sub_alpha.df, sub_alpha.df.mean, by = c(time.var, group.var, strata.var))
      }

      # Set the baseline time point if not specified
      if (is.null(t0.level)) {
        if (is.numeric(meta_tab[, time.var])) {
          t0.level <- sort(unique(meta_tab[, time.var]))[1]
        } else {
          t0.level <- levels(meta_tab[, time.var])[1]
        }
      }

      # Set the subsequent time points if not specified
      if (is.null(ts.levels)) {
        if (is.numeric(meta_tab[, time.var])) {
          ts.levels <- sort(unique(meta_tab[, time.var]))[-1]
        } else {
          ts.levels <- levels(meta_tab[, time.var])[-1]
        }
      }

      # Create a vector of all time points
      time.points <- c(t0.level, ts.levels)

      # Select a subset of time points to display if there are more than 80
      # This improves readability of the plot without affecting the underlying data
      if (length(time.points) > 80) {
        indices <- round(seq(1, length(time.points), length.out = 80))
        message("There are more than 80 time points, so we are selecting a subset of 80 to display on the x-axis. This does not affect any calculations or the resulting spaghetti plot.")
      } else {
        indices <- 1:length(time.points)
      }

      # Create breaks and labels for the x-axis
      breaks <- time.points[indices]
      labels <- time.points[indices]

      # Set the y-axis label, indicating if values are adjusted for covariates
      if (!is.null(adj.vars)) {
        covariates <- paste(adj.vars, collapse = ", ")
        y_label <- paste0(index, " index (adjusted by: ", covariates, ")")
      } else {
        y_label <- paste0(index, " index")
      }

      # Create the spaghetti plot using ggplot2
      plot <- ggplot() +
        geom_point(
          data = sub_alpha.df,
          aes_string(
            x = time.var,
            y = index,
            group = subject.var,
            color = group.var
          ),
          alpha = 0.5,
          size = 3
        ) +
        geom_line(
          data = sub_alpha.df,
          aes_string(
            x = time.var,
            y = "mean_alpha",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        geom_point(
          data = sub_alpha.df,
          aes_string(
            x = time.var,
            y = "mean_alpha",
            group = group.var,
            color = group.var
          ),
          size = 5
        ) +
        labs(x = time.var, y = y_label, color = group.var) +
        scale_color_manual(
          values = col
        ) +
        scale_x_discrete(breaks = breaks, labels = labels) +
        theme_to_use +
        theme(
          strip.text.x = element_text(size = 15, color = "black"),
          axis.title.x = element_text(size = base.size*2),
          axis.title.y = element_text(size = base.size*2),
          axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size*2),
          axis.text.y = element_text(size = base.size*2),
          legend.text = ggplot2::element_text(size = 16 * 2),
          legend.title = ggplot2::element_text(size = 16 * 2),
          legend.key.size = unit(10, "mm"),
          legend.key.spacing = unit(2, "mm"),
          panel.spacing.x = unit(1, "cm"),
          panel.spacing.y = unit(1, "cm"),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm")
        )

      # Add faceting if a strata variable is specified
      if (!is.null(strata.var)) {
        plot <- plot + ggh4x::facet_nested(
          cols = vars(!!sym(strata.var)),
          scale = "free",
          space = "free"
        )
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "alpha_spaghettiplot_long",
          "_",
          index,
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "group_",
          group.var
        )

        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
        }

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }

        pdf_name <- paste0(pdf_name, ".pdf")

        ggsave(
          filename = pdf_name,
          plot = plot,
          width = pdf.wid,
          height = pdf.hei
        )
      }
      return(plot)
    })

    # Name the plots in the list according to the alpha diversity indices
    names(plot_list) <- alpha.name

    return(plot_list)
  }