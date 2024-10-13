#' Compare within-group and between-group beta diversity across time points
#'
#' This function generates boxplots to visualize changes in within-group and between-group
#' beta diversity across multiple time points, for different groups.
#'
#' For each time point, it calculates the within-group and between-group beta diversity (distance)
#' between all pairs of samples at that time point.
#'
#' Boxplots are generated to compare the distribution of within-group and between-group
#' distances at each time point. The boxes are dodge-positioned to allow for easy comparison.
#'
#' A line plot connecting the medians of boxes is overlaid to show the trajectory over time.
#'
#' Adjustment for covariates is supported by using adjusted distances.
#' @name generate_beta_change_boxplot_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param subject.var A string specifying the name of the subject variable
#' @param time.var A string specifying the name of the time variable
#' @param group.var A string specifying the name of the group variable or NULL (default)
#' @param strata.var A string specifying the name of the strata variable or NULL (default)
#' @param adj.vars A string specifying the name of the adjustment variable or NULL (default)
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
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
#' @param pdf A logical value indicating whether to save the plot as a PDF file. Default is TRUE
#' @param file.ann A character string specifying a custom annotation for the PDF file name or NULL (default)
#' @param pdf.wid (Optional) The width of the PDF file if `pdf` is set to `TRUE` (default is 11).
#' @param pdf.hei (Optional) The height of the PDF file if `pdf` is set to `TRUE` (default is 8.5).
#' @param ... Additional parameters passed on to ggsave()
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object.
#'
#' @return A named list of ggplot objects visualizing within-group and between-group beta diversity across time points.
#' The list contains one plot for each distance metric specified in \code{dist.name}.
#' Each plot shows boxplots of the within-group and between-group distances at each time point,
#' with boxes dodge-positioned by the \code{group_var} variable if provided.
#' A line plot connecting the medians of boxes is overlaid to show the trajectory over time.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(vegan)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC")
#' generate_beta_change_boxplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   t0.level = "1",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#' generate_beta_change_boxplot_long(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   t0.level = "Baseline",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_T2D.obj")
#' generate_beta_change_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_gender",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_beta_change_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("ecam.obj")
#' generate_beta_change_boxplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "diet",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_change_boxplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "diet",
#'   strata.var = "delivery",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_beta_change_boxplot_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
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

    # Exit the function if no distance metric is specified
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    # This step ensures we have the necessary distance matrices for the analysis
    if (is.null(dist.obj)) {
      # Process time variable to ensure proper ordering in longitudinal analysis
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Adjust distances for covariates if specified
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # Process time variable if distance object is provided
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
        data.obj <- list(meta.dat = meta_tab)
        meta_tab <- data.obj$meta.dat
      }
    }

    # Add sample names to metadata
    meta_tab <- meta_tab %>% rownames_to_column("sample")

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Set the theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each specified distance metric
    plot_list <- lapply(dist.name, function(dist.name) {

      # Convert distance object to a tibble for easier manipulation
      dist_tibble <- as_tibble(as.matrix(dist.obj[[dist.name]]), rownames = "Sample1")

      # Join metadata with distance information
      if (!is.null(strata.var)){
        dist_meta <- dist_tibble %>%
          tidyr::pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var), Strata = all_of(strata.var)), by = c("Sample1" = "sample")) %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var), Strata = all_of(strata.var)), by = c("Sample2" = "sample"))
      } else {
        dist_meta <- dist_tibble %>%
          tidyr::pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var)), by = c("Sample1" = "sample")) %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var)), by = c("Sample2" = "sample"))
      }

      # Handle different plotting scenarios based on the number of groups
      if (length(unique(dist_meta$Group.x)) <= 2) {
        # For two or fewer groups, calculate within- and between-group distances
        if (!is.null(strata.var)){
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = dplyr::if_else(Group.x == Group.y, "Within", "Between"),
              Group = dplyr::if_else(Type == "Within", Group.x, paste(Group.x, Group.y, sep = "-")),
              Time = Time.x
            ) %>%
            select(Group, Distance, Type, Time, Strata = Strata.x)
        } else {
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = dplyr::if_else(Group.x == Group.y, "Within", "Between"),
              Group = dplyr::if_else(Type == "Within", Group.x, paste(Group.x, Group.y, sep = "-")),
              Time = Time.x
            ) %>%
            select(Group, Distance, Type, Time)
        }

        # Create boxplot for within- and between-group distances
        p <- ggplot(within_between_dist, aes(x = Time, y = Distance, fill = Type)) +
          geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.8),
            width = 0.2
          ) +
          stat_summary(
            fun = mean,
            geom = "point",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 2
          ) +
          stat_summary(
            fun = mean,
            geom = "line",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 1
          ) +
          scale_fill_manual(values = col) +
          scale_color_manual(values = col) +
          labs(x = time.var, y = "Distance") +
          theme_to_use +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                plot.title = element_text(hjust = 0.5))

      } else {
        # For more than two groups, calculate pairwise distances between all groups
        if (!is.null(strata.var)){
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = paste(Group.x, Group.y, sep = "-"),
              Time = Time.x
            ) %>%
            select(Type, Distance, Time, Strata = Strata.x) %>%
            dplyr::mutate(Type = apply(select(., Type), 1, function(x) paste(sort(unlist(strsplit(x, "-"))), collapse = "-"))) %>%
            distinct(Type, Time, Distance, Strata, .keep_all = TRUE)
        } else {
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = paste(Group.x, Group.y, sep = "-"),
              Time = Time.x
            ) %>%
            select(Type, Distance, Time) %>%
            dplyr::mutate(Type = apply(select(., Type), 1, function(x) paste(sort(unlist(strsplit(x, "-"))), collapse = "-"))) %>%
            distinct(Type, Time, Distance, .keep_all = TRUE)
        }

        # Create boxplot for pairwise distances between all groups
        p <- ggplot(within_between_dist, aes(x = Time, y = Distance, fill = Type)) +
          geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.8),
            width = 0.2
          ) +
          stat_summary(
            fun = mean,
            geom = "point",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 3
          ) +
          stat_summary(
            fun = mean,
            geom = "line",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 1
          ) +
          scale_fill_manual(values = col) +
          scale_color_manual(values = col) +
          labs(x = time.var, y = "Distance") +
          theme_to_use +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                plot.title = element_text(hjust = 0.5))
      }

      # Add faceting by strata if a strata variable is specified
      if (!is.null(strata.var)) {
        p <- p + facet_wrap2("Strata", scales = "free", nrow = length(unique(meta_tab[[strata.var]])))
      } else {
        p <- p
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0("beta_change_boxplot_",
                           dist.name,
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "t0_level_", t0.level)

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

    # Name the plots in the list according to the distance metrics
    names(plot_list) <- dist.name
    return(plot_list)
  }