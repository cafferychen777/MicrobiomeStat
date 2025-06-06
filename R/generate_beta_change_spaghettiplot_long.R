#' @title Generate beta diversity change spaghetti plot for longitudinal data
#'
#' @description This function generates a spaghetti plot visualizing the change in beta diversity over time for different groups and strata in longitudinal data. Optionally, a user can specify whether to save the plot as a PDF file. Change is defined as distance between samples from two time points.
#'
#' @note This function requires subjects to have a baseline time point (specified by `t0.level`) to calculate the change in beta diversity. Subjects without a baseline time point will be excluded from the visualization.
#'
#' @name generate_beta_change_spaghettiplot_long
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var (Optional) The variable in the metadata table that represents the grouping factor.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param strata.var (Optional) The variable in the metadata table that represents the stratification factor.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
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
#' @details This function calculates beta diversity distances between all
#'   pairs of samples using the specified distance metric (dist_name).
#'   The change in beta diversity for each subject is then defined as
#'   the distance between the baseline sample (t0_level) and the sample
#'   at each follow-up timepoint (ts_levels). These changes are visualized
#'   using spaghetti lines connecting the distances over time.
#'
#' @return A spaghetti plot displaying the beta diversity change over time, stratified by the specified grouping and/or strata variables (if provided). The plot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_beta_change_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' generate_beta_change_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = NULL,
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_beta_change_spaghettiplot_long <-
  function(data.obj,
           dist.obj = NULL,
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

    # Exit the function if no distance metric is specified
    if (is.null(dist.name)){
      return()
    }

    # Input validation to ensure correct variable types
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

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Set the theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Calculate new sizes based on base.size for consistent plot aesthetics
    title.size = base.size * 1.25
    axis.title.size = base.size * 0.75
    axis.text.size = base.size * 0.5
    legend.title.size = base.size * 1
    legend.text.size = base.size * 0.75

    # Inform users about the requirement for baseline time points
    message("Note: This function requires subjects to have a baseline time point (specified by `t0.level`) to calculate the change in beta diversity. Subjects without a baseline time point will be excluded from the visualization.")

    # Generate plots for each specified distance metric
    plot_list <- lapply(dist.name,function(dist.name){
      # Calculate beta diversity if not provided
      if (is.null(dist.obj)) {
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
        dist.obj <-
          mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        # Adjust distances for covariates if specified
        if (!is.null(adj.vars)){
          dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
        }
      } else {
        # Process time variable if distance object is provided
        if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
          data.obj <-
            mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
        } else {
          meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
          data.obj <- list(meta.dat = meta_tab)
          data.obj <- mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
          meta_tab <- data.obj$meta.dat
        }
      }

      # Check if the specified distance metric exists in the distance object
      if (is.null(dist.obj[[dist.name]])) {
        message(paste("dist.obj does not contain", dist.name))
      } else {
        # Convert distance object to matrix format
        tryCatch({
          dist.df <- as.matrix(dist.obj[[dist.name]])
        }, error = function(e) {
          message(paste("Failed to convert dist.obj[[", dist.name, "]] to a matrix: ", e$message))
        })

        # Check for missing samples in the distance matrix
        missing_rows <- setdiff(rownames(meta_tab), rownames(dist.df))
        missing_cols <- setdiff(rownames(meta_tab), colnames(dist.df))

        if (length(missing_rows) > 0) {
          message(paste("The following rows in meta_tab are not in dist.df:", paste(missing_rows, collapse = ", ")))
        }

        if (length(missing_cols) > 0) {
          message(paste("The following columns in meta_tab are not in dist.df:", paste(missing_cols, collapse = ", ")))
        }

        # Subset distance matrix to match metadata
        if (length(missing_rows) == 0 & length(missing_cols) == 0) {
          dist.df <- dist.df[rownames(meta_tab), rownames(meta_tab)]
        }
      }

      # Convert distance matrix to long format for plotting
      dist.df <- dist.df %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      meta_tab <- meta_tab %>% rownames_to_column("sample")

      # Determine the baseline time point for change calculation
      if (is.factor(meta_tab[, time.var])) {
        change.base <- levels(meta_tab[, time.var])[1]
      } else if (is.numeric(meta_tab[, time.var])) {
        change.base <- min(meta_tab[, time.var], na.rm = TRUE)
      } else {
        stop("The variable is neither factor nor numeric.")
      }

      # Calculate pairwise distances between baseline and follow-up time points for each subject
      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      # Add group and strata information to the long-format data
      if (!is.null(strata.var)&!is.null(group.var)){
        long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(-all_of("sample")) %>% dplyr::distinct(),by = c(subject.var,time.var))
      } else if (is.null(strata.var)&!is.null(group.var)){
        long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(-all_of("sample")) %>% dplyr::distinct(),by = c(subject.var,time.var))
      } else {
        long.df <- long.df
      }

      # Create a dummy group variable if not provided
      if (is.null(group.var)){
        long.df <- long.df %>% dplyr::mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      # Calculate mean distances for each time point and group (and strata if applicable)
      if (is.null(strata.var)) {
        long.df.mean <- long.df %>%
          dplyr::group_by(!!sym(time.var),!!sym(group.var)) %>%
          dplyr::summarize(mean_distance = mean(distance, na.rm = TRUE))
        long.df <-
          dplyr::left_join(long.df, long.df.mean, by = c(time.var, group.var))
      } else {
        long.df.mean <- long.df %>%
          dplyr::group_by(!!sym(time.var),!!sym(group.var),!!sym(strata.var)) %>%
          dplyr::summarize(mean_distance = mean(distance, na.rm = TRUE))
        long.df <-
          dplyr::left_join(long.df, long.df.mean, by = c(time.var, group.var, strata.var))
      }

      long.df <- long.df %>% dplyr::arrange(subject.var,time.var)

      # Set y-axis label based on whether distances are adjusted for covariates
      if (is.null(adj.vars)) {
        y_label <- paste(dist.name, "Distance from Baseline")
      } else {
        y_label <- paste(dist.name, "Distance from Baseline (adjusted by:", paste(adj.vars, collapse = ", "), ")")
      }

      # Create the spaghetti plot
      p <- ggplot() +
        geom_point(
          data = long.df,
          aes_string(
            x = time.var,
            y = "distance",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.3,
          size = 3
        ) +
        geom_line(
          data = long.df,
          aes_string(
            x = time.var,
            y = "mean_distance",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        geom_point(
          data = long.df,
          aes_string(
            x = time.var,
            y = "mean_distance",
            group = group.var,
            color = group.var
          ),
          size = 5
        ) +
        labs(x = time.var, y = y_label, color = group.var) +
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

      # Add faceting by strata if a strata variable is specified
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(
          cols = vars(!!sym(strata.var)),
          scale = "free",
          space = "free"
        )
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0("beta_change_spaghettiplot_",
                           dist.name,
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var)

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