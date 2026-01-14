#' @title Generate longitudinal spaghettiplot plots of taxonomic composition
#'
#' @description This function generates a spaghettiplot plot for the taxonomic composition of microbiome data over time.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string defining the subject variable in meta_tab.
#'                   This should be a column that contains unique subject IDs to
#'                   identify each sample's subject of origin. Required for plotting
#'                   individual subject lines.
#' @param time.var A character string defining the time variable in meta_tab.
#'                This should be a column in meta_tab that represents the time
#'                points for the longitudinal samples. Required to identify
#'                the time axis for plotting.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var A character string defining the grouping variable in meta_tab.
#'                 This variable is used for coloring lines in the plot and faceting by
#'                 group. For example, variable could be treatment group. Optional,
#'                 can be set to NULL.
#' @param strata.var (Optional) A character string defining the stratification variable
#'                    in meta_tab. This variable is used for nested faceting in the plot.
#'                    For example, could be clinical strata like disease status.
#'                    Optional, can be set to NULL. When set, plot will show nested
#'                    facetting by strata within each feature.
#' @param feature.level The column name(s) in the feature annotation matrix (feature.ann) of data.obj
#'                     to use for summarization and plotting. Can be taxonomic levels like "Phylum",
#'                     "Genus", or any other annotation columns like "OTU_ID". Should be a character
#'                     vector specifying one or more column names in feature.ann. The special value
#'                     "original" can also be provided, which will use the original row names from
#'                     feature.ann. Multiple columns can be provided, and data will be plotted
#'                     separately for each column. **\bold{Cannot be NULL, as NULL value will lead to errors.}**
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw count input.
#' @param top.k.plot Integer specifying number of top k features to plot, when `features.plot` is NULL.
#' Default is NULL, in which case all features passing filters will be plotted.
#' @param top.k.func Function to use for selecting top k features, when `features.plot` is NULL.
#' Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
#' case features will be selected by abundance.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size Base font size for the generated plots.
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
#' @param pdf A logical value. If TRUE (default), saves the plot as a PDF file. If FALSE, the plot will be displayed interactively without creating a PDF.
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed.
#'
#' @return A ggplot object containing the longitudinal line plot.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data("subset_T2D.obj")
#' generate_taxa_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "sample_body_site",
#'   feature.level = c("Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "npg",
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' }
#' @export
generate_taxa_spaghettiplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter = 0.01,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Validate the input data object
    mStat_validate_data(data.obj)

    # Check if input variables are of the correct type
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

    # Process time variable in the data object
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

    # Extract relevant metadata
    meta_tab <- data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Calculate sizes for various plot elements based on base.size
    title.size = base.size * 1.25
    axis.title.size = base.size * 0.75
    axis.text.size = base.size * 0.5
    legend.title.size = base.size * 1
    legend.text.size = base.size * 0.75

    # Adjust filters if necessary based on input parameters
    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    # Normalize count data if necessary
    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    # Generate plots for each feature level
    plot_list <- lapply(feature.level, function(feature.level){

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Get the appropriate feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Apply abundance and prevalence filters
      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        rownames_to_column(feature.level)

      # Select top k features if specified
      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Convert counts to numeric type
      otu_tax_agg_numeric <-
        dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      # Reshape data from wide to long format
      df <- otu_tax_agg_numeric %>%
        tidyr::gather(key = "sample", value = "count", -all_of(feature.level)) %>%
        dplyr::left_join(meta_tab %>% rownames_to_column(var = "sample"), by = "sample")

      # Create a dummy group if group variable is not specified
      if (is.null(group.var)) {
        df <- df %>% dplyr::mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      # Calculate mean counts
      if (!is.null(strata.var)) {
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level),
                          !!sym(time.var),
                          !!sym(group.var),
                          !!sym(strata.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df,
                    mean_df,
                    by = c(feature.level, time.var, group.var, strata.var))
      } else {
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level), !!sym(time.var), !!sym(group.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df, mean_df, by = c(feature.level, time.var, group.var))
      }

      # Get unique taxa levels
      taxa.levels <-
        df %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # Select features to plot
      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 5) {
          features.plot <- taxa.levels[1:4]
        } else {
          features.plot <- taxa.levels
        }
      }

      # Subset data for selected features
      sub_df <- df %>% filter(!!sym(feature.level) %in% features.plot)

      # Get number of strata levels if strata variable is specified
      if (!is.null(strata.var)){
        strata.levels <- df %>% select(!!sym(strata.var)) %>% pull() %>% unique() %>% length()
      }

      # Create the spaghetti plot
      lineplot <- ggplot() +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "count",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.5
        ) +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_count",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        geom_point(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_count",
            group = group.var,
            color = group.var
          ),
          size = 3
        ) +
        scale_color_manual(values = col) +
        {
          if (feature.dat.type != "other"){
            labs(
              x = time.var,
              y = "Relative Abundance",
              color = group.var
            )
          } else {
            labs(
              x = time.var,
              y = "Abundance",
              color = group.var
            )
          }
        } +
        {
          if (!is.null(strata.var)) {
            ggh4x::facet_nested_wrap(as.formula(paste('~', feature.level,"+",strata.var)), scales = "free_y", ncol = strata.levels*2)  # Use facet_wrap with strata.var as the faceting variable
          } else {
            ggh4x::facet_nested_wrap(as.formula(paste('~',feature.level)), scales = "free_y", ncol = 3)
          }
        } +
        theme_to_use +
        theme(
          plot.title = element_text(
            size = title.size,
            hjust = 0.5
          ),
          panel.spacing.x = unit(0, "cm"),
          axis.title.x = element_text(size = axis.title.size),
          axis.title.y = element_text(size = axis.title.size),
          axis.text.x = element_text(size = axis.text.size),
          axis.text.y = element_text(size = axis.text.size),
          legend.title = element_text(size = legend.title.size),
          legend.text = element_text(size = legend.text.size)
        )

      # Remove legend for single group case
      if (group.var == "ALL") {
        lineplot <- lineplot + theme(legend.position = "none")
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_spaghettiplot_long",
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "group_",
          group.var,
          "_",
          "strata_",
          strata.var,
          "_",
          "feature_level_",
          feature.level,
          "_",
          "prev_filter_",
          prev.filter,
          "_",
          "abund_filter_",
          abund.filter
        )

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }

        pdf_name <- paste0(pdf_name, ".pdf")

        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Print the plot to the PDF
        print(lineplot)
        # Close the PDF device
        dev.off()
      }

      return(lineplot)
    })

    # Name the plots in the list
    names(plot_list) <- feature.level

    # Return the list of plots
    return(plot_list)
  }