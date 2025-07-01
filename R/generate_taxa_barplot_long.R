#' Generate Taxa Barplots for Longitudinal Data
#'
#' This function generates barplots to visualize the taxonomic composition of samples for longitudinal data.
#' It also provides options for grouping and stratifying data.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var Character string specifying the column name in metadata containing
#'                    unique subject IDs. Required to connect samples from the same subject.
#' @param time.var Character string specifying the column name in metadata containing the
#'                time variable. Required to order and connect samples over time.
#' @param group.var Character string specifying the column name in metadata containing group identifiers. Default is NULL.
#' @param strata.var Character string specifying the column name in metadata containing stratum identifiers. Default is NULL.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param feature.number A numeric value indicating the number of top abundant features to retain in the plot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 20.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param base.size A numeric value indicating the base font size for the plot.
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
#' @param palette A character vector specifying the names of colors to use for the plot.
#' The vector should contain color names that will be matched to the levels of
#' the taxonomic grouping variable (feature.level). Default is NULL, which will
#' use a pre-defined set of colors from the package ggh4x.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF file. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF file. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A ggplot object of the taxa barplot.
#' @details
#' This function generates a barplot of taxa proportions for longitudinal data.
#' The barplot can be stratified by a group variable and/or other variables.
#' It also allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#' library(ggh4x)
#' library(vegan)
#' data(ecam.obj)
#' generate_taxa_barplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Family",
#'   feature.dat.type = "proportion",
#'   feature.number = 30,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_barplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Genus",
#'   feature.dat.type = "proportion",
#'   feature.number = 30,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_barplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Genus",
#'   feature.dat.type = "proportion",
#'   features.plot = unique(ecam.obj$feature.ann[,"Genus"])[1:10],
#'   feature.number = 30,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_barplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.number = 40,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   pdf.wid = 49,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_barplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "sample_body_site",
#'   strata.var = "subject_race",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   feature.number = 20,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   pdf.wid = 49,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_barplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = unique(subset_T2D.obj$feature.ann[,"Genus"][1:6]),
#'   feature.number = 40,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   pdf.wid = 49,
#'   file.ann = NULL
#' )
#' }
#' @export
generate_taxa_barplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           feature.level = "original",
           feature.dat.type = c("count", "proportion", "other"),
           feature.number = 20,
           features.plot = NULL,
           t0.level = NULL,
           ts.levels = NULL,
           base.size = 10,
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

    # Validate input data and variables
    mStat_validate_data(data.obj)

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

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Define color palette for plotting
    if (is.null(palette)){
      # If no palette is provided, use a default color scheme
      pal <- rep(c("#E0E0E0","#E41A1C","#1E90FF","#FF8C00","#4DAF4A","#984EA3","#40E0D0","#FFC0CB",
                   "#00BFFF","#FFDEAD","#90EE90","#EE82EE","#00FFFF","#F0A3FF", "#0075DC",
                   "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00",
                   "#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010",
                   "#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",
                   "#FF6A6A", "#FF8247", "#FFE7BA", "#87CEFA", "#B0E0E6", "#48D1CC", "#5F9EA0", "#66CDAA", "#458B00",
                   "#BCEE68", "#FFF68F", "#EEEE00", "#FFFFE0", "#8B8682", "#FFB6C1", "#9370DB", "#FFDAB9", "#FA8072",
                   "#90EE90", "#00FFFF", "#00BFFF", "#FFA07A", "#F08080", "#FFD700", "#ADFF2F", "#9ACD32", "#9400D3",
                   "#7B68EE", "#BA55D3", "#FFC0CB", "#FAEBD7", "#F0E68C", "#FFFACD", "#D2B48C", "#C0C0C0", "#696969",
                   "#CD5C5C", "#F08080", "#FA8072", "#E9967A", "#FFA07A", "#FF7F50", "#FF6347", "#FF4500", "#FF8C00",
                   "#FFA500", "#FFDAB9", "#FFE4B5", "#FFE4C4", "#FFDEAD", "#EEE8AA", "#F0E68C", "#BDB76B", "#FFD700",
                   "#DAA520", "#808000", "#7CFC00", "#00FF00", "#32CD32", "#00FA9A", "#90EE90", "#98FB98", "#8FBC8F",
                   "#3CB371", "#2E8B57", "#228B22", "#008000", "#006400", "#66CDAA", "#00FFFF", "#40E0D0", "#48D1CC",
                   "#AFEEEE", "#7FFFD4", "#B0C4DE", "#ADD8E6", "#87CEEB", "#87CEFA", "#6495ED", "#00BFFF", "#1E90FF",
                   "#4169E1", "#0000FF", "#00008B", "#000080", "#191970", "#6A5ACD", "#483D8B", "#9400D3", "#8A2BE2",
                   "#4B0082", "#FF00FF", "#FF69B4", "#FF1493", "#C71585", "#DB7093", "#FFC0CB", "#FFB6C1", "#FF69B4",
                   "#FF5F5F", "#DC143C", "#C0C0C0", "#A9A9A9", "#808080", "#696969", "#000000", "#FF1493", "#FF69B4",
                   "#FFB6C1", "#FFC0CB", "#FFDAB9", "#F4A460", "#FFA07A", "#FF7F50", "#FF6347", "#FF4500", "#FF8C00",
                   "#FFA500", "#FFFF00", "#9ACD32", "#32CD32", "#00FF00", "#7FFF00", "#7CFC00", "#00FA9A", "#90EE90",
                   "#98FB98", "#8FBC8F", "#3CB371", "#2E8B57", "#228B22", "#008000", "#006400"
      ),5)
    } else{
      # Check if the provided palette is sufficient for the number of features
      if (feature.number > length(palette)) {
        stop("The number of unique features exceeds the length of the provided palette. Please provide a larger palette.")
      }

      pal = palette
    }

    # Note: For barplot with position="fill", normalization is not required
    # as ggplot2 automatically converts to proportions. We skip normalization
    # for all data types to preserve pre-computed feature.agg.list and improve
    # performance.
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw count format. For barplot visualization, ",
        "normalization is not required as ggplot2's position='fill' will ",
        "automatically convert to proportions."
      )
      # Skip normalization to preserve feature.agg.list
    } else if (feature.dat.type == "proportion") {
      message(
        "Your data is in proportion format. For barplot visualization, ",
        "additional normalization is not required as ggplot2's position='fill' ",
        "will automatically handle the proportions."
      )
      # Skip normalization to preserve feature.agg.list
    } else if (feature.dat.type == "other") {
      message("The 'other' type is suitable for situations where the user has ",
              "analyzed the data using a method not provided in ",
              "'mStat_normalize_data' method, and the 'barplot' is only ",
              "applicable to raw data that has not undergone any processing ",
              "or proportion data that adds up to 1. If you believe your data ",
              "falls into these two categories, please modify 'feature.dat.type'.")
    }

    # Generate plots for each feature level
    plot_list_all <- lapply(feature.level,function(feature.level){

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Select the appropriate feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Subset and normalize features if specified
      if (!is.null(features.plot)){
        otu_tax_agg <- otu_tax_agg[na.omit(features.plot),]
        otu_tax_agg <- apply(otu_tax_agg, 2, function(x) x / sum(x))
      }

      # Prepare the feature table for analysis
      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        rownames_to_column(feature.level)

      # Normalize the feature table
      otu_tab_norm <- apply(t(otu_tax_agg %>% select(-all_of(feature.level))), 1, function(x) x)

      rownames(otu_tab_norm) <- as.matrix(otu_tax_agg[, feature.level])

      meta_tab_sorted <- meta_tab[colnames(otu_tab_norm), ]

      # Calculate the average relative abundance of each taxon
      avg_abund <- rowMeans(otu_tab_norm)

      # Replace taxa with relative abundance lower than the threshold with "Other"
      otu_tab_other <- otu_tab_norm %>%
        as.data.frame() %>%
        rownames_to_column(feature.level)

      # Set the relative abundance threshold after feature.number
      other.abund.cutoff <- sort(avg_abund, decreasing=TRUE)[feature.number]

      if (!is.na(other.abund.cutoff)){
        otu_tab_other[, feature.level][avg_abund < other.abund.cutoff] <- "Other"
      }

      # Convert data frame to long format for plotting
      otu_tab_long <- otu_tab_other %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarize_all(sum) %>%
        tidyr::gather(key = "sample", value = "value", -feature.level)

      # Merge feature data with metadata
      merged_long_df <- otu_tab_long %>%
        dplyr::inner_join(meta_tab_sorted  %>% rownames_to_column("sample"), by = "sample")

      # Sort the merged data
      sorted_merged_long_df <- merged_long_df %>%
        dplyr::arrange(!!sym(subject.var), !!sym(time.var))

      # Identify the last sample for each subject
      last_sample_ids <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::summarize(last_sample_id = dplyr::last(sample))

      # Convert feature level to factor
      sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(feature.level) := as.factor(!!sym(feature.level)))

      # Calculate the average value of each feature and sort them
      df_sorted <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarise(overall_mean = mean(value, na.rm = TRUE)) %>%
        dplyr::mutate(is_other = ifelse(!!sym(feature.level) == "Other", TRUE, FALSE)) %>%
        dplyr::arrange(is_other, overall_mean) %>%
        dplyr::mutate(!!feature.level := factor(!!sym(feature.level), levels = !!sym(feature.level)))

      # Update feature levels
      if (!is.na(other.abund.cutoff)) {
        new_levels <- c("Other", setdiff(levels(df_sorted[[feature.level]]), "Other"))
      } else {
        new_levels <- levels(df_sorted[[feature.level]])
      }

      # Apply new sorting to the data
      sorted_merged_long_df <- sorted_merged_long_df %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels))

      # Create the main dataframe for plotting
      df <- sorted_merged_long_df %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels)) %>%
        dplyr::arrange(match(!!sym(feature.level), new_levels)) %>%
        dplyr::mutate(cumulative_value = (1-cumsum(value))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::mutate(next_cumulative_value = dplyr::if_else(sample %in% last_sample_ids$last_sample_id, NA_real_, dplyr::lead(cumulative_value))) %>%
        dplyr::ungroup()

      # Update color palette
      color_pal <- setNames(pal[1:length(new_levels)], new_levels)

      # Set bar width and spacing for the plot
      bar_width <- 0.6
      bar_spacing <- bar_width / 2

      # Identify the last time point for plotting
      if (is.factor(meta_tab[, time.var])) {
        last_time_ids <- dplyr::last(levels(meta_tab[, time.var]))
      } else if (is.numeric(meta_tab[, time.var])) {
        last_time_ids <- max(meta_tab[, time.var], na.rm = TRUE)
      } else {
        stop("The variable is neither factor nor numeric.")
      }

      # Handle grouping and stratification
      if (!is.null(strata.var)){
        if (!is.null(group.var)){
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var),!!sym(strata.var)))
        } else {
          group.var = ""
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := "")
        }
      } else {
        if (!is.null(group.var)){
        } else {
          group.var = ""
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := "")
        }
      }

      # Calculate average values for plotting
      df_average <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(feature.level),!!sym(group.var),!!sym(time.var)) %>%
        dplyr::summarise(mean_value  = mean(value)) %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels)) %>%
        dplyr::arrange(match(!!sym(feature.level), new_levels),!!sym(group.var),!!sym(time.var)) %>%
        dplyr::group_by(!!sym(group.var),!!sym(time.var)) %>%
        dplyr::mutate(cumulative_mean_value = (1-cumsum(mean_value))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::mutate(next_cumulative_mean_value = dplyr::if_else(!!sym(time.var) %in% last_time_ids, NA_real_, dplyr::lead(cumulative_mean_value))) %>%
        dplyr::ungroup()

      # Handle grouping for plotting
      if (group.var == ""){
        df_average <- df_average %>% dplyr::mutate(!!sym(group.var) := "")
        df_average <- df_average %>%
          dplyr::mutate(joint_factor = interaction(!!sym(time.var), ""))
      } else {
        df_average <- df_average %>%
          dplyr::mutate(joint_factor = interaction(!!sym(time.var), !!sym(group.var)))
      }

      # Calculate x-axis offsets for labels
      df_average <- df_average %>%
        dplyr::mutate(x_offset = ifelse(cumulative_mean_value == 0, (bar_width + bar_spacing) / 2, -(bar_width + bar_spacing) / 2))

      # Drop unused levels from joint_factor
      df_average$joint_factor <- droplevels(df_average$joint_factor)

      # Convert joint_factor to numeric for plotting
      df_average$joint_factor_numeric <- match(df_average$joint_factor, levels(df_average$joint_factor))

      # Extract labels for x-axis
      labels <- sub("\\..*", "", levels(df_average$joint_factor))

      # Separate group and strata variables if necessary
      if(!is.null(strata.var)){
        df_average <- df_average %>%
          tidyr::separate(!!sym(group.var), into = c(group.var, strata.var), sep = "\\.")
      }

      # Create the main stacked barplot
      stack_barplot_average  <- 
        df_average %>%
        ggplot(aes(x = joint_factor_numeric, y = mean_value, fill = !!sym(feature.level))) +
        geom_bar(stat = "identity", position = "fill", width = bar_width) +
        geom_segment(aes(x = joint_factor_numeric + bar_spacing, xend = joint_factor_numeric + 1 - bar_spacing, y = cumulative_mean_value, yend = next_cumulative_mean_value, group = !!sym(feature.level), color = !!sym(feature.level)),linewidth = 0.6) +
        scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
        scale_x_continuous(expand = c(0.01, 0.01), breaks = unique(df_average$joint_factor_numeric), labels = labels) +
        {
          # Add faceting if group variable is present
          if (!is.null(group.var)){
            if (group.var == ""){
            } else {
              if (!is.null(strata.var)){
                ggh4x::facet_nested(as.formula(paste(". ~", group.var, "+", strata.var)), drop = T, scale = "free", space = "free")
              } else {
                ggh4x::facet_nested(as.formula(paste(". ~", group.var)), drop = T, scale = "free", space = "free")
              }
            }
          }
        } +
        labs(fill = feature.level, y = "", x = "") +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        theme_to_use +
        theme(strip.background = element_rect(fill="white",color="black"),
              panel.spacing = unit(0,"lines"),
              strip.text.x = element_text(size= base.size,color="black"),
              axis.text.y=element_text(size= base.size,color="black"),
              axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size),
              axis.title.y = element_text(size= base.size,color="black"),
              legend.key=element_blank(),
              legend.title = element_text(size= base.size), 
              legend.text = element_text(color="black",size= base.size),
              legend.spacing.x=unit(0.1,'cm'),
              legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.4,'cm'),
              legend.key.height=unit(0.4,'cm'),
              legend.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())

      # Save the stacked barplot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0("taxa_barplot_long",
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "feature_level_", feature.level,
                           "_",
                           "feature_number_", feature.number)
        if (!is.null(group.var)) {
          pdf_name <- paste0(pdf_name, "_", "group_", group.var)
        }
        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
        }
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }
        pdf_name <- paste0(pdf_name,"_avergae", ".pdf")
        ggsave(filename = pdf_name, plot = stack_barplot_average, width = pdf.wid, height = pdf.hei)
      }

      # Return the stacked barplot
      return(stack_barplot_average)
    })

    # Name the plots in the list according to feature levels
    names(plot_list_all) <- feature.level

    # Return the list of plots
    return(plot_list_all)
  }