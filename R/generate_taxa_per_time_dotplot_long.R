#' Longitudinal Taxa Dot Plot for Microbiome Data
#'
#' This function generates dot plots for visualizing the results of longitudinal taxa tests in microbiome data. It focuses on displaying how the abundance and significance of microbial taxa change over time across different groups.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param test.list A list of test results from `generate_taxa_per_time_test_long` function, representing the statistical analysis of taxa over time.
#' @param group.var A string specifying the group variable in meta.dat for between-group comparisons in the dot plot.
#' @param time.var A string representing the time variable in the meta.dat. This variable is used for x-axis in the dot plot.
#' @param t0.level Optional; a baseline time level for time series analysis.
#' @param ts.levels Optional; time series levels for analysis.
#' @param feature.level A string or vector of strings indicating the taxonomic level(s) for analysis (e.g., "Phylum", "Class") and for displaying in the dot plot.
#' @param feature.sig.level The significance level cutoff for highlighting taxa
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none"
#' @param filter_significant Boolean; whether to filter and display only significant features. Default is TRUE.
#' @param features.plot Optional; a vector of feature names to be plotted. If provided, overrides the significant feature filter.
#' @param base.size Numeric; base font size for text elements in the plot. Default is 16.
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
#' @param palette Optional; a vector of colors for the gradient scale. If provided, it overrides the default color scale.
#' @param pdf Boolean; whether to save the plot as a PDF file.
#' @param pdf.wid Numeric; width of the saved PDF file.
#' @param pdf.hei Numeric; height of the saved PDF file.
#' @details
#' The function uses ggplot2 to create dot plots, with options for customizing themes and aesthetics. It visualizes the significance and magnitude of changes in microbial taxa, indicated by color and size of the dots. The function allows for detailed customization of plot features, including base font size and theme choice.
#'
#' The function prepares data by merging and transforming the results from `generate_taxa_per_time_test_long` to fit the plotting structure. It handles the display of statistical significance using color coding and sizes dots based on the coefficient value from the statistical test.
#'
#' @return
#' A list of ggplot objects, each representing a dot plot for a specific taxonomic level and group. The plots can be further customized or directly rendered.
#'
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset
#' data("ecam.obj")
#' # Analyzing the impact of delivery method on microbial composition over months
#' result1 <- generate_taxa_per_time_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   adj.vars = "diet",
#'   feature.level = c("Phylum", "Class", "Family"),
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   feature.dat.type = "proportion"
#' )
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_taxa_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = result1,
#'   group.var = "delivery",
#'   time.var = "month_num",
#'   feature.level = c("Phylum", "Class", "Family"),
#'   feature.mt.method = "fdr",
#'   feature.sig.level = 0.1
#' )
#'
#' # Example 2: Analyzing the Type 2 Diabetes dataset
#' data("subset_T2D.obj")
#' # Longitudinal analysis of microbial changes in different racial groups
#' result2 <- generate_taxa_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus", "Family"),
#'   feature.dat.type = "count"
#' )
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_taxa_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result2,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   feature.level = c("Genus", "Family")
#' )
#' }
#' @export
generate_taxa_per_time_dotplot_long <- function(data.obj,
                                                test.list,
                                                group.var,
                                                time.var,
                                                t0.level = NULL,
                                                ts.levels = NULL,
                                                feature.level,
                                                feature.mt.method = "none",
                                                feature.sig.level = 0.05,
                                                features.plot = NULL,
                                                filter_significant = TRUE,
                                                base.size = 16,
                                                theme.choice = "bw",
                                                custom.theme = NULL,
                                                palette = NULL,
                                                pdf = FALSE,
                                                pdf.wid = 7,
                                                pdf.hei = 5
){

  # Process the time variable in the data object
  # This step ensures that the time variable is properly formatted and levels are set correctly
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  # FIXED: Extract time levels from test.list to ensure consistency
  # This ensures time_levels match the actual time points in test.list
  time_levels <- names(test.list)
  
  # Convert to numeric and sort to ensure proper ordering, then back to character
  if (all(grepl("^\\d+(\\.\\d+)?$", time_levels))) {
    time_levels <- as.character(sort(as.numeric(time_levels)))
  } else {
    time_levels <- sort(time_levels)
  }

  # Get the appropriate theme for plotting
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Extract group names from the test list
  group.names <- names(test.list[[1]][[1]])

  # Determine which p-value to use based on the multiple testing correction method
  p_val_var <- if (feature.mt.method == "fdr") "Adjusted.P.Value" else "P.Value"

  # Process test results for each feature level and group
  test.list <- lapply(feature.level, function(feature.level){
    sub.test.list <- lapply(group.names, function(group.names){
      # Function to merge test results across time points
      merge_time_points <- function(time_test_list, feature_level, group_name) {
        # Extract and combine data for each time point
        data_list <- lapply(time_test_list, function(time_point_data) {
          if (feature_level %in% names(time_point_data) &&
              group_name %in% names(time_point_data[[feature_level]])) {
            return(time_point_data[[feature_level]][[group_name]])
          }
          return(NULL)
        })

        # FIXED: Record valid time point names before filtering
        # This preserves the correct mapping to actual time points
        valid_time_names <- names(data_list)[!sapply(data_list, is.null)]
        data_list <- data_list[!sapply(data_list, is.null)]
        
        # FIXED: Manually add time point column instead of relying on .id
        # This ensures correct time point labels instead of sequential indices
        merged_data_list <- mapply(function(data, time_name) {
          data[[time.var]] <- time_name
          return(data)
        }, data_list, valid_time_names, SIMPLIFY = FALSE)
        
        merged_data <- dplyr::bind_rows(merged_data_list)
        return(merged_data)
      }

      # Merge data for the current feature level and group
      merged_data_genus <- merge_time_points(test.list, feature.level, group.names)
      return(merged_data_genus)
    })
    names(sub.test.list) <- group.names
    return(sub.test.list)
  })

  names(test.list) <- feature.level

  # Generate plots for each feature level
  plot.list <- lapply(feature.level, function(feature.level){

    # Generate plots for each group within the current feature level
    sub_plot.list <- lapply(group.names, function(group.names){

      # Prepare data for plotting
      data_for_plot <- test.list[[feature.level]][[group.names]]

      # Ensure time variable is a factor with correct levels
      data_for_plot[[time.var]] <- factor(data_for_plot[[time.var]], levels = time_levels)

      # Add significance labels based on p-values
      data_for_plot$Significance_Label <- ifelse(data_for_plot[[p_val_var]] < feature.sig.level, "*", "")

      # Filter features based on user specifications
      if (!is.null(features.plot)) {
        # If specific features are requested, keep only those
        data_for_plot <- data_for_plot %>%
          dplyr::filter(Variable %in% features.plot)
      } else if (filter_significant) {
        # If filtering for significant features is requested, keep only significant ones
        significant_features <- data_for_plot %>%
          dplyr::group_by(Variable) %>%
          dplyr::summarise(is_significant = any(!!sym(p_val_var) < feature.sig.level)) %>%
          dplyr::filter(is_significant) %>%
          dplyr::pull(Variable)

        data_for_plot <- data_for_plot %>%
          dplyr::filter(Variable %in% significant_features)
      }

      # Create the dot plot
      # Check if time variable contains only numeric values
      if (all(grepl("^\\d+(\\.\\d+)?$", time_levels))) {
        # Convert to numeric for proper spacing
        data_for_plot[[paste0(time.var, "_numeric")]] <- as.numeric(as.character(data_for_plot[[time.var]]))
        
        # Create plot with numeric x-axis
        dotplot <- ggplot(data_for_plot, aes(x = !!sym(paste0(time.var, "_numeric")), 
                                             y = Variable, size = Coefficient)) +
          geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
          geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
          scale_color_gradientn(colors = rev(c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020"))) +
          scale_x_continuous(breaks = as.numeric(time_levels), labels = time_levels) +
          labs(title = group.names,
               x = time.var,
               y = feature.level,
               size = "Coefficient",
               color = p_val_var) +
          scale_radius(range = c(0, 10)) +
          theme_to_use
      } else {
        # Keep current factor-based approach for non-numeric time values
        dotplot <- ggplot(data_for_plot, aes(x = !!sym(time.var), y = Variable, size = Coefficient)) +
          geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
          geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
          scale_color_gradientn(colors = rev(c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020"))) +
          labs(title = group.names,
               x = time.var,
               y = feature.level,
               size = "Coefficient",
               color = p_val_var) +
          scale_radius(range = c(0, 10)) +
          theme_to_use
      }
      
      # Apply common theme settings to the plot
      dotplot <- dotplot +
        theme(
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = base.size),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = base.size),
          axis.text.y = element_text(size = base.size),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = base.size),
          panel.spacing = unit(0, "lines"),
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.text = element_text(size = base.size),
          legend.title = element_text(size = base.size)
        )

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_filename <- paste0("dotplot_", feature.level, "_", group.names, "_", time.var)
        if (!is.null(group.var)) {
          pdf_filename <- paste0(pdf_filename, "_group_", group.var)
        }
        pdf_filename <- paste0(pdf_filename, ".pdf")
        ggsave(pdf_filename, plot = dotplot, width = pdf.wid, height = pdf.hei)
      }

      return(dotplot)
    })

    names(sub_plot.list) <- group.names

    return(sub_plot.list)
  })

  names(plot.list) <- feature.level

  return(plot.list)

}