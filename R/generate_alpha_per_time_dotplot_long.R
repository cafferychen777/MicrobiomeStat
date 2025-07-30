#' Longitudinal Alpha Diversity Dot Plot for Microbiome Data
#'
#' This function generates dot plots for visualizing the results of longitudinal alpha diversity tests in microbiome data. It focuses on displaying how the alpha diversity indices change over time across different groups.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param test.list A list of test results from `generate_alpha_per_time_test_long` or `generate_alpha_change_test_long` function, representing the statistical analysis of alpha diversity over time.
#' @param group.var A string specifying the group variable in meta.dat for between-group comparisons in the dot plot.
#' @param time.var A string representing the time variable in the meta.dat. This variable is used for the x-axis in the dot plot.
#' @param t0.level The baseline time level for time series analysis.
#' @param ts.levels The time series levels for analysis.
#' @param base.size Optional; numeric, base font size for text elements in the plot. Default is 16.
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
#' @param pdf Optional; boolean, whether to save the plot as a PDF file.
#' @param pdf.wid Optional; numeric, width of the saved PDF file.
#' @param pdf.hei Optional; numeric, height of the saved PDF file.
#' @details
#' The function uses ggplot2 to create dot plots, with options for customizing themes and aesthetics. It visualizes the significance and magnitude of changes in alpha diversity indices, indicated by color and size of the dots. The function allows for detailed customization of plot features, including base font size and theme choice.
#'
#' The function prepares data by merging and transforming the results from `generate_alpha_per_time_test_long` to fit the plotting structure. It handles the display of statistical significance using color coding and sizes dots based on the estimate value from the statistical test.
#'
#' @return
#' A list of ggplot objects, each representing a dot plot for a specific group. The plots can be further customized or directly rendered.
#'
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset
#' data("ecam.obj")
#' # Analyzing the impact of delivery method on microbial composition over months
#' result1 <- generate_alpha_per_time_test_long(
#'   data.obj = ecam.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "pielou"),
#'   time.var = "month_num",
#'   t0.level = unique(ecam.obj$meta.dat$month_num)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month_num)[-1],
#'   group.var = "delivery",
#'   adj.vars = "diet"
#' )
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_alpha_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = result1,
#'   group.var = "delivery",
#'   time.var = "month_num",
#'   t0.level = unique(ecam.obj$meta.dat$month_num)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month_num)[-1]
#' )
#'
#' # Example 2: Analyzing the Type 2 Diabetes dataset
#' data("subset_T2D.obj")
#' # Longitudinal analysis of microbial changes in different racial groups
#' result2 <- generate_alpha_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site"
#' )
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_alpha_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result2,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1]
#' )
#' }
#' @export
generate_alpha_per_time_dotplot_long <- function(data.obj,
                                        test.list,
                                        group.var,
                                        time.var,
                                        t0.level,
                                        ts.levels,
                                        base.size = 16,
                                        theme.choice = "bw",
                                        custom.theme = NULL,
                                        palette = rev(c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")),
                                        pdf = FALSE,
                                        pdf.wid = 7,
                                        pdf.hei = 5
){
  # Process the time variable in the data object.
  # This step ensures that the time variable is properly formatted for longitudinal analysis.
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

  # Select the appropriate theme based on user input or default to "bw".
  # This allows for flexible customization of the plot's appearance.
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Filter group names to only include those with "vs".
  # This step focuses on comparisons between groups, which are typically denoted with "vs" in the name.
  group.names <- names(test.list[[1]])
  group.names <- group.names[grep("vs", group.names)]

  # Process the test results for each group comparison.
  # This loop restructures the data for easier plotting.
  test.list <- lapply(group.names, function(group.names){

    # Define a function to merge data across time points for a given group comparison.
    merge_time_points <- function(time_test_list, group_name) {

      # Extract data for the current group comparison from each time point.
      data_list <- lapply(time_test_list, function(time_point_data) {
          return(time_point_data[[group_name]])
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

    # Apply the merging function to the current group comparison.
    merged_data_genus <- merge_time_points(test.list, group.names)

    return(merged_data_genus)
  })

  # Assign group names to the processed test results.
  names(test.list) <- group.names

  # Define the variable name for p-values in the data.
  p_val_var <- "P.Value"

  # Get the color palette for the plot.
  col <- mStat_get_palette(palette)

  # Generate a plot for each group comparison.
  plot.list <- lapply(group.names, function(group.names){

    # Extract data for the current group comparison.
    data_for_plot <- test.list[[group.names]]

    # Ensure that the time variable is properly ordered as a factor.
    data_for_plot[[time.var]] <- factor(data_for_plot[[time.var]], levels = time_levels)

    # Add a significance label based on the p-value.
    # This will be used to highlight statistically significant results in the plot.
    data_for_plot$Significance_Label <- ifelse(data_for_plot[[p_val_var]] < 0.05, "*", "")

    # Create the dot plot using ggplot2.
    dotplot <- ggplot(data_for_plot, aes(x = !!sym(time.var), y = Term, size = Estimate)) +
      geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
      geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
      scale_color_gradientn(colors = col) +
      labs(title = group.names,
           x = time.var,
           y = "Term",
           size = "Coefficient",
           color = p_val_var) +
      theme_to_use +
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
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = base.size),
        legend.title = element_text(size = base.size)
      )

    # Save the plot as a PDF if requested.
    if (pdf) {
      pdf_filename <- paste0("dotplot_", "alpha", "_", group.names, "_", time.var)
      if (!is.null(group.var)) {
        pdf_filename <- paste0(pdf_filename, "_group_", group.var)
      }
      pdf_filename <- paste0(pdf_filename, ".pdf")
      ggsave(pdf_filename, plot = dotplot, width = pdf.wid, height = pdf.hei)
    }

    return(dotplot)
  })

  # Assign group names to the list of plots.
  names(plot.list) <- group.names

  return(plot.list)
}