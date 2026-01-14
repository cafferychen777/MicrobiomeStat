#' Dot Plot for Beta Diversity Change Tests Over Time
#'
#' Visualizes results from longitudinal beta diversity change tests as dot plots.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param test.list A list of test results from \code{\link{generate_beta_change_per_time_test_long}}.
#' @param group.var Character string specifying the group variable for comparisons.
#' @details
#' The function uses ggplot2 to create dot plots, with options for customizing themes and aesthetics. It visualizes the significance and magnitude of changes in beta diversity indices, indicated by color and size of the dots. The function allows for detailed customization of plot features, including base font size and theme choice.
#'
#' The function prepares data by merging and transforming the results from `generate_beta_test_long` to fit the plotting structure. It handles the display of statistical significance using color coding and sizes dots based on the estimate value from the statistical test.
#'
#' @return
#' A list of ggplot objects, each representing a dot plot for a specific group. The plots can be further customized or directly rendered.
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' result1 <- generate_beta_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#'   subject.var = "subject_id",
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_beta_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result1,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number_num)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number_num)[-1],
#' )
#'
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, c('BC', 'Jaccard'))
#' result2 <- generate_beta_change_per_time_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   subject.var = "subject.id",
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_beta_per_time_dotplot_long(
#'   data.obj = ecam.obj,
#'   test.list = result2,
#'   group.var = "delivery",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 15
#' )
#' }
#' @export
generate_beta_per_time_dotplot_long <- function(data.obj,
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
  # Process the time variable in the data object
  # This ensures that the time variable is properly formatted and ordered
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

  # Get the appropriate theme for the plot
  # This allows for consistent styling across different plots
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Extract group names from the test list
  # These names typically represent different comparison groups
  group.names <- names(test.list[[1]])

  # Process the test list to merge data across time points for each group
  test.list <- lapply(group.names, function(group.names){
    # Define a function to merge data from different time points
    merge_time_points <- function(time_test_list, group_name) {
      # Extract data for the current group from each time point
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

    # Apply the merge function to the current group
    merged_data_genus <- merge_time_points(test.list, group.names)
    return(merged_data_genus)
  })

  # Assign group names to the processed test list
  names(test.list) <- group.names

  # Define the variable name for p-values in the data
  p_val_var <- "P.Value"

  # Get the color palette for the plot
  col <- mStat_get_palette(palette)

  # Create a list of plots, one for each group
  plot.list <- lapply(group.names, function(group.names){
    # Extract data for the current group
    data_for_plot <- test.list[[group.names]]

    # Ensure the time variable is a factor with the correct order of levels
    data_for_plot[[time.var]] <- factor(data_for_plot[[time.var]], levels = time_levels)

    # Add a significance label based on the p-value
    # This will be used to add asterisks to significant results in the plot
    data_for_plot$Significance_Label <- ifelse(data_for_plot[[p_val_var]] < 0.05, "*", "")

    # Create the dot plot using ggplot2
    # Check if time variable contains only numeric values
    if (all(grepl("^\\d+(\\.\\d+)?$", time_levels))) {
      # Convert to numeric for proper spacing
      data_for_plot[[paste0(time.var, "_numeric")]] <- as.numeric(as.character(data_for_plot[[time.var]]))
      
      # Create plot with numeric x-axis
      dotplot <- ggplot(data_for_plot, aes(x = !!sym(paste0(time.var, "_numeric")), 
                                           y = Term, size = Estimate)) +
        # Add points, colored by p-value
        geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
        # Add significance labels
        geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
        # Set color scale for p-values
        scale_color_gradientn(colors = col) +
        # Set numeric x-axis with proper breaks and labels
        scale_x_continuous(breaks = as.numeric(time_levels), labels = time_levels) +
        # Set plot labels
        labs(title = group.names,
             x = time.var,
             y = "Term",
             size = "Coefficient",
             color = p_val_var) +
        # Set size scale for estimates
        scale_radius(range = c(0, 10)) +
        # Apply the chosen theme
        theme_to_use
    } else {
      # Keep current factor-based approach for non-numeric time values
      dotplot <- ggplot(data_for_plot, aes(x = !!sym(time.var), y = Term, size = Estimate)) +
        # Add points, colored by p-value
        geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
        # Add significance labels
        geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
        # Set color scale for p-values
        scale_color_gradientn(colors = col) +
        # Set plot labels
        labs(title = group.names,
             x = time.var,
             y = "Term",
             size = "Coefficient",
             color = p_val_var) +
        # Set size scale for estimates
        scale_radius(range = c(0, 10)) +
        # Apply the chosen theme
        theme_to_use
    }
    
    # Apply common theme settings to the plot
    dotplot <- dotplot +
      # Customize theme elements
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
      pdf_filename <- paste0("dotplot_", "beta", "_", group.names, "_", time.var)
      if (!is.null(group.var)) {
        pdf_filename <- paste0(pdf_filename, "_group_", group.var)
      }
      pdf_filename <- paste0(pdf_filename, ".pdf")
      ggsave(pdf_filename, plot = dotplot, width = pdf.wid, height = pdf.hei)
    }

    return(dotplot)
  })

  # Assign group names to the list of plots
  names(plot.list) <- group.names

  # Filter the plot list to include only plots with "vs" in their names
  # This typically selects plots that compare different groups
  plot.list <- plot.list[grep("vs", names(plot.list), value = TRUE)]

  return(plot.list)
}