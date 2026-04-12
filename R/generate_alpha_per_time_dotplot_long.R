#' @title Alpha Diversity Test Results Dot Plot (Longitudinal)
#'
#' @description Generates dot plots visualizing the results of longitudinal alpha
#'   diversity tests, showing effect sizes and p-values across time points.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param test.list A list of test results from `generate_alpha_per_time_test_long`
#'   or `generate_alpha_change_test_long` function.
#'
#' @return A list of ggplot objects, each representing a dot plot for a specific group. The plots can be further customized or directly rendered.
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
  # Order time levels using the shared time semantics helper.
  ordered_test_entries <- mStat_order_named_time_entries(test.list)
  time_levels <- names(ordered_test_entries)
  test.list <- unname(ordered_test_entries)
  names(test.list) <- time_levels

  # Select the appropriate theme based on user input or default to "bw".
  # This allows for flexible customization of the plot's appearance.
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Keep all group-term names from all time points.
  # This includes both categorical contrasts ("... vs ...") and continuous terms.
  group.names <- mStat_collect_time_result_names(test.list)

  # Process the test results for each group comparison.
  # This loop restructures the data for easier plotting.
  test.list <- mStat_bind_time_results(
    test.list = test.list,
    result_names = group.names,
    time.var = time.var
  )

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
    # Check if time variable contains only numeric values
    if (all(grepl("^\\d+(\\.\\d+)?$", time_levels))) {
      # Convert to numeric for proper spacing
      data_for_plot[[paste0(time.var, "_numeric")]] <- as.numeric(as.character(data_for_plot[[time.var]]))
      
      # Create plot with numeric x-axis
      dotplot <- ggplot(data_for_plot, aes(x = !!sym(paste0(time.var, "_numeric")), 
                                           y = Term, size = Estimate)) +
        geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
        geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
        scale_color_gradientn(colors = col) +
        scale_x_continuous(breaks = as.numeric(time_levels), labels = time_levels) +
        labs(title = group.names,
             x = time.var,
             y = "Term",
             size = "Coefficient",
             color = p_val_var) +
        theme_to_use
    } else {
      # Keep current factor-based approach for non-numeric time values
      dotplot <- ggplot(data_for_plot, aes(x = !!sym(time.var), y = Term, size = Estimate)) +
        geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
        geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
        scale_color_gradientn(colors = col) +
        labs(title = group.names,
             x = time.var,
             y = "Term",
             size = "Coefficient",
             color = p_val_var) +
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