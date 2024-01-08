#' Longitudinal Alpha Diversity Dot Plot for Microbiome Data
#'
#' This function generates dot plots for visualizing the results of longitudinal alpha diversity tests in microbiome data. It focuses on displaying how the alpha diversity indices change over time across different groups.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param test.list A list of test results from `generate_alpha_test_long` or `generate_alpha_change_test_long` function, representing the statistical analysis of alpha diversity over time.
#' @param group.var A string specifying the group variable in meta.dat for between-group comparisons in the dot plot.
#' @param time.var A string representing the time variable in the meta.dat. This variable is used for the x-axis in the dot plot.
#' @param t0.level The baseline time level for time series analysis.
#' @param ts.levels The time series levels for analysis.
#' @param base.size Optional; numeric, base font size for text elements in the plot. Default is 16.
#' @param theme.choice Optional; a string indicating the choice of ggplot2 theme for the plot. Supported themes include 'prism', 'classic', 'gray', and 'bw'. Default is 'bw'.
#' @param custom.theme Optional; a custom ggplot2 theme. If provided, it overrides the standard theme choice.
#' @param palette Optional; a vector of colors for the gradient scale. If provided, it overrides the default color scale.
#' @param pdf Optional; boolean, whether to save the plot as a PDF file.
#' @param pdf.wid Optional; numeric, width of the saved PDF file.
#' @param pdf.hei Optional; numeric, height of the saved PDF file.
#' @details
#' The function uses ggplot2 to create dot plots, with options for customizing themes and aesthetics. It visualizes the significance and magnitude of changes in alpha diversity indices, indicated by color and size of the dots. The function allows for detailed customization of plot features, including base font size and theme choice.
#'
#' The function prepares data by merging and transforming the results from `generate_alpha_test_long` to fit the plotting structure. It handles the display of statistical significance using color coding and sizes dots based on the estimate value from the statistical test.
#'
#' @return
#' A list of ggplot objects, each representing a dot plot for a specific group. The plots can be further customized or directly rendered.
#'
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset
#' data("ecam.obj")
#' # Analyzing the impact of delivery method on microbial composition over months
#' result1 <- generate_alpha_test_long(
#'   data.obj = ecam.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "pielou"),
#'   time.var = "month_num",
#'   t0.level = unique(ecam.obj$meta.dat$month_num)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month_num)[-1],
#'   group.var = "delivery",
#'   adj.vars = "diet"
#' )
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_alpha_dotplot_long(
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
#' result2 <- generate_alpha_test_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site"
#' )
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_alpha_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result2,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1]
#' )
#' }
#' @export
generate_alpha_dotplot_long <- function(data.obj,
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
  # Process the time variable in the data
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  # Extract unique time levels
  time_levels <- unique(data.obj$meta.dat[[time.var]])

  # Determine the theme to use for plots
  theme_function <- switch(
    theme.choice,
    prism = ggprism::theme_prism(),
    classic = theme_classic(),
    gray = theme_gray(),
    bw = theme_bw(),
    ggprism::theme_prism()
  )

  theme_to_use <-
    if (!is.null(custom.theme))
      custom.theme else theme_function

  group.names <- names(test.list[[1]])

  test.list <- lapply(group.names, function(group.names){
    merge_time_points <- function(time_test_list, group_name) {

      data_list <- lapply(time_test_list, function(time_point_data) {

          return(time_point_data[[group_name]])

      })

      data_list <- Filter(Negate(is.null), data_list)

      merged_data <- dplyr::bind_rows(data_list, .id = time.var)
      return(merged_data)
    }

    merged_data_genus <- merge_time_points(test.list, group.names)
    return(merged_data_genus)
  })

  names(test.list) <- group.names

  p_val_var <- "P.Value"

  col <- mStat_get_palette(palette)

  plot.list <- lapply(group.names, function(group.names){

    data_for_plot <- test.list[[group.names]]

    data_for_plot[[time.var]] <- factor(data_for_plot[[time.var]], levels = time_levels)

    data_for_plot$Significance_Label <- ifelse(data_for_plot[[p_val_var]] < 0.05, "*", "")

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

  names(plot.list) <- group.names

  return(plot.list)
}
