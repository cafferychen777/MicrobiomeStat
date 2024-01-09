#' Longitudinal Taxa Dot Plot for Microbiome Data
#'
#' This function generates dot plots for visualizing the results of longitudinal taxa tests in microbiome data. It focuses on displaying how the abundance and significance of microbial taxa change over time across different groups.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param time.test.list A list of test results from `generate_taxa_test_long` function, representing the statistical analysis of taxa over time.
#' @param group.var A string specifying the group variable in meta.dat for between-group comparisons in the dot plot.
#' @param time.var A string representing the time variable in the meta.dat. This variable is used for x-axis in the dot plot.
#' @param t0.level Optional; a baseline time level for time series analysis.
#' @param ts.levels Optional; time series levels for analysis.
#' @param feature.level A string or vector of strings indicating the taxonomic level(s) for analysis (e.g., "Phylum", "Class") and for displaying in the dot plot.
#' @param feature.sig.level The significance level cutoff for highlighting taxa
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none"
#' @param base.size Numeric; base font size for text elements in the plot. Default is 16.
#' @param theme.choice A string indicating the choice of ggplot2 theme for the plot. Supported themes include 'prism', 'classic', 'gray', and 'bw'. Default is 'bw'.
#' @param custom.theme Optional; a custom ggplot2 theme. If provided, it overrides the standard theme choice.
#' @param palette Optional; a vector of colors for the gradient scale. If provided, it overrides the default color scale.
#' @param pdf Boolean; whether to save the plot as a PDF file.
#' @param pdf.wid Numeric; width of the saved PDF file.
#' @param pdf.hei Numeric; height of the saved PDF file.
#' @details
#' The function uses ggplot2 to create dot plots, with options for customizing themes and aesthetics. It visualizes the significance and magnitude of changes in microbial taxa, indicated by color and size of the dots. The function allows for detailed customization of plot features, including base font size and theme choice.
#'
#' The function prepares data by merging and transforming the results from `generate_taxa_test_long` to fit the plotting structure. It handles the display of statistical significance using color coding and sizes dots based on the coefficient value from the statistical test.
#'
#' @return
#' A list of ggplot objects, each representing a dot plot for a specific taxonomic level and group. The plots can be further customized or directly rendered.
#'
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset
#' data("ecam.obj")
#' # Analyzing the impact of delivery method on microbial composition over months
#' result1 <- generate_taxa_test_long(
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
#' dotplot_ecam <- generate_taxa_dotplot_long(
#'   data.obj = ecam.obj,
#'   time.test.list = result1,
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
#' result2 <- generate_taxa_test_long(
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
#' dotplot_T2D <- generate_taxa_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   time.test.list = result2,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   feature.level = c("Genus", "Family")
#' )
#' }
#' @export
generate_taxa_dotplot_long <- function(data.obj,
                                       time.test.list,
                                       group.var,
                                       time.var,
                                       t0.level = NULL,
                                       ts.levels = NULL,
                                       feature.level,
                                       feature.mt.method = "none",
                                       feature.sig.level = 0.05,
                                       base.size = 16,
                                       theme.choice = "bw",
                                       custom.theme = NULL,
                                       palette = NULL,
                                       pdf = FALSE,
                                       pdf.wid = 7,
                                       pdf.hei = 5
                                       ){

  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  time_levels <- data.obj$meta.dat %>%
    dplyr::select(!!dplyr::sym(time.var)) %>%
    dplyr::distinct() %>%
    dplyr::pull()

  # Assuming mStat_get_theme function is already defined
  # Replace the existing theme selection code with this:
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  group.names <- names(time.test.list[[1]][[1]])

  p_val_var <- if (feature.mt.method == "fdr") "Adjusted.P.Value" else "P.Value"

  test.list <- lapply(feature.level, function(feature.level){
    sub.test.list <- lapply(group.names, function(group.names){
      merge_time_points <- function(time_test_list, feature_level, group_name) {

        data_list <- lapply(time_test_list, function(time_point_data) {

          if (feature_level %in% names(time_point_data) &&
              group_name %in% names(time_point_data[[feature_level]])) {
            return(time_point_data[[feature_level]][[group_name]])
          }
          return(NULL)
        })

        data_list <- Filter(Negate(is.null), data_list)

        merged_data <- dplyr::bind_rows(data_list, .id = time.var)
        return(merged_data)
      }

      merged_data_genus <- merge_time_points(time.test.list, feature.level, group.names)
      return(merged_data_genus)
    })
    names(sub.test.list) <- group.names
    return(sub.test.list)
  })

  names(test.list) <- feature.level

  plot.list <- lapply(feature.level, function(feature.level){

    sub_plot.list <- lapply(group.names, function(group.names){

      data_for_plot <- test.list[[feature.level]][[group.names]]

      data_for_plot[[time.var]] <- factor(data_for_plot[[time.var]], levels = time_levels)

      data_for_plot$Significance_Label <- ifelse(data_for_plot[[p_val_var]] < feature.sig.level, "*", "")

      dotplot <- ggplot(data_for_plot, aes(x = !!sym(time.var), y = Variable, size = Coefficient)) +
        geom_point(aes(color = !!sym(p_val_var)), alpha = 0.6, shape = 19) +
        geom_text(aes(label = Significance_Label), vjust = 0.8, show.legend = FALSE, color = "white") +
        scale_color_gradientn(colors = rev(c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020"))) +
        labs(title = group.names,
             x = time.var,
             y = feature.level,
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
