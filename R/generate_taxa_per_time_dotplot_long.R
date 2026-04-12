#' Longitudinal Taxa Dot Plot
#'
#' Generates dot plots visualizing longitudinal taxa test results, showing abundance
#' and significance changes over time across groups.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param test.list List of test results from generate_taxa_per_time_test_long.
#' @param feature.sig.level Numeric; significance level cutoff for highlighting taxa.
#' @param feature.mt.method Character; multiple testing method ("fdr" or "none").
#' @param filter_significant Logical; whether to show only significant features. Default TRUE.
#' @param features.plot Character vector of features to plot. Overrides filter_significant.
#'
#' @return A list of ggplot objects for each taxonomic level and group.
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

  # Order time levels using the shared time semantics helper.
  ordered_test_entries <- mStat_order_named_time_entries(test.list)
  time_levels <- names(ordered_test_entries)
  test.list <- unname(ordered_test_entries)
  names(test.list) <- time_levels

  # Get the appropriate theme for plotting
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Determine which p-value to use based on the multiple testing correction method
  p_val_var <- if (feature.mt.method == "fdr") "Adjusted.P.Value" else "P.Value"

  # Process test results for each feature level and group.
  test.list <- lapply(feature.level, function(feature.level) {
    feature_group_names <- mStat_collect_time_result_names(test.list, entry = feature.level)
    mStat_bind_time_results(
      test.list = test.list,
      result_names = feature_group_names,
      time.var = time.var,
      entry = feature.level
    )
  })

  names(test.list) <- feature.level

  # Generate plots for each feature level
  plot.list <- lapply(feature.level, function(feature.level){

    # Generate plots for each group within the current feature level
    feature_group_names <- names(test.list[[feature.level]])

    sub_plot.list <- lapply(feature_group_names, function(group.names){

      # Prepare data for plotting
      data_for_plot <- test.list[[feature.level]][[group.names]]

      # Ensure time variable is a factor with correct levels
      data_for_plot[[time.var]] <- factor(data_for_plot[[time.var]], levels = time_levels)

      # Add significance labels based on p-values
      data_for_plot$Significance_Label <- ifelse(data_for_plot[[p_val_var]] < feature.sig.level, "*", "")

      # Filter features based on user specifications
      data_for_plot <- mStat_filter_test_result_features(data_for_plot, features.plot)
      if (is.null(features.plot) && filter_significant) {
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

    names(sub_plot.list) <- feature_group_names

    return(sub_plot.list)
  })

  names(plot.list) <- feature.level

  return(plot.list)

}