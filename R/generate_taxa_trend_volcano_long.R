#' Generate Volcano Plots for Longitudinal Taxa Trend Test
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param group.var The grouping variable tested, found in metadata.
#' @param time.var The time variable used in the analysis.
#' @param test.list The list of test results returned by generate_taxa_trend_test_long.
#' @param feature.sig.level The significance level cutoff for highlighting taxa.
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none".
#' @param features.plot An optional parameter specifying the features to be plotted. If provided, only these features will be plotted.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf Boolean; whether to save the plot as a PDF file.
#' @param pdf.wid Numeric; width of the saved PDF file.
#' @param pdf.hei Numeric; height of the saved PDF file.
#' @return A list of ggplot objects of volcano plots for each taxonomic level.
#' @details
#' The function generates volcano plots for each taxonomic level based on the test results. It visualizes the longitudinal trends of taxa abundance over time and highlights statistically significant changes.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' test.list <- generate_taxa_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus","Family"),
#'   feature.dat.type = c("count")
#' )
#'
#' plot.list <- generate_taxa_trend_volcano_long(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   time.var = "visit_number_num",
#'   test.list = test.list,
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "none")
#'
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = NULL,
#'   feature.level = c("Phylum","Class"),
#'   feature.dat.type = c("proportion")
#' )
#' }
#'
#' @importFrom dplyr distinct pull
#' @export
generate_taxa_trend_volcano_long <-
  function(data.obj,
           group.var = NULL,
           time.var = NULL,
           test.list,
           feature.sig.level = 0.1,
           feature.mt.method = "fdr",
           features.plot = NULL,
           palette = c("white", "#7FB695", "#006D2C"),
           pdf = FALSE,
           pdf.wid = 7,
           pdf.hei = 5) {

    # Extract relevant metadata and add sample names as a column
    meta_tab <- data.obj$meta.dat %>%
      dplyr::select(all_of(c(group.var))) %>% rownames_to_column("sample")

    # Get the color palette for the plot
    color_palette <- mStat_get_palette(palette)

    # Extract feature levels from the test list
    feature.level <- names(test.list)

    # Determine which p-value to use based on multiple testing correction method
    p_val_var <-
      ifelse(feature.mt.method == "fdr",
             "Adjusted.P.Value",
             "P.Value")

    # Generate plots for each feature level
    plot.list <- lapply(feature.level, function(feature.level) {
      sub_test.list <- test.list[[feature.level]]

      # Handle cases where group variable is provided
      if (!is.null(group.var)) {
        # Extract group levels from metadata
        group_level <-
          meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

        # Set the reference level as the first group level
        reference_level <- group_level[1]

        # Generate plots for each group comparison
        sub_plot.list <-
          lapply(names(sub_test.list), function(group.level) {

            sub_test.result <- sub_test.list[[group.level]]
            
            # Filter features if a specific set is provided
            if (!is.null(features.plot)) {
              sub_test.result <- sub_test.result %>% 
                filter(Variable %in% features.plot)
            }

            # Calculate the maximum absolute coefficient for symmetric x-axis
            max_abs_log2FC <-
              max(abs(sub_test.result$Coefficient), na.rm = TRUE)

            # Create the volcano plot
            p <-
              ggplot(sub_test.result, aes(x = Coefficient, y = -log10(get(p_val_var)),
                                          color = Prevalence, size = Mean.Abundance)) +
              geom_point() +
              geom_vline(aes(xintercept = 0), linetype = "dashed", linewidth = 1.5, color = "grey") +
              geom_hline(aes(yintercept = -log10(feature.sig.level)), linetype = "dashed", linewidth = 1.5, color = "grey") +
              ggrepel::geom_text_repel(
                aes(label = ifelse(get(p_val_var) < feature.sig.level, as.character(Variable), '')),
                color = "black",
                size = 3.5,
                box.padding = 0.35,
                point.padding = 0.5,
                show.legend = FALSE
              ) +
              scale_shape_manual(values = c(16, 17)) +
              labs(title = group.level, x = "Coefficient", y = "-log10(p-value)", color = "Prevalence", size = "Mean Abundance") +
              theme_bw() +
              theme(
                plot.title.position = "plot",
                plot.title = element_text(hjust = 0.5, size = 12),
                panel.grid.major = element_line(color = "grey", linetype = "dashed"),
                panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
                legend.position = "bottom",
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 14)
              ) +
              scale_color_gradientn(colors = color_palette) +
              scale_size_continuous(range = c(3, 7)) +
              coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))

            # Save the plot as PDF if requested
            if (pdf) {
              pdf_filename <- paste0("volcano_", feature.level, "_", group.level, ".pdf")
              ggsave(pdf_filename, plot = p, width = pdf.wid, height = pdf.hei)
            }

            return(p)
          })

        names(sub_plot.list) <-
          names(sub_test.list)
      } else {
        # Handle cases where no group variable is provided (time-based analysis)
        sub_test.result <- sub_test.list[[time.var]]
        
        # Filter features if a specific set is provided
        if (!is.null(features.plot)) {
          sub_test.result <- sub_test.result %>% 
            filter(Variable %in% features.plot)
        }

        # Calculate the maximum absolute coefficient for symmetric x-axis
        max_abs_log2FC <-
          max(abs(sub_test.result$Coefficient), na.rm = TRUE)

        # Generate plot for time-based analysis
        sub_plot.list <- lapply("time", function(time) {
          p <-
            ggplot(sub_test.result, aes(x = Coefficient, y = -log10(get(p_val_var)),
                                        color = Prevalence, size = Mean.Abundance)) +
            geom_point() +
            geom_vline(aes(xintercept = 0), linetype = "dashed", linewidth = 1.5, color = "grey") +
            geom_hline(aes(yintercept = -log10(feature.sig.level)), linetype = "dashed", linewidth = 1.5, color = "grey") +
            ggrepel::geom_text_repel(
              aes(label = ifelse(get(p_val_var) < feature.sig.level, as.character(Variable), '')),
              color = "black",
              size = 3.5,
              box.padding = 0.35,
              point.padding = 0.5,
              show.legend = FALSE
            ) +
            scale_shape_manual(values = c(16, 17)) +
            labs(x = "Coefficient", y = "-log10(p-value)", color = "Mean Prevalence", size = "Mean Abundance") +
            theme_bw() +
            theme(
              plot.title.position = "plot",
              plot.title = element_text(hjust = 0.5, size = 12),
              panel.grid.major = element_line(color = "grey", linetype = "dashed"),
              panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
              legend.position = "bottom",
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14)
            ) +
            scale_color_gradientn(colors = color_palette) +
            scale_size_continuous(range = c(3, 7)) +
            coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))

          return(p)
        })
        names(sub_plot.list) <-
          time.var
      }
      return(sub_plot.list)
    })

    # Assign names to the plot list based on feature levels
    names(plot.list) <- feature.level
    return(plot.list)
  }