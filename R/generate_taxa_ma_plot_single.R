#' Generate MA Plots for Taxa Differential Test for a Single Time Point
#'
#' This function creates MA plots to visualize differential abundance analysis results.
#' MA plots display the log2 fold change (M) against the average abundance (A) for each taxon,
#' helping to identify taxa with significant abundance differences between groups.
#' 
#' @importFrom ggplot2 aes scale_color_manual labs theme element_text scale_size_continuous
#' @importFrom dplyr filter mutate
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var The grouping variable tested, found in metadata.
#' @param test.list The list of test results returned by generate_taxa_test_single.
#' @param feature.sig.level The significance level cutoff for highlighting taxa.
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none".
#' @param features.plot A character vector of taxa to be plotted. If NULL, all taxa will be plotted.
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
#' @return A list of ggplot objects of MA plots for each taxonomic level.
#' @details
#' The function generates MA plots for each taxonomic level based on the test results. 
#' It visualizes the relationship between average abundance and log2 fold change, 
#' highlighting statistically significant changes.
#'
#' MA plots are useful for identifying abundance-dependent biases in the data and 
#' for visualizing which taxa show significant differential abundance between conditions.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' test.list <- generate_taxa_test_single(
#'   data.obj = subset_T2D.obj,
#'   time.var = "visit_number",
#'   time.point = "1",
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   feature.level = c("Genus","Family"),
#'   feature.dat.type = c("count")
#' )
#'
#' plot.list <- generate_taxa_ma_plot_single(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   test.list = test.list,
#'   feature.sig.level = 0.05,
#'   feature.mt.method = "fdr"
#' )
#' }
#' @export
generate_taxa_ma_plot_single <- function(data.obj,
                                        group.var,
                                        test.list,
                                        feature.sig.level = 0.05,
                                        feature.mt.method = "fdr",
                                        features.plot = NULL,
                                        palette = NULL,
                                        pdf = FALSE,
                                        pdf.wid = 11,
                                        pdf.hei = 8.5) {
  
  # Data validation
  if (is.null(data.obj$meta.dat)) {
    stop("Meta data not found in data.obj.")
  }
  
  if (!(group.var %in% colnames(data.obj$meta.dat))) {
    stop(paste("Group variable", group.var, "not found in metadata."))
  }
  
  if (is.null(test.list)) {
    stop("Test list is NULL. Please run generate_taxa_test_single first.")
  }
  
  # Get the feature levels from test.list
  feature.level <- names(test.list)
  
  # Set up color palette
  color_palette <- mStat_get_palette(palette)
  
  # Default theme
  theme_to_use <- mStat_get_theme(theme.choice = "bw")
  
  # Create a list to store plots for each feature level
  plot.list <- lapply(feature.level, function(level) {
    
    # Extract test results for the current feature level
    sub_test.list <- test.list[[level]]
    
    # Create a list to store plots for each group comparison
    sub_plot.list <- lapply(names(sub_test.list), function(group.level) {
      
      # Extract test data for the current group comparison
      df <- sub_test.list[[group.level]]
      
      # Skip if dataframe is empty
      if (nrow(df) == 0) {
        message(paste0("No data available for ", level, " at ", group.level, "."))
        return(NULL)
      }
      
      # Filter features if specified
      if (!is.null(features.plot)) {
        df <- df[df$Variable %in% features.plot, ]
        if (nrow(df) == 0) {
          message(paste0("None of the specified features found for ", level, " at ", group.level, "."))
          return(NULL)
        }
      }
      
      # Determine p-value column based on multiple testing method
      p_val_var <- if (feature.mt.method == "fdr") "Adjusted.P.Value" else "P.Value"
      
      # Add significance indicator
      df$Significant <- ifelse(df[[p_val_var]] < feature.sig.level, "Yes", "No")
      
      # Calculate maximum absolute coefficient for symmetric x-axis limits
      max_abs_log2FC <- max(abs(df$Coefficient), na.rm = TRUE)
      
      # Create MA plot
      p <- ggplot(df, aes(x = Mean.Abundance, 
                          y = Coefficient, 
                          color = Significant,
                          size = Prevalence)) +
        geom_point(alpha = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
        # Add horizontal lines at significance thresholds if needed
        # geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "darkgray") +
        ggrepel::geom_text_repel(
          data = subset(df, df[[p_val_var]] < feature.sig.level),
          aes(label = Variable),
          size = 3.5,
          box.padding = 0.35,
          point.padding = 0.5,
          max.overlaps = 20
        ) +
        # Use colors from the default palette for significant/non-significant points
        scale_color_manual(values = c("Yes" = color_palette[1], "No" = color_palette[2])) +
        labs(
          title = paste("MA Plot -", group.level),
          x = "Average Abundance (A)",
          y = "Log2 Fold Change (M)",
          color = "Significant",
          size = "Prevalence"
        ) +
        theme_to_use +
        theme(
          plot.title.position = "plot",
          plot.title = element_text(hjust = 0.5, size = 12),
          panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
          panel.grid.minor = element_line(color = "grey98", linetype = "dotted"),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)
        ) +
        # Set the size range for prevalence
        scale_size_continuous(range = c(2, 6))
      
      # Save the plot as PDF if requested
      if (pdf) {
        pdf_filename <- paste0("ma_plot_", level, "_", group.level, ".pdf")
        ggsave(pdf_filename, plot = p, width = pdf.wid, height = pdf.hei)
      }
      
      return(p)
    })
    
    # Assign names to the sub-plot list based on group levels
    names(sub_plot.list) <- names(sub_test.list)
    
    # Remove NULL elements (if any)
    sub_plot.list <- sub_plot.list[!sapply(sub_plot.list, is.null)]
    
    return(sub_plot.list)
  })
  
  # Assign names to the plot list based on feature levels
  names(plot.list) <- feature.level
  
  # Remove empty elements (if any)
  plot.list <- plot.list[!sapply(plot.list, function(x) length(x) == 0)]
  
  return(plot.list)
}