#' Generate Volcano Plots for Longitudinal Taxa Abundance Volatility Test
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param group.var The grouping variable tested, found in metadata.
#' @param test.list The list of test results returned by generate_taxa_volatility_test_long.
#' @param feature.sig.level The significance level cutoff for highlighting taxa.
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none".
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
#'
#' @examples
#' # data("subset_T2D.obj")
#' # test_list <- generate_taxa_volatility_test_long(data.obj, ...)
#' # volcano_plots <- generate_taxa_volatility_volcano_long(data.obj,
#'                                                       # group.var,
#'                                                       # test.list,
#'                                                       # feature.sig.level = 0.05,
#'                                                       # feature.mt.method = "fdr")
#'
#' @importFrom dplyr pull
#' @export
generate_taxa_volatility_volcano_long <- function(data.obj,
                                                  group.var,
                                                  test.list,
                                                  feature.sig.level = 0.1,
                                                  feature.mt.method = c("fdr","none"),
                                                  palette = c("#F9F871", "#F4A261", "#FF6347"),
                                                  pdf = FALSE,
                                                  pdf.wid = 7,
                                                  pdf.hei = 5){

  meta_tab <- data.obj$meta.dat %>%
    dplyr::select(all_of(c(
      group.var
    ))) %>% rownames_to_column("sample")

  feature.level <- names(test.list)

  # Define the custom color palette
  color_palette <- mStat_get_palette(palette)

  group_level <- meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

  reference_level <- group_level[1]

  # 使用条件表达式设置要使用的p值变量
  p_val_var <-
    ifelse(feature.mt.method == "fdr",
           "Adjusted.P.Value",
           "P.Value")

  plot.list <- lapply(feature.level, function(feature.level) {
    sub_test.list <- test.list[[feature.level]]

    if (length(group_level) > 2) {
      group.levels <- names(sub_test.list)[-length(names(sub_test.list))]
    } else {
      group.levels <- names(sub_test.list)
    }

    sub_plot.list <-
      lapply(group.levels, function(group.level) {

        sub_test.result <- sub_test.list[[group.level]]

        # Find max absolute log2FoldChange for symmetric x-axis
        max_abs_log2FC <-
          max(abs(sub_test.result$Coefficient), na.rm = TRUE)

        p <-
          ggplot(sub_test.result, aes(x = Coefficient, y = -log10(get(p_val_var)),
                                      color = Prevalence, size = Mean.Abundance)) +
          geom_point() +
          geom_vline(aes(xintercept = 0), linetype = "dashed", linewidth = 1.5, color = "grey") +
          geom_hline(aes(yintercept = -log10(feature.sig.level)), linetype = "dashed", linewidth = 1.5, color = "grey") +
          geom_text(aes(label = ifelse(get(p_val_var) < feature.sig.level, as.character(Variable), '')),
                    vjust = -0.5, hjust = 0.5, size = 3.5) +
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

        # 判断是否保存为 PDF
        if (pdf) {
          pdf_filename <- paste0("volcano_", feature.level, "_", group.level, ".pdf")
          ggsave(pdf_filename, plot = p, width = pdf.wid, height = pdf.hei)
        }

        return(p)
      })

    names(sub_plot.list) <-
      group.levels

    return(sub_plot.list)
  })

  names(plot.list) <- feature.level

  return(plot.list)
}
