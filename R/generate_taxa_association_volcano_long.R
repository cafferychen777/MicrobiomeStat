#' Generate Volcano Plots for Taxa Association Test
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components
#'                such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree,
#'                and feature.agg.list (list). The data.obj can be converted from other formats using
#'                several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj',
#'                'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj',
#'                'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj',
#'                'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and
#'                'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj.
#'                Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var The grouping variable tested, found in metadata.
#' @param test.list The list of test results returned by generate_taxa_trend_test_long.
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
#' @param pdf (optional) Boolean, whether to save the plot as a PDF file.
#' @param pdf.wid (optional) Numeric, width of the saved PDF file.
#' @param pdf.hei (optional) Numeric, height of the saved PDF file.
#'
#' @return A list of ggplot objects of volcano plots for each taxonomic level.
#'
#' @examples
#' # Generate test results
#' # test_list <- generate_taxa_association_test_long(...)

#' # Create volcano plots
#' # volcano_plots <- generate_taxa_association_volcano_long(
#' #   data.obj, group.var, test_list,
#' #   feature.sig.level = 0.05, feature.mt.method = "fdr"
#' # )
#'
#' @importFrom dplyr distinct pull
#' @export
generate_taxa_association_volcano_long <-
  function(data.obj,
           group.var = NULL,
           test.list,
           feature.sig.level = 0.1,
           feature.mt.method = "fdr",
           palette = c("#F9F871", "#F4A261", "#FF6347"),
           pdf = FALSE,
           pdf.wid = 7,
           pdf.hei = 5) {

    meta_tab <- data.obj$meta.dat %>%
      dplyr::select(all_of(c(group.var))) %>% rownames_to_column("sample")

    # Define the custom color palette
    color_palette <- mStat_get_palette(palette)

    feature.level <- names(test.list)

    # 使用条件表达式设置要使用的p值变量
    p_val_var <-
      ifelse(feature.mt.method == "fdr",
             "Adjusted.P.Value",
             "P.Value")

    plot.list <- lapply(feature.level, function(feature.level) {
      sub_test.list <- test.list[[feature.level]]

      group_level <-
        meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

      reference_level <- group_level[1]

      test.result <- sub_test.list %>%
        dplyr::filter(grepl(
          paste0('^', group.var, '.*', NULL, '$'),
          Output.Element
        ))

      sub_plot.list <-
        lapply(group_level[-1], function(group_level) {

          sub_test.result <- test.result %>%
            dplyr::filter(Output.Element == paste0(group.var, group_level))

          # Find max absolute log2FoldChange for symmetric x-axis
          max_abs_log2FC <-
            max(abs(sub_test.result$Log2.Fold.Change), na.rm = TRUE)

          p <-
            ggplot(sub_test.result,
                   aes(
                     x = Log2.Fold.Change,
                     y = -log10(get(p_val_var)),
                     color = -log10(get(p_val_var))
                   )) +
            geom_point(aes(shape = get(p_val_var) < feature.sig.level), size = 7) +
            geom_vline(
              aes(xintercept = 0),
              linetype = "dashed",
              linewidth = 1.5,
              color = "grey"
            ) +
            geom_hline(
              aes(yintercept = -log10(feature.sig.level)),
              linetype = "dashed",
              linewidth = 1.5,
              color = "grey"
            ) +
            geom_text(
              aes(label = ifelse(
                get(p_val_var) < feature.sig.level,
                as.character(Variable),
                ''
              )),
              vjust = -0.5,
              hjust = 0.5,
              size = 3.5
            ) +
            scale_shape_manual(values = c(16, 17)) +
            labs(
              title = paste(group_level, "vs", reference_level, "(Reference)"),
              x = "Log2 Fold Change",
              y = "-log10(p-value)",
              shape = "Significant",
              color = "-log10(p-value)"
            ) +
            theme_minimal() +
            theme(
              plot.title.position = "plot",
              plot.title = element_text(hjust = 0.5),
              panel.grid.major = element_line(color = "grey", linetype = "dashed"),
              panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
              legend.position = "bottom"
            ) +
            scale_color_gradientn(colors = color_palette) +
            coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))

          if (pdf) {
            pdf_filename <- paste0("volcano_", feature.level, "_", group_level, ".pdf")
            ggsave(pdf_filename, plot = p, width = pdf.wid, height = pdf.hei)
          }

          return(p)
        })

      names(sub_plot.list) <-
        paste(group_level[-1], "vs", reference_level, "(Reference)")

      return(sub_plot.list)
    })


    names(plot.list) <- feature.level
    return(plot.list)
  }
