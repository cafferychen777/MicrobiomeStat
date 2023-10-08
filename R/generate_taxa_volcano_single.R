#' Generate volcano plots for taxa differential test for a single time point
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var The grouping variable tested, found in metadata
#' @param test.list The list of test results returned by generate_taxa_trend_test_long
#' @param feature.sig.level The significance level cutoff for highlighting taxa
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none"
#'
#' @return A list of ggplot objects of volcano plots for each taxonomic level
#'
#' @examples
#' data(peerj32.obj)
#' test.list <- generate_taxa_test_single(
#'     data.obj = peerj32.obj,
#'     time.var = "time",
#'     t.level = "2",
#'     group.var = "group",
#'     adj.vars = "sex",
#'     feature.dat.type = "count",
#'     feature.level = c("Phylum","Genus","Family"),
#'     prev.filter = 0.1,
#'     abund.filter = 0.0001,
#'     feature.sig.level = 0.1,
#'     feature.mt.method = "none"
#' )
#'    volcano_plots <- generate_taxa_volcano_single(data.obj = peerj32.obj,
#'                                                  group.var = "group",
#'                                                  test.list = test.list,
#'                                                  feature.sig.level = 0.1,
#'                                                  feature.mt.method = "none")
#'
#' @importFrom dplyr distinct pull
#' @export
generate_taxa_volcano_single <-
  function(data.obj,
           group.var = NULL,
           test.list,
           feature.sig.level = 0.1,
           feature.mt.method = "fdr") {
    meta_tab <- data.obj$meta.dat %>%
      dplyr::select(all_of(c(group.var))) %>% rownames_to_column("sample")

    # Define the custom color palette
    color_palette <- c("#F9F871", "#F4A261", "#FF6347")

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

      sub_plot.list <-
        lapply(names(sub_test.list), function(group.level) {
          sub_test.result <- sub_test.list[[group.level]]

          # Find max absolute log2FoldChange for symmetric x-axis
          max_abs_log2FC <-
            max(abs(sub_test.result$Coefficient), na.rm = TRUE)

          p <-
            ggplot(
              sub_test.result,
              aes(
                x = Coefficient,
                y = -log10(get(p_val_var)),
                color = Prevalence,
                size = Mean.Abundance
              )
            ) +
            geom_point() +
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
              title = group.level,
              x = "Coefficient",
              y = "-log10(p-value)",
              color = "Prevalence",
              size = "Mean Abundance"
            ) +
            theme_bw() +
            theme(
              plot.title.position = "plot",
              plot.title = element_text(hjust = 0.5, size = 12),
              panel.grid.major = element_line(color = "grey", linetype = "dashed"),
              panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
              legend.position = "bottom",
              legend.text = element_text(size = 12),
              # 调整图例文本大小
              legend.title = element_text(size = 14),
              # 调整图例标题大小
              axis.text = element_text(size = 12),
              # 调整轴文本大小
              axis.title = element_text(size = 14)          # 调整轴标题大小
            ) +
            scale_color_gradientn(colors = color_palette) +
            scale_size_continuous(range = c(3, 7)) +
            coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))

          return(p)
        })

      names(sub_plot.list) <-
        names(sub_test.list)

      return(sub_plot.list)
    })


    names(plot.list) <- feature.level
    return(plot.list)
  }
