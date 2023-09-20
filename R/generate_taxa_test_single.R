#' Conduct Differential Abundance Testing Using ZicoSeq Method in MicrobiomeStat Package
#'
#' This function applies a differential abundance analysis using ZicoSeq on a data set. The function filters taxa based on prevalence and abundance, then it aggregates and applies the ZicoSeq method. Finally, it creates a report of significant taxa with relevant statistics.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param time.var Character string specifying the column name in metadata containing time variable.
#'                Used to subset data to a single timepoint if provided. Default NULL does not subset.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var Character string specifying the column name in metadata containing grouping
#'                 categories. This will be used as the predictor in differential abundance testing.
#' @param adj.vars Character vector specifying column names in metadata containing covariates.
#'                These will be used for adjustment in differential abundance testing.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param feature.sig.level A numeric threshold, usually between 0 and 1, for assessing the significance of individual taxa. Default is 0.1.
#' @param feature.mt.method A character string specifying the method employed for multiple testing correction (e.g., "fdr" for False Discovery Rate). Default is "fdr".
#' @param ... Additional arguments to be passed to the ZicoSeq function.
#'
#' @return A list of tibble(s) containing information about significant taxa, including R.Squared, F.Statistic, Estimate, P.Value, Adjusted.P.Value, Mean.Proportion, Mean.Prevalence, SD.Abundance and SD.Prevalence.
#'
#' @details
#' This function performs differential abundance analysis using ZicoSeq method:
#'
#' 1. Subset data to specific time point if time variable and level provided.
#'
#' 2. Extract OTU table, taxonomy table, and sample metadata.
#'
#' 3. Filter OTUs by prevalence and abundance thresholds if counting data.
#'
#' 4. Aggregate OTU table to taxonomic levels specified by \code{feature.level}.
#'
#' 5. Run ZicoSeq on aggregated table, with grouping and adjustment variables.
#'
#' 6. Extract and compile statistics for significant taxa into results tables,
#' including R-squared, F statistic, coefficients, p-values, mean proportions, etc.
#'
#' 7. Return list of tables, with each element corresponding to one taxonomic level.
#'
#' In summary, it applies preprocessing steps, runs ZicoSeq differential abundance testing,
#' and compiles informative results tables for significant taxa. Adjustment for covariates is
#' supported. Customizable taxonomic aggregation allows testing at different levels.
#'
#' @examples
#' \dontrun{
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
#' }
#' @export
generate_taxa_test_single <- function(data.obj,
                                      time.var = NULL,
                                      t.level = NULL,
                                      group.var,
                                      adj.vars,
                                      prev.filter = 0,
                                      abund.filter = 0,
                                      feature.level,
                                      feature.dat.type = c("count", "proportion", "other"),
                                      feature.sig.level = 0.1,
                                      feature.mt.method = "fdr",
                                      ...) {
  # Extract data
  mStat_validate_data(data.obj)

  if (!is.null(time.var)) {
    if (!is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
  }

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% select(all_of(c(time.var, group.var, adj.vars)))

  # 初始化formula为group.var
  formula <- group.var

  # 如果adj.vars不为空，则将其添加到formula中
  if (!is.null(adj.vars)) {
    adj.vars_string <- paste(adj.vars, collapse = " + ")
    formula <- paste(formula, "+", adj.vars_string)
  }

  if (feature.dat.type == "other") {
    prev.filter <- 0
    abund.filter <- 0
  }

  if (feature.dat.type == "count") {
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
    )
    data.obj <-
      mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
  }

  test.list <- lapply(feature.level, function(feature.level) {
    if (is.null(data.obj$feature.agg.list[[feature.level]]) &
        feature.level != "original") {
      data.obj <-
        mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    if (feature.level != "original") {
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- load_data_obj_count(data.obj)
    }

    otu_tax_agg_filter <-  otu_tax_agg %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter)

    linda.obj <- linda(
      feature.dat = otu_tax_agg_filter,
      meta.dat = meta_tab,
      formula = paste("~", formula),
      feature.dat.type = "proportion",
      prev.filter = prev.filter,
      mean.abund.filter = abund.filter,
      ...
    )

    if (!is.null(group.var)) {
      reference_level <- levels(as.factor(meta_tab[, group.var]))[1]
    }

    # 计算每个分组的平均丰度
    prop_prev_data <-
      otu_tax_agg %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::group_by(Var1) %>%  # Var1是taxa
      dplyr::summarise(avg_abundance = mean(Freq),
                       prevalence = sum(Freq > 0) / dplyr::n()) %>% column_to_rownames("Var1") %>%
      rownames_to_column(feature.level)

    extract_data_frames <-
      function(linda_object, group_var = NULL) {
        # 初始化一个空的list来存储提取的数据框
        result_list <- list()

        # 获取所有匹配的数据框名
        matching_dfs <-
          grep(paste0(group_var), names(linda_object$output), value = TRUE)

        # 循环遍历所有匹配的数据框名并提取它们
        for (df_name in matching_dfs) {
          # 从数据框名中提取组值
          group_prefix <- paste0(group_var)

          # 提取group_prefix后面的内容，并在":"之前停止
          group_value <- unlist(strsplit(df_name, split = ":"))[1]
          group_value <-
            gsub(pattern = group_prefix,
                 replacement = "",
                 x = group_value)

          # 将数据框添加到结果列表中
          result_list[[paste0(group_value, " vs ", reference_level, " (Reference)")]] <-
            linda_object$output[[df_name]]
        }

        return(result_list)
      }

    # 使用函数提取数据框
    sub_test.list <-
      extract_data_frames(linda_object = linda.obj, group_var = group.var)

    sub_test.list <- lapply(sub_test.list, function(df) {
      df <- df %>%
        rownames_to_column(feature.level) %>%
        dplyr::left_join(prop_prev_data, by = feature.level) %>%
        dplyr::select(all_of(all_of(
          c(
            feature.level,
            "log2FoldChange",
            "lfcSE",
            "pvalue",
            "padj",
            "avg_abundance",
            "prevalence"
          )
        ))) %>%
        dplyr::rename(
          Variable = feature.level,
          Coefficient = log2FoldChange,
          SE = lfcSE,
          P.Value = pvalue,
          Adjusted.P.Value = padj,
          Mean.Abundance = avg_abundance,
          Prevalence = prevalence
        )

      return(df)
    })

    return(sub_test.list)
  })

  # Assign names to the elements of test.list
  names(test.list) <- feature.level

  # plot.list <-
  #   generate_taxa_volcano_single(
  #     data.obj = data.obj,
  #     group.var = group.var,
  #     test.list = test.list,
  #     feature.sig.level = feature.sig.level,
  #     feature.mt.method = feature.mt.method
  #   )
  #
  # print(plot.list)

  # Return the results table
  return(test.list)
}

#' Generate volcano plots for taxa differential test for a single time point
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var The grouping variable tested, found in metadata
#' @param time.var The time variable used in the analysis
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
#'                                                  feature.sig.level = 0.05,
#'                                                  feature.mt.method = "fdr")
#'
#' @importFrom dplyr distinct pull
#' @export
generate_taxa_volcano_single <-
  function(data.obj,
           group.var = NULL,
           test.list,
           feature.sig.level = 0.1,
           feature.mt.method = "fdr") {
    meta_tab <- load_data_obj_metadata(data.obj) %>%
      dplyr::select(all_of(c(group.var))) %>% rownames_to_column("sample")

    # Define the custom color palette
    color_palette <- c("#2A9D8F", "#F9F871", "#F4A261", "#FF6347")

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
