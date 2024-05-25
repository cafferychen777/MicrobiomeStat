#' Conduct Differential Abundance Testing Using LinDA Method in MicrobiomeStat Package
#'
#' This function applies a differential abundance analysis using LinDA on a data set. The function filters taxa based on prevalence and abundance, then it aggregates and applies the LinDA method. Finally, it creates a report of significant taxa with relevant statistics.
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
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
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
#' @param ... Additional arguments to be passed to the ZicoSeq function.
#'
#' @return A list of tibble(s) containing information about significant taxa, including R.Squared, F.Statistic, Estimate, P.Value, Adjusted.P.Value, Mean.Proportion, Mean.Prevalence, SD.Abundance and SD.Prevalence.
#'
#' @details
#' This function facilitates differential abundance analysis utilizing the LinDA method:
#'
#' 1. If a specific time variable and level are provided, the data is subsetted accordingly.
#'
#' 2. Extracts OTU table and sample metadata.
#'
#' 3. If the feature data type is of "count", it normalizes the data using the "TSS" transformation.
#'
#' 4. If the feature level is not "original", it aggregates the OTU table to the taxonomic levels specified by \code{feature.level}.
#'
#' 5. Executes the LinDA method on the aggregated or original table considering the grouping and adjustment variables.
#'
#' 6. Extracts significant taxa's statistics into results tables, which include coefficients, standard errors, p-values, adjusted p-values, average abundances, and prevalence.
#'
#' 7. Returns a list of result tables where each element corresponds to a particular taxonomic level.
#'
#' In essence, the function streamlines preprocessing, executes LinDA-based differential abundance testing, and assembles tables with pertinent results for significant taxa. It also supports adjusting for covariates and allows taxonomic aggregation at diverse levels for customized analyses.
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
#' )
#' plot.list <- generate_taxa_volcano_single(
#'     data.obj = peerj32.obj, # or other data.obj you want to use
#'     group.var = "group", # or other group.var you want to use
#'     test.list = test.list,
#'     feature.sig.level = 0.1, # or other sig.level you want to use
#'     feature.mt.method = "none" # or other mt.method you want to use
#' )
#' print(plot.list)
#' }
#' @export
generate_taxa_test_single <- function(data.obj,
                                      time.var = NULL,
                                      t.level = NULL,
                                      group.var,
                                      adj.vars = NULL,
                                      prev.filter = 0,
                                      abund.filter = 0,
                                      feature.level,
                                      feature.dat.type = c("count", "proportion", "other"),
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
    data.obj$meta.dat %>% select(all_of(c(time.var, group.var, adj.vars)))

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
      otu_tax_agg <- data.obj$feature.tab
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

  # Return the results table
  return(test.list)
}
