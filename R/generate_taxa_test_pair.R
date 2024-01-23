#' Generate Taxa Test Pair
#'
#' This function takes a MicrobiomeStat data object as input, filters taxa based on
#' prevalence and abundance thresholds, aggregates taxon abundances by sample,
#' fits linear mixed effects models using the linda method to identify significant
#' taxon changes across groups over time, accounting for specified covariates.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A string that specifies the name of the subject variable column in the metadata.
#' @param time.var A string that specifies the name of the time variable column in the metadata. If not provided, it's NULL by default.
#' @param group.var A string that specifies the name of the grouping variable column in the metadata for linear modelling.
#' @param adj.vars A vector of strings that specify the names of additional variables to be used as covariates in the analysis.
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
#' @param ... Additional parameters to be passed to the linda function.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' test.list <- generate_taxa_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.dat.type = "count"
#' )
#' plot.list <-
#' generate_taxa_volcano_single(
#'  data.obj = peerj32.obj,
#'  group.var = "group",
#'  test.list = test.list,
#'  feature.sig.level = 0.1,
#'  feature.mt.method = "none"
#')
#'
#' data("subset_pairs.obj")
#' test.list <- generate_taxa_test_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.dat.type = "count"
#' )
#' plot.list <-
#' generate_taxa_volcano_single(
#'  data.obj = subset_pairs.obj,
#'  group.var = "Sex",
#'  test.list = test.list,
#'  feature.sig.level = 0.1,
#'  feature.mt.method = "none"
#')
#' }
#'
#' @return A named list containing data frames summarizing taxon test results for each taxonomic level.
#' @details Each list element corresponds to a taxonomic level specified in `feature.level`.
#' The data frame contains columns for taxon name, log2 fold change, p-values, adjusted p-values,
#' mean abundance, mean prevalence, and the output element from `linda` where the taxon was found significant.
#' @export
#' @name generate_taxa_test_pair
generate_taxa_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           feature.sig.level = 0.1,
           feature.mt.method = "fdr",
           ...) {
    # Extract data
    mStat_validate_data(data.obj)

    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    # Create the formula
    # Function to generate the formula
    generate_formula <- function(group.var = NULL, adj.vars = NULL, time.var = NULL, subject.var = NULL) {

      # Initialize the fixed_effects and random_effects variables
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine multiple adj.vars into a single string if they exist
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Case where time.var is NULL
      if (is.null(time.var)) {
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var)
          } else {
            fixed_effects <- group.var
          }
        }
        random_effects <- paste("(1 |", subject.var, ")")

        # Case where time.var is NOT NULL
      } else {
        if (is.null(group.var)) {
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        random_effects <- paste("(1 +", time.var, "|", subject.var, ")")
      }

      # Generate the full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    formula <- generate_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = time.var,
                                subject.var = subject.var)

    correct_formula <- function(group.var = NULL, adj.vars = NULL, time.var = NULL, subject.var = NULL) {

      # Initialize the fixed_effects and random_effects variables
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine multiple adj.vars into a single string if they exist
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Case where time.var is NULL
      if (is.null(time.var)) {
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var)
          } else {
            fixed_effects <- group.var
          }
        }
        random_effects <- paste("(1 |", subject.var, ")")

        # Case where time.var is NOT NULL
      } else {
        if (is.null(group.var)) {
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        random_effects <- paste("(1 |", subject.var, ")")
      }

      # Generate the full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    formula_corrected <- correct_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = time.var,
                                subject.var = subject.var)

    test.list <- lapply(feature.level, function(feature.level) {

      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      }

      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      otu_tax_agg_filter <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter)

      if (feature.dat.type == "count"){
        feature.dat.type <- "proportion"
      }

      linda.obj <- tryCatch({
        # 尝试运行 linda 函数
        linda(feature.dat = otu_tax_agg_filter,
              meta.dat = meta_tab,
              formula = paste("~", formula),
              feature.dat.type = feature.dat.type,
              ...)
      }, error = function(e) {
        # 如果出错，打印错误消息
        message("Error in linda: ", e)

        message("Due to the above error, a simpler model will be used for fitting.")

        # 使用修正后的公式重新运行 linda 函数
        linda(feature.dat = otu_tax_agg_filter,
              meta.dat = meta_tab,
              formula = paste("~", formula_corrected),
              feature.dat.type = feature.dat.type)
      })

      if (!is.null(group.var)){
        reference_level <- levels(as.factor(meta_tab[,group.var]))[1]
      }

      # 计算每个分组的平均丰度
      prop_prev_data <-
        otu_tax_agg %>%
        as.matrix() %>%
        as.table() %>%
        as.data.frame() %>%
        dplyr::group_by(Var1) %>%  # Var1是taxa
        dplyr::summarise(
          avg_abundance = mean(Freq),
          prevalence = sum(Freq > 0) / dplyr::n()
        ) %>% column_to_rownames("Var1") %>%
        rownames_to_column(feature.level)

      extract_data_frames <- function(linda_object, group_var = NULL) {

        # 初始化一个空的list来存储提取的数据框
        result_list <- list()

        # 获取所有匹配的数据框名
        matching_dfs <- grep(paste0(group_var), names(linda_object$output), value = TRUE)

        # 循环遍历所有匹配的数据框名并提取它们
        for (df_name in matching_dfs) {
          # 从数据框名中提取组值
          group_prefix <- paste0(group_var)

          # 提取group_prefix后面的内容，并在":"之前停止
          group_value <- unlist(strsplit(df_name, split = ":"))[1]
          group_value <- gsub(pattern = group_prefix, replacement = "", x = group_value)

          if (grepl(pattern = time.var, x = df_name)){
            result_list[[paste0(group_value," vs ", reference_level, " (Reference) [", "Interaction", "]")]] <- linda_object$output[[df_name]]
          } else {
            # 将数据框添加到结果列表中
            result_list[[paste0(group_value," vs ", reference_level, " (Reference) [", "Main Effect", "]")]] <- linda_object$output[[df_name]]
          }
        }

        return(result_list)
      }

      # 使用函数提取数据框
      sub_test.list <- extract_data_frames(linda_object = linda.obj, group_var = group.var)

      sub_test.list <- lapply(sub_test.list, function(df){
        df <- df %>%
          rownames_to_column(feature.level) %>%
          dplyr::left_join(prop_prev_data, by = feature.level) %>%
          dplyr::select(all_of(all_of(c(feature.level,"log2FoldChange","lfcSE","pvalue","padj","avg_abundance","prevalence")))) %>%
          dplyr::rename(Variable = feature.level,
                        Coefficient = log2FoldChange,
                        SE = lfcSE,
                        P.Value = pvalue,
                        Adjusted.P.Value = padj,
                        Mean.Abundance = avg_abundance,
                        Prevalence = prevalence)

        return(df)
      })

      return(sub_test.list)

    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    # plot.list <-
    #   generate_taxa_trend_volcano_long(
    #     data.obj = data.obj,
    #     group.var = group.var,
    #     test.list = test.list,
    #     feature.sig.level = feature.sig.level,
    #     feature.mt.method = feature.mt.method
    #   )
    #
    # print(plot.list)

    return(test.list)

  }
