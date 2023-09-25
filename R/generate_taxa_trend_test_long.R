#' Longitudinal Taxa Trend Test Generation
#'
#' This function is designed to conduct a longitudinal trend test on microbiome data. The primary aim is to discern how the abundance of various microbial taxa changes over time and/or in response to different experimental or observational groups. The function delivers robust statistical insights that enable researchers to draw meaningful conclusions about the dynamics of microbial populations.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string that indicates the column name in the metadata which uniquely identifies each subject or sample.
#' @param time.var A character string representing the time variable column in the metadata. Time points should be numeric. If not, the function will convert it to numeric. Default is NULL.
#' @param group.var A character string specifying the grouping variable column in the metadata. This variable differentiates between different experimental or observational groups.
#' @param adj.vars A vector of character strings. Each string should denote a column name in the metadata that will serve as a covariate in the analysis. These variables might account for potential confounding influences. Default is NULL.
#' @param feature.level A character string indicating the taxonomic resolution for analysis (e.g., "Phylum", "Class"). This choice will determine the granularity of the analysis.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.dat.type A character string, either "count" or "proportion", indicating the nature of the data in the `data.obj`. This helps the function to determine if normalization is required. Default is "count".
#' @param feature.sig.level A numeric threshold, usually between 0 and 1, for assessing the significance of individual taxa. Default is 0.1.
#' @param feature.mt.method A character string specifying the method employed for multiple testing correction (e.g., "fdr" for False Discovery Rate). Default is "fdr".
#' @param ... Additional arguments to cater to any specialized requirements. For now, these are placeholder and not used.
#' @details
#' Based on whether group.var, adj.vars, and time.var are NULL, the formula tests:
#'
#' - When time.var is NULL:
#'   - When group.var is NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars main effects only.
#'     - Adjusted for adj.vars but not for group.var or time.var.
#'   - When group.var is NOT NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars and group.var main effects.
#'     - Adjusted for adj.vars but not for time.var.
#'   - When group.var is NOT NULL and adj.vars is NULL:
#'     - Tests group.var main effects only.
#'     - Unadjusted for adj.vars but adjusted for group.var.
#'   - When both group.var and adj.vars are NULL:
#'     - Tests the intercept only.
#'     - Unadjusted analysis.
#'
#' - When time.var is NOT NULL:
#'   - When group.var is NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars and time.var main effects.
#'     - Adjusted for adj.vars but not for group.var.
#'   - When group.var is NOT NULL and adj.vars is NOT NULL:
#'     - Tests adj.vars main effects.
#'     - Tests group.var and time.var main effects.
#'     - Tests group.var x time.var interaction.
#'     - Adjusted for adj.vars.
#'   - When group.var is NOT NULL and adj.vars is NULL:
#'     - Tests group.var and time.var main effects.
#'     - Tests group.var x time.var interaction.
#'     - Unadjusted analysis.
#'   - When both group.var and adj.vars are NULL:
#'     - Tests time.var main effect only.
#'     - Unadjusted analysis.
#'
#' The formula combines the appropriate terms based on which variables are NULL.
#' Subject variability is accounted for through random effects.
#'
#' When group.var = NULL and adj.vars = NULL and time.var is NOT NULL,
#' the slope of time.var is tested without adjusting for any additional covariates.
#'
#' @return
#' A list of dataframes, with each dataframe representing a specific taxonomic level (as specified in `feature.level`). These dataframes contain essential statistics, including taxa changes, p-values, and other metrics derived from the linear model.
#'
#' @examples
#' \dontrun{
#' # Example 1
#' data("ecam.obj")
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   feature.level = c("Phylum","Class"),
#'   feature.dat.type = c("count")
#' )
#'
#' # Example 2
#' data("subset_T2D.obj")
#' generate_taxa_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus","Family"),
#'   feature.dat.type = c("count")
#' )
#' }
#' @export
generate_taxa_trend_test_long <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var = NULL,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion"),
           feature.sig.level = 0.1,
           feature.mt.method = "fdr",
           ...) {
    # Extract data
    mStat_validate_data(data.obj)

    message(
      "The trend test calculation relies on a numeric time variable.\n",
      "Please check that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the test results.\n",
      "You can ensure the time variable is numeric by mutating it in the metadata."
    )

    if (!is.null(time.var)){
      data.obj$meta.dat <-
        data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))
    }

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

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

      extract_data_frames <- function(linda_object, group_var = NULL, time_var) {

        # 初始化一个空的list来存储提取的数据框
        result_list <- list()

        # 如果group.var不为NULL
        if (!is.null(group_var)) {
          # 获取所有匹配的数据框名
          matching_dfs <- grep(paste0(group_var, ".*:", time_var), names(linda_object$output), value = TRUE)

          # 循环遍历所有匹配的数据框名并提取它们
          for (df_name in matching_dfs) {
            # 从数据框名中提取组值
            group_prefix <- paste0(group_var)

            # 提取group_prefix后面的内容，并在":"之前停止
            group_value <- unlist(strsplit(df_name, split = ":"))[1]
            group_value <- gsub(pattern = group_prefix, replacement = "", x = group_value)

            # 将数据框添加到结果列表中
            result_list[[paste0(group_value," vs ", reference_level, " (Reference)")]] <- linda_object$output[[df_name]]
          }

        } else {
          # 如果group.var为NULL，则直接提取名字为time.var的数据框
          df <- linda_object$output[[time_var]]
          if (!is.null(df)) {
            result_list[[time_var]] <- df
          } else {
            warning(paste("No data frame found with the name:", time_var))
          }
        }

        return(result_list)
      }

      # 使用函数提取数据框
      sub_test.list <- extract_data_frames(linda_object = linda.obj, group_var = group.var, time_var = time.var)

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
    #     time.var = time.var,
    #     test.list = test.list,
    #     feature.sig.level = feature.sig.level,
    #     feature.mt.method = feature.mt.method
    #   )
    #
    # print(plot.list)

    return(test.list)
  }

#' Generate volcano plots for longitudinal taxa trend test
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
#' # data("subset_T2D.obj")
#' # test_list <- generate_taxa_trend_test_long(data.obj, ...)
#' # volcano_plots <- generate_taxa_trend_volcano_long(data.obj,
#'                                                  # group.var,
#'                                                  # time.var,
#'                                                  # test_list,
#'                                                  # feature.sig.level = 0.05,
#'                                                  # feature.mt.method = "fdr")
#'
#' @importFrom dplyr distinct pull
#' @export
generate_taxa_trend_volcano_long <-
  function(data.obj,
           group.var = NULL,
           time.var = NULL,
           test.list,
           feature.sig.level = 0.1,
           feature.mt.method = "fdr") {
    meta_tab <- load_data_obj_metadata(data.obj) %>%
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

      if (!is.null(group.var)) {
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
                legend.text = element_text(size = 12),       # 调整图例文本大小
                legend.title = element_text(size = 14),      # 调整图例标题大小
                axis.text = element_text(size = 12),          # 调整轴文本大小
                axis.title = element_text(size = 14)          # 调整轴标题大小
              ) +
              scale_color_gradientn(colors = color_palette) +
              scale_size_continuous(range = c(3, 7)) +
              coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))

            return(p)
          })

        names(sub_plot.list) <-
          names(sub_test.list)
      } else {
        # 当group.var为NULL时
        sub_test.result <- sub_test.list[[time.var]]

        # Find max absolute log2FoldChange for symmetric x-axis
        max_abs_log2FC <-
          max(abs(sub_test.result$Coefficient), na.rm = TRUE)

        sub_plot.list <- lapply("time", function(time) {
          p <-
            ggplot(sub_test.result, aes(x = Coefficient, y = -log10(get(p_val_var)),
                                        color = Prevalence, size = Mean.Abundance)) +
            geom_point() +
            geom_vline(aes(xintercept = 0), linetype = "dashed", linewidth = 1.5, color = "grey") +
            geom_hline(aes(yintercept = -log10(feature.sig.level)), linetype = "dashed", linewidth = 1.5, color = "grey") +
            geom_text(aes(label = ifelse(get(p_val_var) < feature.sig.level, as.character(Variable), '')),
                      vjust = -0.5, hjust = 0.5, size = 3.5) +
            scale_shape_manual(values = c(16, 17)) +
            labs(x = "Coefficient", y = "-log10(p-value)", color = "Mean Prevalence", size = "Mean Abundance") +
            theme_bw() +
            theme(
              plot.title.position = "plot",
              plot.title = element_text(hjust = 0.5, size = 12),  # 这里的size = 14只是示例，您可以根据需要调整此值
              panel.grid.major = element_line(color = "grey", linetype = "dashed"),
              panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
              legend.position = "bottom",
              legend.text = element_text(size = 12),       # 调整图例文本大小
              legend.title = element_text(size = 14),      # 调整图例标题大小
              axis.text = element_text(size = 12),          # 调整轴文本大小
              axis.title = element_text(size = 14)          # 调整轴标题大小
            ) +
            scale_color_gradientn(colors = color_palette) +
            scale_size_continuous(range = c(3, 7)) +  # 可以根据您的实际需求调整大小范围
            coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))

          return(p)
        })
        names(sub_plot.list) <-
          time.var
      }
      return(sub_plot.list)
    })


    names(plot.list) <- feature.level
    return(plot.list)
  }
