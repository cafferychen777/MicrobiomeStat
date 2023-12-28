#' Longitudinal Differential Abundance Test in Microbiome Data
#'
#' This function performs differential abundance testing across multiple time points in longitudinal microbiome data. It is tailored to analyze how the abundance of microbial taxa varies over time, within different groups or under various conditions.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param subject.var A string specifying the column name in meta.dat that uniquely identifies each subject.
#' @param time.var Optional; a string representing the time variable in the meta.dat. If provided, enables longitudinal analysis.
#' @param group.var Optional; a string specifying the group variable in meta.dat for between-group comparisons.
#' @param adj.vars Optional; a vector of strings representing covariates in meta.dat for adjustment in the analysis.
#' @param feature.level A string or vector of strings indicating the taxonomic level(s) for analysis (e.g., "Phylum", "Class").
#' @param prev.filter Numeric; a minimum prevalence threshold for taxa inclusion in the analysis.
#' @param abund.filter Numeric; a minimum abundance threshold for taxa inclusion in the analysis.
#' @param feature.dat.type Character; "count" or "proportion", indicating the type of feature data.
#' @param feature.sig.level Numeric; significance level for feature testing.
#' @param feature.mt.method Character; method for multiple testing correction, such as "fdr".
#' @param ... Additional arguments passed to other methods.
#' @details
#' The function integrates various data manipulations, normalization procedures, and statistical tests to assess the significance of taxa changes over time or between groups. It allows for the adjustment of covariates and handles both count and proportion data types.
#'
#' The function constructs a mixed-effects model formula based on the provided variables, handling fixed and random effects to account for repeated measures in subjects. It performs filtering based on prevalence and abundance thresholds and applies normalization and aggregation procedures as necessary.
#'
#' Importantly, the function conducts differential abundance analysis separately for each time point in the longitudinal data. This approach allows for the identification of taxa that show significant changes at specific time points, providing insights into the dynamics of the microbiome over time.
#'
#' @return
#' A nested list structure. The top level of the list corresponds to different time points, and each element contains a list of dataframes for each taxonomic level. Each dataframe provides statistical analysis results for taxa at that level and time point.
#'
#' @examples
#' \dontrun{
#' # Example 1: Analyzing the ECAM dataset
#' data("ecam.obj")
#' # Analyzing the impact of delivery method on microbial composition over months
#' result1 <- generate_taxa_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   adj.vars = "diet",
#'   feature.level = c("Phylum", "Class"),
#'   feature.dat.type = "proportion"
#' )
#' # Visualizing the results for the ECAM dataset
#' dotplot_ecam <- generate_taxa_dotplot_long(
#'   data.obj = ecam.obj,
#'   time.test.list = result1,
#'   group.var = "delivery",
#'   time.var = "month_num",
#'   feature.level = c("Phylum", "Class")
#' )
#'
#' # Example 2: Analyzing the Type 2 Diabetes dataset
#' data("subset_T2D.obj")
#' # Longitudinal analysis of microbial changes in different racial groups
#' result2 <- generate_taxa_test_long(
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
#' dotplot_T2D <- generate_taxa_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   time.test.list = result2,
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   feature.level = c("Genus", "Family")
#' )
#' }
#' @export
generate_taxa_test_long <-
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
                                time.var = NULL,
                                subject.var = subject.var)

    time.levels <- data.obj$meta.dat %>% select(all_of(c(time.var))) %>% unique() %>% pull()

    time.test.list <- lapply(time.levels, function(t.level){

      condition <- paste(time.var, "== '", t.level, "'", sep = "")

      subset_data.obj <- mStat_subset_data(data.obj, condition = condition)

      meta_tab <-
        subset_data.obj$meta.dat %>% select(all_of(c(
          time.var, group.var, adj.vars, subject.var
        )))

      test.list <- lapply(feature.level, function(feature.level) {

        if (feature.dat.type == "count"){
          message(
            "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
          )
          subset_data.obj <- mStat_normalize_data(subset_data.obj, method = "TSS")$data.obj.norm
        }

        if (is.null(subset_data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
          subset_data.obj <- mStat_aggregate_by_taxonomy(data.obj = subset_data.obj, feature.level = feature.level)
        }

        if (feature.level != "original"){
          otu_tax_agg <- subset_data.obj$feature.agg.list[[feature.level]]
        } else {
          otu_tax_agg <- subset_data.obj$feature.tab
        }

        otu_tax_agg_filter <-  otu_tax_agg %>%
          as.data.frame() %>%
          mStat_filter(prev.filter = prev.filter,
                       abund.filter = abund.filter)

        if (feature.dat.type == "count"){
          feature.dat.type = "proportion"
        }

        linda.obj <- linda(
          feature.dat = otu_tax_agg_filter,
          meta.dat = meta_tab,
          formula = paste("~", formula),
          feature.dat.type = feature.dat.type,
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

        extract_data_frames <- function(linda_object, group_var = NULL) {

          # 初始化一个空的list来存储提取的数据框
          result_list <- list()

          # 获取所有匹配的数据框名
          matching_dfs <- grep(group_var, names(linda_object$output), value = TRUE)

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

      return(test.list)
    })

    names(time.test.list) <- time.levels

    return(time.test.list)
  }
