#' Determine if a variable is categorical
#'
#' This function checks if a given variable is categorical. A variable is considered categorical if it is a factor, a character, or a logical type. Additionally, numeric variables with a small number of unique values (less than a specified threshold) are also considered categorical.
#'
#' @param x The variable to be tested for being categorical.
#'
#' @return A logical value indicating whether the input variable is categorical (TRUE) or not (FALSE).
#'
#' @examples
#' is_categorical(factor(c("A", "B", "C")))
#' is_categorical(c("A", "B", "C"))
#' is_categorical(c(TRUE, FALSE))
#' is_categorical(c(1, 2, 3, 4))
#' is_categorical(c(1.2, 3.5, 4.6, 7.8))
#'
#' @export
is_categorical <- function(x) {
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    return(TRUE)
  } else if (is.numeric(x) &&
             length(unique(x)) < 10) {
    # 这里的10你可以根据需要调整
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Compute taxa changes between time points and analyze differential abundance between groups
#'
#' This function calculates taxa abundance changes between two time points in the metadata, using the time values specified in `time.var` and baseline `change.base`.
#' It computes changes based on the method in `change.func`.
#' The function then uses the ZicoSeq method to perform differential abundance analysis between the groups in `group.var`, adjusted for `adj.vars`.
#' It returns data frames summarizing the results of the differential abundance tests for each taxonomic level in `feature.level`.
#' If `time.var` is not provided, the first unique value in the metadata will be used as `change.base`.
#' If only one time value exists, the function will exit with a message.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var The name of the subject variable column in the metadata.
#' @param time.var The name of the time variable column in the metadata (optional).
#' @param group.var The name of the grouping variable column for linear modeling in the metadata.
#' @param adj.vars Names of additional variables to be used as covariates in the analysis.
#' @param change.base The baseline time point for detecting changes in taxa. If NULL, the first unique value from the time.var column will be used (optional).
#' @param change.func Function or character string specifying method to compute change between time points. Options are:
#' - 'difference' (default): Computes absolute difference in counts between time points.
#' - 'relative difference': Computes relative difference in counts between time points, calculated as (count_time2 - count_time1) / (count_time2 + count_time1).
#' - 'lfc': Computes log2 fold change between time points. Zero counts are imputed with half the minimum nonzero value before log transform.
#' - Custom function: A user-defined function can also be provided, which should take two vectors of counts (at time 1 and time 2) as input and return the computed change.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param ... Additional parameters to be passed to the ZicoSeq function.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(GUniFrac)
#' library(ape)
#' library(philentropy)
#' library(MicrobiomeStat)
#' data(peerj32.obj)
#'
#' generate_taxa_change_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   change.base = "1",
#'   change.func = "lfc",
#'   feature.level = "original",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   feature.dat.type = "count"
#' )
#' }
#'
#' @return A named list where each element corresponds to a feature level and contains a dataframe with the calculated taxa changes, their corresponding p-values, and other statistics from the linear model.
#' @export
#' @name generate_taxa_change_test_pair
generate_taxa_change_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars,
           change.base,
           change.func,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           ...) {
    # Extract data
    mStat_validate_data(data.obj)

    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      otu_tab <-
        load_data_obj_count(mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm)
    } else{
      otu_tab <- load_data_obj_count(data.obj)
    }

    tax_tab <- load_data_obj_taxonomy(data.obj) %>%
      as.data.frame() %>%
      {
        if ("original" %in% feature.level)
          dplyr::mutate(., original = rownames(.))
        else
          .
      } %>%
      select(all_of(feature.level))

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    if (is.null(change.base)) {
      change.base <- unique(meta_tab %>% select(all_of(c(time.var))))[1,]
      message(
        "The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'meta.dat' data frame: ",
        change.base
      )
    }

    # 将 OTU 表与分类表合并
    otu_tax <-
      cbind(otu_tab, tax_tab)

    if (feature.dat.type == "other") {
      prev.filter <- 0
      abund.filter <- 0
    }

    test.list <- lapply(feature.level, function(feature.level) {
      # Filter taxa based on prevalence and abundance
      otu_tax_filtered <- otu_tax %>%
        tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        dplyr::group_by_at(vars(!!sym(feature.level))) %>%
        dplyr::summarise(total_count = mean(count),
                         prevalence = sum(count > 0) / dplyr::n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count, -prevalence) %>%
        dplyr::left_join(otu_tax, by = feature.level)

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax_filtered %>%
        tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        dplyr::group_by_at(vars(sample, !!sym(feature.level))) %>%
        dplyr::summarise(count = sum(count)) %>%
        tidyr::spread(key = "sample", value = "count")

      # 转换计数为数值类型
      otu_tax_agg_numeric <-
        dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      # 将otu_tax_agg_numeric从宽格式转换为长格式
      otu_tax_long <- otu_tax_agg_numeric %>%
        tidyr::gather(key = "sample", value = "value", -feature.level)

      # 将otu_tax_long和meta_tab按sample列连接
      merged_data <- otu_tax_long %>%
        dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      # 根据time列分组
      grouped_data <- merged_data %>%
        dplyr::group_by(!!sym(time.var))

      change.after <-
        unique(grouped_data %>% select(all_of(c(time.var))))[unique(grouped_data %>% select(all_of(c(time.var)))) != change.base]

      # 拆分成一个列表，每个time值都有一个独立的tibble
      split_data <-
        split(merged_data, f = grouped_data %>% select(all_of(c(time.var))))

      # 提取split_data中的第一个和第二个表
      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # 将这两个表连接在一起，以便计算差值
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # 计算value的差值
      if (is.function(change.func)) {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = change.func(value_time_2, value_time_1))
      } else if (change.func == "lfc") {
        half_nonzero_min_time_2 <- combined_data %>%
          filter(value_time_2 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(value_time_2) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- combined_data %>%
          filter(value_time_1 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(value_time_1) / 2,
                    .groups = "drop")

        combined_data <-
          dplyr::left_join(
            combined_data,
            half_nonzero_min_time_2,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data <-
          dplyr::left_join(
            combined_data,
            half_nonzero_min_time_1,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data$value_time_2[combined_data$value_time_2 == 0] <-
          combined_data$half_nonzero_min_time_2[combined_data$value_time_2 == 0]
        combined_data$value_time_1[combined_data$value_time_1 == 0] <-
          combined_data$half_nonzero_min_time_1[combined_data$value_time_1 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message(
          "Imputation was performed using half the minimum nonzero proportion for each taxon at different time points."
        )

        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (change.func == "relative difference") {
        combined_data <- combined_data %>%
          dplyr::mutate(value_diff = dplyr::case_when(
            value_time_2 == 0 & value_time_1 == 0 ~ 0,
            TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
          ))
      } else {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      }

      value_diff_matrix <- combined_data %>%
        select(feature.level, subject = subject, value_diff) %>%
        tidyr::spread(key = subject, value = value_diff) %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      unique_meta_tab <- meta_tab %>%
        filter(subject %in% colnames(value_diff_matrix)) %>%
        select(all_of(c(subject.var, group.var, adj.vars))) %>%
        dplyr::distinct(subject, .keep_all = TRUE) %>% as_tibble()

      # 获取没有NA的value_diff_matrix的列名
      cols_order <- colnames(na.omit(value_diff_matrix))

      unique_meta_tab <-
        unique_meta_tab %>% column_to_rownames("subject")

      # 使用match()函数生成索引并对unique_meta_tab进行排序，再设回行名为 'subject'
      sorted_unique_meta_tab <- unique_meta_tab[cols_order,]

      # Run ZicoSeq
      zico.obj <- GUniFrac::ZicoSeq(
        meta.dat = sorted_unique_meta_tab,
        feature.dat = na.omit(value_diff_matrix),
        grp.name = group.var,
        adj.name = adj.vars,
        feature.dat.type = "other",
        prev.filter = prev.filter,
        max.abund.filter = abund.filter,
        ...
      )

      # Extract relevant information
      significant_taxa <- names(which(zico.obj$p.adj.fdr <= 1))

      # 计算每个分组的平均丰度
      prop_prev_data <-
        value_diff_matrix %>% as.data.frame() %>% rownames_to_column(feature.level) %>%
        tidyr::gather(-!!sym(feature.level),
               key = !!sym(subject.var),
               value = "count") %>%
        dplyr::inner_join(meta_tab, by = c(subject.var), relationship = "many-to-many") %>%
        dplyr::group_by(!!sym(group.var), !!sym(feature.level)) %>% # Add time.var to dplyr::group_by
        dplyr::summarise(mean_proportion = mean(count),
                  sdev_count = sd(count),
                  prevalence = sum(count > 0) / dplyr::n(),
                  sdev_prevalence = sd(ifelse(count > 0, 1, 0)))

      # Initialize results table
      results <- data.frame()

      for (taxa in significant_taxa) {
        R.Squared <- zico.obj$R2[taxa, 1]
        F.Statistic <- zico.obj$F0[taxa, 1]

        Estimate <-
          zico.obj$coef.list[[1]][startsWith(rownames(zico.obj$coef.list[[1]]), group.var), taxa]  # 选择 group.var 的估计值

        P.Value <- zico.obj$p.raw[taxa]
        Adjusted.P.Value <- zico.obj$p.adj.fdr[taxa]

        if (is_categorical(data.obj$meta.dat[[group.var]])) {
          for (group in unique(data.obj$meta.dat[[group.var]])) {
            group_data <-
              prop_prev_data[which(prop_prev_data[[feature.level]] == taxa &
                                     prop_prev_data[[group.var]] == group), ]
            mean_prop <- group_data$mean_proportion
            mean_prev <- group_data$prevalence
            sd_abundance <- group_data$sdev_count
            sd_prevalence <- group_data$sdev_prevalence

            results <- rbind(
              results,
              data.frame(
                Variable = taxa,
                Group = group,
                R.Squared = R.Squared,
                F.Statistic = F.Statistic,
                Estimate = toString(Estimate),
                P.Value = P.Value,
                Adjusted.P.Value = Adjusted.P.Value,
                Mean.Abundance_Change = mean_prop,
                Mean.Prevalence_Change = mean_prev,
                SD.Abundance_Change = sd_abundance,
                SD.Prevalence_Change = sd_prevalence
              )
            )
          }
        } else {
          total_data <-
            prop_prev_data[which(prop_prev_data[[feature.level]] == taxa), ]
          mean_prop <- total_data$mean_proportion
          mean_prev <- total_data$prevalence
          sd_abundance <- total_data$sdev_count
          sd_prevalence <- total_data$sdev_prevalence

          results <- rbind(
            results,
            data.frame(
              Variable = taxa,
              R.Squared = R.Squared,
              F.Statistic = F.Statistic,
              Estimate = toString(Estimate),
              P.Value = P.Value,
              Adjusted.P.Value = Adjusted.P.Value,
              Mean.Abundance_Change = mean_prop,
              Mean.Prevalence_Change = mean_prev,
              SD.Abundance_Change = sd_abundance,
              SD.Prevalence_Change = sd_prevalence
            )
          )
        }
      }

      return(results)

    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)

  }
