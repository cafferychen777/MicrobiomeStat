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

#' Compute and Analyze Taxa Changes Over Time
#'
#' This function from the MicrobiomeStat package computes the taxa changes over time for different groups. The changes are calculated based on a specified function. The function also performs linear modelling on these changes using a specified group variable and covariates, and returns a report summarizing these analyses.
#' If no time variable is provided, the function will select the first unique value from the metadata as the baseline for comparison.
#'
#' @param data.obj A list object in MicrobiomeStat format.
#' @param subject.var The name of the subject variable column in the metadata.
#' @param time.var The name of the time variable column in the metadata (optional).
#' @param group.var The name of the grouping variable column for linear modeling in the metadata.
#' @param adj.vars Names of additional variables to be used as covariates in the analysis.
#' @param change.base The baseline time point for detecting changes in taxa. If NULL, the first unique value from the time.var column will be used (optional).
#' @param change.func The function to be used for calculating changes. Options include: "lfc" (log fold change), "relative difference", and a user-defined function (default is "lfc").
#' @param feature.level The taxonomic level at which to perform the analysis.
#' @param prev.filter A numeric value indicating the minimum prevalence filter for the taxa (default is 0).
#' @param abund.filter A numeric value indicating the minimum abundance filter for the taxa (default is 0).
#' @param feature.dat.type Type of the feature data. Options include: "count", "proportion", and "other" (default is "count").
#' @param ... Additional parameters to be passed to the ZicoSeq function.
#'
#' @examples
#'
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
#'
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
          mutate(., original = rownames(.))
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
        gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        group_by_at(vars(!!sym(feature.level))) %>%
        dplyr::summarise(total_count = mean(count),
                         prevalence = sum(count > 0) / dplyr::n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count, -prevalence) %>%
        left_join(otu_tax, by = feature.level)

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax_filtered %>%
        gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
        group_by_at(vars(sample, !!sym(feature.level))) %>%
        dplyr::summarise(count = sum(count)) %>%
        spread(key = "sample", value = "count")

      # 转换计数为数值类型
      otu_tax_agg_numeric <-
        mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      # 将otu_tax_agg_numeric从宽格式转换为长格式
      otu_tax_long <- otu_tax_agg_numeric %>%
        gather(key = "sample", value = "value", -feature.level)

      # 将otu_tax_long和meta_tab按sample列连接
      merged_data <- otu_tax_long %>%
        inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      # 根据time列分组
      grouped_data <- merged_data %>%
        group_by(!!sym(time.var))

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
        inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # 计算value的差值
      if (is.function(change.func)) {
        combined_data <-
          combined_data %>% mutate(value_diff = change.func(value_time_2, value_time_1))
      } else if (change.func == "lfc") {
        half_nonzero_min_time_2 <- combined_data %>%
          filter(value_time_2 > 0) %>%
          group_by(!!sym(feature.level)) %>%
          summarize(half_nonzero_min = min(value_time_2) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- combined_data %>%
          filter(value_time_1 > 0) %>%
          group_by(!!sym(feature.level)) %>%
          summarize(half_nonzero_min = min(value_time_1) / 2,
                    .groups = "drop")

        combined_data <-
          left_join(
            combined_data,
            half_nonzero_min_time_2,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data <-
          left_join(
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
          combined_data %>% mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (change.func == "relative difference") {
        combined_data <- combined_data %>%
          mutate(value_diff = case_when(
            value_time_2 == 0 & value_time_1 == 0 ~ 0,
            TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
          ))
      } else {
        combined_data <-
          combined_data %>% mutate(value_diff = value_time_2 - value_time_1)
      }

      value_diff_matrix <- combined_data %>%
        select(feature.level, subject = subject, value_diff) %>%
        spread(key = subject, value = value_diff) %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      unique_meta_tab <- meta_tab %>%
        filter(subject %in% colnames(value_diff_matrix)) %>%
        select(all_of(c(subject.var, group.var, adj.vars))) %>%
        distinct(subject, .keep_all = TRUE) %>% as_tibble()

      # 获取没有NA的value_diff_matrix的列名
      cols_order <- colnames(na.omit(value_diff_matrix))

      unique_meta_tab <-
        unique_meta_tab %>% column_to_rownames("subject")

      # 使用match()函数生成索引并对unique_meta_tab进行排序，再设回行名为 'subject'
      sorted_unique_meta_tab <- unique_meta_tab[cols_order,]

      # Run ZicoSeq
      zico.obj <- ZicoSeq(
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
      significant_taxa <- names(which(zico.obj$p.adj.fdr <= 0.3))

      # If no significant taxa were found, stop the function
      if (length(significant_taxa) == 0) {
        stop("No significant taxa were found. Stopping the function.")
      }

      # 计算每个分组的平均丰度
      prop_prev_data <-
        value_diff_matrix %>% as.data.frame() %>% rownames_to_column(feature.level) %>%
        gather(-!!sym(feature.level),
               key = !!sym(subject.var),
               value = "count") %>%
        inner_join(meta_tab, by = c(subject.var), relationship = "many-to-many") %>%
        group_by(!!sym(group.var), !!sym(feature.level)) %>% # Add time.var to group_by
        summarise(mean_proportion = mean(count),
                  sdev_count = sd(count),
                  prevalence = sum(count > 0) / n(),
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
              # 更新后只包含 group.var 的估计值
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
