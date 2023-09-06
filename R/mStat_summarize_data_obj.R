#' Summarize a MicrobiomeStat Data Object
#'
#' This function takes a MicrobiomeStat data object and provides a comprehensive summary of the key components, including feature abundance matrix (feature.tab), sample metadata (meta.dat), and feature annotations (feature.ann). It also checks for optional components like time variable and phylogenetic tree.
#'
#' The summary aims to give an overview of the input microbiome data object before conducting statistical analysis, allowing users to better understand the basic properties of their data.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param time.var A column name in meta.dat representing the time variable. Optional.
#' @param group.var A column name in meta.dat representing the grouping variable. Optional.
#' @param palette A vector of colors for grouping in histogram. Optional. Default is a preset 10 color palette.
#'
#' @return A tibble containing detailed summaries of:
#' \enumerate{
#'   \item \strong{feature.tab}: Number of features, number of samples, matrix sparsity, singleton features.
#'
#'   \item \strong{meta.dat}: Number of samples, number of metadata fields, missing data, sample distribution over time (if time.var provided).
#'
#'   \item \strong{feature.ann}: Number of features, number of annotation fields, proportion NA for each annotation.
#'
#'   \item \strong{tree}: Whether phylogenetic tree exists.
#' }
#'
#' @details
#' This function checks if each key component of the MicrobiomeStat object exists, and provides a detailed summary if present.
#'
#' For feature.tab, it summarizes number of features, samples, sparsity, and singleton features. For meta.dat, it summarizes sample size, metadata fields, missing values, and temporal distribution (if time variable given). For feature.ann, it summarizes number of features, annotations, and missing values per annotation. It also checks for a phylogenetic tree.
#'
#' If time variable is provided, temporal distribution of samples is visualized using ggplot2 histograms. If group variable is also provided, histograms are grouped by the grouping variable and colored based on the palette.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'data.obj' is your MicrobiomeStat data object
#'   # Summary with time variable
#'   # summary_list <- mStat_summarize_data_obj(data.obj, time.var = "time")
#'
#'   # Summary without time variable
#'   # summary_list <- mStat_summarize_data_obj(data.obj)
#'
#'   # If you have a microbiome data available as a MicrobiomeStat data object
#'   # you can dplyr::summarize it using:
#'   # library(MicrobiomeStat)
#'   # data(data.obj)
#'   # Summary with time variable
#'   # summary_list <- mStat_summarize_data_obj(data.obj, time.var = "time")
#'
#'   # Summary without time variable
#'   # summary_list <- mStat_summarize_data_obj(data.obj)
#'   data(subset_T2D.obj)
#'   summary <- mStat_summarize_data_obj(subset_T2D.obj, "visit_number", "subject_race")
#' }
#'
#' @details
#' The function first checks if each component of the MicrobiomeStat data object is not null. If a component is not null, it is summarized and added to the output list. For the feature.tab, it computes the sparsity and singleton features. For the meta.dat, it computes the number of samples and metadata fields, and the distribution of samples if a time variable is provided. The inclusion of a time variable allows the user to gain insights into how samples are distributed over time. For the feature.ann, it computes the number of features, annotations, and the proportion of NA values for each annotation. It also checks if a phylogenetic tree exists in the data object.
#'
#' @export
mStat_summarize_data_obj <-
  function(data.obj,
           time.var = NULL,
           group.var = NULL,
           palette = NULL) {
    # Verify that feature.tab exists in data.obj
    if (!"feature.tab" %in% names(data.obj)) {
      stop("Data object must contain 'feature.tab'")
    }
    feature_tab <- data.obj$feature.tab

    # If meta.dat exists, compute summary statistics
    if ("meta.dat" %in% names(data.obj)) {
      meta_dat <- data.obj$meta.dat
      num_meta_vars <- ncol(meta_dat)
      meta_var_names <- colnames(meta_dat)
    } else {
      num_meta_vars <- NA
      meta_var_names <- NA
    }

    # If feature.ann exists, compute NA proportions for each column
    if ("feature.ann" %in% names(data.obj)) {
      feature_ann <- data.obj$feature.ann
      NA_props <-
        sapply(colnames(feature_ann), function(x)
          mean(is.na(feature_ann[, x])))
    } else {
      NA_props <- NA
    }

    # Check if tree exists
    tree_exists <-
      ifelse("tree" %in% names(data.obj) &&
               !is.null(data.obj$tree),
             "Yes",
             "No")

    # If feature.agg.list exists, get list of aggregated taxonomies
    if ("feature.agg.list" %in% names(data.obj)) {
      agg_taxonomies <- names(data.obj$feature.agg.list)
    } else {
      agg_taxonomies <- NA
    }

    # Compute summary statistics for feature_tab
    ave_reads_per_sample <-
      sum(colSums(feature_tab)) / nrow(feature_tab)
    min_reads <- min(colSums(feature_tab))
    max_reads <- max(colSums(feature_tab))
    total_reads <- sum(colSums(feature_tab))
    median_reads_per_sample <- median(colSums(feature_tab))
    zero_count_prop <-
      length(which(feature_tab == 0)) / length(feature_tab)
    count_single_occurrence <-
      length(which(rowSums(feature_tab) == 1))

    # Handling the time-series data
    if (!is.null(time.var) && "meta.dat" %in% names(data.obj)) {
      if (time.var %in% colnames(data.obj$meta.dat)) {
        time_var_data <- data.obj$meta.dat[[time.var]]

        # cat("Time-Series Information:\n",
        #     "========================\n")

        # Create a histogram for the time variable
        if (!is.null(group.var) &&
            group.var %in% colnames(data.obj$meta.dat)) {
          # Create a grouped data frame
          grouped_df <- data.obj$meta.dat %>%
            select(all_of(c(time.var, group.var))) %>%
            dplyr::group_by_at(vars(time.var, group.var)) %>%
            dplyr::summarise(SampleCount = dplyr::n(), .groups = "drop")

          # Print the grouped data frame
          #print(grouped_df)

          # Set the palette if it's NULL
          if (is.null(palette)) {
            palette <- c(
              "#E31A1C",
              "#1F78B4",
              "#FB9A99",
              "#33A02C",
              "#FDBF6F",
              "#B2DF8A",
              "#A6CEE3",
              "#BA7A70",
              "#9D4E3F",
              "#829BAB"
            )
          }

          print(
            ggplot(grouped_df, aes(
              x = !!sym(time.var),
              y = SampleCount,
              fill = !!sym(group.var)
            )) +
              geom_bar(stat = "identity", position = "stack") +
              scale_fill_manual(values = palette) +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)) +
              labs(
                title = "Histogram of Sample Counts over Time",
                x = time.var,
                y = "Sample Count"
              ) +
              scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))
          )
        } else {
          # Calculate the table
          time_table <- table(data.obj$meta.dat[[time.var]])

          # Convert to a data frame
          time_df <-
            as.data.frame(table(data.obj$meta.dat[[time.var]]))

          # Rename the columns
          colnames(time_df) <- c("TimePoint", "SampleCount")

          # Print the data frame
          #print(time_df)

          print(
            ggplot(time_df, aes(x = TimePoint, y = SampleCount)) +
              geom_bar(stat = "identity", fill = "steelblue") +
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5)) +
              labs(
                title = "Histogram of Sample Counts over Time",
                x = "Time Point",
                y = "Sample Count"
              ) +
              scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))
          )
        }

        # Explanation for the boxplot
        message("Since the presence of those low-abundance features depends highly on the sequencing depth, the association between the sequence depth and a variable of interest could lead to false associations of the observed abundance of these low-abundance features with the variable of interest. In this case, rarefaction may be needed to reduce the sequencing depth confounding. The sequence depth plot will help diagnose such potential confounding.\n")

        # Compute sequencing depth for each sample
        seq_depth <- colSums(feature_tab)

        # Convert it to a data.frame for plotting
        seq_depth_df <- data.frame(
          TimePoint = data.obj$meta.dat[[time.var]],
          SequencingDepth = seq_depth
        )

        # If group.var is provided, add it to the data.frame
        if (!is.null(group.var) && group.var %in% colnames(data.obj$meta.dat)) {
          seq_depth_df$Group <- data.obj$meta.dat[[group.var]]
        }

        seq_depth_df$TimePoint <- as.factor(seq_depth_df$TimePoint)

        # Plotting the sequencing depth boxplot
        seq_depth_plot <- ggplot(seq_depth_df, aes(x = TimePoint, y = SequencingDepth, fill = Group)) +
          geom_boxplot() +
          theme_minimal() +
          #geom_jitter(position = position_dodge(width = 0.75), alpha = 0.5, size = 2) +
          theme(plot.title = element_text(hjust = 0.5)) +
          labs(
            title = "Boxplot of Sequencing Depth over Time",
            x = time.var,
            y = "Sequencing Depth"
          )
        if (!is.null(group.var) && group.var %in% colnames(data.obj$meta.dat)) {
          seq_depth_plot <- seq_depth_plot + scale_fill_manual(values = palette, name = "Group")
        }
        print(seq_depth_plot)

      } else {
        #cat("The provided time variable does not exist in the metadata.\n")
      }
    }

    # 创建一个 data.frame 来存储统计数据
    table1 <- data.frame(
      Category = character(),
      Variable = character(),
      Value = character(),
      stringsAsFactors = FALSE
    )

    # 基本统计数据
    basic_stats <- data.frame(
      Category = "Basic Statistics",
      Variable = c(
        "Min. reads per sample",
        "Max. reads per sample",
        "Total reads dplyr::across all samples",
        "Average reads per sample",
        "Median reads per sample",
        "Proportion of zero counts",
        "Count of features that only appear once"
      ),
      Value = c(
        min_reads,
        max_reads,
        total_reads,
        ave_reads_per_sample,
        median_reads_per_sample,
        zero_count_prop,
        count_single_occurrence
      )
    )

    table1 <- rbind(table1, basic_stats)

    # 元数据统计
    if (!is.na(num_meta_vars)) {
      # metadata_stats <- data.frame(
      #   Category = "Metadata",
      #   Variable = c("Number of metadata variables", "Metadata variables"),
      #   Value = c(num_meta_vars, paste(meta_var_names, collapse = ", "))
      # )

      metadata_stats <- data.frame(
        Category = "Metadata",
        Variable = c("Number of metadata variables"),
        Value = c(num_meta_vars)
      )

      table1 <- rbind(table1, metadata_stats)
    }

    # feature.ann 缺失注释的统计
    if (length(NA_props) > 0) {
      for (i in 1:length(NA_props)) {
        missing_annotation <- data.frame(
          Category = "Feature Annotations",
          Variable = paste(
            "Proportion of missing annotations in",
            colnames(feature_ann)[i]
          ),
          Value = NA_props[i]
        )
        table1 <- rbind(table1, missing_annotation)
      }
    }

    # 树信息
    tree_info <- data.frame(Category = "Phylogenetic Tree",
                            Variable = "Exists in the dataset",
                            Value = tree_exists)

    table1 <- rbind(table1, tree_info)

    # 聚合分类信息
    if (!is.na(agg_taxonomies)) {
      aggregated_taxonomies <- data.frame(
        Category = "Aggregated Taxonomies",
        Variable = "The taxonomies that have been aggregated",
        Value = paste(agg_taxonomies, collapse = ", ")
      )

      table1 <- rbind(table1, aggregated_taxonomies)
    }

    # 时间序列信息统计
    if (!is.null(time.var) && "meta.dat" %in% names(data.obj)) {
      if (time.var %in% colnames(data.obj$meta.dat)) {
        time_var_data <- data.obj$meta.dat[[time.var]]

        # 数值型时间数据
        if (is.numeric(time_var_data)) {
          time_stats <- data.frame(
            Category = "Time-Series Information",
            Variable = c(
              "Earliest sample time-point",
              "Latest sample time-point"
            ),
            Value = c(min(time_var_data), max(time_var_data))
          )
        }

        # 类别型（字符或因子）时间数据
        if (is.character(time_var_data) || is.factor(time_var_data)) {
          time_stats <- data.frame(
            Category = "Time-Series Information",
            Variable = "Number of unique time-points",
            Value = length(unique(time_var_data))
          )
        }

        table1 <- rbind(table1, time_stats)

        # 分布情况
        time_table <- table(data.obj$meta.dat[[time.var]])
        time_df <- as.data.frame(time_table)
        colnames(time_df) <- c("TimePoint", "SampleCount")

        for (i in 1:nrow(time_df)) {
          distribution <- data.frame(
            Category = "Distribution of sample counts",
            Variable = paste("Sample Count at Time-point:", time_df$TimePoint[i]),
            Value = time_df$SampleCount[i]
          )
          table1 <- rbind(table1, distribution)
        }
      } else {
        #cat("The provided time variable does not exist in the metadata.\n")
      }
    }

    # 添加样本数量和特征数量
    num_samples <- ncol(feature_tab)
    num_features <- nrow(feature_tab)
    sample_feature_stats <- data.frame(
      Category = c("Basic Statistics", "Basic Statistics"),
      Variable = c("Number of samples", "Number of features"),
      Value = c(num_samples, num_features)
    )
    table1 <- rbind(sample_feature_stats, table1)

    # 检测Value列中的哪些值是数值型
    is_numeric <- sapply(table1$Value, function(x) grepl("^-?[0-9.]+$", x))

    # 对于数值型数据，检查是否有超过三位小数
    has_more_than_three_decimals <- sapply(table1$Value[is_numeric], function(x) {
      decimal_part <- sub(".*\\.", "", x)
      nchar(decimal_part) > 3
    })

    # 对于满足条件的数值，进行四舍五入至三位小数
    table1$Value[is_numeric][has_more_than_three_decimals] <- round(as.numeric(table1$Value[is_numeric][has_more_than_three_decimals]), 3)

    # 返回表格形式的数据
    return(as_tibble(table1))
  }
