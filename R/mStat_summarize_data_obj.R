#' Summarize a MicrobiomeStat Data Object
#'
#' This function takes a MicrobiomeStat data object and provides a summary of the various components such as feature.tab, meta.dat, feature.ann, and optionally time.var and tree.
#' The function is useful to give an overview of the input data object before any further microbiome statistical analysis. Providing a time variable can offer additional insights into the temporal dynamics of your samples.
#'
#' @name mStat_summarize_data_obj
#' @param data.obj A MicrobiomeStat data object to be summarized.
#' @param time.var An optional time variable column name from meta.dat. If provided, the function will additionally dplyr::summarize the temporal distribution of samples. Default is NULL.
#' @param group.var An optional grouping variable column name from meta.dat. If provided together with time.var, the histogram of sample counts over time will be grouped by this variable. Default is NULL.
#' @param palette A vector of colors to use for grouping in the histogram. Default is a preset color palette with 10 colors.
#' @return A list containing summaries for:
#' \itemize{
#'   \item feature.tab: The overall summary includes the number of features, number of samples, sparsity, and description of singleton features.
#'   \item meta.dat: The overall summary includes the number of samples, number of metadata fields, missing data, and sample distribution if a time variable is provided. If time.var is provided, a distribution of samples over time is also generated.
#'   \item feature.ann: The overall summary includes number of features, number of annotations, and the proportion of NA values for each annotation.
#'   \item tree: (optional) Whether a phylogenetic tree exists in the data object.
#' }
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

    # Print summary
    cat(
      "Microbiome Feature Table Summary:\n",
      "==============================\n",
      "Min. reads per sample: ",
      min_reads,
      "\n",
      "Max. reads per sample: ",
      max_reads,
      "\n",
      "Total reads dplyr::across all samples: ",
      total_reads,
      "\n",
      "Average reads per sample: ",
      ave_reads_per_sample,
      "\n",
      "Median reads per sample: ",
      median_reads_per_sample,
      "\n",
      "Proportion of zero counts: ",
      zero_count_prop,
      "\n",
      "Count of features that only appear once: ",
      count_single_occurrence,
      "\n",
      "Number of metadata variables: ",
      num_meta_vars,
      "\n",
      "Metadata variables are: ",
      paste(meta_var_names, collapse = ", "),
      "\n\n",

      "Feature Annotation Summary:\n",
      "========================\n",
      paste(
        "Proportion of missing annotations in ",
        colnames(feature_ann),
        " are: ",
        NA_props,
        collapse = "\n"
      ),
      "\n\n",

      "Phylogenetic Tree Information:\n",
      "===========================\n",
      "Exists in the dataset? ",
      tree_exists,
      "\n\n",

      "Aggregated Taxonomies:\n",
      "===================\n",
      "The taxonomies that have been aggregated are: ",
      paste(agg_taxonomies, collapse = ", "),
      "\n\n"
    )

    # Handling the time-series data
    if (!is.null(time.var) && "meta.dat" %in% names(data.obj)) {
      if (time.var %in% colnames(data.obj$meta.dat)) {
        time_var_data <- data.obj$meta.dat[[time.var]]

        cat("Time-Series Information:\n",
            "========================\n")

        # Numeric time data
        if (is.numeric(time_var_data)) {
          cat(
            "Earliest sample time-point: ",
            min(time_var_data),
            "\n",
            "Latest sample time-point: ",
            max(time_var_data),
            "\n\n"
          )
        }

        # Categorical (character or factor) time data
        if (is.character(time_var_data) || is.factor(time_var_data)) {
          cat("Unique time-points: ", length(unique(time_var_data)), "\n\n")
        }

        cat(
          "Distribution of sample counts at each time-point:\n",
          "========================\n"
        )

        # Create a histogram for the time variable
        if (!is.null(group.var) &&
            group.var %in% colnames(data.obj$meta.dat)) {
          # Create a grouped data frame
          grouped_df <- data.obj$meta.dat %>%
            select(all_of(c(time.var, group.var))) %>%
            dplyr::group_by_at(vars(time.var, group.var)) %>%
            dplyr::summarise(SampleCount = dplyr::n(), .groups = "drop")

          # Print the grouped data frame
          print(grouped_df)

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
                x = "Time Point",
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
          print(time_df)

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
      } else {
        cat("The provided time variable does not exist in the metadata.\n")
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
      metadata_stats <- data.frame(
        Category = "Metadata",
        Variable = c("Number of metadata variables", "Metadata variables"),
        Value = c(num_meta_vars, paste(meta_var_names, collapse = ", "))
      )

      table1 <- rbind(table1, metadata_stats)
    }

    # feature.ann 缺失注释的统计
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
            Variable = "Unique time-points",
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
        cat("The provided time variable does not exist in the metadata.\n")
      }
    }

    # 返回表格形式的数据
    return(as_tibble(table1))
  }
