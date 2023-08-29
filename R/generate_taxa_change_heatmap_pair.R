#' @title Generate Taxa Change Heatmap Pair
#'
#' @description This function generates a heatmap showing the pairwise changes in relative abundances of taxa between different time points.
#' The data used in this visualization will be first filtered based on prevalence and abundance thresholds.
#' The plot can either be displayed interactively or saved as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param change.base A numeric value specifying the baseline time point for computing change. This should match one of the time points in the time variable. Default is 1, which assumes the first time point is the baseline.
#' @param feature.change.func A function or character string specifying the method for computing change between time points. Options are:
#' - "difference": Compute simple difference of abundances between time points.
#' - "lfc": Compute log2 fold change between time points.
#' - "relative change": Compute relative change of abundances between time points.
#' - "custom function": Apply a custom function to compute change. The custom function should take two arguments value_time_2 and value_time_1.
#' Default is "difference".
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
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot, when `features.plot` is NULL.
#' Default is NULL, in which case all features passing filters will be plotted.
#' @param top.k.func Function to use for selecting top k features, when `features.plot` is NULL.
#' Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
#' case features will be selected by abundance.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size Base font size for the generated plots.
#' @param palette Color palette used for the plots.
#' @param cluster.rows A logical variable indicating if rows should be clustered. Default is TRUE.
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is NULL.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed to pheatmap function
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the pheatmap::pheatmap plot. If `pdf` is set to FALSE, the function will return the pheatmap plot without creating a PDF file.
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(pheatmap)
#' data(peerj32.obj)
#' generate_taxa_change_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = NULL,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_change_heatmap_pair <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              change.base = NULL,
                                              feature.change.func = "relative change",
                                              feature.level = NULL,
                                              feature.dat.type = c("count", "proportion", "other"),
                                              features.plot = NULL,
                                              top.k.plot = NULL,
                                              top.k.func = NULL,
                                              prev.filter = 0.01,
                                              abund.filter = 0.01,
                                              base.size = 10,
                                              palette = NULL,
                                              cluster.rows = NULL,
                                              cluster.cols = NULL,
                                              pdf = TRUE,
                                              file.ann = NULL,
                                              pdf.wid = 11,
                                              pdf.hei = 8.5,
                                              ...) {
  # Extract data
  mStat_validate_data(data.obj)

  feature.dat.type <- match.arg(feature.dat.type)

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

  meta_tab <-  load_data_obj_metadata(data.obj) %>% select(all_of(c(
    time.var, group.var, strata.var, subject.var
  ))) %>% rownames_to_column("sample")

  if (is.null(cluster.cols)) {
    if (!is.null(group.var)) {
      if (is.numeric(meta_tab %>% select(all_of(group.var)) %>% dplyr::pull()) &&
          !is.integer(meta_tab %>% select(all_of(group.var)) %>% dplyr::pull())) {
        cluster.cols = FALSE
      } else {
        cluster.cols = TRUE
      }
    } else {
      cluster.cols = TRUE
    }
  } else {

  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  # Merge OTU table with taxonomy table
  otu_tax <-
    cbind(otu_tab, tax_tab)

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  plot_list <- lapply(feature.level, function(feature.level) {
    # Filter taxa based on prevalence and abundance
    otu_tax_filtered <- otu_tax %>%
      tidyr::gather(key = "sample", value = "count",-one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(!!sym(feature.level))) %>%
      dplyr::summarise(total_count = mean(count),
                       prevalence = sum(count > 0) / dplyr::n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count,-prevalence) %>%
      dplyr::left_join(otu_tax, by = feature.level)

    # Aggregate OTU table
    otu_tax_agg <- otu_tax_filtered %>%
      tidyr::gather(key = "sample", value = "count",-one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(sample,!!sym(feature.level))) %>%
      dplyr::summarise(count = sum(count)) %>%
      tidyr::spread(key = "sample", value = "count")

    compute_function <- function(top.k.func) {
      if (is.function(top.k.func)) {
        results <-
          top.k.func(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
      } else {
        switch(top.k.func,
               "mean" = {
                 results <-
                   rowMeans(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                            na.rm = TRUE)
               },
               "sd" = {
                 results <-
                   matrixStats::rowSds(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                                       na.rm = TRUE)
                 names(results) <-
                   rownames(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
               },
               stop("Invalid function specified"))
      }

      return(results)
    }

    if (is.null(features.plot) &&
        !is.null(top.k.plot) && !is.null(top.k.func)) {
      features.plot <-
        names(sort(compute_function(top.k.func), decreasing = TRUE)[1:top.k.plot])
    }

    # Convert counts to numeric
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    # 将otu_tax_agg_numeric从宽格式转换为长格式
    otu_tax_long <- otu_tax_agg_numeric %>%
      tidyr::gather(key = "sample", value = "value",-feature.level)

    # 将otu_tax_long和meta_tab按sample列连接
    merged_data <- otu_tax_long %>%
      dplyr::inner_join(meta_tab, by = "sample")

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
    if (is.function(feature.change.func)) {
      combined_data <-
        combined_data %>% dplyr::mutate(value_diff = feature.change.func(value_time_2, value_time_1))
    } else if (feature.change.func == "lfc") {
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
    } else if (feature.change.func == "relative change") {
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
      select(all_of(c(subject.var, group.var, strata.var))) %>%
      dplyr::distinct(subject, .keep_all = TRUE) %>% as_tibble()

    # 获取元素的顺序
    order_index <-
      match(colnames(value_diff_matrix),
            as.matrix(unique_meta_tab %>% select(!!sym(subject.var))))

    suppressWarnings({
      # 根据顺序对 unique_meta_tab 进行排序
      sorted_meta_tab <- unique_meta_tab[order_index,]
      rownames(sorted_meta_tab) <-
        as.matrix(sorted_meta_tab[, subject.var])
    })

    # 如果 group.var 不为空，则根据 group.var 对 sorted_meta_tab 进行排序
    if (!is.null(group.var)) {
      sorted_meta_tab <-
        sorted_meta_tab[order(sorted_meta_tab %>% select(!!sym(group.var)) %>% as.matrix()),]
    }

    # 从 sorted_meta_tab 中提取 subject 列
    subjects <- sorted_meta_tab %>% select(all_of(subject.var))

    # 使用这些索引以正确的顺序重新排列 value_diff_matrix 的列
    value_diff_matrix <- value_diff_matrix[, as.matrix(subjects)]

    # 如果 strata.var 为 NULL，仅使用 group.var 进行注释
    if (is.null(strata.var)) {
      annotation_cols <-
        sorted_meta_tab %>% select(all_of(c(subject.var, group.var))) %>% column_to_rownames(var = subject.var)
      if (is.null(group.var)) {
        annotation_cols <- NULL
      }
    } else {
      annotation_cols <-
        sorted_meta_tab %>% select(all_of(c(group.var, strata.var, subject.var))) %>% column_to_rownames(var = subject.var)
    }

    annotation_col_sorted <-
      annotation_cols[order(annotation_cols[[group.var]]), ]

    if (!is.null(strata.var)) {
      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]], annotation_col_sorted[[group.var]]),]

    }

    value_diff_matrix <-
      value_diff_matrix[, rownames(annotation_col_sorted)]

    if (!is.null(strata.var)) {
      gaps <-
        cumsum(table(sorted_meta_tab[[strata.var]]))[-length(sorted_meta_tab[[strata.var]])]
    } else {
      if (!is.null(group.var)) {
        # 计算分隔线应该出现的位置
        gaps <-
          cumsum(table(sorted_meta_tab[[group.var]]))[-length(sorted_meta_tab[[group.var]])]
      } else {
        gaps <- NULL
      }
    }

    if (!is.null(features.plot)) {
      value_diff_matrix <-
        value_diff_matrix[rownames(value_diff_matrix) %in% features.plot,]
    }

    n_colors <- 100

    col <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
    # 找出数据中的最大绝对值
    max_abs_val <- max(abs(range(na.omit(
      c(value_diff_matrix)
    ))))

    # 计算零值在新色彩向量中的位置
    zero_pos <- round(max_abs_val / (2 * max_abs_val) * n_colors)

    # 创建颜色向量
    my_col <-
      c(
        colorRampPalette(col[1:3])(zero_pos),
        colorRampPalette(col[3:5])(n_colors - zero_pos + 1)
      )

    break_points <-
      seq(-max_abs_val, max_abs_val, length.out = length(my_col) + 1)

    if (is.null(palette)){
      color_vector <- c(
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
    } else {
      color_vector <- palette
    }

    # 为演示目的，假设这些是您的唯一值
    group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
    strata_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()

    # 为 group.var 分配颜色
    group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)

    # 为 strata.var 分配颜色
    strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

    # 创建注释颜色列表
    annotation_colors_list <- setNames(
      list(group_colors, strata_colors),
      c(group.var, strata.var)
    )

    heatmap_plot <- pheatmap::pheatmap(
      value_diff_matrix,
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      annotation_legend = TRUE,
      show_colnames = TRUE,
      show_rownames = TRUE,
      border_color = NA,
      silent = TRUE,
      gaps_col = gaps,
      fontsize = base.size,
      color = my_col,
      breaks = break_points,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    if (is.function(feature.change.func)) {
      feature.change.func = "custom function"
    }

    # Save the heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_pair",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "change_base_",
        change.base,
        "_",
        "feature_level_",
        feature.level,
        "_",
        "prev_filter_",
        prev.filter,
        "_",
        "abund_filter_",
        abund.filter
      )
      if (!is.null(group.var)) {
        pdf_name <- paste0(pdf_name, "_", "group_", group.var)
      }
      if (!is.null(strata.var)) {
        pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
      }
      if (!is.null(file.ann)) {
        pdf_name <- paste0(pdf_name, "_", file.ann)
      }
      pdf_name <- paste0(pdf_name, ".pdf")
      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_heatmap_plot
      )
    }
    return(gg_heatmap_plot)
  })

  return(plot_list)
}
