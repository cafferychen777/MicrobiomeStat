#' Generate Taxa Change Heatmap Pair
#'
#' This function generates a heatmap showing the pairwise changes in relative abundances of taxa between different time points.
#' The data used in this visualization will be first filtered based on prevalence and abundance thresholds.
#' The plot can either be displayed interactively or saved as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string defining subject variable in the `meta_tab`
#' @param time.var A character string defining time variable in the `meta_tab`
#' @param group.var A character string defining group variable in the `meta_tab` used for sorting and facetting
#' @param strata.var A character string defining strata variable in the `meta_tab` used for sorting and facetting
#' @param change.base A numeric value setting base for the change (usually 1)
#' @param change.func A character string specifying the function to apply on the change, default is "log"
#' @param zero.handle A character string specifying how to handle zeros, default is "pseudo"
#' @param feature.level A character string defining the taxonomic level to analyze, e.g., "Phylum", "Family", or "Genus"
#' @param prev.filter A numeric value defining the prevalence threshold to filter taxa, between 0 and 1
#' @param abund.filter A numeric value defining the abundance threshold to filter taxa
#' @param pdf A logical value. If TRUE (default), saves the heatmap as a PDF file. If FALSE, the heatmap will be displayed interactively without creating a PDF
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the pheatmap plot. If `pdf` is set to FALSE, the function will return the pheatmap plot without creating a PDF file.
#' @examples
#' # Load required libraries and example data
#' library(microbiome)
#' library(tidyverse)
#' library(pheatmap)
#' library(circlize)
#' data(peerj32)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#' plot_list <- generate_taxa_change_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   change.base = "1",
#'   change.func = "relative difference",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.1,
#'   abund.filter = 0.1,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = NULL,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' @export
#'
#' @seealso \code{\link{pheatmap}}
generate_taxa_change_heatmap_pair <- function(data.obj,
                                                    subject.var,
                                                    time.var,
                                                    group.var = NULL,
                                                    strata.var = NULL,
                                                    change.base = NULL,
                                                    change.func = "relative difference",
                                                    feature.level = NULL,
                                                    feature.dat.type = c("count", "proportion", "other"),
                                                    features.plot = NULL,
                                                    top.k.plot = NULL,
                                                    top.k.func = NULL,
                                                    prev.filter = 0.1,
                                                    abund.filter = 0.1,
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
    {if("original" %in% feature.level) mutate(., original = rownames(.)) else .} %>%
    select(all_of(feature.level))

  meta_tab <-  load_data_obj_metadata(data.obj) %>% select(all_of(c(
    time.var, group.var, strata.var, subject.var
  ))) %>% rownames_to_column("sample")

  if (is.null(cluster.cols)){
    if (!is.null(group.var)){
      if (is.numeric(meta_tab %>% select(all_of(group.var)) %>% pull()) && !is.integer(meta_tab %>% select(all_of(group.var)) %>% pull())){
        cluster.cols = FALSE
      } else {
        cluster.cols = TRUE
      }
    } else {
      cluster.cols = TRUE
    }
  } else {

  }

  if (is.null(cluster.rows)){
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
      gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      group_by_at(vars(!!sym(feature.level))) %>%
      summarise(total_count = mean(count),
                prevalence = sum(count > 0) / n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count, -prevalence) %>%
      left_join(otu_tax, by = feature.level)

    # Aggregate OTU table
    otu_tax_agg <- otu_tax_filtered %>%
      gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      group_by_at(vars(sample, !!sym(feature.level))) %>%
      summarise(count = sum(count)) %>%
      spread(key = "sample", value = "count")

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
                   rowSds(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                          na.rm = TRUE)
                 names(results) <- rownames(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
               },
               stop("Invalid function specified"))
      }

      return(results)
    }

    if (is.null(features.plot) &&
        !is.null(top.k.plot) && !is.null(top.k.func)) {
      features.plot <- names(sort(compute_function(top.k.func), decreasing = TRUE)[1:top.k.plot])
    }

    # Convert counts to numeric
    otu_tax_agg_numeric <-
      mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    # 将otu_tax_agg_numeric从宽格式转换为长格式
    otu_tax_long <- otu_tax_agg_numeric %>%
      gather(key = "sample", value = "value", -feature.level)

    # 将otu_tax_long和meta_tab按sample列连接
    merged_data <- otu_tax_long %>%
      inner_join(meta_tab, by = "sample")

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
      combined_data <- combined_data %>% mutate(value_diff = change.func(value_time_2, value_time_1))
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

      combined_data <- left_join(combined_data, half_nonzero_min_time_2, by = feature.level, suffix = c("_time_1", "_time_2"))
      combined_data <- left_join(combined_data, half_nonzero_min_time_1, by = feature.level, suffix = c("_time_1", "_time_2"))
      combined_data$value_time_2[combined_data$value_time_2 == 0] <- combined_data$half_nonzero_min_time_2[combined_data$value_time_2 == 0]
      combined_data$value_time_1[combined_data$value_time_1 == 0] <- combined_data$half_nonzero_min_time_1[combined_data$value_time_1 == 0]

      # Add a message to inform users that an imputation operation was performed.
      message("Imputation was performed using half the minimum nonzero proportion for each taxon at different time points.")

      combined_data <- combined_data %>% mutate(value_diff = log2(value_time_2) - log2(value_time_1))
    } else if (change.func == "relative difference"){
      combined_data <- combined_data %>%
        mutate(value_diff = case_when(
          value_time_2 == 0 & value_time_1 == 0 ~ 0,
          TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
        ))
    } else {
      combined_data <- combined_data %>% mutate(value_diff = value_time_2 - value_time_1)
    }

    value_diff_matrix <- combined_data %>%
      select(feature.level, subject = subject, value_diff) %>%
      spread(key = subject, value = value_diff) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    unique_meta_tab <- meta_tab %>%
      filter(subject %in% colnames(value_diff_matrix)) %>%
      select(all_of(c(subject.var, group.var, strata.var))) %>%
      distinct(subject, .keep_all = TRUE) %>% as_tibble()

    # 获取元素的顺序
    order_index <-
      match(colnames(value_diff_matrix),
            as.matrix(unique_meta_tab %>% select(!!sym(subject.var))))

    suppressWarnings({
    # 根据顺序对 unique_meta_tab 进行排序
    sorted_meta_tab <- unique_meta_tab[order_index, ]
    rownames(sorted_meta_tab) <-
      as.matrix(sorted_meta_tab[, subject.var])
    })

    # 如果 group.var 不为空，则根据 group.var 对 sorted_meta_tab 进行排序
    if (!is.null(group.var)) {
      sorted_meta_tab <-
        sorted_meta_tab[order(sorted_meta_tab %>% select(!!sym(group.var)) %>% as.matrix()), ]
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

    if (!is.null(group.var)) {
      # 计算分隔线应该出现的位置
      gaps <-
        cumsum(table(sorted_meta_tab[[group.var]]))[-length(sorted_meta_tab[[group.var]])]
    } else {
      gaps <- NULL
    }

    if (!is.null(features.plot)){
      value_diff_matrix <- value_diff_matrix[rownames(value_diff_matrix) %in% features.plot, ]
    }

    n_colors <- 100

    if (is.null(palette)) {
      palette <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
      # 找出数据中的最大绝对值
      max_abs_val <- max(abs(range(na.omit(c(value_diff_matrix)))))

      # 计算零值在新色彩向量中的位置
      zero_pos <- round(max_abs_val / (2 * max_abs_val) * n_colors)

      # 创建颜色向量
      my_palette <- c(colorRampPalette(palette[1:3])(zero_pos), colorRampPalette(palette[3:5])(n_colors - zero_pos + 1))

      heatmap_plot <- pheatmap(
        value_diff_matrix,
        annotation_col = annotation_cols,
        cluster_rows = cluster.rows,
        cluster_cols = cluster.cols,
        annotation_legend = TRUE,
        show_colnames = TRUE,
        show_rownames = TRUE,
        border_color = NA,
        gaps_col = gaps,
        fontsize = base.size,
        color = my_palette,
        ...
      )
    } else {
      # 创建颜色映射函数
      my_palette <- colorRampPalette(palette)

      # Plot stacked heatmap
      heatmap_plot <- pheatmap(
        value_diff_matrix,
        annotation_col = annotation_cols,
        cluster_rows = cluster.rows,
        cluster_cols = cluster.cols,
        annotation_legend = TRUE,
        show_colnames = TRUE,
        show_rownames = TRUE,
        border_color = NA,
        gaps_col = gaps,
        fontsize = base.size,
        color = my_palette(n_colors),
        ...
      )
    }

    gg_heatmap_plot <- ggplotify::as.ggplot(heatmap_plot)

    if (is.function(change.func)){
      change.func = "custom function"
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
        "change_func_",
        change.func,
        "_",
        "zero_handle_",
        zero.handle,
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
