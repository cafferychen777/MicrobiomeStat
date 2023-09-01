#' @title Generate Taxonomic Change Heatmap Long
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var (Optional) A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param feature.level A character string defining the taxonomic level to analyze ('Phylum', 'Family', or 'Genus').
#' @param feature.change.func A function or character string specifying how to calculate
#' the change from baseline value. This allows flexible options:
#' - If a function is provided, it will be applied to each row to calculate change.
#'   The function should take 2 arguments: value at timepoint t and value at baseline t0.
#' - If a character string is provided, following options are supported:
#'   - 'relative change': (value_t - value_t0) / (value_t + value_t0)
#'   - 'difference': value_t - value_t0
#'   - 'lfc': log2(value_t + 1e-5) - log2(value_t0 + 1e-5)
#' - Default is 'relative change'.
#'
#' If none of the above options are matched, an error will be thrown indicating
#' the acceptable options or prompting the user to provide a custom function.
#' @details This parameter is used to compute the change columns from baseline
#'   (specified by t0.level) for each taxon. The change values are calculated
#'   for each timepoint and appended as new columns in the data frame before
#'   plotting heatmap. This allows flexibly customizing how change is quantified.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
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
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is FALSE.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional arguments passed to pheatmap.
#' @return A list of ggplot heatmap objects, one for each taxonomic level.
#'
#' @details This function generates a separate heatmap for each taxonomic level specified,
#'   with rows clustered and layers arranged by groups over timepoints.
#'   It automatically rarefies raw count data using Rarefy-TSS normalization in MicrobiomeStat.
#'   Annotation columns are generated and ordered properly for visually stacking the layers.
#'   Colormaps are also generated for group and strata variables.
#'
#' @seealso \code{\link{pheatmap}} for heatmap, \code{\link{mStat_normalize_data}} for data normalization.
#'
#' @examples
#' \dontrun{
#' library(pheatmap)
#' data(ecam.obj)
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   feature.level = c("Family","Class"),
#'   feature.change.func = "lfc",
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   palette = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' data(subset_T2D.obj)
#' generate_taxa_change_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Phylum"),
#'   feature.change.func = "lfc",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#'
#' @return An object of class pheatmap, the generated heatmap plot
#' @export
#' @importFrom dplyr distinct pull
#' @seealso \code{\link{pheatmap}}
generate_taxa_change_heatmap_long <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              t0.level = NULL,
                                              ts.levels = NULL,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              feature.level,
                                              feature.change.func = "relative change",
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

  feature.dat.type <- match.arg(feature.dat.type)

  mStat_validate_data(data.obj)

  if (!is.character(subject.var))
    stop("`subject.var` should be a character string.")
  if (!is.character(time.var))
    stop("`time.var` should be a character string.")
  if (!is.null(group.var) &&
      !is.character(group.var))
    stop("`group.var` should be a character string or NULL.")
  if (!is.null(strata.var) &&
      !is.character(strata.var))
    stop("`strata.var` should be a character string or NULL.")

  # 提取数据
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  meta_tab <- load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- "ALL"
  }

  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var), !!sym(strata.var)))
  }

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  if (feature.dat.type == "count"){
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
    )
    data.obj <- mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm
  }

  plot_list <- lapply(feature.level, function(feature.level) {

    if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
      data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    if (feature.level != "original"){
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- load_data_obj_count(data.obj)
    }

    otu_tax_agg <-  otu_tax_agg %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter) %>%
      rownames_to_column(feature.level)

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

    otu_tab_norm <-
      otu_tax_agg %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # Calculate the mean_value for each combination of feature.level, group.var, and time.var
    df_mean_value <- otu_tab_norm %>%
      as.data.frame() %>%
      rownames_to_column(var = feature.level) %>%
      tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
      dplyr::left_join(meta_tab %>%
                  rownames_to_column("sample"), "sample") %>%
      dplyr::group_by(!!sym(feature.level),!!sym(group.var),!!sym(time.var)) %>%
      dplyr::summarise(mean_value = mean(value), .groups = "drop")

    # Spread the data to wide format
    df_wide <- df_mean_value %>%
      tidyr::spread(key = time.var, value = mean_value)

    if (is.null(t0.level)) {
      if (is.numeric(meta_tab[, time.var])) {
        t0.level <- sort(unique(meta_tab[, time.var]))[1]
      } else {
        t0.level <- levels(meta_tab[, time.var])[1]
      }
    }

    if (is.null(ts.levels)) {
      if (is.numeric(meta_tab[, time.var])) {
        ts.levels <- sort(unique(meta_tab[, time.var]))[-1]
      } else {
        ts.levels <- levels(meta_tab[, time.var])[-1]
      }
    }

    # Calculate the changes from t0.level
    for (ts in ts.levels) {
      change_col_name <- paste0("change_", ts)

      if (is.function(feature.change.func)) {
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := feature.change.func(.data[[as.character(ts)]], .data[[as.character(t0.level)]]))
      } else if (feature.change.func == "relative change") {
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := dplyr::case_when(
            (.data[[as.character(ts)]] + .data[[as.character(t0.level)]]) != 0 ~ (.data[[as.character(ts)]] - .data[[as.character(t0.level)]]) / (.data[[as.character(ts)]] + .data[[as.character(t0.level)]]),
            TRUE ~ 0
          ))
      } else if (feature.change.func == "difference") {
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := (.data[[as.character(ts)]] - .data[[as.character(t0.level)]]))
      } else if (feature.change.func == "lfc") {
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := log2(.data[[as.character(ts)]] + 0.00001) - log2(.data[[as.character(t0.level)]] + 0.00001))
      } else {
        stop(
          "`feature.change.func` must be either 'relative change', 'difference', 'lfc' or a function."
        )
      }
    }

    df_wide <-
      df_wide[, c(feature.level, group.var, paste0("change_", ts.levels))]

    # 首先将数据框转化为长格式
    df_long <- df_wide %>%
      tidyr::pivot_longer(
        cols = starts_with("change"),
        names_to = "time",
        values_to = "value"
      )

    # 删除 "time" 列名中的 "change_" 部分
    df_long$time <- gsub("change_", "", df_long$time)

    # 再将数据框转化为宽格式，并在此过程中修改列名
    df_wide_new <- df_long %>%
      tidyr::unite("group_time", c(group.var, "time"), sep = "_") %>%
      tidyr::pivot_wider(names_from = "group_time",
                  values_from = "value")

    wide_data <- df_wide_new %>% column_to_rownames(feature.level)

    df_wide <-
      df_wide[, c(feature.level, group.var, paste0("change_", ts.levels))]

    # 首先将数据框转化为长格式
    df_long <- df_wide %>%
      tidyr::pivot_longer(
        cols = starts_with("change"),
        names_to = "time",
        values_to = "value"
      )

    # 删除 "time" 列名中的 "change_" 部分
    df_long$time <- gsub("change_", "", df_long$time)

    # 再将数据框转化为宽格式，并在此过程中修改列名
    df_wide_new <- df_long %>%
      tidyr::unite("group_time", c(group.var, "time"), sep = "_") %>%
      tidyr::pivot_wider(names_from = "group_time",
                  values_from = "value")

    wide_data <- df_wide_new %>% column_to_rownames(feature.level)

    # Save original column names
    original_colnames <- colnames(wide_data)

    # Remove columns in data frame where all values are NA
    wide_data <-
      wide_data[, colSums(is.na(wide_data)) != nrow(wide_data)]

    # Save new column names
    new_colnames <- colnames(wide_data)

    # Find out removed columns
    removed_colnames <- setdiff(original_colnames, new_colnames)

    # If columns were removed, send a message to the user
    if (length(removed_colnames) > 0) {
      message(
        "The following combinations were all NA, therefore they have been removed: ",
        paste(removed_colnames, collapse = ", ")
      )
    } else {
      message("No columns have been removed.")
    }

    # Sort samples by group.var if not NULL
    annotation_col <- meta_tab %>%
      select(!!sym(time.var), !!sym(group.var)) %>%
      filter(!!sym(time.var) != t0.level) %>%
      as_tibble() %>%
      dplyr::distinct() %>%
      dplyr::mutate(group_time = paste(!!sym(group.var), !!sym(time.var), sep = "_")) %>%
      filter(group_time %in% new_colnames) %>%
      column_to_rownames("group_time")

    annotation_col_sorted <-
      annotation_col[order(annotation_col[[group.var]], annotation_col[[time.var]]),]

    if (!is.null(strata.var)) {
      annotation_col_sorted <- annotation_col_sorted %>%
        tidyr::separate(!!sym(group.var),
                 into = c(group.var, strata.var),
                 sep = "\\.")

      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]], annotation_col_sorted[[group.var]], annotation_col_sorted[[time.var]]), ]

    }

    if (group.var == "ALL") {
      annotation_col_sorted <-
        annotation_col_sorted %>% select(all_of(c(time.var)))
    }

    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted)]

    if (!is.null(features.plot)) {
      wide_data_sorted <-
        wide_data_sorted[rownames(wide_data_sorted) %in% features.plot,]
    }

    #Calculate gaps if group.var is not NULL
    if (!is.null(group.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[group.var]]))[-length(table(annotation_col_sorted[[group.var]]))]
    } else {
      gaps <- NULL
    }

    if (!is.null(strata.var)){
      gaps <-
        cumsum(table(annotation_col_sorted[[strata.var]]))[-length(table(annotation_col_sorted[[strata.var]]))]
    }

    n_colors <- 100

    #if (is.null(palette)) {
      col <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
      # 找出数据中的最大绝对值
      max_abs_val <-
        max(abs(range(na.omit(
          c(as.matrix(wide_data_sorted))
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

      if (!is.null(strata.var) & !is.null(group.var)){
        # 为演示目的，假设这些是您的唯一值
        group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()

        # 为 group.var 分配颜色
        group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)

        strata_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
        # 为 strata.var 分配颜色
        strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

        # 创建注释颜色列表
        annotation_colors_list <- setNames(
          list(group_colors, strata_colors),
          c(group.var, strata.var)
        )
      } else if (!is.null(group.var) & group.var != "ALL"){
        # 为演示目的，假设这些是您的唯一值
        group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
        # 为 group.var 分配颜色
        group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
        # 创建注释颜色列表
        annotation_colors_list <- setNames(
          list(group_colors),
          c(group.var)
        )
      } else {
        annotation_colors_list <- NULL
      }

      # Plot stacked heatmap
      heatmap_plot <- pheatmap::pheatmap(
        wide_data_sorted,
        annotation_col = annotation_col_sorted,
        annotation_colors = annotation_colors_list,
        cluster_rows = cluster.rows,
        cluster_cols = cluster.cols,
        show_colnames = FALSE,
        gaps_col = gaps,
        fontsize = base.size,
        silent = TRUE,
        color = my_col,
        breaks = break_points,
        ...
      )
    # } else {
    #   # 创建颜色映射函数
    #   my_palette <- colorRampPalette(palette)
    #
    #   # Plot stacked heatmap
    #   heatmap_plot <- pheatmap::pheatmap(
    #     wide_data_sorted,
    #     annotation_col = annotation_col_sorted,
    #     cluster_rows = cluster.rows,
    #     cluster_cols = cluster.cols,
    #     show_colnames = FALSE,
    #     gaps_col = gaps,
    #     fontsize = base.size,
    #     silent = TRUE,
    #     color = my_palette(n_colors)
    #   )
    # }

    # 计算颜色的数量
    # 这通常取决于你的数据，你可能需要根据你的实际情况进行调整

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_long",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "group_",
        group.var,
        "_",
        "strata_",
        strata.var,
        "_",
        "taxa_",
        feature.level,
        "_",
        "prev_filter_",
        prev.filter,
        "_",
        "abund_filter_",
        abund.filter
      )

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

  names(plot_list) <- feature.level

  # Return the heatmap plot for display
  return(plot_list)
}
