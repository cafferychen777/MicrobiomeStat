#' @title Generate Taxonomic Heatmap Long
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
#' @param strata.var A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param feature.level A character vector specifying the taxa level(s) to include in the analysis. Default is c('Phylum', 'Family', 'Genus').
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
#' @param pdf A logical value. If TRUE (default), saves the plot as a PDF file. If FALSE, the plot will be displayed interactively without creating a PDF.
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed to the pheatmap() function from the “pheatmap” package.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Family","Phylum","Genus"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 15
#' )
#' }
#' @export
#' @importFrom dplyr distinct pull
#' @seealso \code{\link{pheatmap}}
generate_taxa_heatmap_long <- function(data.obj,
                                       subject.var,
                                       time.var,
                                       t0.level,
                                       ts.levels,
                                       group.var = NULL,
                                       strata.var = NULL,
                                       feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       features.plot = NULL,
                                       top.k.plot = NULL,
                                       top.k.func = NULL,
                                       prev.filter = 0.01,
                                       abund.filter = 0.01,
                                       base.size = 10,
                                       palette = NULL,
                                       cluster.cols = NULL,
                                       cluster.rows = NULL,
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
  data.obj <-
    mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
      subject.var, group.var, time.var, strata.var
    )))

  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- "ALL"
  }

  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(strata.var), !!sym(group.var)))
  }

  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  if (feature.dat.type == "count") {
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
    )
    data.obj <-
      mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm
  }

  plot_list <- lapply(feature.level, function(feature.level) {
    if (is.null(data.obj$feature.agg.list[[feature.level]]) &
        feature.level != "original") {
      data.obj <-
        mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    if (feature.level != "original") {
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

    otu_tab_norm <-
      otu_tax_agg_numeric %>%
      filter(!is.na(!!sym(feature.level))) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    if (!is.null(group.var)) {
      wide_data <- otu_tab_norm %>%
        as.data.frame() %>%
        rownames_to_column(var = feature.level) %>%
        tidyr::gather(key = "sample",
                      value = "value",
                      -one_of(feature.level)) %>%
        dplyr::left_join(meta_tab %>%
                           rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level),
                        !!sym(group.var),
                        !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("group_time", c(group.var, time.var), sep = "_") %>%
        tidyr::spread(key = "group_time", value = "mean_value") %>% column_to_rownames(feature.level)
    } else {
      wide_data <- otu_tab_norm %>%
        as.data.frame() %>%
        rownames_to_column(var = feature.level) %>%
        tidyr::gather(key = "sample",
                      value = "value",
                      -one_of(feature.level)) %>%
        dplyr::left_join(meta_tab %>%
                           rownames_to_column("sample"), "sample") %>%
        dplyr::group_by(!!sym(feature.level),!!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::spread(key = time.var, value = "mean_value") %>% column_to_rownames(feature.level)
    }

    annotation_col <- meta_tab %>%
      select(!!sym(time.var), !!sym(group.var)) %>%
      as_tibble() %>%
      dplyr::distinct() %>%
      dplyr::mutate(group_time = paste(!!sym(group.var), !!sym(time.var), sep = "_")) %>%
      column_to_rownames("group_time")

    annotation_col_sorted <-
      annotation_col[order(annotation_col[[group.var]], annotation_col[[time.var]]),]

    if (!is.null(strata.var)) {
      annotation_col_sorted <- annotation_col_sorted %>%
        tidyr::separate(!!sym(group.var),
                        into = c(strata.var, group.var),
                        sep = "\\.") %>%
        dplyr::select(!!sym(time.var),!!sym(group.var),!!sym(strata.var))

      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]],
                                    annotation_col_sorted[[group.var]],
                                    annotation_col_sorted[[time.var]]),]
    }

    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted)]

    # Calculate gaps if group.var is not NULL
    if (!is.null(group.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[group.var]]))[-length(table(annotation_col_sorted[[group.var]]))]
    } else {
      gaps <- NULL
    }

    if (!is.null(strata.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[strata.var]]))[-length(table(annotation_col_sorted[[strata.var]]))]
    }

    col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")

    # 创建颜色映射函数
    my_col <- colorRampPalette(col)

    # 计算颜色的数量
    # 这通常取决于你的数据，你可能需要根据你的实际情况进行调整
    n_colors <- 100

    if (!is.null(features.plot)) {
      wide_data_sorted <-
        wide_data_sorted[rownames(wide_data_sorted) %in% features.plot,]
    }

    if (is.null(palette)) {
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
    group_levels <-
      annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
    strata_levels <-
      annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()

    # 为 group.var 分配颜色
    group_colors <-
      setNames(color_vector[1:length(group_levels)], group_levels)

    # 为 strata.var 分配颜色
    strata_colors <-
      setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

    # 创建注释颜色列表
    annotation_colors_list <- setNames(list(group_colors, strata_colors),
                                       c(group.var, strata.var))

    # Plot stacked heatmap
    heatmap_plot_average <- pheatmap::pheatmap(
      wide_data_sorted,
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col(n_colors),
      ...
    )

    gg_heatmap_plot_average <- as.ggplot(heatmap_plot_average)

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      )))

    if (!is.null(features.plot)) {
      otu_tab_norm <-
        otu_tab_norm[rownames(otu_tab_norm) %in% features.plot,]
    }

    heatmap_plot_indiv <- pheatmap::pheatmap(
      otu_tab_norm,
      annotation_col = meta_tab,
      annotation_colors = annotation_colors_list,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      gaps_col = gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col(n_colors),
      ...
    )

    gg_heatmap_plot_indiv <- as.ggplot(heatmap_plot_indiv)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_long_average",
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
        "feature_level_",
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
        plot = gg_heatmap_plot_average
      )
    }

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_long_indiv",
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
        "feature_level_",
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
        plot = gg_heatmap_plot_indiv
      )
    }

    sub.plot_list <- list(gg_heatmap_plot_indiv, gg_heatmap_plot_average)

    names(sub.plot_list) <- c("indiv","average")

    return(sub.plot_list)
  })

  names(plot_list) <- feature.level

  # Return the heatmap plot for display
  return(plot_list)
}
