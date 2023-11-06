#' @title Generate Taxonomic Heatmap Pair
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
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
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is FALSE.
#' @param pdf A logical value. If TRUE (default), saves the plot as a PDF file. If FALSE, the plot will be displayed interactively without creating a PDF.
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed to the pheatmap() function from the “pheatmap::pheatmap” package.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(pheatmap)
#' data(peerj32.obj)
#' generate_taxa_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Phylum","Family","Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   cluster.rows = NULL,
#'   cluster.cols = NULL,
#'   base.size = 12,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' subset_T2D.obj2 <- mStat_subset_data(subset_T2D.obj,
#' condition = "visit_number %in% c('   1', '   2')")
#' generate_taxa_heatmap_pair(
#'   data.obj = subset_T2D.obj2,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   feature.level = c("Phylum","Family","Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   cluster.rows = NULL,
#'   cluster.cols = NULL,
#'   base.size = 12,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_heatmap_pair <- function(data.obj,
                                       subject.var,
                                       time.var,
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
                                       cluster.rows = NULL,
                                       cluster.cols = NULL,
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       ...) {
  # Data validation
  mStat_validate_data(data.obj)

  feature.dat.type <- match.arg(feature.dat.type)

  meta_tab <-  data.obj$meta.dat %>% select(all_of(c(
    time.var, group.var, strata.var, subject.var
  )))

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

  if (feature.dat.type == "count"){
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
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
      otu_tax_agg <- data.obj$feature.tab
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

    otu_tab_norm <-
      otu_tax_agg %>% column_to_rownames(var = feature.level) %>% as.matrix()

    # Sort samples by group.var and strata.var
    if (!is.null(group.var) & !is.null(strata.var)) {
      meta_tab_sorted <-
        meta_tab[order(meta_tab[[strata.var]], meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else if (!is.null(group.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else {
      otu_tab_norm_sorted <- otu_tab_norm
    }

    if (!is.null(group.var)){
      # Calculate gaps if group.var is not NULL
      gaps <-
        cumsum(table(meta_tab_sorted[[subject.var]]))[-length(table(meta_tab_sorted[[subject.var]]))]
    } else {
      gaps <- NULL
    }

    # Set up annotation_col based on group.var and strata.var values
    if (!is.null(group.var) & !is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(time.var, group.var, strata.var)))
    } else if (!is.null(group.var) & is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(time.var, group.var)))
    } else {
      annotation_col <- meta_tab %>% select(all_of(c(time.var)))
    }


    if (!is.null(features.plot)) {
      otu_tab_norm_sorted <-
        otu_tab_norm_sorted[rownames(otu_tab_norm_sorted) %in% features.plot,]
    }

    if (is.null(palette)){
      color_vector <- c(
        "#1F78B4",
        "#E31A1C",
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
      group_levels <- annotation_col %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()

      # 为 group.var 分配颜色
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)

      strata_levels <- annotation_col %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      # 为 strata.var 分配颜色
      strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

      # 创建注释颜色列表
      annotation_colors_list <- setNames(
        list(group_colors, strata_colors),
        c(group.var, strata.var)
      )
    } else if (!is.null(group.var)){
      # 为演示目的，假设这些是您的唯一值
      group_levels <- annotation_col %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
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

    col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")

    # 创建颜色映射函数
    my_col <- colorRampPalette(col)

    # 计算颜色的数量
    # 这通常取决于你的数据，你可能需要根据你的实际情况进行调整
    n_colors <- 100

    if (feature.dat.type != "other"){
      quantiles <- c(0,0.1,0.5,1)
      labels <- as.character(round(quantiles^2,4))
    } else {
      quantiles <- NA
      labels <- NA
    }

    if (feature.dat.type != "other") {
      processed_otu_tab <- sqrt(otu_tab_norm_sorted)
    } else {
      processed_otu_tab <- otu_tab_norm_sorted
    }

    # Plot stacked heatmap
    heatmap_plot <- pheatmap::pheatmap(
      processed_otu_tab,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      gaps_col = gaps,
      legend_breaks = quantiles,
      legend_labels = labels,
      color = my_col(n_colors),
      fontsize = base.size,
      silent = TRUE,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    if (!is.null(strata.var) & !is.null(group.var) & !is.null(time.var)){
      wide_data <- otu_tab_norm_sorted %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(cols = -feature.level, names_to = "sample", values_to = "value") %>%
        dplyr::left_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level), !!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("column_name", !!sym(time.var), !!sym(group.var), !!sym(strata.var), sep = "_") %>%
        tidyr::pivot_wider(names_from = column_name, values_from = mean_value) %>%
        dplyr::ungroup()

      annotation_col <- wide_data %>%
        dplyr::select(-all_of(c(feature.level))) %>%
        names() %>%
        tibble(column_name = .) %>%
        dplyr::mutate(
          !!time.var := stringr::str_extract(column_name, "^[^_]+"),
          !!group.var := stringr::str_extract(column_name, "(?<=_).*(?=_)"),
          !!strata.var := stringr::str_extract(column_name, "[^_]+$")
        ) %>%
        column_to_rownames("column_name")

    } else if (!is.null(group.var) & !is.null(time.var)){
      wide_data <- otu_tab_norm_sorted %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(cols = -feature.level, names_to = "sample", values_to = "value") %>%
        dplyr::left_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level), !!sym(group.var), !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("column_name", !!sym(time.var), !!sym(group.var), sep = "_") %>%
        tidyr::pivot_wider(names_from = column_name, values_from = mean_value) %>%
        dplyr::ungroup()

      annotation_col <- wide_data %>%
        dplyr::select(-all_of(c(feature.level))) %>%
        names() %>%
        tibble(column_name = .) %>%
        dplyr::mutate(
          !!time.var := stringr::str_extract(column_name, "^[^_]+"),
          !!group.var := stringr::str_extract(column_name, "(?<=_).*")
        ) %>%
        tidyr::unite("column_name", !!sym(time.var), !!sym(group.var), sep = "_", remove = FALSE) %>%
        column_to_rownames("column_name")

    } else if (!is.null(time.var)){
      wide_data <- otu_tab_norm_sorted %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(cols = -feature.level, names_to = "sample", values_to = "value") %>%
        dplyr::left_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level), !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("column_name", !!sym(time.var), sep = "_") %>%
        tidyr::pivot_wider(names_from = column_name, values_from = mean_value) %>%
        dplyr::ungroup()

      annotation_col <- wide_data %>%
        dplyr::select(-all_of(c(feature.level))) %>%
        names() %>%
        stringr::str_split("_", simplify = TRUE) %>%
        as.data.frame() %>%
        setNames(c(time.var)) %>%
        distinct() %>%
        tidyr::unite("column_name", !!sym(time.var), sep = "_", remove = FALSE) %>%
        column_to_rownames("column_name")
    }

    # Check the feature.dat.type and apply sqrt if not equal to "other"
    if (feature.dat.type != "other") {
      processed_wide_data <- wide_data %>% column_to_rownames(feature.level) %>% sqrt()
    } else {
      processed_wide_data <- wide_data %>% column_to_rownames(feature.level)
    }

    # Plot stacked heatmap
    average_heatmap_plot <- pheatmap::pheatmap(
      processed_wide_data,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = NULL,
      legend_breaks = quantiles,
      legend_labels = labels,
      color = my_col(n_colors),
      fontsize = base.size,
      silent = TRUE,
      ...
    )

    gg_average_heatmap_plot <- as.ggplot(average_heatmap_plot)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_pair_average",
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
        plot = gg_average_heatmap_plot
      )
    }

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_pair_indiv",
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
        plot = gg_heatmap_plot
      )
    }

    sub.plot_list <- list(gg_heatmap_plot, gg_average_heatmap_plot)

    names(sub.plot_list) <- c("indiv", "average")
    # Return the heatmap plot for display
    return(sub.plot_list)
  })

  names(plot_list) <- feature.level
  return(plot_list)
}
