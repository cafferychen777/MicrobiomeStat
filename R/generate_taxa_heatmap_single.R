#' @title Generate Taxonomic Heatmap Single
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var The name of the subject variable in the samples
#' @param time.var The name of the time variable in the samples
#' @param t.level The base level for time points in longitudinal data.
#' @param group.var The name of the grouping variable in the samples
#' @param strata.var The name of the strata variable in the samples
#' @param feature.level The taxonomic level to aggregate, can be “Phylum”, “Family” or “Genus”
#' @param features.plot A character vector specifying the taxa to be plotted. If NULL (default), the top k taxa by mean abundance will be plotted.
#' @param feature.dat.type A character string specifying the type of the data in feature.dat. Options are "count", "proportion", or "other".
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
#' @param prev.filter The prevalence filter to apply
#' @param abund.filter The abundance filter to apply
#' @param base.size Base font size for the generated plots.
#' @param palette Color palette used for the plots.
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is NULL.
#' @param cluster.rows A logical variable indicating if rows should be clustered. Default is NULL.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann The file name annotation (default: NULL)
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional arguments passed to the pheatmap() function from the “pheatmap” package.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#'
#' # Load required libraries and example data
#' library(pheatmap)
#' data(peerj32.obj)
#'
#' # Generate the boxplot pair
#' generate_taxa_heatmap_single(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = NULL,
#'   t.level = NULL,
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' @export
#'
#' @seealso \code{\link{pheatmap}}
generate_taxa_heatmap_single <- function(data.obj,
                                         subject.var,
                                         time.var = NULL,
                                         t.level = NULL,
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

  if (!is.null(time.var)) {
    if (!is.null(t.level)) {
      meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      ))) %>% filter(!!sym(time.var) == t.level)
    } else {
      meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      )))
      if (length(levels(as.factor(meta_tab[, time.var]))) != 1) {
        message <-
          "Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis."
      }
    }
  } else {
    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, group.var, strata.var)))
  }

  otu_tab <-
    otu_tab %>% as.data.frame() %>% select(all_of(c(rownames(meta_tab))))

  if (is.null(palette)) {
    palette <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
  }
  # 创建颜色映射函数
  my_palette <- colorRampPalette(palette)

  # 计算颜色的数量
  # 这通常取决于你的数据，你可能需要根据你的实际情况进行调整
  n_colors <- 100

  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
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
      tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(!!sym(feature.level))) %>%
      dplyr::summarise(total_count = mean(count),
                prevalence = sum(count > 0) / dplyr::n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count, -prevalence) %>%
      dplyr::left_join(otu_tax, by = feature.level)

    # Aggregate OTU table
    otu_tax_agg <- otu_tax_filtered %>%
      tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(sample, !!sym(feature.level))) %>%
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

    otu_tab_norm <-
      otu_tax_agg_numeric %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # Sort samples by group.var if not NULL
    if (!is.null(group.var) & !is.null(strata.var)) {
      meta_tab_sorted <-
        meta_tab[order(meta_tab[[group.var]], meta_tab[[strata.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else if (!is.null(group.var) & is.null(strata.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else {
      otu_tab_norm_sorted <- otu_tab_norm
    }

    # Calculate gaps if group.var is not NULL
    if (!is.null(group.var)) {
      if (!is.numeric(meta_tab[[group.var]])) {
        gaps <-
          cumsum(table(meta_tab_sorted[[group.var]]))[-length(table(meta_tab_sorted[[group.var]]))]
      }
    } else {
      gaps <- NULL
    }

    # Set up annotation_col based on group.var and strata.var values
    if (!is.null(group.var) & !is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(strata.var, group.var)))
    } else if (!is.null(group.var) & is.null(strata.var)) {
      annotation_col <- meta_tab %>% select(all_of(c(group.var)))
    } else {
      annotation_col <- NULL
    }

    # 创建一个空字符串作为标题的默认值
    heatmap_title <- NA

    # 如果time.var不为NULL，并且t.level存在
    if (!is.null(time.var) & !is.null(t.level)) {
      heatmap_title <- paste0("Time = ", t.level)
    }

    if (!is.null(features.plot)) {
      otu_tab_norm_sorted <-
        otu_tab_norm_sorted[rownames(otu_tab_norm_sorted) %in% features.plot,]

    }

    # Plot stacked heatmap
    heatmap_plot <- pheatmap::pheatmap(
      otu_tab_norm_sorted,
      annotation_col = annotation_col,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      gaps_col = gaps,
      color = my_palette(n_colors),
      # 使用自定义颜色
      fontsize = base.size,
      main = heatmap_title,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_single",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "t_level_",
        t.level,
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
    # Return the heatmap plot for display
    return(gg_heatmap_plot)
  })

  return(plot_list)
}
