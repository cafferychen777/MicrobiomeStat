#' Generate Taxonomic Heatmap Long
#'
#' This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var The name of the subject variable in the samples
#' @param time.var The name of the time variable in the samples
#' @param group.var The name of the grouping variable in the samples
#' @param strata.var The name of the strata variable in the samples
#' @param feature.level The taxonomic level to aggregate, can be “Phylum”, “Family” or “Genus”
#' @param prev.filter The prevalence filter to apply
#' @param abund.filter The abundance filter to apply
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann The file name annotation (default: NULL)
#' @param ... Additional arguments to be passed to pheatmap (default: NULL)
#'
#' @examples
#' ecam.obj$meta_tab$new_month <- factor(ecam.obj$meta.dat$new_month, levels =
#' c("Month 0",paste("Month",as.character(sort(as.numeric
#' (as.character(unique(ecam.obj$meta.dat$month)[-1])))))))
#' plot_list <- generate_taxa_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Family",
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' @return An object of class pheatmap, the generated heatmap plot
#' @export
#'
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
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  meta_tab <- load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

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

  if (is.null(group.var)){
    group.var = "ALL"
    meta_tab$ALL <- "ALL"
  }

  if (!is.null(strata.var)){
    meta_tab <- meta_tab %>% mutate(!!sym(group.var) := interaction(!!sym(strata.var),!!sym(group.var)))
  }

  if (is.null(cluster.cols)){
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)){
    cluster.rows = TRUE
  }

  # Merge OTU table with taxonomy table
  otu_tax <-
    cbind(otu_tab, tax_tab %>% select(all_of(feature.level)))

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  plot_list <- lapply(feature.level, function(feature.level) {
  # Filter taxa based on prevalence and abundance
  otu_tax_filtered <- otu_tax %>%
    gather(key = "sample", value = "count", -one_of(feature.level)) %>%
    group_by_at(vars(!!sym(feature.level))) %>%
    summarise(total_count = mean(count),
              prevalence = sum(count > 0) / n()) %>%
    filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
    select(-total_count, -prevalence) %>%
    left_join(otu_tax, by = feature.level)

  # Aggregate OTU table
  otu_tax_agg <- otu_tax_filtered %>%
    gather(key = "sample", value = "count", -one_of(feature.level)) %>%
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

  otu_tab_norm <-
    otu_tax_agg_numeric %>% filter(!is.na(!!sym(feature.level))) %>% column_to_rownames(var = feature.level) %>% as.matrix()

  if (!is.null(group.var)) {
    wide_data <- otu_tab_norm %>%
      as.data.frame() %>%
      rownames_to_column(var = feature.level) %>%
      gather(key = "sample", value = "value", -one_of(feature.level)) %>%
      left_join(meta_tab %>%
                  rownames_to_column("sample"), "sample") %>%
      group_by(!!sym(feature.level), !!sym(group.var), !!sym(time.var)) %>%
      summarise(mean_value = mean(value)) %>%
      unite("group_time", c(group.var, time.var), sep = "_") %>%
      spread(key = "group_time", value = "mean_value") %>% column_to_rownames(feature.level)
  } else {
    wide_data <- otu_tab_norm %>%
      as.data.frame() %>%
      rownames_to_column(var = feature.level) %>%
      gather(key = "sample", value = "value", -one_of(feature.level)) %>%
      left_join(meta_tab %>%
                  rownames_to_column("sample"), "sample") %>%
      group_by(!!sym(feature.level), !!sym(time.var)) %>%
      summarise(mean_value = mean(value)) %>%
      spread(key = time.var, value = "mean_value") %>% column_to_rownames(feature.level)
  }

    annotation_col <- meta_tab %>%
      select(!!sym(time.var),!!sym(group.var)) %>%
      as_tibble() %>%
      distinct() %>%
      mutate(group_time = paste(!!sym(group.var),!!sym(time.var), sep = "_")) %>%
      column_to_rownames("group_time")
    annotation_col_sorted <-
      annotation_col[order(annotation_col[[group.var]], annotation_col[[time.var]]), ]
    if (!is.null(strata.var)){
      annotation_col_sorted <- annotation_col_sorted %>%
        separate(!!sym(group.var), into = c(strata.var, group.var), sep = "\\.")
    }
    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted)]


  # Calculate gaps if group.var is not NULL
  if (!is.null(group.var)) {
    gaps <-
      cumsum(table(annotation_col_sorted[[group.var]]))[-length(table(annotation_col_sorted[[group.var]]))]
  } else {
    gaps <- NULL
  }

  if (is.null(palette)) {
    palette <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
  }

  # 创建颜色映射函数
  my_palette <- colorRampPalette(palette)

  # 计算颜色的数量
  # 这通常取决于你的数据，你可能需要根据你的实际情况进行调整
  n_colors <- 100

  if (!is.null(features.plot)){

    wide_data_sorted <- wide_data_sorted[rownames(wide_data_sorted) %in% features.plot, ]

  }

  # Plot stacked heatmap
  heatmap_plot <- pheatmap(
    wide_data_sorted,
    annotation_col = annotation_col_sorted,
    annotation_colors = NULL,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    gaps_col = gaps,
    fontsize = base.size,
    color = my_palette(n_colors)  # 使用自定义颜色
  )

  gg_heatmap_plot <- as.ggplot(heatmap_plot)

  # Save the stacked heatmap as a PDF file
  if (pdf) {
    pdf_name <- paste0(
      "taxa_heatmap_long",
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
  return(gg_heatmap_plot)
  })

  # Return the heatmap plot for display
  return(plot_list)
}
