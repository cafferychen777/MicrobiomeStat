#' @title Generate Taxonomic Dotplot Single
#'
#' @description This function generates a stacked dotplot of specified taxa level with paired samples. The data used in this
#' visualization will be first filtered based on prevalence and abundance thresholds. The plot can either be displayed
#' interactively or saved as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string defining subject variable in meta_tab
#' @param time.var A character string defining time variable in meta_tab
#' @param t.level The base level for time points in longitudinal data.
#' @param group.var A character string defining group variable in meta_tab used for sorting and facetting
#' @param strata.var (Optional) A character string defining strata variable in meta_tab used for sorting and facetting
#' @param feature.level A character string defining the taxonomic level to analyze ('Phylum', 'Family' or 'Genus')
#' @param features.plot A character vector specifying the taxa to be plotted. If NULL (default), the top k taxa by mean abundance will be plotted.
#' @param feature.dat.type A character string specifying the type of data in feature.dat. Options are "count", "proportion", or "other".
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
#' @param prev.filter A numeric value defining the prevalence threshold to filter taxa, between 0 and 1.
#' @param abund.filter A numeric value defining the abundance threshold to filter taxa.
#' @param base.size Base font size for the generated plots.
#' @param theme.choice Plot theme choice (default: "bw").
#' @param custom.theme Custom ggplot2 theme (optional).
#' @param palette Color palette used for the plots.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann The file name annotation (default: NULL)
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' # Call the function
#' generate_taxa_dotplot_single(
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
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 15,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_dotplot_single <- function(data.obj,
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
                                       base.size = 16,
                                       theme.choice = "bw",
                                       custom.theme = NULL,
                                       palette = NULL,
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       ...) {

  feature.dat.type <- match.arg(feature.dat.type)

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
    {if("original" %in% feature.level) dplyr::mutate(., original = rownames(.)) else .} %>%
    select(all_of(feature.level))

  if (!is.null(time.var)){
    if (!is.null(t.level)){
      meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(
        c(subject.var, time.var, group.var, strata.var))) %>% filter(!!sym(time.var) == t.level)
    } else {
      meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(
        c(subject.var, time.var, group.var, strata.var)))
    }
  } else {
    meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(
      c(subject.var, group.var, strata.var)))
  }

  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- ""
  }

  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var), !!sym(strata.var)))
  }

  # Define the colors
  if (is.null(palette)) {
    colors <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
  } else {
    colors <- palette
  }

  theme_function <- switch(
    theme.choice,
    prism = ggprism::theme_prism(),
    classic = theme_classic(),
    gray = theme_gray(),
    bw = theme_bw(),
    ggprism::theme_prism()
  ) # 根据用户选择设置主题

  # 使用用户自定义主题（如果提供），否则使用默认主题
  theme_to_use <-
    if (!is.null(custom.theme))
      custom.theme else theme_function

  # 将 OTU 表与分类表合并
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

    # 聚合 OTU 表
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

    # 转换计数为数值类型
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    otu_tab_norm <- otu_tax_agg_numeric %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # 计算每个分组的平均丰度和患病率
    otu_tab_norm_agg <- otu_tax_agg_numeric %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level)) %>%
      dplyr::summarise(prevalence = sum(count > 0) / dplyr::n(),
                mean_abundance = mean(count)) %>% dplyr::ungroup()

    # 计算每个分组的平均丰度
    otu_tab_norm_agg <- otu_tax_agg_numeric %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level)) %>% # Add time.var to dplyr::group_by
      dplyr::summarise(mean_abundance = sqrt(mean(count)))

    # 计算所有样本中的prevalence
    prevalence_all <- otu_tax_agg_numeric %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::group_by(!!sym(feature.level)) %>%
      dplyr::summarise(prevalence = sum(count > 0) / dplyr::n())

    # 将两个结果合并
    otu_tab_norm_agg <-
      otu_tab_norm_agg %>% dplyr::left_join(prevalence_all, feature.level)

    # Calculate the midpoint from data
    midpoint <- quantile(otu_tab_norm_agg$mean_abundance, 0.5)

    if (!is.null(strata.var)){
      otu_tab_norm_agg <- otu_tab_norm_agg %>%
        dplyr::mutate(temp = !!sym(group.var)) %>%
        tidyr::separate(temp, into = c(paste0(group.var,"2"), strata.var), sep = "\\.")
    }

    if (!is.null(features.plot)){

      otu_tab_norm_agg <- otu_tab_norm_agg %>% filter(!!sym(feature.level) %in% features.plot)

    }

    # 将患病率添加为点的大小，并将平均丰度作为点的颜色
    dotplot <-
      ggplot(
        otu_tab_norm_agg,
        aes(
          x = !!sym(feature.level),
          y = !!sym(group.var),
          color = mean_abundance,
          size = prevalence
        )
      ) +
      geom_point(aes(group = !!sym(feature.level), fill = mean_abundance), shape = 21, color = "black", position = position_dodge(0.9)) +
      xlab(feature.level) + # Change x-label to "Time"
      ylab(group.var) +
      {
        if(feature.dat.type == "other") {
          quantiles <- quantile(otu_tab_norm_agg$mean_abundance, probs = c(0, 0.25, 0.5, 0.75, 1))
          scale_fill_gradientn(colors = colors,
                               values = scales::rescale(quantiles),
                               name = "Mean Abundance")
        } else {
          scale_fill_gradientn(colors = colors,
                               values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                               name = "Mean Abundance (Sqrt)")
        }
      } +
      scale_size(range = c(4, 10), name = "Prevalence") +
      scale_shape_manual(values = c(19, 1)) +
      {
          if (!is.null(strata.var)){
            ggh4x::facet_nested(rows = vars(!!sym(paste0(group.var,"2")),!!sym(strata.var)), cols = vars(!!sym(feature.level)), scales = "free", switch = "y")
          } else {
            ggh4x::facet_nested(rows = vars(!!sym(group.var)), cols = vars(!!sym(feature.level)), scales = "free", switch = "y")}
      } +
      theme_to_use +
      theme(
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = base.size
        ),
        strip.text.x = element_blank(),
        strip.text.y = if (group.var == "ALL") element_blank() else element_text(size = base.size),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = base.size),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        # 添加主要网格线
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        # 添加次要网格线
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16)
      )

    # Save the stacked dotplot as a PDF file
    if (pdf) {
      dotplot <- as.ggplot(dotplot)
      pdf_name <- paste0(
        "taxa_dotplot_single",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
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
        plot = dotplot,
        width = pdf.wid,
        height = pdf.hei
      )
    }

    return(dotplot)
  })
  return(plot_list)
}
