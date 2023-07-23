#' @title Generate Taxa Stack Dotplot Pair
#'
#' @description This function generates a stacked dotplot of specified taxa level with paired samples. The data used in this
#' visualization will be first filtered based on prevalence and abundance thresholds. The plot can either be displayed
#' interactively or saved as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string defining subject variable in meta_tab
#' @param time.var A character string defining time variable in meta_tab
#' @param group.var A character string defining group variable in meta_tab used for sorting and facetting
#' @param strata.var (Optional) A character string defining strata variable in meta_tab used for sorting and facetting
#' @param change.base A numeric value setting base for the change (usually 1)
#' @param change.func A function or character string specifying the method for computing the change. Default is "difference".
#' @param feature.level A character string defining the taxonomic level to analyze ('Phylum', 'Family', or 'Genus')
#' @param features.plot A character vector specifying the taxa to be plotted. If NULL (default), the top k taxa by mean abundance will be plotted.
#' @param feature.dat.type A character string specifying the type of the data in feature.dat. Options are "count", "proportion", or "other".
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
#' @param prev.filter A numeric value defining the prevalence threshold to filter taxa, between 0 and 1
#' @param abund.filter A numeric value defining the abundance threshold to filter taxa
#' @param base.size Base font size for the generated plots.
#' @param theme.choice The choice of theme for the plot. Default is "bw".
#' @param custom.theme Custom ggplot2 theme (optional).
#' @param palette Color palette used for the plots.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#' @examples
#' # Load required libraries
#' library(tidyverse)
#' library(ggh4x)
#'
#' # Prepare data for the function
#' data(peerj32.obj)
#'
#' # Call the function
#' dotplot_result <- generate_taxa_change_dotplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   change.base = "1",
#'   change.func = "lfc",
#'   feature.level = "Family",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' # View the result
#' print(dotplot_result)
#' @export
generate_taxa_change_dotplot_pair <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              change.base = "1",
                                              change.func = "lfc",
                                              feature.level = NULL,
                                              feature.dat.type = c("count", "proportion", "other"),
                                              features.plot = NULL,
                                              top.k.plot = NULL,
                                              top.k.func = NULL,
                                              prev.filter = 0.001,
                                              abund.filter = 0.001,
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

  # Extract data
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
    {
      if ("original" %in% feature.level)
        mutate(., original = rownames(.))
      else
        .
    } %>%
    select(all_of(feature.level))

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% select(all_of(c(
      time.var, group.var, strata.var, subject.var
    )))

  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- ""
  }

  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% mutate(!!sym(group.var) := interaction(!!sym(group.var),!!sym(strata.var)))
  }

  # Define the colors
  if (is.null(palette)) {
    colors <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
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
      custom.theme else
    theme_function

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
      gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      group_by_at(vars(!!sym(feature.level))) %>%
      dplyr::summarise(total_count = mean(count),
                prevalence = sum(count > 0) / dplyr::n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count, -prevalence) %>%
      left_join(otu_tax, by = feature.level)

    # 聚合 OTU 表
    otu_tax_agg <- otu_tax_filtered %>%
      gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      group_by_at(vars(sample, !!sym(feature.level))) %>%
      dplyr::summarise(count = sum(count)) %>%
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

    # 转换计数为数值类型
    otu_tax_agg_numeric <-
      mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    otu_tab_norm <- otu_tax_agg_numeric %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # 计算每个分组的平均丰度
    otu_tab_norm_agg <- otu_tax_agg_numeric %>%
      gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      group_by(!!sym(group.var), !!sym(feature.level), !!sym(time.var)) %>% # Add time.var to group_by
      summarise(mean_abundance = mean(count))

    change.after <-
      unique(meta_tab %>% select(all_of(c(time.var))))[unique(meta_tab %>% select(all_of(c(time.var)))) != change.base]

    # 将数据从长格式转换为宽格式，将不同的时间点的mean_abundance放到不同的列中
    otu_tab_norm_agg_wide <- otu_tab_norm_agg %>%
      spread(key = !!sym(time.var), value = mean_abundance) %>%
      dplyr::rename(time1_mean_abundance = change.base,
             time2_mean_abundance = change.after)

    # 计算不同时间点的mean_abundance差值
    # 最后，计算新的count值
    if (is.function(change.func)) {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% mutate(abundance_change = change.func(time2_mean_abundance, time2_mean_abundance))
    } else if (change.func == "lfc") {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% mutate(
          abundance_change = log2(time2_mean_abundance + 0.00001) - log2(time1_mean_abundance + 0.00001)
        )
    } else if (change.func == "relative difference") {
      otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
        mutate(
          abundance_change = case_when(
            time2_mean_abundance == 0 & time1_mean_abundance == 0 ~ 0,
            TRUE ~ (time2_mean_abundance - time1_mean_abundance) / (time2_mean_abundance + time1_mean_abundance)
          )
        )
    } else {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% mutate(abundance_change = time2_mean_abundance - time1_mean_abundance)
    }

    # 计算每个taxon在每个时间点的prevalence
    prevalence_time <- otu_tax_agg_numeric %>%
      gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      group_by(!!sym(feature.level), !!sym(time.var)) %>%
      summarise(prevalence = sum(count > 0) / n())

    # 将数据从长格式转换为宽格式，将不同的时间点的prevalence放到不同的列中
    prevalence_time_wide <- prevalence_time %>%
      spread(key = time.var, value = prevalence) %>%
      dplyr::rename(time1_prevalence = change.base,
             time2_prevalence = change.after)

    # 计算不同时间点的prevalence差值
    if (is.function(change.func)) {
      prevalence_time_wide <-
        prevalence_time_wide %>% mutate(prevalence_change = change.func(time2_prevalence, time1_prevalence))
    } else if (change.func == "lfc") {
      prevalence_time_wide <-
        prevalence_time_wide %>% mutate(
          prevalence_change = log2(time2_prevalence + 0.00001) - log2(time1_prevalence + 0.00001)
        )
    } else if (change.func == "relative difference") {
      prevalence_time_wide <- prevalence_time_wide %>%
        mutate(prevalence_change = case_when(
          time2_prevalence == 0 & time1_prevalence == 0 ~ 0,
          TRUE ~ (time2_prevalence - time1_prevalence) / (time2_prevalence + time1_prevalence)
        ))
    } else {
      prevalence_time_wide <-
        prevalence_time_wide %>% mutate(prevalence_change = time2_prevalence - time1_prevalence)
    }

    # 将两个结果合并
    otu_tab_norm_agg_wide <-
      otu_tab_norm_agg_wide %>% left_join(prevalence_time_wide, feature.level)

    # 找到 abundance_change 的最小值和最大值
    abund_change_min <- min(otu_tab_norm_agg_wide$abundance_change)
    abund_change_max <- max(otu_tab_norm_agg_wide$abundance_change)

    # 归一化的函数
    normalize <- function(x, min, max) {
      return((x - min) / (max - min))
    }

    # 归一化的值
    abund_change_min_norm <-
      normalize(abund_change_min, abund_change_min, abund_change_max)
    abund_change_max_norm <-
      normalize(abund_change_max, abund_change_min, abund_change_max)
    abund_change_mid_norm <-
      normalize(0, abund_change_min, abund_change_max)  # abundance_change 的中点为0

    # 计算其他颜色的归一化值
    first_color_norm <-
      abund_change_min_norm + (abund_change_mid_norm - abund_change_min_norm) / 2
    second_color_norm <-
      abund_change_mid_norm + (abund_change_max_norm - abund_change_mid_norm) / 2

    if (!is.null(strata.var)) {
      otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
        mutate(temp = !!sym(group.var)) %>%
        separate(temp,
                 into = c(paste0(group.var, "2"), strata.var),
                 sep = "\\.")
    }

    if (!is.null(features.plot)){
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% filter(!!sym(feature.level) %in% features.plot)
    }

    # 将患病率添加为点的大小，并将平均丰度作为点的颜色
    dotplot <-
      ggplot(
        otu_tab_norm_agg_wide,
        aes(
          x = !!sym(feature.level),
          y = !!sym(group.var),
          color = abundance_change,
          size = prevalence_change
        )
      ) + # Change x to time.var
      geom_point(aes(group = !!sym(feature.level), fill = abundance_change), shape = 21) +
      xlab(feature.level) + # Change x-label to "Time"
      ylab(group.var) +
      scale_fill_gradientn(
        colors = colors,
        values = c(
          abund_change_min_norm,
          first_color_norm,
          abund_change_mid_norm,
          second_color_norm,
          abund_change_max_norm
        ),
        name = "Abundance Change"
      ) +
      scale_size(range = c(4, 10), name = "Prevalence Change") +
      {
        if (!is.null(strata.var)) {
          ggh4x::facet_nested(
            rows = vars(!!sym(paste0(
              group.var, "2"
            )), !!sym(strata.var)),
            cols = vars(!!sym(feature.level)),
            scales = "free",
            switch = "y"
          )
        } else {
          ggh4x::facet_nested(
            rows = vars(!!sym(group.var)),
            cols = vars(!!sym(feature.level)),
            scales = "free",
            switch = "y"
          )
        }
      } +
      theme_to_use +
      theme(
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = base.size
        ),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = base.size),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        strip.text.x = element_blank(),
        strip.text.y = if (group.var == "ALL")
          element_blank() else
          element_text(size = base.size),
        panel.spacing = unit(0, "lines"),
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        # 添加主要网格线
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        # 添加次要网格线
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16)
      ) +
      guides(color = FALSE) # Add this line to remove the color legend

    # Save the stacked dotplot as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_dotplot_pair",
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
