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
#' @param feature.change.func A method or function specifying how to compute the change in feature abundance or prevalence between two time points.
#' The following options are available:
#'
#' - A custom function: If you provide a user-defined function, it should take two numeric arguments corresponding to the values at the two time points and return the computed change. This function can be applied to compute changes both in abundance (`time1_mean_abundance` and `time2_mean_abundance`) and prevalence (`time1_prevalence` and `time2_prevalence`).
#'
#' - "log fold change": Computes the log2 fold change between the two time points. To handle zeros, a small offset (0.00001) is added before taking the logarithm. This method can be applied for both abundance and prevalence changes.
#'
#' - "relative change": Computes the relative change as `(time2_value - time1_value) / (time2_value + time1_value)`. If both time points have a value of 0, the change is defined as 0. This method can be applied for both abundance and prevalence changes.
#'
#' - "absolute change": Computes the difference between the values at the two time points. This method can be applied for both abundance and prevalence changes.
#'
#' - Any other value (or if the parameter is omitted): By default, the function will compute the absolute change as described above, regardless of whether it is abundance or prevalence data.
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
#' @param top.k.plot Integer specifying number of top abundant features to plot, when `features.plot` is NULL.
#' Default is NULL, in which case all features passing filters will be plotted.
#' @param top.k.func Function to use for selecting top abundant features, when `features.plot` is NULL.
#' Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
#' case features will be selected by mean abundance.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size Base font size for the generated plots.
#' @param theme.choice Plot theme choice. Can be one of:
#'   - "prism": ggprism::theme_prism()
#'   - "classic": theme_classic()
#'   - "gray": theme_gray()
#'   - "bw": theme_bw()
#' Default is "bw".
#' @param custom.theme A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. To use a custom theme, first create a theme object with ggplot2::theme(), then pass it to this argument. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=16, color="red"),
#'   legend.position = "none"
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. Default is NULL, which will use the default theme based on `theme.choice`.
#' @param palette Color palette used for the plots.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#' @examples
#' \dontrun{
#'
#' # Note: In the RStudio viewer, the plot might appear cluttered if there are many taxa.
#' # It's recommended to view the generated PDF for better clarity. If it still feels
#' # overcrowded in the PDF, consider increasing the 'pdf.wid' value to adjust the width of the plot.
#'
#' data(peerj32.obj)
#' generate_taxa_change_dotplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = "Family",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 20,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 1e-4,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 10
#' )
#'
#' data("subset_pairs.obj")
#' generate_taxa_change_dotplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = "Family",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 20,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 1e-4,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 10
#' )
#' }
#' # View the result
#' @export
generate_taxa_change_dotplot_pair <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              change.base = "1",
                                              feature.change.func = "log fold change",
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

  meta_tab <-
    data.obj$meta.dat %>% select(all_of(c(
      time.var, group.var, strata.var, subject.var
    )))

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
    colors <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
  } else {
    colors <- palette
  }

  # Assuming mStat_get_theme function is already defined
  # Replace the existing theme selection code with this:
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  plot_list <- lapply(feature.level, function(feature.level) {
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <-
        mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm
    }

    if (is.null(data.obj$feature.agg.list[[feature.level]]) &
        feature.level != "original") {
      data.obj <-
        mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    if (feature.level != "original") {
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- data.obj$feature.tab
    }

    otu_tax_agg <-  otu_tax_agg %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter) %>%
      tibble::rownames_to_column(feature.level)

    if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
    }

    # 计算每个分组的平均丰度
    otu_tab_norm_agg <- otu_tax_agg %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level),!!sym(time.var)) %>% # Add time.var to dplyr::group_by
      dplyr::summarise(mean_abundance = mean(count))

    change.after <-
      unique(meta_tab %>% select(all_of(c(time.var))))[unique(meta_tab %>% select(all_of(c(time.var)))) != change.base]

    # 将数据从长格式转换为宽格式，将不同的时间点的mean_abundance放到不同的列中
    otu_tab_norm_agg_wide <- otu_tab_norm_agg %>%
      tidyr::spread(key = !!sym(time.var), value = mean_abundance) %>%
      dplyr::rename(
        time1_mean_abundance = all_of(change.base),
        time2_mean_abundance = all_of(change.after)
      )

    if (is.function(feature.change.func)) {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% dplyr::mutate(abundance_change = feature.change.func(time2_mean_abundance, time2_mean_abundance))
    } else if (feature.change.func == "log fold change") {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% dplyr::mutate(
          abundance_change = log2(time2_mean_abundance + 0.00001) - log2(time1_mean_abundance + 0.00001)
        )
    } else if (feature.change.func == "relative change") {
      otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
        dplyr::mutate(
          abundance_change = dplyr::case_when(
            time2_mean_abundance == 0 & time1_mean_abundance == 0 ~ 0,
            TRUE ~ (time2_mean_abundance - time1_mean_abundance) / (time2_mean_abundance + time1_mean_abundance)
          )
        )
    } else if (feature.change.func == "absolute change") {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% dplyr::mutate(abundance_change = time2_mean_abundance - time1_mean_abundance)
    } else {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% dplyr::mutate(abundance_change = time2_mean_abundance - time1_mean_abundance)
    }

    # 计算每个taxon在每个时间点的prevalence
    prevalence_time <- otu_tax_agg %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level),!!sym(time.var)) %>%
      dplyr::summarise(prevalence = sum(count > 0) / dplyr::n())

    prevalence_time_wide <- prevalence_time %>%
      tidyr::spread(key = time.var, value = prevalence) %>%
      dplyr::rename(time1_prevalence = change.base,
                    time2_prevalence = change.after)

    # 计算不同时间点的prevalence差值
    if (is.function(feature.change.func)) {
      prevalence_time_wide <-
        prevalence_time_wide %>% dplyr::mutate(prevalence_change = feature.change.func(time2_prevalence, time1_prevalence))
    } else if (feature.change.func == "log fold change") {
      prevalence_time_wide <-
        prevalence_time_wide %>% dplyr::mutate(
          prevalence_change = log2(time2_prevalence + 0.00001) - log2(time1_prevalence + 0.00001)
        )
    } else if (feature.change.func == "relative change") {
      prevalence_time_wide <- prevalence_time_wide %>%
        dplyr::mutate(
          prevalence_change = dplyr::case_when(
            time2_prevalence == 0 & time1_prevalence == 0 ~ 0,
            TRUE ~ (time2_prevalence - time1_prevalence) / (time2_prevalence + time1_prevalence)
          )
        )
    } else if (feature.change.func == "absolute change") {
      prevalence_time_wide <-
        prevalence_time_wide %>% dplyr::mutate(prevalence_change = time2_prevalence - time1_prevalence)
    } else {
      prevalence_time_wide <-
        prevalence_time_wide %>% dplyr::mutate(prevalence_change = time2_prevalence - time1_prevalence)
    }

    # 将两个结果合并
    otu_tab_norm_agg_wide <-
      otu_tab_norm_agg_wide %>% dplyr::left_join(prevalence_time_wide, by = c(feature.level, group.var))

    if (!is.null(strata.var)) {
      otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
        dplyr::mutate(temp = !!sym(group.var)) %>%
        tidyr::separate(temp,
                        into = c(paste0(group.var, "2"), strata.var),
                        sep = "\\.")
    }

    if (!is.null(features.plot)) {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% filter(!!sym(feature.level) %in% features.plot)
    }

    prop_prev_data <-
      otu_tax_agg %>%
      column_to_rownames(feature.level) %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::group_by(Var1) %>%  # Var1是taxa
      dplyr::summarise(avg_abundance = mean(Freq),
                       prevalence = sum(Freq > 0) / dplyr::n()) %>% column_to_rownames("Var1") %>%
      rownames_to_column(feature.level)

    otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
      tidyr::gather(key = "Type",
                    value = "change",
                    abundance_change,
                    prevalence_change) %>%
      select(-all_of(
        c(
          "time1_mean_abundance",
          "time2_mean_abundance",
          "time1_prevalence",
          "time2_prevalence"
        )
      )) %>%
      dplyr::left_join(prop_prev_data, by = feature.level) %>%
      dplyr::mutate(
        Base = dplyr::case_when(
          Type == "abundance_change" ~ avg_abundance,
          Type == "prevalence_change" ~ prevalence,
          TRUE ~ NA_real_
        )
      ) %>%
      select(-all_of(c("avg_abundance", "prevalence")))

    adjust_size_range <- function(taxa.levels) {
      if (taxa.levels <= 2) {
        return(c(40, 57))
      } else if (taxa.levels <= 4) {
        return(c(35, 42))
      } else if (taxa.levels <= 6) {
        return(c(30, 37))
      } else if (taxa.levels <= 8) {
        return(c(25, 33))
      } else if (taxa.levels < 10) {
        return(c(20, 17))
      } else if (taxa.levels < 20) {
        return(c(10, 15))
      } else if (taxa.levels < 30) {
        return(c(8, 13))
      } else if (taxa.levels < 40) {
        return(c(6, 10))
      } else if (taxa.levels < 50) {
        return(c(4, 8))
      } else {
        return(c(1, 4))
      }
    }

    # 找到 abundance_change 的最小值和最大值
    change_min <- min(otu_tab_norm_agg_wide$change)
    change_max <- max(otu_tab_norm_agg_wide$change)

    # 归一化的函数
    normalize <- function(x, min, max) {
      return((x - min) / (max - min))
    }

    # 归一化的值
    change_min_norm <-
      normalize(change_min, change_min, change_max)
    change_max_norm <-
      normalize(change_max, change_min, change_max)
    change_mid_norm <-
      normalize(0, change_min, change_max)  # abundance_change 的中点为0

    # 计算其他颜色的归一化值
    first_color_norm <-
      change_min_norm + (change_mid_norm - change_min_norm) / 2
    second_color_norm <-
      change_mid_norm + (change_max_norm - change_mid_norm) / 2

    taxa.levels <-
      otu_tab_norm_agg_wide %>% dplyr::ungroup() %>% select(all_of(c(feature.level))) %>% pull() %>% unique() %>% length()

    # 将患病率添加为点的大小，并将平均丰度作为点的颜色
    dotplot <-
      ggplot(
        otu_tab_norm_agg_wide,
        aes(
          x = !!sym(feature.level),
          y = !!sym(group.var),
          size = Base,
          shape = Type,
          color = Type
        )
      ) + # Change x to time.var
      geom_point(
        aes(group = interaction(Type, !!sym(feature.level)), fill = change),
        shape = 21,
        position = position_dodge(0.9)
      ) +
      xlab(feature.level) +
      ylab(group.var) +
      scale_colour_manual(values = c("transparent", "black")) +
      scale_size_continuous(range = adjust_size_range(taxa.levels)) +
      scale_fill_gradientn(
        colors = colors,
        values = c(
          change_min_norm,
          first_color_norm,
          change_mid_norm,
          second_color_norm,
          change_max_norm
        ),
        name = "Change"
      ) +
      {
        if (!is.null(strata.var)) {
          ggh4x::facet_nested(
            rows = vars(!!sym(strata.var),!!sym(paste0(
              group.var, "2"
            ))),
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
          element_blank()
        else
          element_text(size = base.size),
        panel.spacing = unit(0, "lines"),
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        legend.text = ggplot2::element_text(size = base.size),
        legend.title = ggplot2::element_text(size = base.size),
      ) + guides(color = guide_legend(override.aes = list(
        size = 5, fill = "#92c5de"
      )),
      shape = guide_legend(override.aes = list(size = 5)))

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

  names(plot_list) <- feature.level
  return(plot_list)
}
