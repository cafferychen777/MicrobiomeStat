#' @title Generate longitudinal spaghettiplot plots of taxonomic composition
#'
#' @description This function generates a spaghettiplot plot for the taxonomic composition of microbiome data over time.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string defining the subject variable in meta_tab.
#'                   This should be a column that contains unique subject IDs to
#'                   identify each sample's subject of origin. Required for plotting
#'                   individual subject lines.
#' @param time.var A character string defining the time variable in meta_tab.
#'                This should be a column in meta_tab that represents the time
#'                points for the longitudinal samples. Required to identify
#'                the time axis for plotting.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var A character string defining the grouping variable in meta_tab.
#'                 This variable is used for coloring lines in the plot and faceting by
#'                 group. For example, variable could be treatment group. Optional,
#'                 can be set to NULL.
#' @param strata.var (Optional) A character string defining the stratification variable
#'                    in meta_tab. This variable is used for nested faceting in the plot.
#'                    For example, could be clinical strata like disease status.
#'                    Optional, can be set to NULL. When set, plot will show nested
#'                    facetting by strata within each feature.
#' @param feature.level The column name(s) in the feature annotation matrix (feature.ann) of data.obj
#'                     to use for summarization and plotting. Can be taxonomic levels like "Phylum",
#'                     "Genus", or any other annotation columns like "OTU_ID". Should be a character
#'                     vector specifying one or more column names in feature.ann. The special value
#'                     "original" can also be provided, which will use the original row names from
#'                     feature.ann. Multiple columns can be provided, and data will be plotted
#'                     separately for each column. **\bold{Cannot be NULL, as NULL value will lead to errors.}**
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
#' @param palette Color palette to use for mapping groups to colors in the plot. Must be a
#'               character vector of colors. Default palette will be used if NULL.
#'               If specified, should be same length as number of groups/levels in
#'               group.var. Colors will be mapped to groups in order of levels.
#' @param pdf A logical value. If TRUE (default), saves the plot as a PDF file. If FALSE, the plot will be displayed interactively without creating a PDF.
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed.
#'
#' @return A ggplot object containing the longitudinal line plot.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Phylum"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   top.k.plot = 3,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' }
#' @export
generate_taxa_spaghettiplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
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
      {
        if ("original" %in% feature.level)
          dplyr::mutate(., original = rownames(.))
        else
          .
      } %>%
      select(all_of(feature.level))

    if (is.null(palette)) {
      col <-
        c(
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
    } else{
      col <- palette
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

    # Calculate new sizes based on base.size
    title.size = base.size * 1.25
    axis.title.size = base.size * 0.75
    axis.text.size = base.size * 0.5
    legend.title.size = base.size * 1
    legend.text.size = base.size * 0.75

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list <- lapply(feature.level, function(feature.level){
      # 将 OTU 表与分类表合并
      otu_tax <-
        cbind(otu_tab, tax_tab %>% select(all_of(feature.level)))

      # Filter taxa based on prevalence and abundance
      otu_tax_filtered <- otu_tax %>%
        tidyr::gather(key = "sample", value = "count", -one_of(feature.level)) %>%
        dplyr::group_by_at(vars(!!sym(feature.level))) %>%
        dplyr::summarise(total_count = mean(count),
                  prevalence = sum(count > 0) / dplyr::n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count, -prevalence) %>%
        dplyr::left_join(otu_tax, by = feature.level)

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax_filtered %>%
        tidyr::gather(key = "sample",  value = "count", -one_of(feature.level)) %>%
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

      df <- otu_tax_agg_numeric %>%
        tidyr::gather(key = "sample", value = "count", -one_of(feature.level)) %>%
        dplyr::left_join(meta_tab %>% rownames_to_column(var = "sample"), by = "sample")

      if (is.null(group.var)) {
        df <- df %>% dplyr::mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      if (!is.null(strata.var)) {
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level),
                          !!sym(time.var),
                          !!sym(group.var),
                          !!sym(strata.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df,
                    mean_df,
                    by = c(feature.level, time.var, group.var, strata.var))
      } else {
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level), !!sym(time.var), !!sym(group.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df, mean_df, by = c(feature.level, time.var, group.var))
      }

      taxa.levels <-
        df %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()

      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 5) {
          features.plot <- taxa.levels[1:4]
        } else {
          features.plot <- taxa.levels
        }
      }

      sub_df <- df %>% filter(!!sym(feature.level) %in% features.plot)

      lineplot <- ggplot() +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "count",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.5
        ) +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_count",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        scale_color_manual(values = col) +
        {
          if (feature.dat.type != "other"){
            labs(
              x = time.var,
              y = "Relative Abundance",
              color = group.var
            )
          } else {
            labs(
              x = time.var,
              y = "Abundance",
              color = group.var
            )
          }
        } +
        {
          if (!is.null(strata.var)) {
            ggh4x::facet_nested_wrap(as.formula(paste('~', feature.level,"+",strata.var)))  # Use facet_wrap with strata.var as the faceting variable
          } else {
            ggh4x::facet_nested_wrap(as.formula(paste('~',feature.level)))
          }
        } +
        theme_to_use +
        theme(
          plot.title = element_text(
            size = title.size,
            face = "bold",
            hjust = 0.5
          ),
          panel.spacing.x = unit(0, "cm"),
          axis.title.x = element_text(size = axis.title.size, face = "bold"),
          axis.title.y = element_text(size = axis.title.size, face = "bold"),
          axis.text.x = element_text(size = axis.text.size),
          axis.text.y = element_text(size = axis.text.size),
          legend.title = element_text(size = legend.title.size, face = "bold"),
          legend.text = element_text(size = legend.text.size)
        )

      if (group.var == "ALL") {
        lineplot <- lineplot + theme(legend.position = "none")
      }

      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "taxa_spaghettiplot_long",
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
          abund.filter,
          "_",
          "base_size_",
          base.size,
          "_",
          "theme_choice_",
          theme.choice,
          "_",
          "pdf_wid_",
          pdf.wid,
          "_",
          "pdf_hei_",
          pdf.hei
        )

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }

        pdf_name <- paste0(pdf_name, ".pdf")

        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Use lapply to print each ggplot object in the list to a new PDF page
        print(lineplot)
        # Close the PDF device
        dev.off()
      }

      return(lineplot)
    })

    return(plot_list)
  }
