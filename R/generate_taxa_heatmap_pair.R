#' @title Generate Taxonomic Heatmap Pair
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param feature.level A character vector specifying the taxa level(s) to include in the analysis. Default is c('Phylum', 'Family', 'Genus').
#' @param features.plot A character vector specifying the taxa to be plotted. If NULL (default), the top k taxa by mean abundance will be plotted.
#' @param feature.dat.type A character string specifying the type of data in feature.dat. Options are "count", "proportion", or "other".
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
#' @param prev.filter A numeric value defining the prevalence threshold to filter taxa, between 0 and 1.
#' @param abund.filter A numeric value defining the abundance threshold to filter taxa.
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
#'
#' # Load required libraries and example data
#' library(tidyverse)
#' library(pheatmap)
#' data(peerj32.obj)
#' plot_list <- generate_taxa_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.1,
#'   abund.filter = 0.1,
#'   base.size = 11,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' @export
#'
#' @seealso \code{\link{pheatmap}}
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

  # Extract data
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
  )))

  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  if (is.null(palette)) {
    palette <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
  }
  # 创建颜色映射函数
  my_palette <- colorRampPalette(palette)

  # 计算颜色的数量
  # 这通常取决于你的数据，你可能需要根据你的实际情况进行调整
  n_colors <- 100

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
      gather(key = "sample", value = "count",-one_of(colnames(tax_tab))) %>%
      group_by_at(vars(!!sym(feature.level))) %>%
      summarise(total_count = mean(count),
                prevalence = sum(count > 0) / n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count,-prevalence) %>%
      left_join(otu_tax, by = feature.level)

    # Aggregate OTU table
    otu_tax_agg <- otu_tax_filtered %>%
      gather(key = "sample", value = "count",-one_of(colnames(tax_tab))) %>%
      group_by_at(vars(sample,!!sym(feature.level))) %>%
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
      otu_tax_agg_numeric %>% column_to_rownames(var = feature.level) %>% as.matrix()

    # Sort samples by group.var and strata.var
    if (!is.null(group.var) & !is.null(strata.var)) {
      meta_tab_sorted <-
        meta_tab[order(meta_tab[[group.var]], meta_tab[[strata.var]]),]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else if (!is.null(group.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[group.var]]),]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else {
      otu_tab_norm_sorted <- otu_tab_norm
    }

    # Calculate gaps if group.var is not NULL
    gaps <-
      cumsum(table(meta_tab_sorted[[subject.var]]))[-length(table(meta_tab_sorted[[subject.var]]))]

    # Set up annotation_col based on group.var and strata.var values
    if (!is.null(group.var) & !is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(time.var, strata.var, group.var)))
    } else if (!is.null(group.var) & is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(time.var, group.var)))
    } else {
      annotation_col <- meta_tab %>% select(all_of(c(time.var)))
    }


    if (!is.null(features.plot)){

      otu_tab_norm_sorted <- otu_tab_norm_sorted[rownames(otu_tab_norm_sorted) %in% features.plot, ]

    }

    # Plot stacked heatmap
    heatmap_plot <- pheatmap(
      otu_tab_norm_sorted,
      annotation_col = annotation_col,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      gaps_col = gaps,
      color = my_palette(n_colors),
      # 使用自定义颜色
      fontsize = base.size,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_pair",
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
    # Return the heatmap plot for display
    return(gg_heatmap_plot)
  })

  return(plot_list)
}
