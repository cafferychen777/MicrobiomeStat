#' @title Generate Taxonomic Change Heatmap Long
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param t0.level The base level for time points in longitudinal data.
#' @param ts.levels The levels for time points in longitudinal data.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var (Optional) A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param change.func A function or character string specifying the method for computing the change. Default is "difference".
#' @param feature.level A character string defining the taxonomic level to analyze ('Phylum', 'Family', or 'Genus').
#' @param features.plot A character vector specifying the taxa to be plotted. If NULL (default), the top k taxa by mean abundance will be plotted.
#' @param feature.dat.type A character string specifying the type of the data in feature.dat. Options are "count", "proportion", or "other".
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
#' @param prev.filter A numeric value defining the prevalence threshold to filter taxa, between 0 and 1.
#' @param abund.filter A numeric value defining the abundance threshold to filter taxa.
#' @param base.size Base font size for the generated plots.
#' @param palette Color palette used for the plots.
#' @param cluster.rows A logical variable indicating if rows should be clustered. Default is TRUE.
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is FALSE.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional arguments to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the pheatmap plot. If `pdf` is set to FALSE, the function will return the pheatmap plot without creating a PDF file.
#'
#' @examples
#' library(tidyverse)
#' library(pheatmap)
#' data(ecam.obj)
#'
#' plot_list <- generate_taxa_change_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   feature.level = c("Family"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   change.func = "lfc",
#'   palette = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = "test"
#' )
#'
#'
#' @return An object of class pheatmap, the generated heatmap plot
#' @export
#'
#' @seealso \code{\link{pheatmap}}
generate_taxa_change_heatmap_long <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              t0.level = NULL,
                                              ts.levels = NULL,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              feature.level,
                                              feature.dat.type = c("count", "proportion", "other"),
                                              features.plot = NULL,
                                              top.k.plot = NULL,
                                              top.k.func = NULL,
                                              change.func = "relative difference",
                                              prev.filter = 0.00000001,
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
      meta_tab %>% mutate(!!sym(group.var) := interaction(!!sym(strata.var), !!sym(group.var)))
  }

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  plot_list <- lapply(feature.level, function(feature.level) {
    # Merge OTU table with taxonomy table
    otu_tax <-
      cbind(otu_tab, tax_tab %>% select(all_of(feature.level)))

    # Filter taxa based on prevalence and abundance
    otu_tax_filtered <- otu_tax %>%
      gather(key = "sample", value = "count",-one_of(feature.level)) %>%
      group_by_at(vars(!!sym(feature.level))) %>%
      summarise(total_count = mean(count),
                prevalence = sum(count > 0) / n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count,-prevalence) %>%
      left_join(otu_tax, by = feature.level)

    # Aggregate OTU table
    otu_tax_agg <- otu_tax_filtered %>%
      gather(key = "sample", value = "count",-one_of(feature.level)) %>%
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
      otu_tax_agg_numeric %>% filter(!is.na(!!sym(feature.level))) %>% column_to_rownames(var = feature.level) %>% as.matrix()

    # Calculate the mean_value for each combination of feature.level, group.var, and time.var
    df_mean_value <- otu_tab_norm %>%
      as.data.frame() %>%
      rownames_to_column(var = feature.level) %>%
      gather(key = "sample", value = "value",-one_of(feature.level)) %>%
      left_join(meta_tab %>%
                  rownames_to_column("sample"), "sample") %>%
      group_by(!!sym(feature.level),!!sym(group.var),!!sym(time.var)) %>%
      summarise(mean_value = mean(value), .groups = "drop")

    # Spread the data to wide format
    df_wide <- df_mean_value %>%
      spread(key = time.var, value = mean_value)
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

      if (is.function(change.func)) {
        df_wide <- df_wide %>%
          rowwise() %>%
          mutate("{change_col_name}" := change.func(.data[[as.character(ts)]], .data[[as.character(t0.level)]]))
      } else if (change.func == "relative difference") {
        df_wide <- df_wide %>%
          rowwise() %>%
          mutate("{change_col_name}" := case_when(
            (.data[[as.character(ts)]] + .data[[as.character(t0.level)]]) != 0 ~ (.data[[as.character(ts)]] - .data[[as.character(t0.level)]]) / (.data[[ts]] + .data[[as.character(t0.level)]]),
            TRUE ~ 0
          ))
      } else if (change.func == "difference") {
        df_wide <- df_wide %>%
          rowwise() %>%
          mutate("{change_col_name}" := (.data[[as.character(ts)]] - .data[[as.character(t0.level)]]))
      } else if (change.func == "lfc") {
        df_wide <- df_wide %>%
          rowwise() %>%
          mutate("{change_col_name}" := log2(.data[[as.character(ts)]] + 0.00001) - log2(.data[[as.character(t0.level)]] + 0.00001))
      } else {
        stop(
          "`change.func` must be either 'relative difference', 'difference', 'lfc' or a function."
        )
      }
    }

    df_wide <-
      df_wide[, c(feature.level, group.var, paste0("change_", ts.levels))]

    # 首先将数据框转化为长格式
    df_long <- df_wide %>%
      pivot_longer(
        cols = starts_with("change"),
        names_to = "time",
        values_to = "value"
      )

    # 删除 "time" 列名中的 "change_" 部分
    df_long$time <- gsub("change_", "", df_long$time)

    # 再将数据框转化为宽格式，并在此过程中修改列名
    df_wide_new <- df_long %>%
      unite("group_time", c(group.var, "time"), sep = "_") %>%
      pivot_wider(names_from = "group_time",
                  values_from = "value")

    wide_data <- df_wide_new %>% column_to_rownames(feature.level)

    df_wide <-
      df_wide[, c(feature.level, group.var, paste0("change_", ts.levels))]

    # 首先将数据框转化为长格式
    df_long <- df_wide %>%
      pivot_longer(
        cols = starts_with("change"),
        names_to = "time",
        values_to = "value"
      )

    # 删除 "time" 列名中的 "change_" 部分
    df_long$time <- gsub("change_", "", df_long$time)

    # 再将数据框转化为宽格式，并在此过程中修改列名
    df_wide_new <- df_long %>%
      unite("group_time", c(group.var, "time"), sep = "_") %>%
      pivot_wider(names_from = "group_time",
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
      distinct() %>%
      mutate(group_time = paste(!!sym(group.var), !!sym(time.var), sep = "_")) %>%
      filter(group_time %in% new_colnames) %>%
      column_to_rownames("group_time")
    annotation_col_sorted <-
      annotation_col[order(annotation_col[[group.var]], annotation_col[[time.var]]),]
    if (!is.null(strata.var)) {
      annotation_col_sorted <- annotation_col_sorted %>%
        separate(!!sym(group.var),
                 into = c(strata.var, group.var),
                 sep = "\\.")
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

    n_colors <- 100

    if (is.null(palette)) {
      palette <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
      # 找出数据中的最大绝对值
      max_abs_val <-
        max(abs(range(na.omit(
          c(as.matrix(wide_data_sorted))
        ))))

      # 计算零值在新色彩向量中的位置
      zero_pos <- round(max_abs_val / (2 * max_abs_val) * n_colors)

      # 创建颜色向量
      my_palette <-
        c(
          colorRampPalette(palette[1:3])(zero_pos),
          colorRampPalette(palette[3:5])(n_colors - zero_pos + 1)
        )

      # Plot stacked heatmap
      heatmap_plot <- pheatmap::pheatmap(
        wide_data_sorted,
        annotation_col = annotation_col_sorted,
        cluster_rows = cluster.rows,
        cluster_cols = cluster.cols,
        show_colnames = FALSE,
        gaps_col = gaps,
        fontsize = base.size,
        color = my_palette
      )
    } else {
      # 创建颜色映射函数
      my_palette <- colorRampPalette(palette)

      # Plot stacked heatmap
      heatmap_plot <- pheatmap(
        wide_data_sorted,
        annotation_col = annotation_col_sorted,
        cluster_rows = cluster.rows,
        cluster_cols = cluster.cols,
        show_colnames = FALSE,
        gaps_col = gaps,
        fontsize = base.size,
        color = my_palette(n_colors)
      )
    }

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


  # Return the heatmap plot for display
  return(plot_list)
}
