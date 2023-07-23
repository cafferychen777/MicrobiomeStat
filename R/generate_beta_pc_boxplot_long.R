#' Create Boxplots for Beta Diversity Metrics
#'
#' This function generates boxplots of beta diversity metrics for different groups over time.
#' It allows the use of Principal Coordinates Analysis (PCoA), non-metric multidimensional scaling (NMDS), t-SNE, or UMAP for dimension reduction.
#' Default method of dimension reduction is PCoA.
#'
#' @param data.obj A list object containing all data.
#' @param dist.obj List object containing distance matrices.
#' @param pc.obj List object containing Principal Coordinates matrices.
#' @param pc.ind Numeric vector indicating which Principal Coordinates to plot.
#' @param subject.var Character string indicating the variable for subject identifiers.
#' @param time.var Character string indicating the variable for time points.
#' @param t0.level Character string indicating the baseline time point.
#' @param ts.levels Vector of character strings indicating the time points to include in the plot.
#' @param group.var Character string indicating the variable for group identifiers.
#' @param strata.var Character string indicating the variable for stratum identifiers.
#' @param dist.name Vector of character strings indicating the distance measures to include in the plot.
#' @param base.size Numeric value indicating the base font size for the plot.
#' @param theme.choice Character string indicating the ggplot theme to use for the plot.
#' @param custom.theme Custom ggplot theme to use.
#' @param palette Vector of colors for the plot.
#' @param pdf Logical indicating whether to save the plot as a PDF file.
#' @param file.ann Optional character string to append to the file name.
#' @param pdf.wid Numeric value indicating the width of the PDF file.
#' @param pdf.hei Numeric value indicating the height of the PDF file.
#' @param ... Additional arguments to be passed to the plotting function.
#'
#' @return A list of ggplot objects for each distance measure and Principal Coordinate.
#' @examples
#' # Load required libraries and data
#' library(vegan)
#'
#'   data(ecam.obj)
#'   generate_beta_pc_boxplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2],
#'   group.var = "diet",
#'   strata.var = "delivery",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#'   generate_beta_pc_boxplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_beta_pc_boxplot_long <- function(data.obj = NULL,
                                          dist.obj = NULL,
                                          pc.obj = NULL,
                                          pc.ind = c(1, 2),
                                          subject.var,
                                          time.var,
                                          t0.level = NULL,
                                          ts.levels = NULL,
                                          group.var = NULL,
                                          strata.var = NULL,
                                          dist.name = c("BC"),
                                          base.size = 16,
                                          theme.choice = "prism",
                                          custom.theme = NULL,
                                          palette = NULL,
                                          pdf = TRUE,
                                          file.ann = NULL,
                                          pdf.wid = 11,
                                          pdf.hei = 8.5,
                                          ...) {
  if (is.null(dist.obj)) {
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
    dist.obj <-
      mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
  } else {
    if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
    } else {
      meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
      data.obj <- list(meta.dat = meta_tab)
      data.obj <- mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
      meta_tab <- load_data_obj_metadata(data.obj)
      dist.obj <- mStat_subset_dist(dist.obj, colnames(meta_tab))
    }
  }

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

  # 根据 strata.var 和 group.var 的值调整 ggplot() 函数
  aes_function <- if (!is.null(group.var)) {
    aes(
      x = !!sym(time.var),
      y = value,
      fill = !!sym(group.var)
    )
  } else {
    aes(
      x = !!sym(time.var),
      y = value,
      fill = !!sym(time.var)
    )
  }

  line_aes_function <- if (!is.null(group.var)) {
    aes(
      x = !!sym(time.var),
      y = value,
      group = !!sym(subject.var),
      color = !!sym(group.var)
    )
  } else {
    aes(
      x = !!sym(time.var),
      y = value,
      group = !!sym(subject.var)
    )
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

  plot_list <- lapply(dist.name, function(dist.name) {
    if (is.null(pc.obj)) {
      message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
      message("If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method.")
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj[dist.name],
          method = "mds",
          k = max(pc.ind),
          dist.name = dist.name
        )
    }

      pc.mat <- pc.obj[[dist.name]]$points

      colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

      pc.mat <- pc.mat %>% as_tibble()

      df <-
        cbind(pc.mat[, paste0("PC", pc.ind)], meta_tab[, c(subject.var, time.var, group.var, strata.var)])

      df <-
        df %>%
        as_tibble() %>%
        gather(
          key = "PC",
          value = "value",-one_of(subject.var, group.var, time.var, strata.var)
        )

      n_subjects <- length(unique(df[[subject.var]]))
      n_times <- length(unique(df[[time.var]]))

      lapply(unique(df$PC), function(pc.index) {
        sub_df <- df %>% filter(PC == pc.index)

        # 在数据处理部分创建一个新的数据框
        average_sub_df <- NULL
        if (n_times > 10 || n_subjects > 25) {
          if (!is.null(strata.var) & !is.null(group.var)){
            average_sub_df <- sub_df %>%
              dplyr::group_by(!!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
              summarise(dplyr::across(value, mean, na.rm = TRUE), .groups = "drop") %>%
              dplyr::ungroup() %>%
              dplyr::mutate(!!sym(subject.var) := "ALL")
          } else if (!is.null(group.var)) {
            average_sub_df <- sub_df %>%
              dplyr::group_by(!!sym(group.var), !!sym(time.var)) %>%
              summarise(dplyr::across(value, mean, na.rm = TRUE), .groups = "drop") %>%
              dplyr::ungroup() %>%
              dplyr::mutate(!!sym(subject.var) := "ALL")
          } else {
            average_sub_df <- sub_df %>%
              dplyr::group_by(!!sym(time.var)) %>%
              summarise(dplyr::across(value, mean, na.rm = TRUE), .groups = "drop") %>%
              dplyr::ungroup() %>%
              dplyr::mutate(!!sym(subject.var) := "ALL")
          }
        }

        boxplot <- ggplot(sub_df,
                          aes_function) +
          geom_violin(trim = FALSE, alpha = 0.8) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.1
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.1,
            fill = "white"
          ) +
          geom_line(
            line_aes_function,
            alpha = 0.8,
            linewidth = 0.6,
            color = "black",
            linetype = "dashed", # 更改线条类型为虚线
            data = if (!is.null(average_sub_df)) average_sub_df else sub_df
          ) +
          scale_fill_manual(values = col) +
          labs(
            x = time.var,
            y = paste(
              "Distance:",
              dist.name,
              " - Axis",
              gsub("PC", "", pc.index)
            )
          ) +
        theme_to_use +
          theme(
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0, "cm"),
            strip.text.x = element_text(size = 12, color = "black"),
            axis.text.x = element_text(color = "black", size = base.size),
            axis.text.y = element_text(color = "black", size = base.size),
            axis.title.x = element_text(size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            legend.text = ggplot2::element_text(size = base.size),
            legend.title = ggplot2::element_text(size = base.size)
          )

        if (!is.null(group.var)) {
          if (is.null(strata.var)) {
            boxplot <-
              boxplot + ggh4x::facet_nested(
                cols = vars(!!sym(group.var)),
                scales = "free",
                space = "free"
              )
          } else {
            boxplot <-
              boxplot + ggh4x::facet_nested(
                cols = vars(!!sym(group.var),!!sym(strata.var)),
                scales = "free",
                space = "free"
              )
          }
        }

        # Add geom_jitter() if the number of unique time points or subjects is greater than 10
        if (n_subjects > 10 || n_times > 10) {
          boxplot <- boxplot + geom_jitter(width = 0.1, alpha = 0.1, size = 1)
        }

        # Save the plots as a PDF file
        if (pdf) {
          pdf_name <- paste0(
            "beta_pc_boxplot_long_",
            "pc.ind_",
            pc.index,
            "_",
            dist.name,
            "_",
            "subject_",
            subject.var,
            "_",
            "time_",
            time.var
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
            plot = boxplot,
            width = pdf.wid,
            height = pdf.hei,
            dpi = 300
          )
        }
        return(boxplot)
    })
  })


  return(plot_list)
}
