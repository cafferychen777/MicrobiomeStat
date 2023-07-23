#' Beta Diversity Change Boxplot Pairs
#'
#' This function generates boxplots to visualize the changes in beta diversity principal components (PCs) over time.
#' It allows the use of Principal Coordinates Analysis (PCoA), non-metric multidimensional scaling (NMDS), t-SNE, or UMAP for dimension reduction.
#'
#' @param data.obj A list containing the input data. Default is NULL.
#' @param dist.obj A list containing the distance object. Default is NULL.
#' @param pc.obj A list containing the Principal Component Analysis object. Default is NULL.
#' @param pc.ind A numeric vector indicating the PC indexes used. Default is c(1, 2).
#' @param subject.var A string specifying the variable for subjects.
#' @param time.var A string specifying the variable for time.
#' @param group.var A string specifying the variable for groups. Default is NULL.
#' @param strata.var A string specifying the variable for strata. Default is NULL.
#' @param change.base The baseline for calculating changes in beta diversity. Default is NULL.
#' @param change.func A function or string specifying how to calculate changes. Default is "difference".
#' @param dist.name A character vector indicating the distance metrics used. Default is c("BC", "Jaccard").
#' @param base.size A numeric value for the base size of the plot. Default is 16.
#' @param theme.choice A string specifying the theme of the plot. Default is "prism".
#' @param custom.theme A ggplot2 theme object for user-defined theme. Default is NULL.
#' @param palette A character vector specifying the color palette. Default is NULL.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#' @return A list of ggplot objects for each PC index and distance metric.
#' @details
#' This function generates a boxplot of changes in beta diversity based on PCoA coordinates for longitudinal data.
#' The boxplot can be stratified by a group variable and/or other variables. It also allows for different
#' distance metrics and principal component indexes to be used.
#' The function can handle a large number of time points or subjects by averaging the data and adding jitter to the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' # Load required libraries and example data
#' library(tidyverse)
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC'))
#' # Add metadata information
#' attr(dist.obj[["BC"]], "labels") <- peerj32_obj$meta.dat
#'
#' # Generate the boxplot pair
#' generate_beta_pc_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   change.base = "1",
#'   change.func = "difference",
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
#'
#' @export
generate_beta_pc_change_boxplot_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           pc.ind = c(1, 2),
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           change.func = "difference",
           dist.name = c('BC', 'Jaccard'),
           base.size = 16,
           theme.choice = "prism",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    if (is.null(data.obj) & is.null(dist.obj)) {
      stop("Both data.obj and dist.obj cannot be NULL. Please provide at least one.")
    }

    if (is.null(dist.obj)) {
      message("No dist.obj provided, calculating beta diversity.")
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      metadata <- load_data_obj_metadata(data.obj)
      if (is.null(metadata)) {
        stop("No metadata could be loaded from data.obj. Please ensure it contains the necessary metadata.")
      }
    } else {
      if (!all(dist.name %in% names(dist.obj))) {
        stop(paste0("The requested dist.name(s) ", paste(dist.name[!dist.name %in% names(dist.obj)], collapse = ", "),
                    " are not available in the provided dist.obj. Please check again."))
      }
      metadata <- attr(dist.obj[[dist.name[1]]], "labels")
      if (is.null(metadata)) {
        message("No metadata found in dist.obj. Attempting to load metadata from data.obj.")
        if (is.null(data.obj)) {
          stop("No data.obj provided to load metadata from. Please ensure either dist.obj or data.obj contain the necessary metadata.")
        }
        metadata <- load_data_obj_metadata(data.obj)
        if (is.null(metadata)) {
          stop("No metadata could be loaded from data.obj. Please ensure it contains the necessary metadata.")
        }
      }
    }

    if (is.null(pc.obj)) {
      message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
      message("If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method.")
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = max(pc.ind),
          dist.name = dist.name
        )
    }

    if (is.null(palette)){
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
        custom.theme else theme_function

    plot_list <- lapply(dist.name, function(dist.name) {
        pc.mat <- pc.obj[[dist.name]]$points

        # 将PC和元数据组合成一个数据框
        colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

        pc.mat <- pc.mat %>% as_tibble()

        df <-
          cbind(pc.mat[, paste0("PC", pc.ind)], metadata[, c(subject.var, time.var, group.var, strata.var)])

        change.after <-
          unique(df %>% select(all_of(c(time.var))))[unique(df %>% select(all_of(c(time.var)))) != change.base]

        df <-
          df %>%
          as_tibble() %>%
          gather(
            key = "PC",
            value = "value",-one_of(subject.var, group.var, time.var, strata.var)
          )

        # 拆分成一个列表，每个time值都有一个独立的tibble
        split_data <-
          split(df, f = df %>%
                  group_by(!!sym(time.var)) %>% select(!!sym(time.var)))

        # 提取split_data中的第一个和第二个表
        data_time_1 <- split_data[[change.base]]
        data_time_2 <- split_data[[change.after]]

        combined_data <- data_time_1 %>%
          inner_join(
            data_time_2,
            by = c("PC", subject.var),
            suffix = c("_time_1", "_time_2")
          )

        combined_data <- combined_data %>%
          mutate(value_diff = if (is.function(change.func)) {
            change.func(value_time_2, value_time_1)
          } else if (change.func == "difference"){
            value_time_2 - value_time_1
          } else {
            value_time_2 - value_time_1
          })

        combined_data <-
          combined_data %>% left_join(metadata %>% select(all_of(
            c(subject.var, time.var, group.var, strata.var)
          )) %>% filter(!!sym(time.var) == change.after),
          by = subject.var)

        if (is.null(group.var)) {
          group.var = "ALL"
          combined_data$ALL <- "ALL"
        }

        lapply(paste0("PC", pc.ind), function(pc.index) {
          boxplot <- ggplot(
            combined_data %>% filter(PC == pc.index),
            aes(
              x = !!sym(group.var),
              y = value_diff,
              fill = !!sym(group.var)
            )
          ) +
            geom_violin(trim = F, alpha = 0.8) +
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
            geom_jitter(width = 0.1,
                        alpha = 0.5,
                        size = 1.7) +
            scale_alpha_manual(values = c(0.5, 0.5)) +
            scale_fill_manual(values = col) +
            labs(x = group.var, y = paste("Change in ", "Axis ", gsub("PC", "", pc.index), " - ",
                                          if(is.function(change.func)){
                                            "custom function"
                                          } else {
                                            change.func
                                          })) +
            theme_to_use +
            theme(
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0, "cm"),
              strip.text.x = element_text(size = 12, color = "black"),
              axis.text.x = element_blank(),
              axis.text.y = element_text(color = "black", size = base.size),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = base.size),
              axis.ticks.x = element_blank(),
              plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
              legend.text = ggplot2::element_text(size = base.size),
              legend.title = ggplot2::element_text(size = base.size)
            )

          if (is.null(strata.var)) {
            boxplot <- boxplot
          } else {
            boxplot <- boxplot +
              ggh4x::facet_nested(cols = vars(!!sym(strata.var)),
                         scales = "fixed",
                         space = "free")
          }


            # 不显示图例的条件
            if (group.var == "ALL") {
              boxplot <- boxplot  + theme(
                axis.text.x = element_blank(),  # 隐藏x轴text
                axis.title.x = element_blank(),  # 隐藏x轴title
                legend.position = "none",  # 隐藏图例
                strip.text.x = element_blank()  # 隐藏分面title
              )
            }


          # Save the plots as a PDF file
          if (pdf) {
              pdf_name <- paste0(
                "beta_pc_change_boxplot_pair_",
                "pc.ind_",
                pc.index,
                "_",
                dist.name,
                "_",
                "method_",
                method,
                "_",
                "subject_",
                subject.var,
                "_",
                "time_",
                time.var,
                "_",
                "change_base_",
                change.base
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
