#' @title Generate box plot of beta diversity change
#'
#' @description This function generates a box plot that visualizes the change in beta diversity
#' between two time points for different groups and strata. It takes a \code{dist.obj},
#' and can also handle a matrix with coordinates for different groups and strata.
#' Optionally, a user can specify whether to save the plot as a PDF file.
#' @name generate_beta_change_boxplot_pair
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj A dissimilarity object created with the vegdist() function or NULL. If NULL, the function will automatically compute the beta diversity based on the input \code{data.obj}.
#' @param subject.var A string specifying the name of the subject variable
#' @param time.var A string specifying the name of the time variable
#' @param group.var A string specifying the name of the group variable or NULL (default)
#' @param strata.var A string specifying the name of the strata variable or NULL (default)
#' @param change.base A string or numeric value specifying the first time point used for computing the beta diversity change
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#' @param base.size (Optional) Base font size for the plot (default is 16).
#' @param theme.choice (Optional) Name of the theme for the plot. Default is "prism". Other options include "plain", "classic", and any other themes compatible with ggplot2.
#' @param custom.theme (Optional) A custom ggplot2 theme.
#' @param palette (Optional) A palette function or character vector with the colors for the plot.
#' @param pdf A logical value indicating whether to save the plot as a PDF file. Default is TRUE
#' @param file.ann A character string specifying a custom annotation for the PDF file name or NULL (default)
#' @param pdf.wid (Optional) The width of the PDF file if `pdf` is set to `TRUE` (default is 11).
#' @param pdf.hei (Optional) The height of the PDF file if `pdf` is set to `TRUE` (default is 8.5).
#' @param ... Additional parameters passed on to ggsave()
#'
#' @return A ggplot object with the box plot of beta diversity change.
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object. This package was created by Caffery Yang. For more information, you can contact them at cafferychen7850@gmail.com.
#'
#' @return A list of ggplot objects with the box plot of beta diversity change for each diversity index specified.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(vegan)
#' data(peerj32.obj)
#'
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC")
#' generate_beta_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = dist.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
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
#' }
#' @export
generate_beta_change_boxplot_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base,
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

    if (is.null(dist.obj)&!is.null(data.obj)) {
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,time.var,strata.var))) %>% rownames_to_column("sample")
    } else {
      metadata <- attr(dist.obj[[dist.name[1]]], "labels")  %>% select(all_of(c(subject.var,group.var,time.var,strata.var))) %>% rownames_to_column("sample")
    }

    if (is.null(change.base)){
      change.base <- unique(metadata %>% select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'alpha_df' data frame: ", change.base)
    }

    change.after <-
      unique(metadata %>% select(all_of(c(time.var))))[unique(metadata %>% select(all_of(c(time.var)))) != change.base]

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
    } else {
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

    plot_list <- lapply(dist.name,function(dist.name){
      dist.df <- as.matrix(dist.obj[[dist.name]]) %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(metadata, by = "sample") %>%
        dplyr::left_join(metadata, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      long.df <- long.df %>% dplyr::left_join(metadata %>% select(-time.var) %>% dplyr::distinct(), by = subject.var)

      facet_formula <-
        if (!is.null(strata.var)) {
          paste(". ~", strata.var)
        } else {
          ". ~ 1"
        }

      long.df$x_alternative <- "ALL"
      aes_function <- if (!is.null(group.var)){
        aes(
          x = !!sym(group.var),
          y = distance,
          fill = !!sym(group.var)
        )
      }else{
        aes(
          x = x_alternative,
          y = distance,
          fill = x_alternative
        )
      }

      p <-
        ggplot(long.df, aes_function) +
        geom_violin(trim = F,alpha = 0.8) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
          fill = "white"
        ) +
        geom_jitter(width = 0.1, alpha = 0.5, size = 1.7) +
        scale_fill_manual(values = col) +
        facet_wrap(as.formula(facet_formula), scales = "fixed") +
        xlab(group.var) +
        ylab(paste0("Distance from ", change.base," to ", change.after)) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0.2, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
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

      if(is.null(group.var) && is.null(strata.var)) {
        p <- p + theme(
          axis.text.x = element_blank(),  # 隐藏x轴text
          axis.title.x = element_blank(),  # 隐藏x轴title
          legend.position = "none",  # 隐藏图例
          strip.text.x = element_blank()  # 隐藏分面title
        )
      }

      if (!is.null(group.var) & is.null(strata.var)){
        p <- p + theme(
          strip.text.x = element_blank()  # 隐藏分面title
        )
      }


      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0("beta_change_boxplot_pair_",
                           dist.name,
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "change_base_", change.base)

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
          plot = p,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      }
      return(p)
    })
    return(plot_list)
  }
