#' Compare beta diversity between time points
#'
#' This function generates boxplots to visualize changes in beta diversity
#' between two time points, for different groups and strata.
#'
#' For each subject, it calculates the beta diversity (distance)
#' between the sample at baseline time point and the sample at later time point.
#' This within-subject change in beta diversity over time is used as the outcome.
#'
#' Boxplots are generated to compare the distribution of beta diversity changes
#' across groups and strata. Violin plots are also included to show the underlying distribution.
#'
#' Adjustment for covariates is supported by using adjusted distances.
#'
#' @name generate_beta_change_boxplot_pair
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param subject.var A string specifying the name of the subject variable
#' @param time.var A string specifying the name of the time variable
#' @param group.var A string specifying the name of the group variable or NULL (default)
#' @param strata.var A string specifying the name of the strata variable or NULL (default)
#' @param adj.vars A string specifying the name of the adjustment variable or NULL (default)
#' @param change.base A string or numeric value specifying the first time point used for computing the beta diversity change
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param base.size (Optional) Base font size for the plot (default is 16).
#' @param theme.choice Character string specifying ggplot2 theme to use.
#' Options include:
#' - "prism": ggprism::theme_prism()
#' - "classic": theme_classic()
#' - "gray": theme_gray()
#' - "bw": theme_bw()
#'
#' Can also take a custom ggplot2 theme object.
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
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf A logical value indicating whether to save the plot as a PDF file. Default is TRUE
#' @param file.ann A character string specifying a custom annotation for the PDF file name or NULL (default)
#' @param pdf.wid (Optional) The width of the PDF file if `pdf` is set to `TRUE` (default is 11).
#' @param pdf.hei (Optional) The height of the PDF file if `pdf` is set to `TRUE` (default is 8.5).
#' @param ... Additional parameters passed on to ggsave()
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object.
#'
#' @return A named list of ggplot objects visualizing beta diversity change.
#' The list contains one plot for each distance metric specified in \code{dist.name}.
#' Each plot shows boxplots of the beta diversity changes, with samples faceted by
#' the \code{group_var} and \code{strata_var} variables if provided.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(vegan)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC")
#' generate_beta_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#' generate_beta_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   change.base = "Baseline",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
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
           adj.vars = NULL,
           change.base,
           dist.name = c('BC', 'Jaccard'),
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    if (is.null(dist.name)){
      return()
    }

    if (is.null(dist.obj)) {
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
        data.obj <- list(meta.dat = meta_tab)
        meta_tab <- data.obj$meta.dat
      }
    }

    meta_tab <- meta_tab %>% rownames_to_column("sample")

    if (is.null(change.base)){
      change.base <- unique(meta_tab %>% dplyr::select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'alpha_df' data frame: ", change.base)
    }

    change.after <-
      unique(meta_tab %>% dplyr::select(all_of(c(time.var))))[unique(meta_tab %>% dplyr::select(all_of(c(time.var)))) != change.base]

    col <- mStat_get_palette(palette)

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

      if (is.null(adj.vars)) {
        y_label <- paste0("Distance from ", change.base, " to ", change.after)
      } else {
        y_label <- paste0("Distance from ", change.base, " to ", change.after, " (adjusted by: ", paste(adj.vars, collapse = ", "), ")")
      }

      dist.df <- as.matrix(dist.obj[[dist.name]]) %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(-time.var) %>% dplyr::distinct(), by = subject.var)

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
        ylab(y_label) +
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
    names(plot_list) <- dist.name
    return(plot_list)
  }
