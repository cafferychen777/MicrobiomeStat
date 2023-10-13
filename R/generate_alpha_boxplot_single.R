#' @title Generate boxplot for alpha diversity index at a single time point
#'
#' @description This function generates a boxplot of specified alpha diversity indices at a single time point dplyr::across different groupings, with optional stratification. The output can be saved as a PDF.
#' @name generate_alpha_boxplot_single
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var An optional variable in the metadata table that represents the grouping factor.
#' @param strata.var An optional variable in the metadata table that represents the stratification factor.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param base.size The base font size for the plot. Default is 16.
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
#' @param palette An optional palette to use for the plot. If NULL (default), a pre-defined palette will be used.
#' @param pdf A boolean indicating whether to save the output as a PDF file. Default is TRUE.
#' @param file.ann A string for annotating the output file name.
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @return A list of boxplots displaying the specified alpha diversity indices at the specified time point dplyr::across different groupings, stratified by the specified stratification variable (if provided). Each boxplot in the list corresponds to one of the alpha diversity indices specified in `alpha.name`. The boxplots will be saved as PDF files if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(ggh4x)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # First example with peerj32.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("simpson"),
#'   subject.var  = "subject",
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' # Load another dataset
#' data("subset_T2D.obj")
#'
#' # Second example with subset_T2D.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = subset_T2D.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("shannon"),
#'   subject.var  = "subject_id",
#'   time.var     = "visit_number",
#'   t.level      = "   3",
#'   group.var    = "subject_race",
#'   strata.var   = "subject_gender",
#'   adj.vars     = "sample_body_site",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 20,
#'   pdf.hei      = 8.5
#' )
#'
#' }
#' library(vegan)
#' library(ggh4x)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # First example with peerj32.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("simpson"),
#'   subject.var  = "subject",
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = FALSE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' # Load another dataset
#' data("subset_T2D.obj")
#'
#' # Second example with subset_T2D.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = subset_T2D.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("shannon"),
#'   subject.var  = "subject_id",
#'   time.var     = "visit_number",
#'   t.level      = "   3",
#'   group.var    = "subject_race",
#'   strata.var   = "subject_gender",
#'   adj.vars     = "sample_body_site",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = FALSE,
#'   file.ann     = NULL,
#'   pdf.wid      = 20,
#'   pdf.hei      = 8.5
#' )
#' @export
generate_alpha_boxplot_single <- function (data.obj,
                                           alpha.obj = NULL,
                                           alpha.name = c("shannon",
                                                          "observed_species"),
                                           depth = NULL,
                                           subject.var,
                                           time.var = NULL,
                                           t.level = NULL,
                                           group.var = NULL,
                                           strata.var = NULL,
                                           adj.vars = NULL,
                                           base.size = 16,
                                           theme.choice = "bw",
                                           custom.theme = NULL,
                                           palette = NULL,
                                           pdf = TRUE,
                                           file.ann = NULL,
                                           pdf.wid = 11,
                                           pdf.hei = 8.5,
                                           ...) {

  if (is.null(alpha.name)){
    return()
  }

  if (is.null(alpha.obj)) {
    if (!is_rarefied(data.obj)) {
      message(
        "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj, depth = depth)
    }

    if (!is.null(time.var) & !is.null(t.level)){
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    otu_tab <- data.obj$feature.tab
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
  }

  meta_tab <- data.obj$meta.dat

  # Convert the alpha.obj list to a data frame
  alpha_df <-
    dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    dplyr::inner_join(meta_tab %>% rownames_to_column(var = "sample"),
                      by = c("sample"))

  if (is.null(group.var)) {
    alpha_df <- alpha_df %>% dplyr::mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  theme_function <- switch(
    theme.choice,
    prism = ggprism::theme_prism(),
    classic = theme_classic(),
    gray = theme_gray(),
    bw = theme_bw(),
    ggprism::theme_prism()
  )

  theme_to_use <-
    if (!is.null(custom.theme))
      custom.theme
  else
    theme_function

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

  # Create a plot for each alpha diversity index
  plot_list <- lapply(alpha.name, function(index) {
    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(group.var),
        y = !!sym(index),
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(group.var),
        y = !!sym(index),
        fill = !!sym(time.var)
      )
    }

    if (!is.null(adj.vars)) {
      # 对非数值型协变量进行因子转换
      data_subset <- alpha_df %>%
        dplyr::select(all_of(adj.vars)) %>%
        dplyr::mutate(dplyr::across(where(is.character) & !is.factor, factor))

      # 创建模型矩阵，并为非数值型协变量设定对比度
      M <- model.matrix(
        ~ 0 + .,
        data = data_subset,
        contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
      )

      # Center the covariates (no scaling)
      M_centered <- scale(M, scale = FALSE)

      # Fit the regression model
      fit <- lm(alpha_df[[index]] ~ M_centered)

      # 计算调整后的alpha多样性值
      adjusted_value <- fit$coefficients[1] + residuals(fit)

      # 在alpha_df中更新alpha多样性值
      alpha_df[[index]] <- adjusted_value

      # 显示消息，表示已经为特定的协变量调整了alpha多样性
      message(
        "Alpha diversity has been adjusted for the following covariates: ",
        paste(adj.vars, collapse = ", "),
        "."
      )
    }

    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      y_label <-
        paste0(index, " index (adjusted by: ", covariates, ")")
    } else {
      y_label <- paste0(index, " index")
    }

    boxplot <- ggplot(alpha_df,
                      aes_function) +
      geom_violin(trim = FALSE, alpha = 0.8) +
      geom_jitter(width = 0.1,
                  alpha = 0.5,
                  size = 1) +
      stat_boxplot(geom = "errorbar",
                   position = position_dodge(width = 0.2),
                   width = 0.1) +
      geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.1,
        fill = "white"
      ) +
      scale_fill_manual(values = col) +
      {
        if (!is.null(strata.var) & !is.null(group.var)) {
          ggh4x::facet_nested(
            as.formula(paste(". ~", strata.var, "+", group.var)),
            drop = T,
            scale = "free",
            space = "free"
          )
        } else {
          if (group.var != "ALL") {
            ggh4x::facet_nested(
              as.formula(paste(". ~", group.var)),
              drop = T,
              scale = "free",
              space = "free"
            )
          }
        }
      } +
      labs(y = y_label,
           title = dplyr::if_else(
             !is.null(time.var) &
               !is.null(t.level),
             paste0(time.var, " = ", t.level),
             ""
           ))  +
      theme_to_use +
      theme(
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0, "cm"),
        strip.text = element_text(size = base.size),
        axis.text.y = element_text(color = "black", size = base.size),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = base.size),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 20)
      ) + {
        if (group.var == "ALL") {
          guides(fill = "none")
        }
      }

    return(boxplot)
  })



  # Save the plots as a PDF file
  if (pdf) {
    for (plot_index in seq_along(plot_list)) {
      plot <- plot_list[[plot_index]]
      current_alpha_name <- alpha.name[plot_index]

      pdf_name <- paste0(
        "alpha_boxplot_single_",
        current_alpha_name,
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

      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      print(plot)
      dev.off()
    }
  }

  names(plot_list) <- alpha.name

  return(plot_list)
}
