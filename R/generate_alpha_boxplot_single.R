#' @title Generate boxplot for alpha diversity index at a single time point
#'
#' @description This function generates a boxplot of specified alpha diversity indices at a single time point across different groupings, with optional stratification. The output can be saved as a PDF.
#' @name generate_alpha_boxplot_single
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t.level The specific time point to be analyzed.
#' @param group.var An optional variable in the metadata table that represents the grouping factor.
#' @param strata.var An optional variable in the metadata table that represents the stratification factor.
#' @param base.size The base font size for the plot. Default is 16.
#' @param theme.choice A character string indicating the choice of pre-defined ggplot2 theme for the plot. Supported choices include "prism" (default), "classic", "gray", and "bw".
#' @param custom.theme An optional custom ggplot2 theme. If provided, this theme will be used instead of the pre-defined themes.
#' @param palette An optional palette to use for the plot. If NULL (default), a pre-defined palette will be used.
#' @param pdf A boolean indicating whether to save the output as a PDF file. Default is TRUE.
#' @param file.ann A string for annotating the output file name.
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @return A list of boxplots displaying the specified alpha diversity indices at the specified time point across different groupings, stratified by the specified stratification variable (if provided). Each boxplot in the list corresponds to one of the alpha diversity indices specified in `alpha.name`. The boxplots will be saved as PDF files if `pdf` is set to `TRUE`.
#'
#' @examples
#' library(microbiome)
#' library(vegan)
#' library(ggh4x)
#' data(peerj32)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#'
#' plot_alpha_list <- generate_alpha_boxplot_single(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("simpson"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   t.level = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   base.size = 16,
#'   theme.choice = "classic",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5)
#'
#'
#' # Display the boxplot
#' print(plot_alpha_list)
#'
#' @export
generate_alpha_boxplot_single <- function (data.obj,
                                         alpha.obj = NULL,
                                         alpha.name = c("shannon",
                                                        "simpson",
                                                        "observed_species",
                                                        "chao1",
                                                        "ace",
                                                        "pielou"),
                                         subject.var,
                                         time.var =NULL,
                                         t.level = NULL,
                                         group.var = NULL,
                                         strata.var = NULL,
                                         base.size = 16,
                                         theme.choice = "prism",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         ...) {

  if (is.null(alpha.obj)) {
    if (!is_rarefied(data.obj)){
      message(
        "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj)
    }
    otu_tab <- as.data.frame(load_data_obj_count(data.obj))

    if (!is.null(time.var)){
      if (!is.null(t.level)){
        meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(
          c(subject.var, time.var, group.var, strata.var))) %>% filter(!!sym(time.var) == t.level)
      } else {
        meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(
          c(subject.var, time.var, group.var, strata.var)))
        if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
          message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
        }
      }
    } else {
      meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(
        c(subject.var, group.var, strata.var)))
    }

    otu_tab <- otu_tab %>% select(all_of(c(rownames(meta_tab))))
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
  }

  # Convert the alpha.obj list to a data frame
  alpha_df <-
    bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    inner_join(
      meta_tab %>% rownames_to_column(var = "sample"),
      by = c("sample")
    )

  if (is.null(group.var)){
    alpha_df <- alpha_df %>% mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  theme_function <- switch(theme.choice,
                           prism = ggprism::theme_prism(),
                           classic = theme_classic(),
                           gray = theme_gray(),
                           bw = theme_bw(),
                           ggprism::theme_prism())

  theme_to_use <- if (!is.null(custom.theme)) custom.theme else theme_function

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

    boxplot <- ggplot(alpha_df,
                      aes_function) +
      geom_violin(trim = FALSE, alpha = 0.8) +
      geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
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
        if (!is.null(strata.var) & !is.null(group.var)){
          facet_nested(as.formula(paste(". ~", group.var, "+", strata.var)), drop = T, scale = "free", space = "free")
        } else {
          if (group.var != "ALL"){
            facet_nested(as.formula(paste(". ~", group.var)), drop = T, scale = "free", space = "free")
          }
        }
      } +
      labs(y = index,
           title = if_else(!is.null(time.var) & !is.null(t.level),paste0(time.var," = ", t.level), ""))  +
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
        plot.title = element_text(hjust = 0.5,size = 20)
      ) + {
        if (group.var == "ALL"){
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

      pdf_name <- paste0("alpha_diversity_boxplot_single_",
                         current_alpha_name,
                         "_",
                         "subject_", subject.var,
                         "_",
                         "time_", time.var)

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

  return(plot_list)
}
