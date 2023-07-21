#' @title Generate beta diversity change spaghetti plot for longitudinal data
#'
#' @description This function generates a spaghetti plot visualizing the change in beta diversity over time. The plot is produced for different groups and strata in the given longitudinal data. The function allows users to customize the appearance of the plot and the data processing parameters. If desired, the generated plot can be saved as a PDF file.
#' @name generate_beta_ordination_single
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj (Optional) An object storing the pre-calculated distance matrices. If provided, the function will not compute the distances again. This can save a lot of computational time when working with large datasets.
#' @param subject.var A character string specifying the name of the variable that indicates the subject in the metadata table.
#' @param time.var A character string specifying the name of the variable that indicates the time in the metadata table.
#' @param t0.level (Optional) The baseline time point level. If not provided, the function will infer it from the data.
#' @param ts.levels (Optional) A vector specifying the levels of the time points to be considered. If not provided, the function will infer it from the data.
#' @param group.var (Optional) A character string specifying the name of the variable that indicates the grouping factor in the metadata table.
#' @param strata.var (Optional) A character string specifying the name of the variable that indicates the stratification factor in the metadata table.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#' @param base.size (Optional) A numeric value specifying the base size of the plot elements. The default is 16.
#' @param theme.choice (Optional) A character string specifying the name of the ggplot2 theme to be used for the plot. The default is "bw".
#' @param custom.theme (Optional) A custom ggplot2 theme. If provided, this theme will be used instead of the theme specified by `theme.choice`.
#' @param palette (Optional) A character vector specifying the color palette to be used for the plot. If not provided, a default color palette will be used.
#' @param pdf (Optional) A logical value indicating whether to save the plot as a PDF file. The default is TRUE.
#' @param file.ann (Optional) A character string to be added as an annotation to the file name of the output PDF file.
#' @param pdf.wid (Optional) A numeric value specifying the width of the output PDF file. The default is 11.
#' @param pdf.hei (Optional) A numeric value specifying the height of the output PDF file. The default is 8.5.
#' @param ... (Optional) Additional arguments to pass to the `ggplot()` function.
#'
#' @details The function is flexible and allows for various modifications, including the choice of distance measure and stratification factor, providing a comprehensive tool for microbiome beta diversity exploration. It integrates well with other MicrobiomeStat functions and takes their output as input.
#'
#' @return A PCoA plot displaying the beta diversity ordination, stratified by the specified grouping and/or strata variables (if provided). The plot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object, \code{\link[MicrobiomeStat]{mStat_calculate_PC}} for computing the principal coordinates, and \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{geom_boxplot}} for the underlying plot functions used, and \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_DESeqDataSet_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_phyloseq_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_SummarizedExperiment_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_qiime2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_mothur_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_dada2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_biom_as_data_obj}} for data conversion.
#'
#' @author Caffery Yang \email{cafferychen7850@@gmail.com}
#'
#' @examples
#' # Load required libraries and example data
#' library(microbiome)
#' library(tidyverse)
#' library(vegan)
#' data(peerj32)
#' peerj32.obj <- list()
#' peerj32.phy <- peerj32$phyloseq
#' peerj32.obj <- mStat_convert_phyloseq_to_data_obj(peerj32.phy)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC','Jaccard'))
#' # Generate the boxplot pair
#' generate_beta_ordination_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t.level = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' @export
generate_beta_ordination_single <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           subject.var,
           time.var = NULL,
           t.level = NULL,
           group.var = NULL,
           strata.var = NULL,
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

    if (is.null(dist.obj)) {
      if (!is.null(time.var)){
        if (!is.null(t.level)){
          metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,strata.var,time.var))) %>% filter(!!sym(time.var) == t.level)
          data.obj <- update_data_obj_count(data.obj,load_data_obj_count(data.obj)[,rownames(metadata)])
          dist.obj <-
            mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        } else {
          if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
            message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
          }
          metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,strata.var,time.var)))
          dist.obj <-
            mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        }
      } else {
        metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,strata.var,time.var)))
        dist.obj <-
          mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      }
    } else {
      if (!is.null(data.obj)){
        if (!is.null(time.var)){
          if (!is.null(t.level)){
            metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,strata.var,time.var))) %>% filter(!!sym(time.var) == t.level)
            data.obj <- update_data_obj_count(data.obj,load_data_obj_count(data.obj)[,rownames(metadata)])
          } else {
            if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
              message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
            }
            metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,strata.var,time.var)))
          }
        } else {
          metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,strata.var,time.var)))
        }
      }
      if (!is.null(attr(dist.obj[[dist.name[1]]], "labels"))){
        metadata <- attr(dist.obj[[dist.name[1]]], "labels")  %>% select(all_of(c(subject.var,group.var,strata.var,time.var)))
      }
    }

    if (is.null(pc.obj)) {
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = 2,
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

    if (is.null(group.var)){
      group.var = "ALL"
      metadata$ALL <- "ALL"
    }

    aes_function <- if (!is.null(strata.var)) {
      aes(color = !!sym(group.var),
          shape = !!sym(strata.var))
    } else {
      aes(color = !!sym(group.var))
    }

    theme_function <- switch(theme.choice,
                             prism = ggprism::theme_prism(),
                             classic = theme_classic(),
                             gray = theme_gray(),
                             bw = theme_bw(),
                             ggprism::theme_prism()) # 根据用户选择设置主题

    # 使用用户自定义主题（如果提供），否则使用默认主题
    theme_to_use <- if (!is.null(custom.theme)) custom.theme else theme_function

    plot_list <- lapply(dist.name, function(dist.name) {
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]
      df <-
        cbind(pc.mat, metadata[, c(subject.var, time.var, group.var, strata.var)])
      colnames(df)[1:2] <- c("PC1", "PC2")

      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(size = 12, aes_function, show.legend = T) +
        ggplot2::labs(
          x = ifelse(!is.null(pc.obj[[dist.name]]$eig),paste0("Axis 1 (", round(pc.obj[[dist.name]]$eig[1]/sum(pc.obj[["BC"]]$eig)*100,2),"%)"),"Axis 1"),
          y = ifelse(!is.null(pc.obj[[dist.name]]$eig),paste0("Axis 2 (", round(pc.obj[[dist.name]]$eig[2]/sum(pc.obj[["BC"]]$eig)*100,2),"%)"),"Axis 2")
        ) +
        ggplot2::stat_ellipse(ggplot2::aes(color = !!sym(group.var)),fill="white",geom = "polygon",
                              level=0.95,alpha = 0.01,show.legend = F) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed",
          color = "black"
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "black"
        ) +
        theme_to_use  +
        ggplot2::theme(
          axis.line.x = ggplot2::element_line(size = 1, colour = "black"),
          axis.line.y = ggplot2::element_line(size = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title = ggplot2::element_text(color = "black", size = 20),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(
            color = "black",
            size = base.size
          ),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = ggplot2::element_text(size = base.size),
          legend.title = ggplot2::element_text(size = base.size)
        )

      if (group.var == "ALL") {
        p <- p + scale_color_manual(values = col, guide = "none")
      } else {
        p <- p + scale_color_manual(values = col)
      }

      # Create a ggplot object for the bar plot of PC1
      Fig1a.taxa.pc1.boxplot <-
        ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x=!!sym(group.var), y=PC1, fill=!!sym(group.var)), color="black", alpha=0.5, show.legend = F) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand = c(0,0.001)) +
        ggplot2::labs(x=NULL, y=NULL) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.text.y =ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())+
        ggplot2::coord_flip()

      # Create a ggplot object for the bar plot of PC2
      Fig1a.taxa.pc2.boxplot <-
        ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x=!!sym(group.var), y=PC2, fill=!!sym(group.var)), color="black", alpha=0.5, show.legend = F) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand = c(0,0.001)) +
        ggplot2::labs(x=NULL, y=NULL) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.text.y =ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())

      # Combine the two plots into a single plot with the PC1 plot on top and the PC2 plot on the right
      p <- p %>%
        aplot::insert_top(Fig1a.taxa.pc1.boxplot, height = 0.2) %>%
        aplot::insert_right(Fig1a.taxa.pc2.boxplot, width=0.2) %>%
        ggplotify::as.ggplot()

      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "beta_ordination_single_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "dist.name_",
          dist.name
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
