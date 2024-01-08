#' Visualize Beta Diversity Changes Over Time
#'
#' This function generates a series of boxplots visualizing changes in beta diversity PCoA/PCA coordinates over time.
#' It is designed to work with longitudinal microbiome composition data across multiple time points.
#' The primary output is boxplots showing distributions of coordinate changes between two user-specified time points.
#' These changes can be stratified by groups. Violin plots are overlaid to show density distributions.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param pc.obj A list containing the results of dimension reduction/Principal Component Analysis.
#' This should be the output from functions like \code{\link[MicrobiomeStat]{mStat_calculate_PC}}, containing the PC coordinates and other metadata.
#' If NULL (default), dimension reduction will be automatically performed using metric multidimensional scaling (MDS) via \code{\link[MicrobiomeStat]{mStat_calculate_PC}}.
#' The pc.obj list structure should contain:
#' \itemize{
#'  \item{$points:}{A matrix with samples as rows and PCs as columns containing the coordinates.}
#'  \item{$eig:}{Eigenvalues for each PC dimension.}
#'  \item{$vectors:}{Loadings vectors for features onto each PC.}
#'  \item{Other metadata like $method, $dist.name, etc.}
#' }
#' See \code{\link[MicrobiomeStat]{mStat_calculate_PC}} function for details on output format.
#' @param pc.ind A numeric vector specifying which principal coordinate (PC) axes to use for plotting.
#' This refers to the PC axes from the dimension reduction method specified in pc.obj or calculated by default.
#' For example, c(1,2) will generate plots for PC1 and PC2.
#' Default is c(1,2) to plot the first two PCs.
#' @param subject.var Character string specifying the column name in metadata containing unique subject IDs.
#'                    Required to connect samples from the same subject across timepoints.
#' @param time.var Character string specifying the column name in metadata containing time values for each
#'                sample. Required to identify pairs of timepoints to calculate changes between.
#' @param group.var Character string specifying the column in metadata containing grouping categories.
#'                 Used for stratification in plots. Optional, can be NULL.
#' @param strata.var Character string specifying the column in metadata containing stratification categories.
#'                  Used for nested faceting in plots. Optional, can be NULL.
#' @param adj.vars Character vector specifying columns in metadata containing covariates
#'                to adjust for in distance matrix calculation. Optional, can be NULL.
#' @param change.base The baseline time point value in the time variable to be used
#'                    as the reference for calculating changes. Required if time.var
#'                    contains multiple time points. Changes will be calculated from
#'                    this baseline time point to the other later time point(s). If
#'                    NULL, the first time point will be used as change.base automatically.
#' @param change.func A function or string specifying how to calculate changes between time points.
#' If a function is provided, it should take two arguments "value_time_2" and "value_time_1" representing the PC values at the two time points. The function should return the change value.
#' If a string, currently only "absolute change" is supported, which calculates simple differences between time points.
#' More options could be added in the future as needed. Default is "absolute change".
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param base.size A numeric value for the base size of the plot. Default is 16.
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
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#' @return A named list of ggplot objects, with one element per combination of pc.ind and dist.name.
#' Each element contains the plot for that PC index and distance metric.
#' @details
#' This function generates boxplots visualizing changes in beta diversity PCoA coordinates over time.
#' It is designed for longitudinal microbiome data with multiple time points.
#' The primary output is a boxplot showing distributions of changes in PCoA coordinates between two
#' user-specified time points, optionally stratified by groups. Violin plots are also overlaid to show density.
#' The function offers flexibility to control:
#' - Distance metrics used (via dist.name argument)
#' - PCoA axes to plot (via pc.ind argument)
#' - Subject variable for pairing time points (subject.var)
#' - Time points to compare (time.var and change.base)
#' - Stratification variable(s) (group.var and strata.var)
#' - Calculation of change between time points (change.func argument)
#' - Plot aesthetics like theme, color, file saving, etc.
#' For large datasets, the data are subset to the two time points of interest before plotting.
#' Jitter is added to handle overlapping points. These steps help in generating readable plots.
#' The output plot list allows downstream iteration through multiple PC axes and/or distance metrics.
#' Plots can be accessed via e.g. plotlist$BC_PC1 for BC dissimilarity PC1 coordinates.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' generate_beta_pc_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   change.func = "absolute change",
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
generate_beta_pc_change_boxplot_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           pc.ind = c(1, 2),
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           change.base = NULL,
           change.func = "absolute change",
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

    if (is.null(data.obj) & is.null(dist.obj)) {
      stop("Both data.obj and dist.obj cannot be NULL. Please provide at least one.")
    }

    if (is.null(dist.obj)) {
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      metadata <- data.obj$meta.dat
      if (is.null(metadata)) {
        stop("No metadata could be loaded from data.obj. Please ensure it contains the necessary metadata.")
      }
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
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
        metadata <- data.obj$meta.dat
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

    col <- mStat_get_palette(palette)

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
        custom.theme else theme_function

    plot_list <- lapply(dist.name, function(dist.name) {
        pc.mat <- pc.obj[[dist.name]]$points

        colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

        pc.mat <- pc.mat %>% as_tibble()

        df <-
          cbind(pc.mat[, paste0("PC", pc.ind)], metadata[, c(subject.var, time.var, group.var, strata.var)])

        change.after <-
          unique(df %>% select(all_of(c(time.var))))[unique(df %>% select(all_of(c(time.var)))) != change.base]

        df <-
          df %>%
          as_tibble() %>%
          tidyr::gather(
            key = "PC",
            value = "value",-one_of(subject.var, group.var, time.var, strata.var)
          )

        split_data <-
          split(df, f = df %>%
                  dplyr::group_by(!!sym(time.var)) %>% select(!!sym(time.var)))

        data_time_1 <- split_data[[change.base]]
        data_time_2 <- split_data[[change.after]]

        combined_data <- data_time_1 %>%
          dplyr::inner_join(
            data_time_2,
            by = c("PC", subject.var),
            suffix = c("_time_1", "_time_2")
          )

        combined_data <- combined_data %>%
          dplyr::mutate(value_diff = if (is.function(change.func)) {
            change.func(value_time_2, value_time_1)
          } else if (change.func == "absolute change"){
            value_time_2 - value_time_1
          } else {
            value_time_2 - value_time_1
          })

        combined_data <-
          combined_data %>% dplyr::left_join(metadata %>% select(all_of(
            c(subject.var, time.var, group.var, strata.var)
          )) %>% filter(!!sym(time.var) == change.after),
          by = subject.var)

        if (is.null(group.var)) {
          group.var = "ALL"
          combined_data$ALL <- "ALL"
        }

        sub.plot_list <- lapply(paste0("PC", pc.ind), function(pc.index) {
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


            if (group.var == "ALL") {
              boxplot <- boxplot  + theme(
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                legend.position = "none",
                strip.text.x = element_blank()
              )
            }


          # Save the plots as a PDF file
          if (pdf) {
              pdf_name <- paste0(
                "beta_pc_change_boxplot_pair_",
                dist.name,
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
        names(sub.plot_list) <- paste0("PC", pc.ind)
        return(sub.plot_list)
    })
    names(plot_list) <- dist.name
    return(plot_list)
  }
