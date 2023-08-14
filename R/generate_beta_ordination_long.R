#' Generate Beta Ordination Plots
#'
#' This function generates beta ordination plots using the specified distance measure.
#' The plots can be stratified by a given variable and the color and shape aesthetics
#' can be mapped to groups and time variables, respectively. The function also supports
#' saving the plots as PDF files.
#'
#' @param data.obj a list with the data to be used for the plots. If NULL, the distance object will be used instead.
#' @param dist.obj a list with distance matrices. If NULL, it will be calculated from the data object.
#' @param pc.obj a list with principal coordinates. If NULL, it will be calculated from the distance object.
#' @param subject.var a character string specifying the subject variable.
#' @param time.var a character string specifying the time variable.
#' @param t0.level a character string specifying the baseline time level.
#' @param ts.levels a vector of character strings specifying the time series levels.
#' @param group.var a character string specifying the grouping variable. If NULL, only color will be used in the plots.
#' @param strata.var a character string specifying the stratification variable. If NULL, no stratification will be done.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param dist.name a character vector specifying the distance measures to use. Defaults to c('BC', 'Jaccard').
#' @param base.size a numeric value specifying the base size for the plot text.
#' @param theme.choice a character string specifying the theme to use. Defaults to "prism".
#' @param custom.theme a ggplot2 theme object to use instead of the default theme.
#' @param palette a character vector specifying the colors to use for the plots. If NULL, a default palette will be used.
#' @param pdf a logical value indicating whether to save the plots as PDF files. Defaults to TRUE.
#' @param file.ann a character string specifying an annotation to add to the file names of the saved plots.
#' @param pdf.wid a numeric value specifying the width of the saved PDF files.
#' @param pdf.hei a numeric value specifying the height of the saved PDF files.
#' @param ... further arguments to be passed to the underlying functions.
#'
#' @return A list of ggplot2 objects representing the beta ordination plots.
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object, \code{\link[MicrobiomeStat]{mStat_calculate_PC}} for computing the principal coordinates, and \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{geom_boxplot}} for the underlying plot functions used, and \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_DESeqDataSet_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_phyloseq_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_SummarizedExperiment_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_qiime2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_mothur_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_dada2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_biom_as_data_obj}} for data conversion.
#'
#' @author Caffery Yang \email{cafferychen7850@@gmail.com}
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = sort(unique(subset_T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = sort(unique(subset_T2D.obj$meta.dat$visit_number))[2:6],
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   adj.vars = "sample_body_site",
#'   dist.name = c("BC","Jaccard"),
#'   base.size = 16,
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
generate_beta_ordination_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           subject.var,
           time.var,
           t0.level,
           ts.levels,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
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
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
      } else {
        metadata <- attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
        data.obj <- list(meta.dat = metadata)
        data.obj <- mStat_process_time_variable(metadata, time.var, t0.level, ts.levels)
        metadata <- load_data_obj_metadata(data.obj)
      }
    }

    if (is.null(pc.obj)) {
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "nmds",
          k = 2,
          dist.name = dist.name
        )
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

    aes_function <- if (!is.null(group.var)) {
      aes(color = !!sym(group.var),
          shape = !!sym(time.var))
    } else {
      aes(color = !!sym(time.var))
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

    plot_list <- lapply(dist.name, function(dist.name) {
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]
      df <- as.data.frame(pc.mat) %>%
        setNames(c("PC1", "PC2")) %>%
        dplyr::bind_cols(metadata[, c(subject.var, time.var, group.var, strata.var)]) %>%
        dplyr::mutate(x_start = PC1,
               y_start = PC2,
               x_end = NA,
               y_end = NA)
      Time_choices <-
        df %>% select(all_of(time.var)) %>% dplyr::pull() %>% unique()
      df <- df %>%
        dplyr::group_by(.data[[subject.var]]) %>%
        dplyr::mutate(x_end = dplyr::if_else(.data[[time.var]] == max(levels(Time_choices)), NA_real_, dplyr::lead(PC1)),
               y_end = dplyr::if_else(.data[[time.var]] == max(levels(Time_choices)), NA_real_, dplyr::lead(PC2))) %>%
        dplyr::ungroup()
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(
          size = 15,
          aes_function,
          show.legend = T,
          alpha = 0.8
        ) +
        ggplot2::geom_segment(
          aes(
            x = x_start,
            y = y_start,
            xend = x_end,
            yend = y_end
          ),
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
          size = 1,
          color = "black"
        ) +
        ggplot2::labs(
          x = ifelse(
            !is.null(pc.obj[[dist.name]]$eig),
            paste0("Axis 1 (", round(
              pc.obj[[dist.name]]$eig[1] / sum(pc.obj[["BC"]]$eig) * 100, 2
            ), "%)"),
            "Axis 1"
          ),
          y = ifelse(
            !is.null(pc.obj[[dist.name]]$eig),
            paste0("Axis 2 (", round(
              pc.obj[[dist.name]]$eig[2] / sum(pc.obj[["BC"]]$eig) * 100, 2
            ), "%)"),
            "Axis 2"
          )
        ) +
        scale_color_manual(values = col) +
        ggplot2::geom_vline(xintercept = 0,
                            linetype = "dashed",
                            color = "black") +
        ggplot2::geom_hline(yintercept = 0,
                            linetype = "dashed",
                            color = "black") +
        theme_to_use  +
        ggplot2::theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title = ggplot2::element_text(color = "black"),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color = "black",
                                            size = base.size),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(as.formula(paste(".~", strata.var)))
      }
      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "beta_ordination_long_",
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
