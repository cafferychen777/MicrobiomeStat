#' Generate boxplot for specified alpha diversity index
#'
#' This function generates a boxplot of a specified alpha diversity index across different groupings and time points, with optional stratification. The output can be saved as a PDF.
#' @name generate_alpha_boxplot_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t0.level The level in the metadata table that represents the initial time point. This parameter is used to determine the starting point for time series data. If NULL, the first unique time point in the data will be used.
#' @param ts.levels The levels in the metadata table that represent the subsequent time points.
#' These parameters are used to determine the time points that follow the initial time point (t0.level) for time series data. If NULL, all unique time points in the data (except for the initial time point) will be used.
#' @param palette An optional color palette for the plot. If not provided, a default color palette will be used. The palette should be a vector of color codes in a format accepted by ggplot2 (e.g., hexadecimal color codes). The number of colors in the palette should be at least as large as the number of groups being plotted.
#' @param group.var An optional variable in the metadata table that represents the grouping factor.
#' @param strata.var An optional variable in the metadata table that represents the stratification factor.
#' @param base.size The base font size for the plot.
#' @param theme.choice A character string indicating the choice of pre-defined ggplot2 theme for the plot. Supported choices include "prism" (default), "classic", "gray", and "bw".
#' @param custom.theme An optional custom ggplot2 theme. If provided, this theme will be used instead of the pre-defined themes.
#' @param pdf A boolean indicating whether to save the output as a PDF file.
#' @param file.ann A string for annotating the output file name.
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @return A boxplot displaying the specified alpha diversity index across different groupings and time points, stratified by the specified stratification variable (if provided). The boxplot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @examples
#'
#' data("subset_T2D.obj")
#'
#' generate_alpha_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon"),
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = sort(unique(subset_T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = sort(unique(subset_T2D.obj$meta.dat$visit_number))[2:6],
#'   group.var = "subject_race",
#'   strata.var = NULL,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#'
#' data(peerj32.obj)
#' generate_alpha_boxplot_long(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("simpson"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = NULL,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5)
#'
#' @export
#'
generate_alpha_boxplot_long <- function (data.obj,
                                         alpha.obj = NULL,
                                         alpha.name = c("shannon",
                                                        "simpson",
                                                        "observed_species",
                                                        "chao1",
                                                        "ace",
                                                        "pielou"),
                                         subject.var,
                                         time.var,
                                         t0.level = NULL,
                                         ts.levels = NULL,
                                         group.var = NULL,
                                         strata.var = NULL,
                                         base.size = 12,
                                         theme.choice = "prism",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         ...) {
  # Data validation
  if (!is(data.obj, "list"))
    stop("`data.obj` should be a list.")
  if (!is.null(alpha.obj) &&
      !is(alpha.obj, "list"))
    stop("`alpha.obj` should be a list or NULL.")
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

  if (is.null(alpha.obj)) {
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    if (!is_rarefied(data.obj)) {
      message(
        "Diversity analysis needs rarefaction! Call \"mStat_rarefy_data\" to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj)
    }
    otu_tab <- load_data_obj_count(data.obj)
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
  } else {
    # Verify that all alpha.name are present in alpha.obj
    if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x)
      colnames(x))))) {
      missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
      stop(
        "The following alpha diversity indices are not available in alpha.obj: ",
        paste(missing_alphas, collapse = ", "),
        call. = FALSE
      )
    }
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
  }

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(
      subject.var, group.var, time.var, strata.var
    )))

  # Convert the alpha.obj list to a data frame
  alpha_df <-
    bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    inner_join(meta_tab %>% rownames_to_column(var = "sample"),
               by = c("sample"))

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

  if (is.null(group.var)) {
    alpha_df <- alpha_df %>% mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  # Create a plot for each alpha diversity index
  plot_list <- lapply(alpha.name, function(index) {
    line_aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        group = !!sym(subject.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        group = !!sym(subject.var)
      )
    }

    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        fill = !!sym(time.var)
      )
    }

    average_alpha_df <- NULL
    if (length(unique(alpha_df[[time.var]])) > 10 ||
        length(unique(alpha_df[[subject.var]])) > 25) {
      if (!is.null(group.var)) {
        average_alpha_df <- alpha_df %>%
          group_by(!!sym(group.var),!!sym(time.var)) %>%
          dplyr::summarise(across(!!sym(index), mean, na.rm = TRUE), .groups = "drop") %>%
          ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL")
      } else {
        average_alpha_df <- alpha_df %>%
          group_by(!!sym(time.var)) %>%
          dplyr::summarise(across(!!sym(index), mean, na.rm = TRUE), .groups = "drop") %>%
          ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL")
      }
    }

    boxplot <- ggplot(alpha_df,
                      aes_function) +
      geom_violin(trim = FALSE, alpha = 0.8) +
      stat_boxplot(geom = "errorbar",
                   position = position_dodge(width = 0.2),
                   width = 0.1) +
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
        linetype = "dashed",
        data = if (!is.null(average_alpha_df))
          average_alpha_df
        else
          alpha_df
      ) +
      scale_fill_manual(values = col) +
      {
        if (!is.null(strata.var) & !is.null(group.var)) {
          ggh4x::facet_nested(
            as.formula(paste(". ~", group.var, "+", strata.var)),
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
      labs(x = time.var,
           y = paste0(index, " index"))  +
      theme_to_use +
      theme(
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0, "cm"),
        strip.text.x = element_text(size = 15, color = "black"),
        axis.text.x = element_text(
          angle = 90,
          color = "black",
          vjust = 0.5,
          size = base.size * 0.75
        ),
        axis.text.y = element_text(color = "black", size = base.size),
        axis.title.x = element_text(size = base.size),
        axis.title.y = element_text(size = base.size),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16)
      ) + {
        if (group.var == "ALL") {
          guides(fill = "none")
        }
      }

    # Add geom_jitter() if the number of unique time points or subjects is greater than 10
    if (length(unique(alpha_df[[time.var]])) > 10 ||
        length(unique(alpha_df[[subject.var]])) > 10) {
      boxplot <- boxplot + geom_jitter(width = 0.1,
                                       alpha = 0.1,
                                       size = 1)
    }

    return(boxplot)
  })



  # Save the plots as a PDF file
  if (pdf) {
    for (plot_index in seq_along(plot_list)) {
      plot <- plot_list[[plot_index]]
      current_alpha_name <- alpha.name[plot_index]

      pdf_name <- paste0(
        "alpha_boxplot_long_",
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

  return(plot_list)
}
