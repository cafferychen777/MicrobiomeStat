#' Generate an alpha diversity line plot for longitudinal data
#'
#' This function creates a ggplot object of alpha diversity (e.g., Shannon index) line plot for longitudinal data,
#' showing individual subject trajectories and the mean trajectory for each group.
#' @name generate_alpha_spaghettiplot_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be calculated (e.g., "Shannon").
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param t0.level The level in the metadata table that represents the initial time point. This is used to determine the starting point for time series data.
#' @param ts.levels The levels in the metadata table that represent the subsequent time points. These are used to determine the time points that follow the initial time point for time series data.
#' @param base.size The base font size for the plot.
#' @param palette An optional color palette for the plot. If not provided, a default color palette will be used. The palette should be a vector of color codes in a format accepted by ggplot2 (e.g., hexadecimal color codes). The number of colors in the palette should be at least as large as the number of groups being plotted.
#' @param theme.choice A character string indicating the choice of pre-defined ggplot2 theme for the plot. Supported choices include "prism" (default), "classic", "gray", and "bw".
#' @param custom.theme An optional custom ggplot2 theme. If provided, this theme will be used instead of the pre-defined themes.
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param subject.var The name of the subject variable.
#' @param time.var The name of the time variable.
#' @param group.var The name of the group variable.
#' @param strata.var The name of the strata variable (default is NULL).
#' @param pdf Logical, whether to save the plot as a PDF (default is TRUE).
#' @param file.ann The annotation to be added to the PDF file name (default is NULL).
#' @param ... Additional arguments passed to ggplot().
#'
#' @return A ggplot object of the alpha diversity line plot.
#'
#' @examples
#' \dontrun{
#' library("HMP2Data")
#' T2D <- T2D16S()
#'
#' T2D.obj <- mStat_convert_phyloseq_to_data_obj(T2D)
#'
#' subset_T2D.obj <- mStat_subset_data(T2D.obj,colnames(T2D.obj$feature.tab
#' [,colSums(T2D.obj$feature.tab) >= 2000]))
#' T2D.alpha.obj <- mStat_calculate_alpha_diversity(subset_T2D.obj$feature.tab,"shannon")
#'
#' generate_alpha_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = T2D.alpha.obj,
#'   alpha.name = c("shannon"),
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = sort(unique(T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = sort(unique(T2D.obj$meta.dat$visit_number))[-1],
#'   group.var = "subject_race",
#'   strata.var = NULL,
#'   theme.choice = "bw",
#'   palette = ggsci::pal_npg()(9),
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' print(plot)
#' }
#' @export
generate_alpha_spaghettiplot_long <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = c("shannon",
                          "simpson",
                          "observed_species",
                          "chao1",
                          "ace",
                          "pielou"),
           subject.var,
           time.var,
           t0.level,
           ts.levels,
           group.var = NULL,
           strata.var = NULL,
           base.size = 16,
           palette = NULL,
           theme.choice = "bw",
           custom.theme = NULL,
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
          "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj)
      }

      otu_tab <- load_data_obj_count(data.obj)

      alpha.obj <-
        mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
    } else {
      # Verify that all alpha.name are present in alpha.obj
      if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x) colnames(x))))) {
        missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
        stop("The following alpha diversity indices are not available in alpha.obj: ",
             paste(missing_alphas, collapse = ", "), call. = FALSE)
      }
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    }

    meta_tab <-
      load_data_obj_metadata(data.obj)

    # Convert the alpha.obj list to a data frame
    alpha.df <-
      bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
      inner_join(
        meta_tab %>% select(all_of(
          c(subject.var, time.var, group.var, strata.var)
        )) %>% rownames_to_column(var = "sample"),
        by = c("sample")
      )

    if (is.null(group.var)){
      alpha.df <- alpha.df %>% mutate("ALL" = "ALL")
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

    # Calculate new sizes based on base.size
    title.size = base.size * 1.25
    axis.title.size = base.size * 1
    axis.text.size = base.size * 0.5
    legend.title.size = base.size * 1
    legend.text.size = base.size * 0.75

    # Create a plot for each alpha diversity index
    plot_list <- lapply(alpha.name, function(index) {
      sub_alpha.df <- alpha.df %>% select(all_of(c(index,subject.var, time.var, group.var, strata.var)))

      if (is.null(strata.var)) {
        sub_alpha.df.mean <- sub_alpha.df %>%
          group_by(!!sym(time.var), !!sym(group.var)) %>%
          summarize(mean_alpha = mean(!!sym(index), na.rm = TRUE))
        sub_alpha.df <-
          left_join(sub_alpha.df, sub_alpha.df.mean, by = c(time.var, group.var))
      } else {
        sub_alpha.df.mean <- sub_alpha.df %>%
          group_by(!!sym(time.var), !!sym(group.var), !!sym(strata.var)) %>%
          summarize(mean_alpha = mean(!!sym(index), na.rm = TRUE))
        sub_alpha.df <-
          left_join(sub_alpha.df, sub_alpha.df.mean, by = c(time.var, group.var, strata.var))
      }

      # create a vector of time points
      time.points <- c(t0.level, ts.levels)

      # create a sequence of indices for time points to be shown
      if (length(time.points) > 80) {
        indices <- round(seq(1, length(time.points), length.out = 80))
        message("There are more than 80 time points, so we are selecting a subset of 80 to display on the x-axis. This does not affect any calculations or the resulting spaghetti plot.")
      } else {
        indices <- 1:length(time.points)
      }

      # create breaks and labels using the selected indices
      breaks <- time.points[indices]
      labels <- time.points[indices]

      plot <- ggplot() +
        geom_line(
          data = sub_alpha.df,
          aes_string(
            x = time.var,
            y = index,
            group = subject.var,
            color = group.var
          ),
          alpha = 0.5
        ) +
        geom_line(
          data = sub_alpha.df,
          aes_string(
            x = time.var,
            y = "mean_alpha",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        labs(x = time.var, y = index, color = group.var) +
        scale_color_manual(
          values = col
        ) +
        scale_x_discrete(breaks = breaks, labels = labels) +
        theme_to_use +
        theme(
          plot.title = element_text(
            size = title.size,
            face = "bold",
            hjust = 0.5
          ),
          axis.title.x = element_text(size = axis.title.size, face = "bold"),
          axis.title.y = element_text(size = axis.title.size, face = "bold"),
          axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = axis.text.size),
          axis.text.y = element_text(size = axis.text.size),
          legend.title = element_text(size = legend.title.size, face = "bold"),
          legend.text = element_text(size = legend.text.size)
        )

      if (!is.null(strata.var)) {
        plot <- plot + facet_wrap(as.formula(paste0("~", strata.var)))
      }

      if (pdf) {
        pdf_name <- paste0(
          "alpha_spaghettiplot_long",
          index,
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "group_",
          group.var
        )

        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
        }

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }

        pdf_name <- paste0(pdf_name, ".pdf")

        ggsave(
          filename = pdf_name,
          plot = plot,
          width = pdf.wid,
          height = pdf.hei
        )
      }
      return(plot)
    })

    return(plot_list)
  }
