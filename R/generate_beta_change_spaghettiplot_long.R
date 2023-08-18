#' @title Generate beta diversity change spaghetti plot for longitudinal data
#'
#' @description This function generates a spaghetti plot visualizing the change in beta diversity over time for different groups and strata in longitudinal data. Optionally, a user can specify whether to save the plot as a PDF file.
#' @name generate_beta_change_spaghettiplot_long
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj (Optional) A distance object generated from distance matrices using 'mStat_calculate_beta_diversity' function on data.obj. If data.obj is not provided, metadata can be retrieved from dist.obj.
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param time.var The variable in the metadata table that represents the time.
#' @param t0.level (Optional) The baseline time point level.
#' @param ts.levels (Optional) Time point levels.
#' @param group.var (Optional) The variable in the metadata table that represents the grouping factor.
#' @param adj.vars A character vector containing the names of the columns in data.obj$meta.dat to include as covariates in the PERMANOVA analysis. If no covariates are needed, use NULL (default).
#' @param strata.var (Optional) The variable in the metadata table that represents the stratification factor.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned.
#' @param base.size (Optional) Base font size for the plot (default is 16).
#' @param theme.choice (Optional) Name of the theme for the plot. Default is "bw". Other options include "plain", "classic", and any other themes compatible with ggplot2.
#' @param custom.theme (Optional) A custom ggplot2 theme.
#' @param palette (Optional) A palette function or character vector with the colors for the plot.
#' @param pdf (Optional) A boolean indicating whether to save the output as a PDF file (default is TRUE).
#' @param file.ann (Optional) A string for annotating the output file name.
#' @param pdf.wid (Optional) The width of the PDF file if `pdf` is set to `TRUE` (default is 11).
#' @param pdf.hei (Optional) The height of the PDF file if `pdf` is set to `TRUE` (default is 8.5).
#' @param ... (Optional) Additional arguments to pass to the plotting function.
#'
#' @return A spaghetti plot displaying the beta diversity change over time, stratified by the specified grouping and/or strata variables (if provided). The plot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#'
#' generate_beta_change_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_beta_change_spaghettiplot_long <-
  function(data.obj,
           dist.obj = NULL,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           dist.name = c("BC"),
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    # Data validation
    mStat_validate_data(data.obj)
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

    theme_function <- switch(theme.choice,
                             prism = ggprism::theme_prism(),
                             classic = theme_classic(),
                             gray = theme_gray(),
                             bw = theme_bw(),
                             ggprism::theme_prism())

    theme_to_use <- if (!is.null(custom.theme)) custom.theme else theme_function

    # Calculate new sizes based on base.size
    title.size = base.size * 1.25
    axis.title.size = base.size * 0.75
    axis.text.size = base.size * 0.5
    legend.title.size = base.size * 1
    legend.text.size = base.size * 0.75

    plot_list <- lapply(dist.name,function(dist.name){
      if (is.null(dist.obj)) {
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
        dist.obj <-
          mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        if (!is.null(adj.vars)){
          dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
        }
      } else {
        if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
          data.obj <-
            mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
          meta_tab <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
        } else {
          meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(subject.var, time.var, group.var, strata.var)))
          data.obj <- list(meta.dat = meta_tab)
          data.obj <- mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
          meta_tab <- load_data_obj_metadata(data.obj)
        }
      }

      if (is.null(dist.obj[[dist.name]])) {
        message(paste("dist.obj does not contain", dist.name))
      } else {
        tryCatch({
          dist.df <- as.matrix(dist.obj[[dist.name]])
        }, error = function(e) {
          message(paste("Failed to convert dist.obj[[", dist.name, "]] to a matrix: ", e$message))
        })

        missing_rows <- setdiff(rownames(meta_tab), rownames(dist.df))
        missing_cols <- setdiff(rownames(meta_tab), colnames(dist.df))

        if (length(missing_rows) > 0) {
          message(paste("The following rows in meta_tab are not in dist.df:", paste(missing_rows, collapse = ", ")))
        }

        if (length(missing_cols) > 0) {
          message(paste("The following columns in meta_tab are not in dist.df:", paste(missing_cols, collapse = ", ")))
        }

        # 获取在 meta_tab 的 rownames 中的行和列
        if (length(missing_rows) == 0 & length(missing_cols) == 0) {
          dist.df <- dist.df[rownames(meta_tab), rownames(meta_tab)]
        }
      }

      dist.df <- dist.df %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      meta_tab <- meta_tab %>% rownames_to_column("sample")

      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == levels(meta_tab[,time.var])[1]) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      if (!is.null(strata.var)&!is.null(group.var)){
        long.df <- long.df %>% dplyr::left_join(meta_tab %>% select(-all_of("sample")) %>% dplyr::distinct(),by = c(subject.var,time.var))
      } else if (is.null(strata.var)&!is.null(group.var)){
        long.df <- long.df %>% dplyr::left_join(meta_tab %>% select(-all_of("sample")) %>% dplyr::distinct(),by = c(subject.var,time.var))
      } else {
        long.df <- long.df
      }

      if (is.null(group.var)){
        long.df <- long.df %>% dplyr::mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      if (is.null(strata.var)) {
        long.df.mean <- long.df %>%
          dplyr::group_by(!!sym(time.var),!!sym(group.var)) %>%
          dplyr::summarize(mean_distance = mean(distance, na.rm = TRUE))
        long.df <-
          dplyr::left_join(long.df, long.df.mean, by = c(time.var, group.var))
      } else {
        long.df.mean <- long.df %>%
          dplyr::group_by(!!sym(time.var),!!sym(group.var),!!sym(strata.var)) %>%
          dplyr::summarize(mean_distance = mean(distance, na.rm = TRUE))
        long.df <-
          dplyr::left_join(long.df, long.df.mean, by = c(time.var, group.var, strata.var))
      }

      long.df <- long.df %>% dplyr::arrange(subject.var,time.var)

      # summary.df <- long.df %>%
      #   dplyr::group_by(!!sym(time.var), !!sym(group.var)) %>%
      #   dplyr::summarise(mean_distance = mean(distance),
      #             sd_distance = sd(distance),
      #             lower = mean_distance - sd_distance,
      #             upper = mean_distance + sd_distance)

      if (is.null(adj.vars)) {
        y_label <- paste(dist.name, "Distance from Baseline")
      } else {
        y_label <- paste(dist.name, "Distance from Baseline (adjusted by:", paste(adj.vars, collapse = ", "), ")")
      }

      p <- ggplot() +
        geom_line(
          data = long.df,
          aes_string(
            x = time.var,
            y = "distance",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.3
        ) +
        geom_line(
          data = long.df,
          aes_string(
            x = time.var,
            y = "mean_distance",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        # geom_errorbar(
        #   data = summary.df,
        #   aes_string(
        #     x = time.var,
        #     ymin = "lower",
        #     ymax = "upper",
        #     color = group.var
        #   ),
        #   width = 0.2
        # ) +
        labs(x = time.var, y = y_label, color = group.var) +
        scale_color_manual(
          values = col
        ) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = axis.title.size),
          axis.title.y = element_text(size = axis.title.size),
          axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size * 0.75),
          axis.text.y = element_text(size = base.size),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )
      if (group.var == "ALL"){
        p <- p + theme(legend.position = "none")
      }

      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(
          cols = vars(!!sym(strata.var)),
          scale = "free",
          space = "free"
        )
      }

      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0("beta_change_spaghettiplot_",
                           dist.name,
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
