#' @title Generate Taxa Heatmap for Single Time Point
#'
#' @description Generates hierarchical clustered heatmaps for microbiome data at a single time point
#' using the pheatmap package, with options for group and strata annotations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param t.level Character string specifying the time level/value to subset data to.
#'   Default NULL does not subset data.
#' @param other.vars A character vector specifying additional metadata variables to include
#'   in the heatmap annotation. Default is NULL.
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param cluster.cols Logical indicating if columns should be clustered. Default is NULL.
#' @param cluster.rows Logical indicating if rows should be clustered. Default is NULL.
#' @param ... Additional arguments passed to the pheatmap() function.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(pheatmap)
#'
#' data(peerj32.obj)
#' generate_taxa_heatmap_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   other.vars = NULL,
#'   feature.level = c("Phylum", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_heatmap_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   other.vars = NULL,
#'   feature.level = c("Phylum", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   cluster.rows = FALSE,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_heatmap_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = NULL,
#'   other.vars = NULL,
#'   feature.level = c("Phylum", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_heatmap_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = NULL,
#'   strata.var = NULL,
#'   other.vars = NULL,
#'   feature.level = c("Phylum", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(ecam.obj)
#' generate_taxa_heatmap_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "0",
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   other.vars = "delivery",
#'   feature.level = c("Order", "Family", "Genus"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_heatmap_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "0",
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   other.vars = "delivery",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "proportion",
#'   features.plot = unique(ecam.obj$feature.ann[,"Genus"])[-c(1,9)],
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.001,
#'   abund.filter = 0.01,
#'   base.size = 10,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
#'
#' @seealso \code{\link{pheatmap}}
generate_taxa_heatmap_single <- function(data.obj,
                                         time.var = NULL,
                                         t.level = NULL,
                                         group.var = NULL,
                                         strata.var = NULL,
                                         other.vars = NULL,
                                         feature.level,
                                         feature.dat.type = c("count", "proportion", "other"),
                                         features.plot = NULL,
                                         top.k.plot = NULL,
                                         top.k.func = NULL,
                                         prev.filter = 0.01,
                                         abund.filter = 0.01,
                                         base.size = 10,
                                         palette = NULL,
                                         cluster.cols = NULL,
                                         cluster.rows = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         ...) {
  # Validate the input data object
  data.obj <- mStat_validate_data(data.obj)

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Subset the data if time variable and level are provided
  if (!is.null(time.var) && !is.null(t.level)) {
    data.obj <- mStat_subset_by_meta_values(data.obj, time.var, t.level)
  }

  # Select relevant variables from metadata
  meta_tab <- data.obj$meta.dat %>% select(all_of(
    c(time.var, group.var, strata.var, other.vars)))

  col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
  my_col <- colorRampPalette(col)
  n_colors <- 100

  # Set clustering options for columns and rows
  cluster.cols <- mStat_resolve_optional_flag(cluster.cols, FALSE, "cluster.cols")
  cluster.rows <- mStat_resolve_optional_flag(cluster.rows, TRUE, "cluster.rows")

  # Disable filtering if specific conditions are met
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  abund.filter <- mStat_adjust_other_abundance_filter(
    data.obj = data.obj,
    feature.dat.type = feature.dat.type,
    abund.filter = abund.filter
  )

  # Normalize count data if necessary.
  analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

  # Generate plots for each feature level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    otu_tax_agg <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter)

    current_features_plot <- mStat_resolve_selected_features(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      features.plot = features.plot,
      top.k.plot = top.k.plot,
      top.k.func = top.k.func
    )

    # Convert counts to numeric
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    # Prepare the normalized feature table
    otu_tab_norm <-
      otu_tax_agg_numeric %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      mStat_as_taxa_feature_matrix(
        feature.level = feature.level,
        feature_in_column = TRUE
      )

    if (!is.null(group.var) && !is.null(strata.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[strata.var]], meta_tab[[group.var]]), , drop = FALSE]
      otu_tab_norm_sorted <- otu_tab_norm[, rownames(meta_tab_sorted), drop = FALSE]
    } else if (!is.null(group.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[group.var]]), , drop = FALSE]
      otu_tab_norm_sorted <- otu_tab_norm[, rownames(meta_tab_sorted), drop = FALSE]
    } else if (!is.null(strata.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[strata.var]]), , drop = FALSE]
      otu_tab_norm_sorted <- otu_tab_norm[, rownames(meta_tab_sorted), drop = FALSE]
    } else {
      meta_tab_sorted <- meta_tab
      otu_tab_norm_sorted <- otu_tab_norm
    }

    if (!is.null(strata.var) && !is.numeric(meta_tab_sorted[[strata.var]])) {
      gaps <- cumsum(table(meta_tab_sorted[[strata.var]]))[-length(table(meta_tab_sorted[[strata.var]]))]
    } else if (!is.null(group.var) && !is.numeric(meta_tab_sorted[[group.var]])) {
      gaps <- cumsum(table(meta_tab_sorted[[group.var]]))[-length(table(meta_tab_sorted[[group.var]]))]
    } else {
      gaps <- NULL
    }

    annotation_col <- meta_tab %>% select(all_of(c(other.vars, group.var, strata.var)))
    if (ncol(annotation_col) == 0) {
      annotation_col <- NULL
    }

    heatmap_title <- if (!is.null(time.var) && !is.null(t.level)) paste0("Time = ", t.level) else NA

    # Filter features to plot if specified
    if (!is.null(current_features_plot)) {
      otu_tab_norm_sorted <-
        otu_tab_norm_sorted[rownames(otu_tab_norm_sorted) %in% current_features_plot, , drop = FALSE]
    }

    color_vector <- mStat_get_palette(palette)

    if (!is.null(strata.var) && !is.null(group.var)) {
      group_levels <- if (is.factor(annotation_col[[group.var]])) levels(annotation_col[[group.var]]) else annotation_col %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      strata_levels <- if (is.factor(annotation_col[[strata.var]])) levels(annotation_col[[strata.var]]) else annotation_col %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      group_colors <- setNames(color_vector[seq_along(group_levels)], group_levels)
      strata_colors <- setNames(rev(color_vector)[seq_along(strata_levels)], strata_levels)
      annotation_colors_list <- setNames(list(group_colors, strata_colors), c(group.var, strata.var))
    } else if (!is.null(group.var)) {
      group_levels <- if (is.factor(annotation_col[[group.var]])) levels(annotation_col[[group.var]]) else annotation_col %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      group_colors <- setNames(color_vector[seq_along(group_levels)], group_levels)
      annotation_colors_list <- setNames(list(group_colors), c(group.var))
    } else if (!is.null(strata.var)) {
      strata_levels <- if (is.factor(annotation_col[[strata.var]])) levels(annotation_col[[strata.var]]) else annotation_col %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      strata_colors <- setNames(rev(color_vector)[seq_along(strata_levels)], strata_levels)
      annotation_colors_list <- setNames(list(strata_colors), c(strata.var))
    } else {
      annotation_colors_list <- NULL
    }

    heatmap_gaps <- if (!is.null(gaps) && length(gaps) > 0 && ncol(otu_tab_norm_sorted) > 1) {
      gaps
    } else {
      NULL
    }

    # Generate the heatmap
    heatmap_plot <- pheatmap::pheatmap(
      mat = otu_tab_norm_sorted[order(rowMeans(otu_tab_norm_sorted, na.rm = TRUE), decreasing = TRUE), , drop = FALSE],
      annotation_col = annotation_col,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      gaps_col = heatmap_gaps,
      color = my_col(n_colors),
      fontsize = base.size,
      main = heatmap_title,
      silent = TRUE,
      ...
    )

    # Convert pheatmap to ggplot object
    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Save the heatmap as a PDF if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_single",
        "_",
        "feature_level_",
        feature.level,
        "_",
        "prev_filter_",
        prev.filter,
        "_",
        "abund_filter_",
        abund.filter
      )
      if (!is.null(time.var)) {
        pdf_name <- paste0(pdf_name, "_", "time_", time.var)
      }
      if (!is.null(t.level)) {
        pdf_name <- paste0(pdf_name, "_", "t_level_", t.level)
      }
      pdf_name <- mStat_append_pdf_group_suffixes(
        pdf_name = pdf_name,
        group.var = group.var,
        strata.var = strata.var
      )
      if (!is.null(file.ann)) {
        pdf_name <- paste0(pdf_name, "_", file.ann)
      }
      pdf_name <- paste0(pdf_name, ".pdf")
      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_heatmap_plot
      )
    }
    # Return the heatmap plot
    return(gg_heatmap_plot)
  })

  # Name the plots in the list
  names(plot_list) <- feature.level
  return(plot_list)
}
