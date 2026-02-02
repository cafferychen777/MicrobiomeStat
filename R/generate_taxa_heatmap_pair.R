#' @title Generate Taxa Heatmap for Paired Samples
#'
#' @description Generates hierarchical clustered heatmaps for paired microbiome data
#' using the pheatmap package, with options for group and strata annotations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param cluster.rows Logical indicating if rows should be clustered. Default is TRUE.
#' @param cluster.cols Logical indicating if columns should be clustered. Default is FALSE.
#' @param ... Additional parameters to be passed to the pheatmap() function.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(pheatmap)
#' data(peerj32.obj)
#' generate_taxa_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Phylum","Family","Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   cluster.rows = NULL,
#'   cluster.cols = NULL,
#'   base.size = 12,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' subset_T2D.obj2 <- mStat_subset_data(subset_T2D.obj,
#' condition = "visit_number %in% c('   1', '   2')")
#' generate_taxa_heatmap_pair(
#'   data.obj = subset_T2D.obj2,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   feature.level = c("Phylum","Family","Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   cluster.rows = NULL,
#'   cluster.cols = NULL,
#'   base.size = 12,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#' generate_taxa_heatmap_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   feature.level = c("Phylum","Family","Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   cluster.rows = NULL,
#'   cluster.cols = NULL,
#'   base.size = 12,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_heatmap_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   feature.level = c("Phylum","Family","Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   cluster.rows = FALSE,
#'   cluster.cols = NULL,
#'   base.size = 12,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_heatmap_pair <- function(data.obj,
                                       subject.var,
                                       time.var,
                                       group.var = NULL,
                                       strata.var = NULL,
                                       feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       features.plot = NULL,
                                       top.k.plot = NULL,
                                       top.k.func = NULL,
                                       prev.filter = 0.01,
                                       abund.filter = 0.01,
                                       base.size = 10,
                                       palette = NULL,
                                       cluster.rows = NULL,
                                       cluster.cols = NULL,
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       ...) {
  # Validate the input data object
  mStat_validate_data(data.obj)

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Extract relevant variables from the metadata
  meta_tab <-  data.obj$meta.dat %>% select(all_of(c(
    time.var, group.var, strata.var, subject.var
  )))

  # Set default clustering options if not specified
  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  # Determine whether to apply abundance and prevalence filters
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  # Normalize count data if necessary
  if (feature.dat.type == "count"){
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
    )
    data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
  }

  # Create a list to store plots for each taxonomic level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

    # Select top k features if specified
    if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
    }

    # Normalize the feature table
    otu_tab_norm <-
      otu_tax_agg %>% column_to_rownames(var = feature.level) %>% as.matrix()

    # Sort samples by group.var and strata.var
    if (!is.null(group.var) & !is.null(strata.var)) {
      meta_tab_sorted <-
        meta_tab[order(meta_tab[[strata.var]], meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else if (!is.null(group.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else {
      otu_tab_norm_sorted <- otu_tab_norm
    }

    # Calculate gaps for the heatmap
    if (!is.null(group.var)){
      gaps <-
        cumsum(table(meta_tab_sorted[[subject.var]]))[-length(table(meta_tab_sorted[[subject.var]]))]
    } else {
      gaps <- NULL
    }

    # Set up annotation_col based on group.var and strata.var values
    if (!is.null(group.var) & !is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(time.var, group.var, strata.var)))
    } else if (!is.null(group.var) & is.null(strata.var)) {
      annotation_col <-
        meta_tab %>% select(all_of(c(time.var, group.var)))
    } else {
      annotation_col <- meta_tab %>% select(all_of(c(time.var)))
    }

    # Filter features if specified
    if (!is.null(features.plot)) {
      otu_tab_norm_sorted <-
        otu_tab_norm_sorted[rownames(otu_tab_norm_sorted) %in% features.plot,]
    }

    # Get color palette
    color_vector <- mStat_get_palette(palette)

    # Set up annotation colors
    if (!is.null(strata.var) & !is.null(group.var)){
      # Check if variables are already factors to preserve their order
      if (is.factor(annotation_col[[group.var]])) {
        group_levels <- levels(annotation_col[[group.var]])
      } else {
        group_levels <- annotation_col %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      }
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)

      if (is.factor(annotation_col[[strata.var]])) {
        strata_levels <- levels(annotation_col[[strata.var]])
      } else {
        strata_levels <- annotation_col %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      }
      strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

      annotation_colors_list <- setNames(
        list(group_colors, strata_colors),
        c(group.var, strata.var)
      )
    } else if (!is.null(group.var)){
      # Check if variable is already a factor to preserve its order
      if (is.factor(annotation_col[[group.var]])) {
        group_levels <- levels(annotation_col[[group.var]])
      } else {
        group_levels <- annotation_col %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      }
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      annotation_colors_list <- setNames(
        list(group_colors),
        c(group.var)
      )
    } else {
      annotation_colors_list <- NULL
    }

    # Set up color scheme for heatmap
    col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
    my_col <- colorRampPalette(col)
    n_colors <- 100

    # Set up quantiles and labels for legend
    if (feature.dat.type != "other"){
      quantiles <- c(0,0.1,0.5,1)
      labels <- as.character(round(quantiles^2,4))
    } else {
      quantiles <- NA
      labels <- NA
    }

    # Apply square root transformation if necessary
    if (feature.dat.type != "other") {
      processed_otu_tab <- sqrt(otu_tab_norm_sorted)
    } else {
      processed_otu_tab <- otu_tab_norm_sorted
    }

    # Generate individual sample heatmap
    heatmap_plot <- pheatmap::pheatmap(
      mat = processed_otu_tab[order(rowMeans(processed_otu_tab, na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_col,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      gaps_col = gaps,
      legend_breaks = quantiles,
      legend_labels = labels,
      color = my_col(n_colors),
      fontsize = base.size,
      silent = TRUE,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Generate average heatmap
    if (!is.null(strata.var) & !is.null(group.var) & !is.null(time.var)){
      wide_data <- otu_tab_norm_sorted %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(cols = -feature.level, names_to = "sample", values_to = "value") %>%
        dplyr::left_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level), !!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("column_name", !!sym(time.var), !!sym(group.var), !!sym(strata.var), sep = "_") %>%
        tidyr::pivot_wider(names_from = column_name, values_from = mean_value) %>%
        dplyr::ungroup()

      annotation_col <- wide_data %>%
        dplyr::select(-all_of(c(feature.level))) %>%
        names() %>%
        tibble(column_name = .) %>%
        dplyr::mutate(
          !!time.var := stringr::str_extract(column_name, "^[^_]+"),
          !!group.var := stringr::str_extract(column_name, "(?<=_).*(?=_)"),
          !!strata.var := stringr::str_extract(column_name, "[^_]+$")
        ) %>%
        column_to_rownames("column_name")

    } else if (!is.null(group.var) & !is.null(time.var)){
      wide_data <- otu_tab_norm_sorted %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(cols = -feature.level, names_to = "sample", values_to = "value") %>%
        dplyr::left_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level), !!sym(group.var), !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("column_name", !!sym(time.var), !!sym(group.var), sep = "_") %>%
        tidyr::pivot_wider(names_from = column_name, values_from = mean_value) %>%
        dplyr::ungroup()

      annotation_col <- wide_data %>%
        dplyr::select(-all_of(c(feature.level))) %>%
        names() %>%
        tibble(column_name = .) %>%
        dplyr::mutate(
          !!time.var := stringr::str_extract(column_name, "^[^_]+"),
          !!group.var := stringr::str_extract(column_name, "(?<=_).*")
        ) %>%
        tidyr::unite("column_name", !!sym(time.var), !!sym(group.var), sep = "_", remove = FALSE) %>%
        column_to_rownames("column_name")

    } else if (!is.null(time.var)){
      wide_data <- otu_tab_norm_sorted %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(cols = -feature.level, names_to = "sample", values_to = "value") %>%
        dplyr::left_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::group_by(!!sym(feature.level), !!sym(time.var)) %>%
        dplyr::summarise(mean_value = mean(value)) %>%
        tidyr::unite("column_name", !!sym(time.var), sep = "_") %>%
        tidyr::pivot_wider(names_from = column_name, values_from = mean_value) %>%
        dplyr::ungroup()

      annotation_col <- wide_data %>%
        dplyr::select(-all_of(c(feature.level))) %>%
        names() %>%
        stringr::str_split("_", simplify = TRUE) %>%
        as.data.frame() %>%
        setNames(c(time.var)) %>%
        distinct() %>%
        tidyr::unite("column_name", !!sym(time.var), sep = "_", remove = FALSE) %>%
        column_to_rownames("column_name")
    }

    # Apply square root transformation if necessary
    if (feature.dat.type != "other") {
      processed_wide_data <- wide_data %>% column_to_rownames(feature.level) %>% sqrt()
    } else {
      processed_wide_data <- wide_data %>% column_to_rownames(feature.level)
    }

    # Generate average heatmap
    average_heatmap_plot <- pheatmap::pheatmap(
      mat = processed_wide_data[order(rowMeans(processed_wide_data, na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_col,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = NULL,
      legend_breaks = quantiles,
      legend_labels = labels,
      color = my_col(n_colors),
      fontsize = base.size,
      silent = TRUE,
      ...
    )

    gg_average_heatmap_plot <- as.ggplot(average_heatmap_plot)

    # Save the average heatmap as a PDF file if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_pair_average",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "group_",
        group.var,
        "_",
        "strata_",
        strata.var,
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
      if (!is.null(file.ann)) {
        pdf_name <- paste0(pdf_name, "_", file.ann)
      }
      pdf_name <- paste0(pdf_name, ".pdf")
      ggsave(
        filename = pdf_name,
        width = pdf.wid,
        height = pdf.hei,
        plot = gg_average_heatmap_plot
      )
    }

    # Save the individual sample heatmap as a PDF file if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_pair_indiv",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "group_",
        group.var,
        "_",
        "strata_",
        strata.var,
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

    sub.plot_list <- list(gg_heatmap_plot, gg_average_heatmap_plot)

    names(sub.plot_list) <- c("indiv", "average")
    # Return the heatmap plot for display
    return(sub.plot_list)
  })

  names(plot_list) <- feature.level
  return(plot_list)
}
