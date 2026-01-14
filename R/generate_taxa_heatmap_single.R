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
                                         feature.level = NULL,
                                         feature.dat.type = c("count", "proportion", "other"),
                                         features.plot = NULL,
                                         top.k.plot = NULL,
                                         top.k.func = NULL,
                                         prev.filter = 0.01,
                                         abund.filter = 0.0001,
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
  mStat_validate_data(data.obj)

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Subset the data if time variable and level are provided
  if (!is.null(time.var)) {
    if (!is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
  }

  # Select relevant variables from metadata
  meta_tab <- data.obj$meta.dat %>% select(all_of(
    c(time.var, group.var, strata.var, other.vars)))

  # Define color palette for the heatmap
  col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")

  # Create color mapping function
  my_col <- colorRampPalette(col)

  # Set the number of colors for the heatmap
  n_colors <- 100

  # Set clustering options for columns and rows
  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  # Disable filtering if specific conditions are met
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  # For "other" data type, check if data contains negative values
  # If so, adjust abundance filter to handle negative values appropriately
  if (feature.dat.type == "other") {
    # Check if any feature table contains negative values
    has_negative <- FALSE
    if (!is.null(data.obj$feature.tab)) {
      has_negative <- any(data.obj$feature.tab < 0, na.rm = TRUE)
    }
    if (!has_negative && !is.null(data.obj$feature.agg.list)) {
      for (agg_table in data.obj$feature.agg.list) {
        if (any(agg_table < 0, na.rm = TRUE)) {
          has_negative <- TRUE
          break
        }
      }
    }

    if (has_negative) {
      message("Note: Negative values detected in 'other' data type. Abundance filtering is disabled to preserve all features.")
      abund.filter <- -Inf  # Set to negative infinity to include all features regardless of abundance
    }
  }

  # Normalize count data if necessary
  if (feature.dat.type == "count"){
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
    )
    data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
  }

  # Generate plots for each feature level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
      data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    # Select appropriate feature table
    if (feature.level != "original"){
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- data.obj$feature.tab
    }

    # Filter and prepare the feature table
    otu_tax_agg <-  otu_tax_agg %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter) %>%
      rownames_to_column(feature.level)

    # Select top k features if specified
    if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
    }

    # Convert counts to numeric
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    # Prepare the normalized feature table
    otu_tab_norm <-
      otu_tax_agg_numeric %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # Sort samples by group and strata variables if provided
    if (!is.null(group.var) & !is.null(strata.var)) {
      meta_tab_sorted <-
        meta_tab[order(meta_tab[[strata.var]], meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else if (!is.null(group.var) & is.null(strata.var)) {
      meta_tab_sorted <- meta_tab[order(meta_tab[[group.var]]), ]
      otu_tab_norm_sorted <-
        otu_tab_norm[, rownames(meta_tab_sorted)]
    } else {
      otu_tab_norm_sorted <- otu_tab_norm
    }

    # Calculate gaps for the heatmap
    if (!is.null(strata.var)){
      if (!is.numeric(meta_tab[[strata.var]])){
        gaps <-
          cumsum(table(meta_tab_sorted[[strata.var]]))[-length(table(meta_tab_sorted[[strata.var]]))]
      }
    } else if (!is.null(group.var)){
      if (!is.numeric(meta_tab[[group.var]])) {
        gaps <-
          cumsum(table(meta_tab_sorted[[group.var]]))[-length(table(meta_tab_sorted[[group.var]]))]
      }
    } else {
      gaps <- NULL
    }

    # Prepare annotation for columns
    annotation_col <-
      meta_tab %>% select(all_of(c(other.vars, group.var, strata.var)))

    if (ncol(annotation_col) == 0){
      annotation_col <- NULL
    }

    # Create title for the heatmap
    heatmap_title <- NA
    if (!is.null(time.var) & !is.null(t.level)) {
      heatmap_title <- paste0("Time = ", t.level)
    }

    # Filter features to plot if specified
    if (!is.null(features.plot)) {
      otu_tab_norm_sorted <-
        otu_tab_norm_sorted[rownames(otu_tab_norm_sorted) %in% features.plot,]
    }

    # Get color palette
    color_vector <- mStat_get_palette(palette)

    # Prepare annotation colors
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

    # Generate the heatmap
    heatmap_plot <- pheatmap::pheatmap(
      mat = otu_tab_norm_sorted[order(rowMeans(otu_tab_norm_sorted, na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_col,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      gaps_col = gaps,
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