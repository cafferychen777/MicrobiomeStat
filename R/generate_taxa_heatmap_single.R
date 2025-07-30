#' @title Generate Taxonomic Heatmap Single
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var The name of the subject variable in the samples
#' @param time.var The name of the time variable in the samples
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var The name of the grouping variable in the samples
#' @param strata.var The name of the strata variable in the samples
#' @param other.vars A character vector specifying additional variables from the metadata to include
#'   in the heatmap annotation. These variables will be added to the annotation columns alongside
#'   `group.var` and `strata.var`. This allows for the visualization of additional metadata
#'   information in the heatmap. Default is NULL, which means no additional variables are included.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param top.k.plot Integer specifying number of top k features to plot, when `features.plot` is NULL.
#' Default is NULL, in which case all features passing filters will be plotted.
#' @param top.k.func Function to use for selecting top k features, when `features.plot` is NULL.
#' Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
#' case features will be selected by abundance.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size Base font size for the generated plots.
#' @param palette The color palette to be used for annotating the plots.
#'                This parameter can be specified in several ways:
#'                - As a character string representing a predefined palette name.
#'                  Available predefined palettes include 'npg', 'aaas', 'nejm',
#'                  'lancet', 'jama', 'jco', and 'ucscgb'.
#'                - As a vector of color codes in a format accepted by ggplot2
#'                  (e.g., hexadecimal color codes).
#'                The function uses `mStat_get_palette` to retrieve or generate
#'                the color palette. If `palette` is NULL or an unrecognized string,
#'                a default color palette will be used. The colors are applied to
#'                the specified grouping variables (`group.var`, `strata.var`) in the
#'                heatmap, ensuring each level of these variables is associated with a
#'                unique color. If both `group.var` and `strata.var` are specified,
#'                the function assigns colors to `group.var` from the start of the
#'                palette and to `strata.var` from the end, ensuring distinct color
#'                representations for each annotation layer.
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is NULL.
#' @param cluster.rows A logical variable indicating if rows should be clustered. Default is NULL.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann The file name annotation (default: NULL)
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional arguments passed to the pheatmap() function from the “pheatmap” package.
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
#'   subject.var = "subject",
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
#'   subject.var = "subject",
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
#'   subject.var = "subject",
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
#'   subject.var = "subject",
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
#'   subject.var = "subject.id",
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
#'   subject.var = "subject.id",
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
                                         subject.var,
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

  # If subject variable is not provided, create a default one
  if (is.null(subject.var)){
    data.obj$meta.dat$subject.id <- rownames(data.obj$meta.dat)
    subject.var <- "subject.id"
  }

  # Subset the data if time variable and level are provided
  if (!is.null(time.var)) {
    if (!is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
  }

  # Select relevant variables from metadata
  meta_tab <- data.obj$meta.dat %>% select(all_of(
    c(subject.var, time.var, group.var, strata.var, other.vars)))

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

    # Prepare the normalized OTU table
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
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "t_level_",
        t.level,
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
    # Return the heatmap plot
    return(gg_heatmap_plot)
  })

  # Name the plots in the list
  names(plot_list) <- feature.level
  return(plot_list)
}