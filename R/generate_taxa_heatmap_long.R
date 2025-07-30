#' @title Generate Taxonomic Heatmap Long
#'
#' @description This function performs hierarchical clustering on microbiome data based on grouping
#' variables and strata variables in sample metadata and generates stacked heatmaps
#' using the “pheatmap” package. It can also save the resulting heatmap as a PDF file.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
#' @param subject.var A character string specifying the subject variable in the metadata.
#' @param time.var A character string specifying the time variable in the metadata.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var A character string specifying the grouping variable in the metadata. Default is NULL.
#' @param strata.var A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param feature.level A character vector specifying the taxa level(s) to include in the analysis. Default is c('Phylum', 'Family', 'Genus').
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param top.k.plot A numeric value specifying the number of top taxa to be plotted if features.plot is NULL. If NULL (default), all taxa will be plotted.
#' @param top.k.func A function to compute the top k taxa if features.plot is NULL. If NULL (default), the mean function will be used.
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
#' @param cluster.rows A logical variable indicating if rows should be clustered. Default is TRUE.
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is FALSE.
#' @param pdf A logical value. If TRUE (default), saves the plot as a PDF file. If FALSE, the plot will be displayed interactively without creating a PDF.
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional parameters to be passed to the pheatmap() function from the “pheatmap” package.
#'
#' @return An object of class pheatmap, the generated heatmap plot
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[2:18],
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Family","Phylum","Genus", "Class"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[2:18],
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Family","Phylum","Genus", "Class"),
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   cluster.rows = FALSE,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "sample_body_site",
#'   strata.var = "subject_gender",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 15
#' )
#' }
#' @export
#' @importFrom dplyr distinct pull
#' @seealso \code{\link{pheatmap}}
generate_taxa_heatmap_long <- function(data.obj,
                                       subject.var,
                                       time.var,
                                       t0.level,
                                       ts.levels,
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
                                       cluster.cols = NULL,
                                       cluster.rows = NULL,
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       ...) {
  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Validate the input data object
  mStat_validate_data(data.obj)

  # Input validation for key variables
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

  # Process time variable in the data object
  data.obj <-
    mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  # Extract relevant metadata
  meta_tab <-
    data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
      subject.var, group.var, time.var, strata.var
    )))

  # If no group variable is provided, create a dummy "ALL" group
  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- "ALL"
  }

  # If a strata variable is provided, create an interaction term with the group variable
  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(strata.var), !!sym(group.var)))
  }

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
  if (feature.dat.type == "count") {
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
    )
    data.obj <-
      mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
  }

  # Create a list to store plots for each taxonomic level
  plot_list <- lapply(feature.level, function(feature.level) {
    # Aggregate data by taxonomy if necessary
    if (is.null(data.obj$feature.agg.list[[feature.level]]) &
        feature.level != "original") {
      data.obj <-
        mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    # Get the appropriate feature table
    if (feature.level != "original") {
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- data.obj$feature.tab
    }

    # Apply abundance and prevalence filters
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

    # Prepare the data for heatmap visualization
    otu_tab_norm <-
      otu_tax_agg_numeric %>%
      filter(!is.na(!!sym(feature.level))) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # Calculate mean values for each group and time point
    wide_data <- otu_tab_norm %>%
      as.data.frame() %>%
      rownames_to_column(var = feature.level) %>%
      tidyr::gather(key = "sample",
                    value = "value",
                    -one_of(feature.level)) %>%
      dplyr::left_join(meta_tab %>%
                         rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(feature.level),
                      !!sym(group.var),
                      !!sym(time.var)) %>%
      dplyr::summarise(mean_value = mean(value)) %>%
      tidyr::unite("group_time", c(group.var, time.var), sep = "_") %>%
      tidyr::spread(key = "group_time", value = "mean_value") %>% column_to_rownames(feature.level)

    # Prepare column annotations for the heatmap
    annotation_col <- meta_tab %>%
      select(!!sym(time.var), !!sym(group.var)) %>%
      as_tibble() %>%
      dplyr::distinct() %>%
      dplyr::mutate(group_time = paste(!!sym(group.var), !!sym(time.var), sep = "_")) %>%
      column_to_rownames("group_time")

    # Sort the annotation columns
    annotation_col_sorted <-
      annotation_col[order(annotation_col[[group.var]], annotation_col[[time.var]]),]

    # Handle strata variable if present
    if (!is.null(strata.var)) {
      annotation_col_sorted <- annotation_col_sorted %>%
        tidyr::separate(!!sym(group.var),
                        into = c(strata.var, group.var),
                        sep = "\\.") %>%
        dplyr::select(!!sym(time.var),!!sym(group.var),!!sym(strata.var))

      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]],
                                    annotation_col_sorted[[group.var]],
                                    annotation_col_sorted[[time.var]]),]
    }

    # Sort the data according to the annotation
    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted)]

    # Calculate gaps for the heatmap
    if (!is.null(group.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[group.var]]))[-length(table(annotation_col_sorted[[group.var]]))]
    } else {
      gaps <- NULL
    }

    if (!is.null(strata.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[strata.var]]))[-length(table(annotation_col_sorted[[strata.var]]))]
    }

    # Define color palette for the heatmap
    col <- c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020")
    my_col <- colorRampPalette(col)
    n_colors <- 100

    # Filter features if specified
    if (!is.null(features.plot)) {
      wide_data_sorted <-
        wide_data_sorted[rownames(wide_data_sorted) %in% features.plot,]
    }

    # Define colors for group and strata variables
    color_vector <- mStat_get_palette(palette)
    
    # Check if group.var is already a factor to preserve its order
    if (is.factor(annotation_col_sorted[[group.var]])) {
      group_levels <- levels(annotation_col_sorted[[group.var]])
    } else {
      group_levels <-
        annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
    }
    group_colors <-
      setNames(color_vector[1:length(group_levels)], group_levels)

    if (!is.null(strata.var)){
      # Check if strata.var is already a factor to preserve its order
      if (is.factor(annotation_col_sorted[[strata.var]])) {
        strata_levels <- levels(annotation_col_sorted[[strata.var]])
      } else {
        strata_levels <-
          annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      }
      strata_colors <-
        setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)
    }

    # Create annotation color list
    if (!is.null(strata.var)){
      annotation_colors_list <- setNames(list(group_colors, strata_colors),
                                         c(group.var, strata.var))
    } else {
      annotation_colors_list <- setNames(list(group_colors),
                                         c(group.var))
    }

    if (group.var == "ALL"){
      annotation_col_sorted <- annotation_col_sorted %>% select(-all_of("ALL"))
      annotation_colors_list <- NULL
    }

    # Generate heatmap for average values
    heatmap_plot_average <- pheatmap::pheatmap(
      mat = wide_data_sorted[order(rowMeans(wide_data_sorted, na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col(n_colors),
      ...
    )

    gg_heatmap_plot_average <- as.ggplot(heatmap_plot_average)

    # Filter features for individual heatmap if specified
    if (!is.null(features.plot)) {
      otu_tab_norm <-
        otu_tab_norm[rownames(otu_tab_norm) %in% features.plot,]
    }

    # Prepare metadata for individual heatmap
    if (group.var == "ALL"){
      meta_tab <- meta_tab %>% select(-all_of("ALL"))
      annotation_colors_list <- NULL
    } else {
      meta_tab <-
        data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
          subject.var, group.var, time.var, strata.var
        )))
    }

    # Generate heatmap for individual samples
    heatmap_plot_indiv <- pheatmap::pheatmap(
      mat = otu_tab_norm[order(rowMeans(otu_tab_norm, na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = meta_tab,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = TRUE,
      show_colnames = FALSE,
      gaps_col = gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col(n_colors),
      ...
    )

    gg_heatmap_plot_indiv <- as.ggplot(heatmap_plot_indiv)

    # Save the average heatmap as a PDF file if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_long_average",
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
        plot = gg_heatmap_plot_average
      )
    }

    # Save the individual heatmap as a PDF file if requested
    if (pdf) {
      pdf_name <- paste0(
        "taxa_heatmap_long_indiv",
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
        plot = gg_heatmap_plot_indiv
      )
    }

    # Combine individual and average heatmaps into a list
    sub.plot_list <- list(gg_heatmap_plot_indiv, gg_heatmap_plot_average)
    names(sub.plot_list) <- c("indiv","average")

    return(sub.plot_list)
  })

  # Name the plots in the list by feature level
  names(plot_list) <- feature.level

  # Return the list of heatmap plots
  return(plot_list)
}