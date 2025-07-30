#' @title Generate Taxonomic Change Heatmap Long
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
#' @param strata.var (Optional) A character string specifying the stratification variable in the metadata. Default is NULL.
#' @param feature.level A character string defining the taxonomic level to analyze ('Phylum', 'Family', or 'Genus').
#' @param feature.change.func A function or character string specifying how to calculate
#' the change from baseline value. This allows flexible options:
#' - If a function is provided, it will be applied to each row to calculate change.
#'   The function should take 2 arguments: value at timepoint t and value at baseline t0.
#' - If a character string is provided, following options are supported:
#'   - 'relative change': (value_t - value_t0) / (value_t + value_t0)
#'   - 'absolute change': value_t - value_t0
#'   - 'log fold change': log2(value_t + 1e-5) - log2(value_t0 + 1e-5)
#' - Default is 'relative change'.
#'
#' If none of the above options are matched, an error will be thrown indicating
#' the acceptable options or prompting the user to provide a custom function.
#' @details This parameter is used to compute the change columns from baseline
#'   (specified by t0.level) for each taxon. The change values are calculated
#'   for each timepoint and appended as new columns in the data frame before
#'   plotting heatmap. This allows flexibly customizing how change is quantified.
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
#' @param palette Specifies the color palette to be used for annotating groups and strata in the heatmap.
#'                The parameter can be provided in several ways:
#'                - As a character string denoting a predefined palette name.
#'                  Available predefined palettes include 'npg', 'aaas', 'nejm',
#'                  'lancet', 'jama', 'jco', and 'ucscgb', sourced from the `mStat_get_palette` function.
#'                - As a vector of color codes in a format accepted by ggplot2
#'                  (e.g., hexadecimal color codes).
#'                If `palette` is NULL or an unrecognized string, a default color palette will be used.
#'                The function assigns colors from this palette to the unique levels of
#'                `group.var` and, if provided, `strata.var`. When both `group.var` and
#'                `strata.var` are present, `group.var` levels are colored using the
#'                beginning of the palette, while `strata.var` levels are colored using
#'                the reversed palette, ensuring a distinct color representation for each.
#'                If only `group.var` is provided, its levels are assigned colors from the
#'                palette sequentially. If neither `group.var` nor `strata.var` is provided,
#'                no annotation colors are applied.
#' @param cluster.rows A logical variable indicating if rows should be clustered. Default is TRUE.
#' @param cluster.cols A logical variable indicating if columns should be clustered. Default is FALSE.
#' @param pdf If TRUE, save the plot as a PDF file (default: TRUE)
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param ... Additional arguments passed to pheatmap.
#' @return A list of ggplot heatmap objects, one for each taxonomic level.
#'
#' @details This function generates a separate heatmap for each taxonomic level specified,
#'   with rows clustered and layers arranged by groups over timepoints.
#'   It automatically rarefies raw count data using Rarefy-TSS normalization in MicrobiomeStat.
#'   Annotation columns are generated and ordered properly for visually stacking the layers.
#'   Colormaps are also generated for group and strata variables.
#'
#' @seealso \code{\link{pheatmap}} for heatmap, \code{\link{mStat_normalize_data}} for data normalization.
#'
#' @examples
#' \dontrun{
#' library(pheatmap)
#' data(ecam.obj)
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   feature.level = c("Family","Class"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   palette = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   feature.level = c("Family","Class"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "proportion",
#'   features.plot = NULL,
#'   cluster.rows = FALSE,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   palette = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_change_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Phylum"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 50,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_change_heatmap_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Phylum"),
#'   feature.change.func = "log fold change",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   cluster.rows = FALSE,
#'   top.k.plot = 50,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#'
#' @return An object of class pheatmap, the generated heatmap plot
#' @export
#' @importFrom dplyr distinct pull
#' @seealso \code{\link{pheatmap}}
generate_taxa_change_heatmap_long <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              t0.level = NULL,
                                              ts.levels = NULL,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              feature.level,
                                              feature.change.func = "relative change",
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

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Validate the input data object
  mStat_validate_data(data.obj)

  # Check if the input variables are of the correct type
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

  # Process the time variable in the data object
  data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

  # Extract relevant columns from the metadata
  meta_tab <- data.obj$meta.dat %>%
    as.data.frame() %>%
    select(all_of(c(subject.var,group.var,time.var,strata.var)))

  # If group.var is not provided, create a dummy group variable
  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- "ALL"
  }

  # Set default values for clustering if not provided
  if (is.null(cluster.cols)) {
    cluster.cols = FALSE
  } else {
    cluster.cols = TRUE
  }

  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
  }

  # If strata.var is provided, create an interaction term with group.var
  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var), !!sym(strata.var)))
  }

  # Adjust filtering parameters based on input conditions
  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  # Normalize count data if necessary
  if (feature.dat.type == "count"){
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
    )
    data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
  }

  # Generate plots for each taxonomic level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
      data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    # Extract the appropriate feature table
    if (feature.level != "original"){
      otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
    } else {
      otu_tax_agg <- data.obj$feature.tab
    }

    # Filter the feature table based on prevalence and abundance
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

    # Calculate mean values for each combination of feature, group, and time
    df_mean_value <- otu_tax_agg %>%
      tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
      dplyr::left_join(meta_tab %>%
                  rownames_to_column("sample"), "sample") %>%
      dplyr::group_by(!!sym(feature.level),!!sym(group.var),!!sym(time.var)) %>%
      dplyr::summarise(mean_value = mean(value), .groups = "drop")

    # Spread the data to wide format
    df_wide <- df_mean_value %>%
      tidyr::spread(key = time.var, value = mean_value)

    # Determine baseline time point if not provided
    if (is.null(t0.level)) {
      if (is.numeric(meta_tab[, time.var])) {
        t0.level <- sort(unique(meta_tab[, time.var]))[1]
      } else {
        t0.level <- levels(meta_tab[, time.var])[1]
      }
    }

    # Determine follow-up time points if not provided
    if (is.null(ts.levels)) {
      if (is.numeric(meta_tab[, time.var])) {
        ts.levels <- sort(unique(meta_tab[, time.var]))[-1]
      } else {
        ts.levels <- levels(meta_tab[, time.var])[-1]
      }
    }

    # Calculate changes from baseline for each follow-up time point
    for (ts in ts.levels) {
      change_col_name <- paste0("change_", ts)

      # Apply the specified change function
      if (is.function(feature.change.func)) {
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := feature.change.func(.data[[as.character(ts)]], .data[[as.character(t0.level)]]))
      } else if (feature.change.func == "relative change") {
        # Calculate relative change: (value_t - value_t0) / (value_t + value_t0)
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := dplyr::case_when(
            (.data[[as.character(ts)]] + .data[[as.character(t0.level)]]) != 0 ~ (.data[[as.character(ts)]] - .data[[as.character(t0.level)]]) / (.data[[as.character(ts)]] + .data[[as.character(t0.level)]]),
            TRUE ~ 0
          ))
      } else if (feature.change.func == "absolute change") {
        # Calculate absolute change: value_t - value_t0
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := (.data[[as.character(ts)]] - .data[[as.character(t0.level)]]))
      } else if (feature.change.func == "log fold change") {
        # Calculate log fold change: log2(value_t + pseudocount) - log2(value_t0 + pseudocount)
        df_wide <- df_wide %>%
          dplyr::rowwise() %>%
          dplyr::mutate(!!sym(change_col_name) := log2(.data[[as.character(ts)]] + 0.00001) - log2(.data[[as.character(t0.level)]] + 0.00001))
      } else {
        stop(
          "`feature.change.func` must be either 'relative change', 'absolute change', 'log fold change' or a function."
        )
      }
    }

    # Select relevant columns
    df_wide <-
      df_wide[, c(feature.level, group.var, paste0("change_", ts.levels))]

    # Reshape data from wide to long format
    df_long <- df_wide %>%
      tidyr::pivot_longer(
        cols = starts_with("change"),
        names_to = "time",
        values_to = "value"
      )

    # Remove the "change_" prefix from the time column
    df_long$time <- gsub("change_", "", df_long$time)

    # Reshape data back to wide format, combining group and time
    df_wide_new <- df_long %>%
      tidyr::unite("group_time", c(group.var, "time"), sep = "_") %>%
      tidyr::pivot_wider(names_from = "group_time",
                  values_from = "value")

    # Set row names to feature names
    wide_data <- df_wide_new %>% column_to_rownames(feature.level)

    # Save original column names
    original_colnames <- colnames(wide_data)

    # Remove columns where all values are NA
    wide_data <-
      wide_data[, colSums(is.na(wide_data)) != nrow(wide_data)]

    # Save new column names
    new_colnames <- colnames(wide_data)

    # Find out removed columns
    removed_colnames <- setdiff(original_colnames, new_colnames)

    # Notify user about removed columns
    if (length(removed_colnames) > 0) {
      message(
        "The following combinations were all NA, therefore they have been removed: ",
        paste(removed_colnames, collapse = ", ")
      )
    } else {
      message("No columns have been removed.")
    }

    # Prepare annotation for columns
    annotation_col <- meta_tab %>%
      select(!!sym(time.var), !!sym(group.var)) %>%
      filter(!!sym(time.var) != t0.level) %>%
      as_tibble() %>%
      dplyr::distinct() %>%
      dplyr::mutate(group_time = paste(!!sym(group.var), !!sym(time.var), sep = "_")) %>%
      filter(group_time %in% new_colnames) %>%
      column_to_rownames("group_time")

    # Sort annotation by group and time
    annotation_col_sorted <-
      annotation_col[order(annotation_col[[group.var]], annotation_col[[time.var]]),]

    # Handle strata variable if present
    if (!is.null(strata.var)) {
      annotation_col_sorted <- annotation_col_sorted %>%
        tidyr::separate(!!sym(group.var),
                 into = c(group.var, strata.var),
                 sep = "\\.")

      annotation_col_sorted <-
        annotation_col_sorted[order(annotation_col_sorted[[strata.var]], annotation_col_sorted[[group.var]], annotation_col_sorted[[time.var]]), ]

    }

    # Handle the case when group.var is "ALL"
    if (group.var == "ALL") {
      annotation_col_sorted <-
        annotation_col_sorted %>% select(all_of(c(time.var)))
    }

    # Sort wide data according to annotation
    wide_data_sorted <- wide_data[, rownames(annotation_col_sorted)]

    # Filter features if specified
    if (!is.null(features.plot)) {
      wide_data_sorted <-
        wide_data_sorted[rownames(wide_data_sorted) %in% features.plot,]
    }

    # Calculate gaps for heatmap
    if (!is.null(group.var)) {
      gaps <-
        cumsum(table(annotation_col_sorted[[group.var]]))[-length(table(annotation_col_sorted[[group.var]]))]
    } else {
      gaps <- NULL
    }

    if (!is.null(strata.var)){
      gaps <-
        cumsum(table(annotation_col_sorted[[strata.var]]))[-length(table(annotation_col_sorted[[strata.var]]))]
    }

    # Set up color scheme for heatmap
    n_colors <- 100
    col <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
    
    # Find the maximum absolute value in the data
    max_abs_val <-
      max(abs(range(na.omit(
        c(as.matrix(wide_data_sorted))
      ))))

    # Calculate the position of the zero value in the new color vector
    zero_pos <- round(max_abs_val / (2 * max_abs_val) * n_colors)

    # Create color vectors
    my_col <-
      c(
        colorRampPalette(col[1:3])(zero_pos),
        colorRampPalette(col[3:5])(n_colors - zero_pos + 1)
      )

    # Create break points for color scale
    break_points <-
      seq(-max_abs_val, max_abs_val, length.out = length(my_col) + 1)

    # Get color palette
    color_vector <- mStat_get_palette(palette)

    # Set up annotation colors
    if (!is.null(strata.var) & !is.null(group.var)){
      # Get unique levels for group and strata variables
      # Check if variables are already factors to preserve their order
      if (is.factor(annotation_col_sorted[[group.var]])) {
        group_levels <- levels(annotation_col_sorted[[group.var]])
      } else {
        group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      }
      
      if (is.factor(annotation_col_sorted[[strata.var]])) {
        strata_levels <- levels(annotation_col_sorted[[strata.var]])
      } else {
        strata_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      }

      # Assign colors to group and strata levels
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

      # Create a list of annotation colors
      annotation_colors_list <- setNames(
        list(group_colors, strata_colors),
        c(group.var, strata.var)
      )
    } else if (!is.null(group.var) & group.var != "ALL"){
      # Get unique levels for group variable
      # Check if variable is already a factor to preserve its order
      if (is.factor(annotation_col_sorted[[group.var]])) {
        group_levels <- levels(annotation_col_sorted[[group.var]])
      } else {
        group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      }
      
      # Assign colors to group levels
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      
      # Create annotation color list
      annotation_colors_list <- setNames(
        list(group_colors),
        c(group.var)
      )
    } else {
      annotation_colors_list <- NULL
    }

    # Plot stacked heatmap
    heatmap_plot <- pheatmap::pheatmap(
      mat = wide_data_sorted[order(rowMeans(abs(wide_data_sorted), na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      show_colnames = FALSE,
      gaps_col = gaps,
      fontsize = base.size,
      silent = TRUE,
      color = my_col,
      breaks = break_points,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Save the stacked heatmap as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_long",
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
        "taxa_",
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
    return(gg_heatmap_plot)
  })

  names(plot_list) <- feature.level

  # Return the heatmap plot for display
  return(plot_list)
}