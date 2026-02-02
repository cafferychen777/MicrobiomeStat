#' @title Generate Taxa Change Heatmap for Paired Data
#'
#' @description Creates heatmaps showing taxa abundance changes between paired time points.
#' Generates both individual and group-averaged views using pheatmap.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param features.plot Character vector of specific feature IDs to plot.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#' @param cluster.rows Logical, whether to cluster rows. Default TRUE.
#' @param cluster.cols Logical, whether to cluster columns. Default NULL.
#'
#' @return A list containing individual and average heatmap objects.
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(pheatmap)
#' data(peerj32.obj)
#' generate_taxa_change_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = NULL,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_change_heatmap_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = FALSE,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(subset_pairs.obj)
#' generate_taxa_change_heatmap_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = NULL,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_taxa_change_heatmap_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "relative change",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "sd",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   base.size = 10,
#'   palette = NULL,
#'   cluster.rows = FALSE,
#'   cluster.cols = FALSE,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_change_heatmap_pair <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              change.base = NULL,
                                              feature.change.func = "relative change",
                                              feature.level,
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

  # Extract relevant metadata
  meta_tab <-  data.obj$meta.dat %>%
    select(all_of(c(
    time.var, group.var, strata.var, subject.var
  ))) %>%
    rownames_to_column("sample")

  # Determine whether to cluster columns based on group variable
  if (is.null(cluster.cols)) {
    if (!is.null(group.var)) {
      if (is.numeric(meta_tab %>% select(all_of(group.var)) %>% dplyr::pull()) &&
          !is.integer(meta_tab %>% select(all_of(group.var)) %>% dplyr::pull())) {
        cluster.cols = FALSE
      } else {
        cluster.cols = TRUE
      }
    } else {
      cluster.cols = TRUE
    }
  }

  # Set default for clustering rows if not specified
  if (is.null(cluster.rows)) {
    cluster.rows = TRUE
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

  # Generate plots for each feature level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

    # Select top k features if specified
    if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
    }

    # Reshape data to long format
    otu_tax_long <- otu_tax_agg %>%
      tidyr::gather(key = "sample", value = "value",-feature.level)

    # Merge with metadata
    merged_data <- otu_tax_long %>%
      dplyr::inner_join(meta_tab, by = "sample")

    # Group data by time variable
    grouped_data <- merged_data %>%
      dplyr::group_by(!!sym(time.var))

    # Identify the time point to compare against the baseline
    change.after <-
      unique(grouped_data %>% select(all_of(c(time.var))))[unique(grouped_data %>% select(all_of(c(time.var)))) != change.base]

    # Split data into baseline and comparison time points
    split_data <-
      split(merged_data, f = grouped_data %>% select(all_of(c(time.var))))

    data_time_1 <- split_data[[change.base]]
    data_time_2 <- split_data[[change.after]]

    # Combine data from both time points
    combined_data <- data_time_1 %>%
      dplyr::inner_join(
        data_time_2,
        by = c(feature.level, subject.var),
        suffix = c("_time_1", "_time_2")
      )

    # Calculate the change in feature values between time points
    combined_data <- combined_data %>%
      dplyr::mutate(value_diff = compute_taxa_change(
        value_after  = value_time_2,
        value_before = value_time_1,
        method       = feature.change.func,
        feature_id   = .data[[feature.level]]
      ))

    # Create a matrix of value differences
    value_diff_matrix <- combined_data %>%
      select(feature.level, !!sym(subject.var), value_diff) %>%
      tidyr::spread(key = !!sym(subject.var), value = value_diff) %>%
      column_to_rownames(var = feature.level) %>%
      as.matrix()

    # Prepare metadata for annotation
    unique_meta_tab <- meta_tab %>%
      filter(!!sym(subject.var) %in% colnames(value_diff_matrix)) %>%
      select(all_of(c(subject.var, group.var, strata.var))) %>%
      dplyr::distinct(!!sym(subject.var), .keep_all = TRUE) %>% as_tibble()

    # Match the order of metadata to the matrix columns
    order_index <-
      match(colnames(value_diff_matrix),
            as.matrix(unique_meta_tab %>% select(!!sym(subject.var))))

    suppressWarnings({
      sorted_meta_tab <- unique_meta_tab[order_index,]
      rownames(sorted_meta_tab) <-
        as.matrix(sorted_meta_tab[, subject.var])
    })

    # Sort metadata by group variable if present
    if (!is.null(group.var)) {
      sorted_meta_tab <-
        sorted_meta_tab[order(sorted_meta_tab %>% select(!!sym(group.var)) %>% as.matrix()),]
    }

    # Extract subject column from sorted metadata
    subjects <- sorted_meta_tab %>% select(all_of(subject.var))

    # Rearrange the columns of the value difference matrix
    value_diff_matrix <- value_diff_matrix[, as.matrix(subjects)]

    # Prepare annotation columns for the heatmap
    if (!is.null(strata.var) & !is.null(group.var)){
      annotation_cols <-
        sorted_meta_tab %>%
        select(all_of(c(group.var, strata.var, subject.var))) %>%
        column_to_rownames(var = subject.var)
      annotation_col_sorted <-
        annotation_cols[order(annotation_cols[[strata.var]], annotation_cols[[group.var]]),]
    } else if (!is.null(group.var)){
      annotation_cols <-
        sorted_meta_tab %>%
        select(all_of(c(subject.var, group.var))) %>%
        column_to_rownames(var = subject.var)
      annotation_col_sorted <- annotation_cols
    } else {
      annotation_col_sorted <- NULL
    }

   # Rearrange value difference matrix columns if annotation is present
   if (!is.null(group.var) | !is.null(strata.var)){
     value_diff_matrix <-
       value_diff_matrix[, rownames(annotation_col_sorted)]
   }

    # Determine gaps for the heatmap based on strata or group variables
    if (!is.null(strata.var)) {
      gaps <-
        cumsum(table(sorted_meta_tab[[strata.var]]))[-length(sorted_meta_tab[[strata.var]])]
    } else {
      if (!is.null(group.var)) {
        gaps <-
          cumsum(table(sorted_meta_tab[[group.var]]))[-length(sorted_meta_tab[[group.var]])]
      } else {
        gaps <- NULL
      }
    }

    # Filter features to plot if specified
    if (!is.null(features.plot)) {
      value_diff_matrix <-
        value_diff_matrix[rownames(value_diff_matrix) %in% features.plot,]
    }

    # Set up color scheme for the heatmap
    n_colors <- 100
    col <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
    max_abs_val <- max(abs(range(na.omit(
      c(value_diff_matrix)
    ))))
    zero_pos <- round(max_abs_val / (2 * max_abs_val) * n_colors)
    my_col <-
      c(
        colorRampPalette(col[1:3])(zero_pos),
        colorRampPalette(col[3:5])(n_colors - zero_pos + 1)
      )
    break_points <-
      seq(-max_abs_val, max_abs_val, length.out = length(my_col) + 1)

    # Set up color palette for annotations
    color_vector <- mStat_get_palette(palette)

    if (!is.null(strata.var) & !is.null(group.var)){
      group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)
      strata_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(strata.var))) %>% distinct() %>% pull()
      strata_colors <- setNames(rev(color_vector)[1:length(strata_levels)], strata_levels)

      annotation_colors_list <- setNames(
        list(group_colors, strata_colors),
        c(group.var, strata.var)
      )
    } else if(!is.null(group.var)){
      group_levels <- annotation_col_sorted %>% dplyr::select(all_of(c(group.var))) %>% distinct() %>% pull()
      group_colors <- setNames(color_vector[1:length(group_levels)], group_levels)

      annotation_colors_list <- setNames(
        list(group_colors),
        c(group.var)
      )
    } else {
      annotation_colors_list <- NULL
    }

    # Generate the heatmap for individual data
    heatmap_plot <- pheatmap::pheatmap(
      mat = value_diff_matrix[order(rowMeans(abs(value_diff_matrix), na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_col_sorted,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      annotation_legend = TRUE,
      show_colnames = TRUE,
      show_rownames = TRUE,
      border_color = NA,
      silent = TRUE,
      gaps_col = gaps,
      fontsize = base.size,
      color = my_col,
      breaks = break_points,
      ...
    )

    gg_heatmap_plot <- as.ggplot(heatmap_plot)

    # Calculate average value differences for group-level heatmap
    if (!is.null(strata.var)){
      average_value_diff_matrix <- combined_data %>%
        select(feature.level, !!sym(subject.var), value_diff) %>%
        dplyr::left_join(sorted_meta_tab, by = subject.var) %>%
        dplyr::group_by(!!sym(feature.level), !!sym(group.var), !!sym(strata.var)) %>%
        dplyr::summarize(mean_value_diff = mean(value_diff), .groups = 'drop') %>%
        tidyr::pivot_wider(names_from = all_of(c(group.var, strata.var)), values_from = mean_value_diff) %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      # Create annotation dataframe for the average heatmap
      create_annotation_df <- function(colnames_vec, delimiter = "_") {
        # Extract group and strata variables from column names
        split_names <- strsplit(colnames_vec, delimiter)

        # Ensure all column names have the same number of delimiters
        if (length(unique(sapply(split_names, length))) != 1) {
          stop("All column names must have the same number of delimiters.")
        }

        # Create data frame from split names
        annotation_df <- do.call(rbind.data.frame, split_names)
        colnames(annotation_df) <- c(group.var, strata.var)

        # Convert to factors with unique levels
        annotation_df[] <- lapply(annotation_df, function(x) factor(x, levels = unique(x)))

        # Convert to Matrix
        annotation_matrix <- as.data.frame(annotation_df)
        rownames(annotation_matrix) <- colnames_vec

        return(annotation_matrix)
      }

      annotation_df <- create_annotation_df(colnames(average_value_diff_matrix))

      # Sort annotation dataframe by strata variable
      annotation_df <- annotation_df[order(annotation_df[,strata.var]),]
    } else if (!is.null(group.var)){
      # Calculate average value differences for group-level heatmap without strata
      average_value_diff_matrix <- combined_data %>%
        select(feature.level, !!sym(subject.var), value_diff) %>%
        dplyr::left_join(sorted_meta_tab, by = subject.var) %>%
        dplyr::group_by(!!sym(feature.level), !!sym(group.var)) %>%
        dplyr::summarize(mean_value_diff = mean(value_diff), .groups = 'drop') %>%
        tidyr::pivot_wider(names_from = all_of(c(group.var)), values_from = mean_value_diff) %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      # Create simple annotation dataframe for group variable
      annotation_df <- data.frame(Var1 = colnames(average_value_diff_matrix)) %>%
        dplyr::mutate(Var2 = Var1) %>%
        column_to_rownames("Var2") %>%
        dplyr::rename(!!sym(group.var) := Var1)
    } else {
      # Calculate overall average value differences when no group or strata variables are provided
      average_value_diff_matrix <- combined_data %>%
        select(feature.level, !!sym(subject.var), value_diff) %>%
        dplyr::left_join(sorted_meta_tab, by = subject.var) %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarize(mean_value_diff = mean(value_diff), .groups = 'drop') %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      annotation_df <- NULL
    }

    # Filter features for plotting if specified
    if (!is.null(features.plot)) {
      average_value_diff_matrix <-
        average_value_diff_matrix[rownames(average_value_diff_matrix) %in% features.plot,]
    }

    # Reorder columns of average_value_diff_matrix to match annotation_df if applicable
    if (!is.null(group.var) | !is.null(strata.var)){
      average_value_diff_matrix <- average_value_diff_matrix[, rownames(annotation_df)]
    }

    # Generate the heatmap for average data
    average_heatmap_plot <- pheatmap::pheatmap(
      mat = average_value_diff_matrix[order(rowMeans(abs(average_value_diff_matrix), na.rm = TRUE), decreasing = TRUE), ],
      annotation_col = annotation_df,
      annotation_colors = annotation_colors_list,
      cluster_rows = cluster.rows,
      cluster_cols = cluster.cols,
      annotation_legend = TRUE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      border_color = NA,
      silent = TRUE,
      gaps_col = NULL,
      fontsize = base.size,
      color = my_col,
      breaks = break_points
    )

    gg_average_heatmap_plot <- as.ggplot(average_heatmap_plot)

    # Convert feature.change.func to string if it's a custom function
    if (is.function(feature.change.func)) {
      feature.change.func = "custom function"
    }

    # Save individual heatmap as PDF if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_pair_indiv",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "change_base_",
        change.base,
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

    # Save average heatmap as PDF if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_heatmap_pair_average",
        "_",
        "subject_",
        subject.var,
        "_",
        "time_",
        time.var,
        "_",
        "change_base_",
        change.base,
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
        plot = gg_average_heatmap_plot
      )
    }
    
    # Create a list of both heatmap plots
    sub.plot_list <- list(gg_average_heatmap_plot, gg_heatmap_plot)
    names(sub.plot_list) <- c("average", "indiv")

    return(sub.plot_list)
  })

  # Name the plot list with feature levels and return
  names(plot_list) <- feature.level
  return(plot_list)
}
