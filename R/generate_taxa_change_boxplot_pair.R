#' Generate Taxa Change Boxplots for Paired Data
#'
#' Creates boxplots showing changes in taxa abundance between paired time points.
#' Supports multiple change metrics and grouping/stratification.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param features.plot Character vector of specific feature IDs to plot.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#'
#' @return A list of ggplot objects showing abundance changes.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and data
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' # Generate the boxplot pair
#' generate_taxa_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   feature.level = "original",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 6,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#'
#' # Generate the boxplot pair
#' generate_taxa_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "relative change",
#'   feature.level = "original",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 6,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_change_boxplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           feature.change.func = "relative change",
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter = 0,
           abund.filter = 0,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    # Validate input data
    mStat_validate_data(data.obj)

    # Match and validate the feature data type
    feature.dat.type <- match.arg(feature.dat.type)

    # Validate input variables
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

    # Extract metadata
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      ))) %>% rownames_to_column("sample")

    # Get plot theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Get color palette
    col <- mStat_get_palette(palette)

    # Set y-axis label based on feature data type and change function
    ylab_label <- if (feature.dat.type != "other") {
      if (is.function(feature.change.func)) {
        paste0("Change in Relative Abundance", " (custom function)")
      } else {
        paste0("Change in Relative Abundance", " (", feature.change.func, ")")
      }
    }
    else {
      if (is.function(feature.change.func)) {
        paste0("Change in Abundance", " (custom function)")
      } else {
        paste0("Change in Abundance", " (", feature.change.func, ")")
      }
    }

    # Adjust filtering parameters if necessary
    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    # Generate plots for each taxonomic level
    plot_list <- lapply(feature.level, function(feature.level) {

      # Normalize count data if necessary
      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      }

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Extract aggregated feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Filter feature table based on prevalence and abundance
      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        tibble::rownames_to_column(feature.level)

      # Select top k features if specified
      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Reshape data for plotting
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value", -all_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Merge OTU data with metadata
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab, by = "sample") %>%
        select(all_of(
          c(
            "sample",
            feature.level,
            subject.var,
            time.var,
            group.var,
            strata.var,
            "value"
          )
        ))

      # Identify the time point after the change base
      change.after <-
        unique(otu_tax_agg_merged %>% select(all_of(c(time.var))))[unique(otu_tax_agg_merged %>% select(all_of(c(time.var)))) != change.base]

      # Split data into separate tibbles for each time point
      split_data <-
        split(otu_tax_agg_merged,
              f = otu_tax_agg_merged %>%
                dplyr::group_by(!!sym(time.var)) %>% select(all_of(c(time.var))))

      # Extract data for the base time point and the change time point
      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # Combine data from both time points
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # Calculate the change in abundance based on the specified function
      combined_data <- combined_data %>%
        dplyr::mutate(value_diff = compute_taxa_change(
          value_after  = value_time_2,
          value_before = value_time_1,
          method       = feature.change.func,
          feature_id   = .data[[feature.level]]
        ))

      # Add metadata for the change time point
      combined_data <-
        combined_data %>% dplyr::left_join(meta_tab %>% filter(!!sym(time.var) == change.after), by = subject.var)

      # Get unique taxa levels
      taxa.levels <-
        combined_data %>% select(all_of(c(feature.level))) %>% dplyr::distinct() %>% dplyr::pull()

      # Set default group if not provided
      if (is.null(group.var)) {
        group.var = "group"
        combined_data$group <- "ALL"
      }

      # Select features to plot
      if (!is.null(features.plot)) {
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      } else {
        if (length(taxa.levels) >= 5) {
          taxa.levels <- taxa.levels[1:4]
        } else {
          taxa.levels <- taxa.levels
        }
      }

      # Create the boxplot
      boxplot <-
        ggplot(
          combined_data %>% filter(!!sym(feature.level) %in% taxa.levels),
          aes(
            x = !!sym(group.var),
            y = value_diff,
            fill = !!sym(group.var)
          )
        ) +
        geom_jitter(width = 0.1,
                    alpha = 0.3,
                    size = 1.5) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
        ) +
        scale_alpha_manual(values = c(0.5, 0.5)) +
        scale_fill_manual(values = col) +
        labs(x = group.var,
             y = ylab_label) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          plot.title = element_text(hjust = 0.5, size = 20),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16),
          ...
        )

      # Add facets if strata variable is provided
      if (!is.null(strata.var)) {
        boxplot <- boxplot +
          ggh4x::facet_nested_wrap(
            as.formula(paste("~",feature.level,"+",strata.var)),
            scales = "fixed",
            nrow = ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 2, (length(taxa.levels) + 1) / 2)
          )
      } else {
        boxplot <- boxplot +
          ggh4x::facet_nested_wrap(
            as.formula(paste("~",feature.level)),
            scales = "fixed",
            nrow = ifelse(length(taxa.levels) %% 2 == 0, length(taxa.levels) / 2, (length(taxa.levels) + 1) / 2)
          )
      }

      # Adjust theme if there's only one group
      if (group.var == "group" &&
          unique(combined_data$group)[1] == "ALL") {
        boxplot <- boxplot +
          theme(
            legend.position = "none",
            axis.title.x = element_blank(),
            axis.text.x = element_blank()
          )
      }

      return(boxplot)

    })

    # Name the plot list with feature levels
    names(plot_list) <- feature.level

    # Save the plots as a PDF file if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_boxplot_pair",
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
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Use lapply to print each ggplot object in the list to a new PDF page
      lapply(plot_list, print)
      # Close the PDF device
      dev.off()
    }

    # Return the boxplot for display
    return(plot_list)
  }