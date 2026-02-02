#' @title Generate Individual Change Boxplots for Paired Samples
#'
#' @description Creates boxplots showing the change in taxonomic composition between two time points
#' in a longitudinal study, with options for grouping and stratification.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param change.base A string indicating the base time point for change computation.
#' @param feature.change.func Method for computing change: "absolute change", "log fold change",
#'   "relative change", or a custom function.
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of ggplot objects, one for each taxonomic level.
#' @details
#' This function generates a boxplot of the change in taxa abundances between two time points in a longitudinal study.
#' The boxplot can be stratified by a group variable and/or other variables.
#' It allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and data
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' # Generate the boxplot pair
#' generate_taxa_indiv_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 20,
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
#' generate_taxa_indiv_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 20,
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
generate_taxa_indiv_change_boxplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           feature.change.func = "relative change",
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter = 0.01,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    # Validate the input data object
    mStat_validate_data(data.obj)

    # Match the feature data type argument
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

    # Extract metadata and add sample column
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      ))) %>% rownames_to_column("sample")

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Get the color palette for plotting
    col <- mStat_get_palette(palette)

    # Determine the y-axis label based on the feature data type and change function
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

    # Adjust filters if necessary
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

    # Process each feature level
    plot_list_all <- lapply(feature.level, function(feature.level) {

      # Aggregate data by taxonomy if necessary
      otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

      # Select top k features if specified
      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
        computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
        features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Reshape the data for plotting
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value",-all_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Merge feature data with metadata
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab, by = "sample") %>%
        select(all_of(c("sample",
                        feature.level,
                        subject.var,
                        time.var,
                        group.var,
                        strata.var,
                        "value")))

      # Determine the time point after the base time point
      change.after <-
        unique(otu_tax_agg_merged %>% select(all_of(c(time.var))))[unique(otu_tax_agg_merged %>% select(all_of(c(time.var)))) != change.base]

      # Split data by time points
      split_data <-
        split(otu_tax_agg_merged,
              f = otu_tax_agg_merged %>%
                dplyr::group_by(!!sym(time.var)) %>% select(all_of(c(time.var))))

      # Extract data for the two time points
      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # Join data from the two time points
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # Calculate the change in abundance
      combined_data <- combined_data %>%
        dplyr::mutate(value_diff = compute_taxa_change(
          value_after  = value_time_2,
          value_before = value_time_1,
          method       = feature.change.func,
          feature_id   = .data[[feature.level]]
        ))

      # Merge with metadata for the second time point
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

      # Filter taxa levels if specific features are requested
      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      # Create a plot for each taxon
      plot_list <- lapply(taxa.levels, function(tax) {
        # Create the boxplot
        boxplot <-
          ggplot(
            combined_data %>% filter(!!sym(feature.level) == tax),
            aes(
              x = !!sym(group.var),
              y = value_diff,
              fill = !!sym(group.var)
            )
          ) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.3
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.3,
          ) +
          geom_jitter(width = 0.3,
                      alpha = 0.3,
                      size = 1.7) +
          scale_alpha_manual(values = c(0.5, 0.5)) +
          scale_fill_manual(values = col) +
          labs(
            x = group.var,
            y = ylab_label,
            title = tax
          ) +
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
            legend.title = ggplot2::element_text(size = 16)
          )

        # Add faceting if strata variable is provided
        if (!is.null(strata.var)) {
          boxplot <- boxplot +
            ggh4x::facet_nested(cols = vars(!!sym(strata.var)),
                                scales = "fixed",
                                space = "free")
        }

        # Adjust theme for single group case
        if (group.var == "group" && unique(combined_data$group)[1] == "ALL") {
          boxplot <- boxplot +
            theme(
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_blank()
            )
        }

        return(boxplot)
      })

      # Save plots as PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_change_boxplot_pair",
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
        # Print each plot to a new page in the PDF
        lapply(plot_list, print)
        # Close the PDF device
        dev.off()
      }

      # Name the plots in the list
      names(plot_list) <- taxa.levels

      return(plot_list)
    })
    # Name the list of plot lists
    names(plot_list_all) <- feature.level
    # Return the complete list of plots
    return(plot_list_all)
  }
