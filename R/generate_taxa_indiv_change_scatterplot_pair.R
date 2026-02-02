#' Check if a Variable is Continuous Numeric
#'
#' This function checks if a given variable is continuous numeric by checking if it is numeric and has at least 10 unique values.
#'
#' @param x A variable that you want to check.
#'
#' @return Logical. Returns TRUE if the variable is continuous numeric, FALSE otherwise.
#'
#' @examples
#' is_continuous_numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) # TRUE
#' is_continuous_numeric(c(1, 2, 3, 4, 5)) # FALSE
#' is_continuous_numeric(c('a', 'b', 'c')) # FALSE
#'
#' @export
is_continuous_numeric <- function(x) {
  if (is.numeric(x) && length(unique(x)) >= 10) {

    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @title Generate Individual Change Scatterplots for Paired Samples
#'
#' @description Creates scatterplots showing the change in taxonomic composition between two time points
#' in a longitudinal study, useful for visualizing associations with continuous variables.
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
#' This function generates a scatterplot of the change in taxa abundances between two time points in a longitudinal study.
#' The scatterplot can be stratified by a group variable and/or other variables.
#' It also allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#'  library(vegan)
#' data(peerj32.obj)
#' peerj32.obj$meta.dat <- peerj32.obj$meta.dat %>%
#' dplyr::select(all_of("subject")) %>% dplyr::distinct() %>%
#' dplyr::mutate(cons = runif(dplyr::n(),0,5)) %>%
#' dplyr::left_join(peerj32.obj$meta.dat %>% rownames_to_column("sample"),by = "subject") %>%
#' tibble::column_to_rownames("sample")
#'
#'  # Generate the scatterplot pairs
#'  generate_taxa_indiv_change_scatterplot_pair(
#'    data.obj = peerj32.obj,
#'    subject.var = "subject",
#'    time.var = "time",
#'    group.var = "cons",
#'    strata.var = "sex",
#'    change.base = "1",
#'    feature.change.func = "log fold change",
#'    feature.level = "Genus",
#'    top.k.plot = NULL,
#'    top.k.func = NULL,
#'    prev.filter = 0.01,
#'    abund.filter = 0.01
#'  )
#'
#' data("subset_pairs.obj")
#' subset_pairs.obj$meta.dat <- subset_pairs.obj$meta.dat %>%
#' dplyr::select(all_of("MouseID")) %>% dplyr::distinct() %>%
#' dplyr::mutate(cons = runif(dplyr::n(),0,5)) %>%
#' dplyr::left_join(subset_pairs.obj$meta.dat %>% rownames_to_column("sample"),by = "MouseID") %>%
#' tibble::column_to_rownames("sample")
#'
#'  # Generate the scatterplot pairs
#'  generate_taxa_indiv_change_scatterplot_pair(
#'    data.obj = subset_pairs.obj,
#'    subject.var = "MouseID",
#'    time.var = "Antibiotic",
#'    group.var = "cons",
#'    strata.var = NULL,
#'    change.base = "Baseline",
#'    feature.change.func = "log fold change",
#'    feature.level = "Genus",
#'    top.k.plot = NULL,
#'    top.k.func = NULL,
#'    prev.filter = 0.01,
#'    abund.filter = 0.01
#'  )
#' }
#' @export
generate_taxa_indiv_change_scatterplot_pair <-
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

    # Extract metadata
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        time.var, group.var, strata.var, subject.var
      )))

    # Determine the time point after the base time point
    change.after <-
      unique(meta_tab %>% select(all_of(c(time.var))))[unique(meta_tab %>% select(all_of(c(time.var)))) != change.base]

    # If group variable is not provided, create a dummy group
    if (is.null(group.var)) {
      group.var = "ALL"
      meta_tab$ALL <- ""
    }

    # If strata variable is not provided, create a dummy strata
    if (is.null(strata.var)) {
      strata.var = "ALL2"
      meta_tab$ALL2 <- ""
    }

    # Get color palette
    colors <- mStat_get_palette(palette)

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Define aesthetic function based on whether strata variable is provided
    aes_function <- if (!is.null(strata.var)){
      aes(shape = !!sym(strata.var), color = !!sym(strata.var))
    } else {
      aes(color = !!sym(time.var))
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

      # Convert counts to numeric type
      otu_tax_agg_numeric <-
        dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      # Create normalized feature table
      otu_tab_norm <- otu_tax_agg_numeric %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      # Reshape data for plotting
      otu_tab_norm_agg <- otu_tax_agg_numeric %>%
        tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
        dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      # Get unique taxa levels
      taxa.levels <-
        otu_tab_norm_agg %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # Split data into two subsets: one for change.base, one for change.after
      df_t0 <- otu_tab_norm_agg %>% filter(!!sym(time.var) == change.base)
      df_ts <- otu_tab_norm_agg %>% filter(!!sym(time.var) == change.after)

      # Join the two subsets
      df <- dplyr::inner_join(df_ts, df_t0, by = c(feature.level, subject.var), suffix = c("_ts", "_t0"), relationship = "many-to-many")

      # Calculate the change in abundance based on the specified function
      df <- df %>%
        dplyr::mutate(new_count = compute_taxa_change(
          value_after  = count_ts,
          value_before = count_t0,
          method       = feature.change.func,
          feature_id   = .data[[feature.level]]
        ))

      # Join with metadata
      df <- df %>% dplyr::left_join(meta_tab %>% select(-all_of(time.var)) %>% dplyr::distinct(), by = c(subject.var))

      # Rename columns
      df <- df %>% setNames(ifelse(names(.) == paste0(time.var,"_ts"), time.var, names(.)))

      # Define y-axis label based on feature data type and change function
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

      # Filter taxa levels if specific features are requested
      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      # Create a plot for each taxon
      plot_list <- lapply(taxa.levels, function(tax) {
        scatterplot <-
          ggplot(df %>% filter(!!sym(feature.level) == tax),
                 aes(
                   x = !!sym(group.var),
                   y = new_count,
                   fill = !!sym(strata.var),
                   color = !!sym(strata.var),
                   group = !!sym(strata.var)
                 )) +
          geom_smooth(se = TRUE, method = 'lm') +
          geom_point(aes_function,data = df %>% filter(!!sym(feature.level) == tax),
                     size = 4) +
          scale_shape_manual(values = c(21, 22, 24, 25)) +
          scale_fill_manual(values = colors) +
          ylab(ylab_label) +
          ggtitle(tax) +
          scale_color_manual(values = colors, guide = guide_legend(override.aes = list(size = 0))) +
          theme_to_use +
          theme(
            axis.text.x = element_text(
              vjust = 0.5,
              hjust = 1,
              size = base.size
            ),
            axis.text.y = element_text(color = "black", size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.title.x = element_text(size = base.size),
            legend.position = "right",
            legend.direction = "vertical",
            legend.box = "vertical",
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            axis.text = ggplot2::element_text(color = "black",
                                              size = base.size),
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_text(size = 16),
            plot.title = element_text(hjust = 0.5, size = 20)
          )

        # Adjust theme for single group case
        if (group.var == "ALL"){
          scatterplot <- scatterplot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
        }

        # Adjust theme for single strata case
        if (strata.var == "ALL2") {
          scatterplot <- scatterplot + theme(
            legend.title = element_blank()
          )
        }

        return(scatterplot)
      })

      # Save the plots as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_change_scatterplot",
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
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
