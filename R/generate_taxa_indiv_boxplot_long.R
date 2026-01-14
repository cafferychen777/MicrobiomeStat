#' @title Generate Individual Taxa Boxplots for Longitudinal Data
#'
#' @description Creates boxplots showing the abundance distribution of individual taxa
#' at specified taxonomic levels over time from longitudinal data.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to compute the top k taxa. Default is NULL (uses mean).
#' @param transform Transformation to apply: "identity" (default), "sqrt", or "log".
#' @param ... Additional arguments passed to ggplot2 functions.
#'
#' @return A ggplot object showing the abundance distribution of taxa over time.
#'
#' @examples
#' \dontrun{
#' # Generate the boxplot pair
#' data(ecam.obj)
#' generate_taxa_indiv_boxplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = NULL,
#'   feature.level = c("Phylum"),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#'   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8",
#'   "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2",
#'   "#c7c7c7", "#dbdb8d", "#9edae5", "#f0f0f0", "#3182bd"),
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 25,
#'   pdf.hei = 8.5
#' )
#'
#' data(peerj32.obj)
#'
#' generate_taxa_indiv_boxplot_long(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Family"),
#'   features.plot = NULL,
#'   feature.dat.type = "other",
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
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
generate_taxa_indiv_boxplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           top.k.plot = NULL,
           top.k.func = NULL,
           transform = c("sqrt", "identity", "log"),
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

    # Match and validate input arguments
    feature.dat.type <- match.arg(feature.dat.type)
    transform <- match.arg(transform)

    # Validate the input data object
    mStat_validate_data(data.obj)

    # Check if input variables are properly specified
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

    # Process time variable and extract relevant data
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    meta_tab <- data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

    # Define aesthetic functions for plotting
    # These functions determine how the data will be mapped to visual properties in the plot
    line_aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = value,
        group = !!sym(subject.var),
        color = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = value,
        group = !!sym(subject.var)
      )
    }

    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = value,
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = value,
        fill = !!sym(time.var)
      )
    }

    # Get color palette for the plot
    col <- mStat_get_palette(palette)

    # Get the appropriate theme for the plot
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Adjust filtering parameters based on input conditions
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

    # Generate plots for each taxonomic level
    plot_list_all <- lapply(feature.level, function(feature.level) {

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

      # Reshape the data for plotting
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value", -all_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Merge feature data with metadata
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric, meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
        select(all_of(c("sample",
                        feature.level,
                        subject.var,
                        time.var,
                        group.var,
                        strata.var,
                        "value")))

      # Apply data transformation if specified
      if (feature.dat.type %in% c("count","proportion")){
        if (transform %in% c("identity", "sqrt", "log")) {
          if (transform == "identity") {
            # No transformation needed
          } else if (transform == "sqrt") {
            otu_tax_agg_merged$value <- sqrt(otu_tax_agg_merged$value)
          } else if (transform == "log") {
            # For log transformation, we need to handle zeros
            # We replace zeros with half of the minimum non-zero value for each taxon
            min_half_nonzero <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::summarise(min_half_value = min(value[value > 0]) / 2) %>%
              dplyr::ungroup()
            otu_tax_agg_merged <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::ungroup() %>%
              dplyr::left_join(min_half_nonzero, by = feature.level) %>%
              dplyr::mutate(value = ifelse(value == 0, log10(min_half_value), log10(value))) %>%
              select(-min_half_value)
          }
        }
      }

      # Get unique taxa levels
      taxa.levels <-
        otu_tax_agg_merged %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()

      # Count number of subjects and time points
      n_subjects <- length(unique(otu_tax_agg_merged[[subject.var]]))
      n_times <- length(unique(otu_tax_agg_merged[[time.var]]))

      # Filter taxa levels if specified
      if (!is.null(features.plot)){
        taxa.levels <- taxa.levels[taxa.levels %in% features.plot]
      }

      # Generate individual plots for each taxon
      plot_list <- lapply(taxa.levels, function(tax) {

        sub_otu_tax_agg_merged <- otu_tax_agg_merged %>% filter(!!sym(feature.level) == tax)

        # Calculate average values if there are many subjects or time points
        average_sub_otu_tax_agg_merged <- NULL
        if (n_times > 10 || n_subjects > 25) {
          if (!is.null(group.var) && !is.null(strata.var)) {
            average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
              dplyr::summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
              dplyr::ungroup() %>%
              dplyr::mutate(!!sym(subject.var) := "ALL")
          } else if (!is.null(group.var)) {
            average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(group.var), !!sym(time.var)) %>%
              dplyr::summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
              dplyr::ungroup() %>%
              dplyr::mutate(!!sym(subject.var) := "ALL")
          } else {
            average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(time.var)) %>%
              dplyr::summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
              dplyr::ungroup() %>%
              dplyr::mutate(!!sym(subject.var) := "ALL")
          }
        }

        # Create the boxplot
        boxplot <-
          ggplot(sub_otu_tax_agg_merged  %>%
                   dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var))),
                 aes_function) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.3
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.3,
          ) +
          geom_line(
            line_aes_function,
            alpha = 1,
            linewidth = 0.6,
            color = "black",
            linetype = "dashed",
            data = if (!is.null(average_sub_otu_tax_agg_merged)) average_sub_otu_tax_agg_merged %>%
              dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var))) else sub_otu_tax_agg_merged %>%
              dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
          ) +
          scale_fill_manual(values = col) +
          {
            if (feature.dat.type == "other"){
              labs(
                x = time.var,
                y = "Abundance",
                title = tax
              )
            } else {
              labs(
                x = time.var,
                y = paste("Relative Abundance(", transform, ")"),
                title = tax
              )
            }
          } +
          theme_to_use +
          theme(
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0, "cm"),
            plot.title = element_text(hjust = 0.5, size = 20),
            strip.text.x = element_text(size = 12, color = "black"),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(color = "black", size = base.size),
            axis.text.y = element_text(color = "black", size = (base.size-2)),
            axis.title.x = element_text(size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_text(size = 16)
          )

        # Add faceting if group or strata variables are provided
        if (!is.null(group.var)) {
          if (is.null(strata.var)) {
            boxplot <-
              boxplot + ggh4x::facet_nested(as.formula(paste("~", group.var)), scales = "fixed")
          } else {
            boxplot <-
              boxplot + ggh4x::facet_nested(as.formula(paste("~", strata.var, "+", group.var)), scales = "free", space = "free") + theme(panel.spacing = unit(0,"lines"))
          }
        }

        # Add jitter points if there are many subjects or time points
        if (n_subjects > 20 || n_times > 10) {
          boxplot <- boxplot + geom_jitter(width = 0.1, alpha = 0.1, size = 1)
        }

        # Modify y-axis scale based on the transformation
        if (feature.dat.type != "other"){
          if (transform == "sqrt") {
            boxplot <- boxplot + scale_y_continuous(
              labels = function(x) sapply(x, function(i) as.expression(substitute(a^b, list(a = i, b = 2))))
            )
          } else if (transform == "log") {
            boxplot <- boxplot + scale_y_continuous(
              labels = function(x) sapply(x, function(i) as.expression(substitute(10^a, list(a = i))))
            )
          }
        }

        return(boxplot)
      })

      # Save the plots as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_boxplot_long",
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
          "transform_",
          transform,
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
        pdf_name <- paste0(pdf_name,"_", feature.level, ".pdf")
        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Use lapply to print each ggplot object in the list to a new PDF page
        lapply(plot_list, print)
        # Close the PDF device
        dev.off()
      }
      names(plot_list) <- taxa.levels
      return(plot_list)
    })

    names(plot_list_all) <- feature.level
    return(plot_list_all)
  }