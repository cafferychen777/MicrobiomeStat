#' Generate Taxa Boxplots for Longitudinal Data
#'
#' Creates boxplots showing taxa abundance distributions across time points.
#' Supports grouping, stratification, and various transformations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param features.plot Character vector of specific feature IDs to plot.
#'   If NULL, features are selected based on top.k.plot and top.k.func.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#' @param transform Transformation to apply: "identity", "sqrt", or "log".
#'
#' @return A list of ggplot objects for each taxonomic level.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' data(ecam.obj)
#' # Generate the boxplot pair
#' generate_taxa_boxplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = NULL,
#'   feature.level = c("Class"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   top.k.plot = 1,
#'   top.k.func = "sd",
#'   transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 1e-7,
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 25,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_boxplot_long(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = "Family",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   transform = "sqrt",
#'   prev.filter = 0.001,
#'   abund.filter = 0.001,
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 13
#' )
#' }
#' @export
generate_taxa_boxplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
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

    # Match and validate the feature data type
    feature.dat.type <- match.arg(feature.dat.type)

    # Match and validate the transformation method
    transform <- match.arg(transform)

    # Validate the input data object
    mStat_validate_data(data.obj)

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

    # Process time variable in the data object
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

    # Extract metadata and select relevant variables
    meta_tab <- data.obj$meta.dat %>%
      as.data.frame() %>%
      select(all_of(c(subject.var,group.var,time.var,strata.var)))

    # Define aesthetic mapping function for line plots
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

    # Define aesthetic mapping function for box plots
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

    # Get color palette
    col <- mStat_get_palette(palette)

    # Get plot theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Adjust filtering parameters if necessary
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
    plot_list <- lapply(feature.level, function(feature.level) {

      # Aggregate data by taxonomy if necessary
      otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

      # Select top k features if specified
      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # Reshape data for plotting
      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value",-all_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      # Merge OTU data with metadata
      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric,
                  meta_tab %>% rownames_to_column("sample"),
                  by = "sample") %>%
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

      # Apply data transformation if necessary
      if (feature.dat.type %in% c("count","proportion")){
        if (transform %in% c("identity", "sqrt", "log")) {
          if (transform == "identity") {
            # No transformation needed
          } else if (transform == "sqrt") {
            otu_tax_agg_merged$value <- sqrt(otu_tax_agg_merged$value)
          } else if (transform == "log") {
            # Find the half of the minimum non-zero proportion for each taxon
            min_half_nonzero <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::summarise(min_half_value = min(value[value > 0]) / 2) %>%
              dplyr::ungroup()
            # Replace zeros with the log of the half minimum non-zero proportion
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
        otu_tax_agg_merged %>% select(all_of(c(feature.level))) %>% dplyr::distinct() %>% dplyr::pull()

      # Count number of subjects and time points
      n_subjects <-
        length(unique(otu_tax_agg_merged[[subject.var]]))
      n_times <- length(unique(otu_tax_agg_merged[[time.var]]))

      sub_otu_tax_agg_merged <- otu_tax_agg_merged

      # Calculate average values if number of subjects or time points is large
      average_sub_otu_tax_agg_merged <- NULL
      if (n_times > 10 || n_subjects > 25) {
        if (!is.null(group.var) & !is.null(strata.var)) {
          average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
            dplyr::group_by(
              !!sym(strata.var),
              !!sym(group.var),
              !!sym(time.var),
              !!sym(feature.level)
            ) %>%
            dplyr::summarise(dplyr::across(value, \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
            dplyr::ungroup() %>%
            dplyr::mutate(!!sym(subject.var) := "ALL")
        } else if (!is.null(group.var)) {
          average_sub_otu_tax_agg_merged <- sub_otu_tax_agg_merged %>%
            dplyr::group_by(!!sym(group.var),
                     !!sym(time.var),
                     !!sym(feature.level)) %>%
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

      # Select features to plot
      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 6) {
          features.plot <- taxa.levels[1:6]
        } else {
          features.plot <- taxa.levels
        }
      }

      # Get number of group levels if group variable is specified
      if (!is.null(group.var)){
        group.levels <- sub_otu_tax_agg_merged %>% select(!!sym(group.var)) %>% pull() %>% unique() %>% length()
      }

      # Create box plot
      boxplot <-
        ggplot(sub_otu_tax_agg_merged  %>%
                 dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var))) %>%
                 filter(!!sym(feature.level) %in% features.plot),
               aes_function) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
        ) +
        geom_line(
          line_aes_function,
          alpha = 1,
          linewidth = 0.6,
          color = "black",
          linetype = "dashed",
          data = if (!is.null(average_sub_otu_tax_agg_merged))
            average_sub_otu_tax_agg_merged   %>% dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var))) %>% filter(!!sym(feature.level) %in% features.plot)
          else
            sub_otu_tax_agg_merged  %>% dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var))) %>% filter(!!sym(feature.level) %in% features.plot)
        ) +
        scale_fill_manual(values = col) +
        {
          if (feature.dat.type == "other"){
            labs(x = time.var,
                 y = "Abundance")
          } else {
            labs(x = time.var,
                 y = paste("Relative Abundance(", transform, ")"))
          }
        } +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          plot.title = element_text(hjust = 0.5, size = 20),
          strip.text.x = element_text(size = base.size*0.8, color = "black"),
          strip.text.y = element_text(size = base.size*0.8, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = (base.size -
                                                              2)),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      # Add facets if group or strata variables are specified
      if (!is.null(group.var)) {
        if (is.null(strata.var)) {
          boxplot <-
            boxplot + ggh4x::facet_nested_wrap(as.formula(paste(
              "~", feature.level, "+", group.var
            )), scales = "free_x", ncol = group.levels*4)
        } else {
          boxplot <- boxplot + ggh4x::facet_nested_wrap(
            as.formula(paste('~', feature.level, '+', strata.var, '+', group.var)),
            scales = "free_x", ncol = group.levels*4
          )
        }
      }

      # Add jitter points if number of subjects or time points is large
      if (n_subjects > 20 || n_times > 10) {
        boxplot <- boxplot + geom_jitter(width = 0.1,
                                         alpha = 0.1,
                                         size = 1)
      }

      # Modify y-axis scale based on transformation
      if (feature.dat.type != "other"){
        if (transform == "sqrt") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  a ^ b, list(a = i, b = 2)
                )))
          )
        } else if (transform == "log") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  10 ^ a, list(a = i)
                )))
          )
        }
      }

      return(boxplot)
    })

    # Save plots as PDF if specified
    if (pdf) {
      pdf_name <- paste0(
        "taxa_boxplot_long",
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
      pdf_name <- paste0(pdf_name, "_", feature.level, ".pdf")
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Use lapply to print each ggplot object in the list to a new PDF page
      lapply(plot_list, print)
      # Close the PDF device
      dev.off()
    }

    names(plot_list) <- feature.level

    return(plot_list)
  }
