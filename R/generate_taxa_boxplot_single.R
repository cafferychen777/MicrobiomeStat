#' Generate Taxa Boxplots for Single Time Point
#'
#' Creates boxplots showing taxa abundance distributions at a single time point.
#' Supports grouping, stratification, and various transformations.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param t.level Character string specifying the time level to subset data to.
#'   Default NULL uses all data.
#' @param features.plot Character vector of specific feature IDs to plot.
#'   If NULL, features are selected based on top.k.plot and top.k.func.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#' @param transform Transformation to apply: "identity", "sqrt", or "log".
#' @param point.alpha Numeric value (0-1) for jitter point transparency. Default 0.6.
#'
#' @return A list of ggplot objects for each taxonomic level.
#'
#' @examples
#' \dontrun{
#' # Generate the boxplot pair
#' data(ecam.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = "diet",
#'   strata.var = NULL,
#'   feature.level = c("Phylum"),
#'   features.plot = sample(unique(ecam.obj$feature.ann[,"Phylum"]),3),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   point.alpha = 0.6,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Phylum"),
#'   features.plot = sample(unique(ecam.obj$feature.ann[,"Phylum"]),3),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   point.alpha = 0.6,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = NULL,
#'   strata.var = NULL,
#'   feature.level = c("Order", "Phylum", "Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   point.alpha = 0.6,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = NULL,
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   base.size = 12,
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
generate_taxa_boxplot_single <-
  function(data.obj,
           time.var = NULL,
           t.level = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           transform = c("sqrt", "identity", "log"),
           prev.filter = 0.05,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           point.alpha = 0.6,
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
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    # Create a subject identifier from sample names (rownames)
    data.obj$meta.dat$subject.id <- rownames(data.obj$meta.dat)
    subject.var <- "subject.id"

    # Subset data if time variable and level are provided
    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <-
        mStat_subset_data(data.obj, condition = condition)
    }

    # If group variable is not provided, create a default one
    if (is.null(group.var)) {
      data.obj$meta.dat <- data.obj$meta.dat %>% dplyr::mutate("ALL" = "ALL")
      group.var <- "ALL"
    }

    # Extract metadata
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        subject.var, group.var, time.var, strata.var
      )))

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
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <-
        mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    # Generate plots for each taxonomic level
    plot_list <- lapply(feature.level, function(feature.level) {

      # Define aesthetic mapping for the plot
      aes_function <- aes(
        x = !!sym(feature.level),
        y = value,
        fill = !!sym(group.var)
      )

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) &
          feature.level != "original") {
        data.obj <-
          mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Extract aggregated feature table
      if (feature.level != "original") {
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
      if (feature.dat.type %in% c("count", "proportion")) {
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
        otu_tax_agg_merged %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # Count number of subjects
      n_subjects <-
        length(unique(otu_tax_agg_merged[[subject.var]]))

      sub_otu_tax_agg_merged <- otu_tax_agg_merged

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
        ggplot(sub_otu_tax_agg_merged %>% filter(!!sym(feature.level) %in% features.plot),
               aes_function) +
        # Add error bars
        stat_boxplot(
          geom = "errorbar",
          position = position_dodge(width = 0.8),
          width = 0.2
        ) +
        # Add boxplots without outliers
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.2,
          outlier.shape = NA,  # Hide outliers since we show all points
          alpha = 0.8  # Make boxplot slightly transparent so points are visible
        ) +
        # Draw jitter points last so they appear on top of boxplots
        geom_jitter(
          aes(fill = !!sym(group.var)),  # Need fill aesthetic for dodging
          position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
          alpha = point.alpha,
          size = 2,
          color = "gray30"
        ) +
        scale_fill_manual(values = col) +
        {
          if (feature.dat.type == "other") {
            labs(y = "Abundance")
          } else {
            labs(y = paste("Relative Abundance(", transform, ")"))
          }
        } +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          strip.text.x = element_text(size = base.size, color = "black"),
          strip.text.y = element_text(size = base.size, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = (base.size -
                                                                                   2)),
          axis.text.y = element_text(color = "black", size = (base.size -
                                                                2)),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      # Add facets if strata variable is specified
      if (!is.null(group.var)) {
        if (is.null(strata.var)) {
        } else {
          boxplot <- boxplot + ggh4x::facet_nested_wrap(
            as.formula(paste('~', strata.var)),
            scales = "free_x"
          )
        }
      }

      # Modify y-axis scale based on transformation
      if (feature.dat.type != "other") {
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
        "taxa_boxplot_single",
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

    # Name the plot list with feature levels
    names(plot_list) <- feature.level

    # Return the list of plots
    return(plot_list)
  }