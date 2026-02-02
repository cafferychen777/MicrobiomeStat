#' @title Generate Taxa Dotplot for Paired Samples
#'
#' @description Generates a stacked dotplot of taxa abundance for paired samples, showing mean
#' abundance and prevalence across groups and time points.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#' @examples
#' \dontrun{
#' # Load required libraries
#' data(peerj32.obj)
#'
#' # Call the function
#' generate_taxa_dotplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 11
#' )
#'
#' data("subset_pairs.obj")
#'
#' # Call the function
#' generate_taxa_dotplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.001,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 11
#' )
#' }
#' @export
generate_taxa_dotplot_pair <- function(data.obj,
                                       subject.var,
                                       time.var,
                                       group.var = NULL,
                                       strata.var = NULL,
                                       feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       features.plot = NULL,
                                       top.k.plot = NULL,
                                       top.k.func = NULL,
                                       prev.filter = 0.01,
                                       abund.filter = 0.01,
                                       base.size = 16,
                                       theme.choice = "bw",
                                       custom.theme = NULL,
                                       palette = c("white", "#92c5de", "#0571b0", "#f4a582", "#ca0020"),
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       ...) {

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Validate the input data object
  mStat_validate_data(data.obj)

  # Capture original factor levels before any data extraction
  fl <- mStat_capture_factor_levels(data.obj, group.var, strata.var)
  data.obj <- fl$data.obj

  # Extract relevant variables from the metadata
  meta_tab <-
    data.obj$meta.dat %>% select(all_of(c(
      time.var, group.var, strata.var, subject.var
    )))

  # If no group variable is provided, create a dummy "ALL" group
  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- ""
  }

  # If a strata variable is provided, create an interaction term with the group variable
  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var), !!sym(strata.var), sep = .STRATA_SEP))
  }

  # Get the color palette
  colors <- mStat_get_palette(palette)

  # Get the appropriate theme
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Determine whether to apply abundance and prevalence filters
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

  # Create a list to store plots for each taxonomic level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
      data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
    }

    # Get the appropriate feature table
    if (feature.level != "original"){
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

    # Convert counts to numeric type
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    # Calculate mean abundance for each group and time point
    otu_tab_norm_agg <- otu_tax_agg_numeric %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level),!!sym(time.var)) %>%
      dplyr::summarise(mean_abundance = sqrt(mean(count)))

    # Calculate prevalence for each group and time point
    prevalence_all <- otu_tax_agg_numeric %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level),!!sym(time.var)) %>%
      dplyr::summarise(prevalence = sum(count > 0) / dplyr::n()) %>% as.data.frame()

    # Merge mean abundance and prevalence data
    otu_tab_norm_agg <-
      otu_tab_norm_agg %>% dplyr::left_join(prevalence_all, c(feature.level,group.var,time.var))

    # Calculate the midpoint of mean abundance for color scaling
    midpoint <- quantile(otu_tab_norm_agg$mean_abundance, 0.5)

    # Handle strata variable if present
    if (!is.null(strata.var)){
      otu_tab_norm_agg <- otu_tab_norm_agg %>%
        dplyr::mutate(temp = !!sym(group.var)) %>%
        tidyr::separate(temp, into = c(paste0(group.var,"2"), strata.var), sep = .STRATA_SEP)
      otu_tab_norm_agg <- mStat_restore_factor_levels(
        otu_tab_norm_agg, fl$levels, paste0(group.var, "2"), strata.var)
    }

    # Filter features if specified
    if (!is.null(features.plot)){
      otu_tab_norm_agg <- otu_tab_norm_agg %>% filter(!!sym(feature.level) %in% features.plot)
    }

    # Count unique taxa levels
    taxa.levels <- otu_tab_norm_agg %>% select(all_of(c(feature.level))) %>% pull() %>% unique() %>% length()

    # Function to adjust point size range based on number of taxa levels
    adjust_size_range <- function(taxa.levels) {
      if (taxa.levels <= 2) {
        return(c(40, 57))
      } else if (taxa.levels <= 4) {
        return(c(35, 42))
      } else if (taxa.levels <= 6) {
        return(c(30, 37))
      } else if (taxa.levels <= 8) {
        return(c(25, 33))
      } else if (taxa.levels < 10) {
        return(c(20, 17))
      } else if (taxa.levels < 20) {
        return(c(10, 15))
      } else if (taxa.levels < 30) {
        return(c(8, 13))
      } else if (taxa.levels < 40) {
        return(c(6, 10))
      } else if (taxa.levels < 50) {
        return(c(4, 8))
      } else {
        return(c(1, 4))
      }
    }

    # Create the dot plot
    dotplot <-
      ggplot(
        otu_tab_norm_agg,
        aes(
          x = !!sym(feature.level),
          y = !!sym(group.var),
          size = prevalence,
          shape = !!sym(time.var),
          color = !!sym(time.var)
        )
      ) +
      geom_point(aes(group = interaction(!!sym(time.var),!!sym(feature.level)), fill = mean_abundance),
                 shape = 21,
                 position = position_dodge(0.9)) +
      xlab(feature.level) +
      ylab(group.var) +
      scale_colour_manual(values = c("transparent","black")) +
      scale_size_continuous(range = adjust_size_range(taxa.levels)) +
      {
        # Set up color scale based on feature data type
        if(feature.dat.type == "other") {
          quantiles <- quantile(otu_tab_norm_agg$mean_abundance, probs = c(0, 0.25, 0.5, 0.75, 1))
          scale_fill_gradientn(colors = colors,
                               values = scales::rescale(quantiles),
                               name = "Mean Abundance")
        } else {
          quantiles <- quantile(otu_tab_norm_agg$mean_abundance, probs = c(0, 1/3, 2/3, 1))
          labels <- sapply(quantiles, function(x) {
            as.expression(bquote(.(round(x, 2))^2))
          })

          scale_fill_gradientn(colors = colors,
                               values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                               name = "Mean Abundance",
                               breaks = quantiles,
                               labels = labels)
        }
      } +
      {
        # Set up faceting based on presence of strata variable
        if (!is.null(strata.var)){
          ggh4x::facet_nested(rows = vars(!!sym(strata.var),!!sym(paste0(group.var,"2"))), cols = vars(!!sym(feature.level)), scales = "free", switch = "y")
        } else {
          ggh4x::facet_nested(rows = vars(!!sym(group.var)), cols = vars(!!sym(feature.level)), scales = "free", switch = "y")}
      } +
      theme_to_use +
      theme(
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = base.size
        ),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        strip.text.y = if (group.var == "ALL") element_blank() else element_text(size = base.size),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        axis.text = ggplot2::element_text(color = "black",
                                          size = base.size),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16)
      ) + guides(
        color = guide_legend(override.aes = list(size = 5, fill = "#92c5de")),
        shape = guide_legend(override.aes = list(size = 5))
      )

    # Save the plot as a PDF if requested
    if (pdf) {
      dotplot <- as.ggplot(dotplot)
      pdf_name <- paste0(
        "taxa_dotplot_pair",
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
      ggsave(
        filename = pdf_name,
        plot = dotplot,
        width = pdf.wid,
        height = pdf.hei
      )
    }

    return(dotplot)
  })

  # Name the plots in the list by feature level
  names(plot_list) <- feature.level
  return(plot_list)
}