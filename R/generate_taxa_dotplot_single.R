#' @title Generate Taxa Dotplot for Single Time Point
#'
#' @description Generates a dotplot of taxa abundance for a single time point, showing mean
#' abundance and prevalence across groups.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param t.level Character string specifying the time level/value to subset data to.
#'   Default NULL does not subset data.
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#'
#' # Call the function
#' generate_taxa_dotplot_single(
#'   data.obj = peerj32.obj,
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 20,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 15,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_dotplot_single <- function(data.obj,
                                       time.var = NULL,
                                       t.level = NULL,
                                       group.var = NULL,
                                       strata.var = NULL,
                                       feature.level,
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
  data.obj <- mStat_validate_data(data.obj)

  # Capture original factor levels before any data extraction
  fl <- mStat_capture_factor_levels(data.obj, group.var, strata.var)
  data.obj <- fl$data.obj

  context <- mStat_prepare_taxa_single_context(
    data.obj = data.obj,
    time.var = time.var,
    t.level = t.level,
    group.var = group.var,
    strata.var = strata.var
  )
  data.obj <- context$data.obj
  meta_tab <- context$meta_tab

  placeholder_group <- mStat_ensure_group_placeholder(
    meta_tab,
    group.var = group.var,
    value = "ALL",
    column_name = "ALL"
  )
  meta_tab <- placeholder_group$df
  resolved_group_var <- placeholder_group$group.var

  # If a strata variable is provided, create an interaction term with the group variable
  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(resolved_group_var) := interaction(!!sym(resolved_group_var), !!sym(strata.var), sep = .STRATA_SEP))
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

  # Normalize count data if necessary.
  analysis_data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

  # Create a list to store plots for each taxonomic level
  plot_list <- lapply(feature.level, function(feature.level) {

    # Aggregate data by taxonomy if necessary
    otu_tax_agg <- get_taxa_data(analysis_data.obj, feature.level, prev.filter, abund.filter)

    current_features_plot <- mStat_resolve_selected_features(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      features.plot = features.plot,
      top.k.plot = top.k.plot,
      top.k.func = top.k.func
    )

    otu_tax_long <- mStat_prepare_taxa_long_data(
      feature.dat = otu_tax_agg,
      feature.level = feature.level,
      value_col = "count",
      meta.dat = meta_tab,
      join = "inner"
    )

    current_features_plot <- mStat_resolve_selected_features(
      feature.level = feature.level,
      features.plot = current_features_plot,
      taxa.levels = otu_tax_long[[feature.level]] %>% unique(),
      fallback_n = 4
    )

    otu_tab_norm_agg <- mStat_summarize_grouped_taxa_long(
      long.df = otu_tax_long,
      feature.level = feature.level,
      group_vars = resolved_group_var,
      value_col = "count",
      mean_col = "mean_abundance",
      prevalence_col = "group_prevalence",
      mean_transform = sqrt
    ) %>%
      dplyr::select(-group_prevalence)

    prevalence_all <- mStat_summarize_grouped_taxa_long(
      long.df = otu_tax_long,
      feature.level = feature.level,
      group_vars = NULL,
      value_col = "count",
      mean_col = "avg_abundance",
      prevalence_col = "prevalence"
    ) %>%
      dplyr::select(all_of(c(feature.level, "prevalence")))

    # Merge the abundance and prevalence data
    otu_tab_norm_agg <-
      otu_tab_norm_agg %>% dplyr::left_join(prevalence_all, feature.level)

    # Calculate the midpoint of mean abundance for color scaling
    midpoint <- quantile(otu_tab_norm_agg$mean_abundance, 0.5)

    # Handle strata variable if present
    if (!is.null(strata.var)){
      otu_tab_norm_agg <- otu_tab_norm_agg %>%
        dplyr::mutate(group_key = !!sym(resolved_group_var)) %>%
        tidyr::separate(group_key, into = c(paste0(group.var,"2"), strata.var), sep = .STRATA_SEP)
      otu_tab_norm_agg <- mStat_restore_factor_levels(
        otu_tab_norm_agg, fl$levels, paste0(group.var, "2"), strata.var)
    }

    # Filter features if specified
    if (!is.null(current_features_plot)){
      otu_tab_norm_agg <- otu_tab_norm_agg %>% filter(!!sym(feature.level) %in% current_features_plot)
    }

    # Create the dot plot
    # The size of the points represents prevalence, while the color represents mean abundance
    dotplot <-
      ggplot(
        otu_tab_norm_agg,
        aes(
          x = !!sym(feature.level),
          y = !!sym(resolved_group_var),
          size = prevalence
        )
      ) +
      geom_point(aes(group = !!sym(feature.level), fill = mean_abundance), shape = 21, color = "black", position = position_dodge(0.9)) +
      {
        # Set up color scale based on feature data type
        if(feature.dat.type == "other") {
          quantiles <- quantile(otu_tab_norm_agg$mean_abundance, probs = c(0, 0.25, 0.5, 0.75, 1))
          scale_fill_gradientn(colors = colors,
                               values = scales::rescale(quantiles),
                               name = "Mean Abundance")
        } else {
          scale_fill_gradientn(colors = colors,
                               values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                               name = "Mean Abundance (Sqrt)")
        }
      } +
      scale_size(range = c(4, 10), name = "Prevalence") +
      {
          # Set up faceting based on presence of strata variable
          if (!is.null(strata.var)){
            ggh4x::facet_nested(rows = vars(!!sym(strata.var), !!sym(paste0(group.var,"2"))), cols = vars(!!sym(feature.level)), scales = "free", switch = "y")
          } else {
            ggh4x::facet_nested(rows = vars(!!sym(resolved_group_var)), cols = vars(!!sym(feature.level)), scales = "free", switch = "y")}
      } +
      theme_to_use +
      theme(
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1,
          size = base.size
        ),
        strip.text.x = element_blank(),
        strip.text.y = if (is.null(group.var)) element_blank() else element_text(size = base.size),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16)
      )

    # Save the dot plot as a PDF file if requested
    if (pdf) {
      dotplot <- as.ggplot(dotplot)
      pdf_name <- paste0(
        "taxa_dotplot_single",
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
      pdf_name <- mStat_append_pdf_group_suffixes(
        pdf_name = pdf_name,
        group.var = group.var,
        strata.var = strata.var
      )
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
