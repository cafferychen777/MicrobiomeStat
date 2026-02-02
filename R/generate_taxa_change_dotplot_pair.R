#' @title Generate Taxa Change Dotplot for Paired Data
#'
#' @description Creates stacked dotplots showing abundance and prevalence changes between
#' paired time points. Point size represents baseline values, color represents change magnitude.
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
#' @return A list of ggplot dotplot objects.
#' @examples
#' \dontrun{
#'
#' # Note: In the RStudio viewer, the plot might appear cluttered if there are many taxa.
#' # It's recommended to view the generated PDF for better clarity. If it still feels
#' # overcrowded in the PDF, consider increasing the 'pdf.wid' value to adjust the width of the plot.
#'
#' data(peerj32.obj)
#' generate_taxa_change_dotplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = "Family",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 20,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 1e-4,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 10
#' )
#'
#' data("subset_pairs.obj")
#' generate_taxa_change_dotplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = "Family",
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = 20,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 1e-4,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 30,
#'   pdf.hei = 10
#' )
#' }
#' # View the result
#' @export
generate_taxa_change_dotplot_pair <- function(data.obj,
                                              subject.var,
                                              time.var,
                                              group.var = NULL,
                                              strata.var = NULL,
                                              change.base = "1",
                                              feature.change.func = "relative change",
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
                                              palette = NULL,
                                              pdf = TRUE,
                                              file.ann = NULL,
                                              pdf.wid = 11,
                                              pdf.hei = 6,
                                              ...) {

  feature.dat.type <- match.arg(feature.dat.type)

  # Extract data
  mStat_validate_data(data.obj)

  # Capture original factor levels before any data extraction
  fl <- mStat_capture_factor_levels(data.obj, group.var, strata.var)
  data.obj <- fl$data.obj

  meta_tab <-
    data.obj$meta.dat %>% select(all_of(c(
      time.var, group.var, strata.var, subject.var
    )))

  if (is.null(group.var)) {
    group.var = "ALL"
    meta_tab$ALL <- ""
  }

  if (!is.null(strata.var)) {
    meta_tab <-
      meta_tab %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var), !!sym(strata.var), sep = .STRATA_SEP))
  }

  # Define the colors
  if (is.null(palette)) {
    colors <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
  } else {
    colors <- palette
  }

  # Assuming mStat_get_theme function is already defined
  # Replace the existing theme selection code with this:
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  if (feature.dat.type == "other" || !is.null(features.plot) ||
      (!is.null(top.k.func) && !is.null(top.k.plot))) {
    prev.filter <- 0
    abund.filter <- 0
  }

  plot_list <- lapply(feature.level, function(feature.level) {
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <-
        mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    otu_tax_agg <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

    if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
    }

    # Calculate the average abundance for each group
    otu_tab_norm_agg <- otu_tax_agg %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var), !!sym(feature.level), !!sym(time.var), !!sym(subject.var)) %>% # Add time.var to dplyr::group_by
      dplyr::summarise(mean_abundance = mean(count))

    change.after <-
      unique(meta_tab %>% select(all_of(c(time.var))))[unique(meta_tab %>% select(all_of(c(time.var)))) != change.base]

    # Convert data from long format to wide format, placing the mean_abundance values at different time points in different columns.
    otu_tab_norm_agg_wide <- otu_tab_norm_agg %>%
      tidyr::spread(key = !!sym(time.var), value = mean_abundance) %>%
      dplyr::rename(
        time1_mean_abundance = all_of(change.base),
        time2_mean_abundance = all_of(change.after)
      )

    # Calculate abundance change (fix: was passing time2 twice for custom func,
    # and using fixed +0.00001 instead of per-taxon pseudocount for log fold change)
    otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
      dplyr::mutate(abundance_change = compute_taxa_change(
        value_after  = time2_mean_abundance,
        value_before = time1_mean_abundance,
        method       = feature.change.func,
        feature_id   = .data[[feature.level]]
      ))

    # Compute the prevalence of each taxon at each time point
    prevalence_time <- otu_tax_agg %>%
      tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample") %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level),!!sym(time.var)) %>%
      dplyr::summarise(prevalence = sum(count > 0) / dplyr::n())

    prevalence_time_wide <- prevalence_time %>%
      tidyr::spread(key = time.var, value = prevalence) %>%
      dplyr::rename(time1_prevalence = change.base,
                    time2_prevalence = change.after)

    # Compute the difference in prevalence at different time points
    prevalence_time_wide <- prevalence_time_wide %>%
      dplyr::mutate(prevalence_change = compute_taxa_change(
        value_after  = time2_prevalence,
        value_before = time1_prevalence,
        method       = feature.change.func,
        feature_id   = .data[[feature.level]],
        verbose      = FALSE
      ))

    otu_tab_norm_agg_wide <-
      otu_tab_norm_agg_wide %>% dplyr::left_join(prevalence_time_wide, by = c(feature.level, group.var))

    if (!is.null(strata.var)) {
      otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
        dplyr::mutate(temp = !!sym(group.var)) %>%
        tidyr::separate(temp,
                        into = c(paste0(group.var, "2"), strata.var),
                        sep = .STRATA_SEP)
      otu_tab_norm_agg_wide <- mStat_restore_factor_levels(
        otu_tab_norm_agg_wide, fl$levels, paste0(group.var, "2"), strata.var)
    }

    if (!is.null(features.plot)) {
      otu_tab_norm_agg_wide <-
        otu_tab_norm_agg_wide %>% filter(!!sym(feature.level) %in% features.plot)
    }

    prop_prev_data <-
      otu_tax_agg %>%
      column_to_rownames(feature.level) %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::group_by(Var1) %>%
      dplyr::summarise(avg_abundance = mean(Freq),
                       prevalence = sum(Freq > 0) / dplyr::n()) %>% column_to_rownames("Var1") %>%
      rownames_to_column(feature.level)

    otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide %>%
      tidyr::gather(key = "Type",
                    value = "change",
                    abundance_change,
                    prevalence_change) %>%
      select(-all_of(
        c(
          "time1_mean_abundance",
          "time2_mean_abundance",
          "time1_prevalence",
          "time2_prevalence"
        )
      )) %>%
      dplyr::left_join(prop_prev_data, by = feature.level) %>%
      dplyr::mutate(
        Base = dplyr::case_when(
          Type == "abundance_change" ~ avg_abundance,
          Type == "prevalence_change" ~ prevalence,
          TRUE ~ NA_real_
        )
      ) %>%
      select(-all_of(c("avg_abundance", "prevalence")))

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

    # Find the minimum and maximum values of abundance_change
    change_min <- min(otu_tab_norm_agg_wide$change)
    change_max <- max(otu_tab_norm_agg_wide$change)

    # Normalized function
    normalize <- function(x, min, max) {
      return((x - min) / (max - min))
    }

    # Normalized value
    change_min_norm <-
      normalize(change_min, change_min, change_max)
    change_max_norm <-
      normalize(change_max, change_min, change_max)
    change_mid_norm <-
      normalize(0, change_min, change_max)  # The midpoint of abundance_change is 0

    # Calculate the normalized value of other colors
    first_color_norm <-
      change_min_norm + (change_mid_norm - change_min_norm) / 2
    second_color_norm <-
      change_mid_norm + (change_max_norm - change_mid_norm) / 2

    taxa.levels <-
      otu_tab_norm_agg_wide %>% dplyr::ungroup() %>% select(all_of(c(feature.level))) %>% pull() %>% unique() %>% length()

    otu_tab_norm_agg_wide <- otu_tab_norm_agg_wide  %>%
      select(-all_of(subject.var)) %>%
      dplyr::group_by(!!sym(group.var),!!sym(feature.level)) %>%
      dplyr::mutate(
        change = mean(change)
      ) %>%
      dplyr::distinct()

    # Add disease prevalence as the size of the points, and use the average abundance as the color of the points
    dotplot <-
      ggplot(
        otu_tab_norm_agg_wide,
        aes(
          x = !!sym(feature.level),
          y = !!sym(group.var),
          size = Base,
          shape = Type,
          color = Type
        )
      ) + # Change x to time.var
      geom_point(
        aes(group = interaction(Type, !!sym(feature.level)), fill = change),
        shape = 21,
        position = position_dodge(0.9)
      ) +
      xlab(feature.level) +
      ylab(group.var) +
      scale_colour_manual(values = c("transparent", "black")) +
      scale_size_continuous(range = adjust_size_range(taxa.levels)) +
      scale_fill_gradientn(
        colors = colors,
        values = c(
          change_min_norm,
          first_color_norm,
          change_mid_norm,
          second_color_norm,
          change_max_norm
        ),
        name = "Change"
      ) +
      {
        if (!is.null(strata.var)) {
          ggh4x::facet_nested(
            rows = vars(!!sym(strata.var),!!sym(paste0(
              group.var, "2"
            ))),
            cols = vars(!!sym(feature.level)),
            scales = "free",
            switch = "y"
          )
        } else {
          ggh4x::facet_nested(
            rows = vars(!!sym(group.var)),
            cols = vars(!!sym(feature.level)),
            scales = "free",
            switch = "y"
          )
        }
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
        axis.title.x = element_text(size = base.size),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        strip.text.x = element_blank(),
        strip.text.y = if (group.var == "ALL")
          element_blank()
        else
          element_text(size = base.size),
        panel.spacing = unit(0, "lines"),
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
        legend.text = ggplot2::element_text(size = base.size),
        legend.title = ggplot2::element_text(size = base.size),
      ) + guides(color = guide_legend(override.aes = list(
        size = 5, fill = "#92c5de"
      )),
      shape = guide_legend(override.aes = list(size = 5)))

    # Save the stacked dotplot as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_change_dotplot_pair",
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

  names(plot_list) <- feature.level
  return(plot_list)
}
