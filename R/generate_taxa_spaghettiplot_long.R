#' @title Generate Longitudinal Spaghetti Plots of Taxonomic Composition
#'
#' @description Generates spaghetti plots showing taxonomic composition trends over time,
#' with individual subject trajectories and group-level mean lines.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param features.plot A character vector specifying which feature IDs to plot.
#'   Default is NULL, in which case features are selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot. Default is NULL.
#' @param top.k.func Function to use for selecting top k features (e.g., "mean", "sd"). Default is NULL.
#' @param ... Additional parameters to be passed.
#'
#' @return A ggplot object containing the longitudinal line plot.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data("subset_T2D.obj")
#' generate_taxa_spaghettiplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "sample_body_site",
#'   feature.level = c("Genus"),
#'   features.plot = NULL,
#'   feature.dat.type = "count",
#'   top.k.plot = 10,
#'   top.k.func = "mean",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "npg",
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' }
#' @export
generate_taxa_spaghettiplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
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
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Validate the input data object
    mStat_validate_data(data.obj)

    # Check if input variables are of the correct type
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

    # Extract relevant metadata
    meta_tab <- data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Calculate sizes for various plot elements based on base.size
    title.size = base.size * 1.25
    axis.title.size = base.size * 0.75
    axis.text.size = base.size * 0.5
    legend.title.size = base.size * 1
    legend.text.size = base.size * 0.75

    # Adjust filters if necessary based on input parameters
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

    # Generate plots for each feature level
    plot_list <- lapply(feature.level, function(feature.level){

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

      # Reshape data from wide to long format
      df <- otu_tax_agg_numeric %>%
        tidyr::gather(key = "sample", value = "count", -all_of(feature.level)) %>%
        dplyr::left_join(meta_tab %>% rownames_to_column(var = "sample"), by = "sample")

      # Create a dummy group if group variable is not specified
      if (is.null(group.var)) {
        df <- df %>% dplyr::mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      # Calculate mean counts
      if (!is.null(strata.var)) {
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level),
                          !!sym(time.var),
                          !!sym(group.var),
                          !!sym(strata.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df,
                    mean_df,
                    by = c(feature.level, time.var, group.var, strata.var))
      } else {
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level), !!sym(time.var), !!sym(group.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df, mean_df, by = c(feature.level, time.var, group.var))
      }

      # Get unique taxa levels
      taxa.levels <-
        df %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # Select features to plot
      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 5) {
          features.plot <- taxa.levels[1:4]
        } else {
          features.plot <- taxa.levels
        }
      }

      # Subset data for selected features
      sub_df <- df %>% filter(!!sym(feature.level) %in% features.plot)

      # Get number of strata levels if strata variable is specified
      if (!is.null(strata.var)){
        strata.levels <- df %>% select(!!sym(strata.var)) %>% pull() %>% unique() %>% length()
      }

      # Create the spaghetti plot
      lineplot <- ggplot() +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "count",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.5
        ) +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_count",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        geom_point(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_count",
            group = group.var,
            color = group.var
          ),
          size = 3
        ) +
        scale_color_manual(values = col) +
        {
          if (feature.dat.type != "other"){
            labs(
              x = time.var,
              y = "Relative Abundance",
              color = group.var
            )
          } else {
            labs(
              x = time.var,
              y = "Abundance",
              color = group.var
            )
          }
        } +
        {
          if (!is.null(strata.var)) {
            ggh4x::facet_nested_wrap(as.formula(paste('~', feature.level,"+",strata.var)), scales = "free_y", ncol = strata.levels*2)  # Use facet_wrap with strata.var as the faceting variable
          } else {
            ggh4x::facet_nested_wrap(as.formula(paste('~',feature.level)), scales = "free_y", ncol = 3)
          }
        } +
        theme_to_use +
        theme(
          plot.title = element_text(
            size = title.size,
            hjust = 0.5
          ),
          panel.spacing.x = unit(0, "cm"),
          axis.title.x = element_text(size = axis.title.size),
          axis.title.y = element_text(size = axis.title.size),
          axis.text.x = element_text(size = axis.text.size),
          axis.text.y = element_text(size = axis.text.size),
          legend.title = element_text(size = legend.title.size),
          legend.text = element_text(size = legend.text.size)
        )

      # Remove legend for single group case
      if (group.var == "ALL") {
        lineplot <- lineplot + theme(legend.position = "none")
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_spaghettiplot_long",
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "group_",
          group.var,
          "_",
          "strata_",
          strata.var,
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

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }

        pdf_name <- paste0(pdf_name, ".pdf")

        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Print the plot to the PDF
        print(lineplot)
        # Close the PDF device
        dev.off()
      }

      return(lineplot)
    })

    # Name the plots in the list
    names(plot_list) <- feature.level

    # Return the list of plots
    return(plot_list)
  }