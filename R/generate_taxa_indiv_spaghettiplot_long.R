#' @title Generate Individual Taxa Spaghetti Plots for Longitudinal Data
#'
#' @description Creates longitudinal line plots (spaghetti plots) for individual taxa showing
#' abundance trajectories over time with group means overlaid.
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
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_taxa_indiv_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = c("Phylum"),
#'   features.plot = NULL,
#'   feature.dat.type = "proportion",
#'   top.k.plot = 5,
#'   top.k.func = "sd",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' }
#' @export
generate_taxa_indiv_spaghettiplot_long <-
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

    # Process time variable and extract relevant data
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

    # Extract metadata
    meta_tab <- data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

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

      # Reshape data for plotting
      df <- otu_tax_agg_numeric %>%
        tidyr::gather(key = "sample", value = "count", -all_of(feature.level)) %>%
        dplyr::left_join(meta_tab %>% rownames_to_column(var = "sample"), by = "sample")

      # Create a dummy group if group variable is not provided
      if (is.null(group.var)) {
        df <- df %>% dplyr::mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      # Calculate mean counts
      if (!is.null(strata.var)) {
        # Calculate mean counts for each combination of feature, time, group, and strata
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level),!!sym(time.var),!!sym(group.var),!!sym(strata.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df,
                           mean_df,
                           by = c(feature.level, time.var, group.var, strata.var))
      } else {
        # Calculate mean counts for each combination of feature, time, and group
        mean_df <-
          df %>% dplyr::group_by(!!sym(feature.level),
                                 !!sym(time.var),
                                 !!sym(group.var)) %>%
          dplyr::summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          dplyr::left_join(df, mean_df, by = c(feature.level, time.var, group.var))
      }

      # Get color palette
      col <- mStat_get_palette(palette)

      # Get the appropriate theme for plotting
      theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

      # Calculate new sizes based on base.size
      title.size = base.size * 1.25
      axis.title.size = base.size * 0.75
      axis.text.size = base.size * 0.5
      legend.title.size = base.size * 1
      legend.text.size = base.size * 0.75

      # Get taxa levels to plot
      if (is.null(features.plot)){
        taxa.levels <- df %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()
      } else {
        taxa.levels <- df %>% filter(!!sym(feature.level) %in% features.plot) %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()
      }

      # Create a plot for each taxon
      plot_list <- lapply(taxa.levels, function(tax) {
        sub_df <- df %>% filter(!!sym(feature.level) == tax)
        lineplot <- ggplot() +
          # Add individual subject points
          geom_point(
            data = sub_df,
            aes_string(
              x = time.var,
              y = "count",
              group = subject.var,
              color = group.var
            ),
            alpha = 0.5
          ) +
          # Add mean abundance line
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
          # Add mean abundance points
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
          # Set appropriate labels based on feature data type
          {
            if (feature.dat.type != "other") {
              labs(
                x = time.var,
                y = "Relative Abundance",
                color = group.var,
                title = tax
              )
            } else {
              labs(
                x = time.var,
                y = "Abundance",
                color = group.var,
                title = tax
              )
            }
          } +
          # Add faceting if strata variable is provided
          {
            if (!is.null(strata.var)) {
              facet_wrap(as.formula(paste('~', strata.var)))
            }
          } +
          theme_to_use +
          theme(
            plot.title = element_text(
              size = title.size,
              hjust = 0.5
            ),
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
        return(lineplot)
      })

      # Save the plots as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_spaghettiplot_long",
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
