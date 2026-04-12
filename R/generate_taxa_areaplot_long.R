#' Generate Taxa Area Plots for Longitudinal Data
#'
#' Creates stacked area plots showing relative abundance of taxa over time.
#' Supports grouping and stratification for comparative visualization.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param features.plot Character vector of specific feature IDs to plot.
#'   If NULL, top features by mean abundance are displayed.
#' @param feature.number Integer specifying number of top features to display.
#'   Lower-ranked features are grouped into "Other". Default is 20.
#'
#' @return A list of ggplot objects, each representing a taxa area plot for the specified feature level.
#'
#' @examples
#' \dontrun{
#' library(ggh4x)
#' library(vegan)
#' data(ecam.obj)
#' generate_taxa_areaplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "proportion",
#'   feature.number = 40,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' generate_taxa_areaplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "proportion",
#'   feature.number = 20,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' generate_taxa_areaplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "proportion",
#'   feature.number = 20,
#'   features.plot = unique(ecam.obj$feature.ann[,"Genus"])[1:15],
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#' data(subset_T2D.obj)
#' generate_taxa_areaplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_gender",
#'   strata.var = "subject_race",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.number = 40,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   pdf.wid = 49,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_areaplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_id",
#'   strata.var = "subject_gender",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.number = 40,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   pdf.wid = 49,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_areaplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "sample_body_site",
#'   strata.var = "subject_race",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.number = 40,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   pdf.wid = 49,
#'   file.ann = NULL
#' )
#' }
#' @import rlang
#' @import tibble
#' @export
generate_taxa_areaplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           feature.level,
           feature.dat.type = c("count", "proportion", "other"),
           feature.number = 20,
           features.plot = NULL,
           t0.level = NULL,
           ts.levels = NULL,
           base.size = 10,
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
    data.obj <- mStat_validate_data(data.obj)

    # Check if the input variables are of the correct type
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

    # Process the time variable in the data object
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

    # Capture original factor levels before any data extraction
    fl <- mStat_capture_factor_levels(data.obj, group.var, strata.var)
    data.obj <- fl$data.obj

    # Extract relevant metadata
    meta_tab <- data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Define a color palette if not provided
    if (is.null(palette)){
      pal <- rep(c("#E0E0E0","#E41A1C","#1E90FF","#FF8C00","#4DAF4A","#984EA3","#40E0D0","#FFC0CB",
                   "#00BFFF","#FFDEAD","#90EE90","#EE82EE","#00FFFF","#F0A3FF", "#0075DC",
                   "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00",
                   "#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010",
                   "#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",
                   "#FF6A6A", "#FF8247", "#FFE7BA", "#87CEFA", "#B0E0E6", "#48D1CC", "#5F9EA0", "#66CDAA", "#458B00",
                   "#BCEE68", "#FFF68F", "#EEEE00", "#FFFFE0", "#8B8682", "#FFB6C1", "#9370DB", "#FFDAB9", "#FA8072",
                   "#90EE90", "#00FFFF", "#00BFFF", "#FFA07A", "#F08080", "#FFD700", "#ADFF2F", "#9ACD32", "#9400D3",
                   "#7B68EE", "#BA55D3", "#FFC0CB", "#FAEBD7", "#F0E68C", "#FFFACD", "#D2B48C", "#C0C0C0", "#696969",
                   "#CD5C5C", "#F08080", "#FA8072", "#E9967A", "#FFA07A", "#FF7F50", "#FF6347", "#FF4500", "#FF8C00",
                   "#FFA500", "#FFDAB9", "#FFE4B5", "#FFE4C4", "#FFDEAD", "#EEE8AA", "#F0E68C", "#BDB76B", "#FFD700",
                   "#DAA520", "#808000", "#7CFC00", "#00FF00", "#32CD32", "#00FA9A", "#90EE90", "#98FB98", "#8FBC8F",
                   "#3CB371", "#2E8B57", "#228B22", "#008000", "#006400", "#66CDAA", "#00FFFF", "#40E0D0", "#48D1CC",
                   "#AFEEEE", "#7FFFD4", "#B0C4DE", "#ADD8E6", "#87CEEB", "#87CEFA", "#6495ED", "#00BFFF", "#1E90FF",
                   "#4169E1", "#0000FF", "#00008B", "#000080", "#191970", "#6A5ACD", "#483D8B", "#9400D3", "#8A2BE2",
                   "#4B0082", "#FF00FF", "#FF69B4", "#FF1493", "#C71585", "#DB7093", "#FFC0CB", "#FFB6C1", "#FF69B4",
                   "#FF5F5F", "#DC143C", "#C0C0C0", "#A9A9A9", "#808080", "#696969", "#000000", "#FF1493", "#FF69B4",
                   "#FFB6C1", "#FFC0CB", "#FFDAB9", "#F4A460", "#FFA07A", "#FF7F50", "#FF6347", "#FF4500", "#FF8C00",
                   "#FFA500", "#FFFF00", "#9ACD32", "#32CD32", "#00FF00", "#7FFF00", "#7CFC00", "#00FA9A", "#90EE90",
                   "#98FB98", "#8FBC8F", "#3CB371", "#2E8B57", "#228B22", "#008000", "#006400"
      ),5)
    } else{
      # Check if the provided palette is sufficient
      if (feature.number > length(palette)) {
        stop("The number of unique features exceeds the length of the provided palette. Please provide a larger palette.")
      }

      pal <- palette
    }

    # Normalize the data if it's in count format.
    data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)
    if (feature.dat.type == "other"){
      stop("The 'other' type is suitable for situations where the user has analyzed the data using a method not provided in 'mStat_normalize_data' method, and the 'areaplot' is only applicable to raw data that has not undergone any processing or proportion data that adds up to 1. If you believe your data falls into these two categories, please modify 'feature.dat.type'.")
    }

    # Generate plots for each feature level
    plot_list_all <- lapply(feature.level,function(feature.level){

      otu_tax_agg <- get_taxa_data(data.obj, feature.level)
      selected_features <- mStat_resolve_selected_features(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        features.plot = features.plot
      )
      otu_tax_agg <- mStat_filter_taxa_features(
        feature.dat = otu_tax_agg,
        feature.level = feature.level,
        features.plot = selected_features
      )

      otu_tab_norm <- mStat_as_taxa_composition_matrix(
        feature.dat = otu_tax_agg,
        feature.level = feature.level
      )

      meta_tab_sorted <- meta_tab[colnames(otu_tab_norm), ]

      stack_levels <- mStat_prepare_stack_levels(
        composition.mat = otu_tab_norm,
        feature.level = feature.level,
        feature.number = feature.number,
        other_first = TRUE,
        other_inclusive = FALSE
      )
      otu_tab_long <- stack_levels$long.df
      new_levels <- stack_levels$new_levels
      other.abund.cutoff <- stack_levels$other.abund.cutoff

      # Merge feature data with metadata
      merged_long_df <- otu_tab_long %>%
        dplyr::inner_join(meta_tab_sorted  %>%
                            tibble::rownames_to_column("sample"), by = "sample")

      # Sort the merged data
      sorted_merged_long_df <- merged_long_df %>%
        dplyr::arrange(!!sym(subject.var), !!sym(time.var))

      last_sample_ids <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::summarize(last_sample_id = dplyr::last(sample))

      sorted_merged_long_df <- sorted_merged_long_df %>%
        dplyr::mutate(!!sym(feature.level) := as.factor(!!sym(feature.level)))

      # Calculate the average value of each feature and sort them.
      df_sorted <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarise(overall_mean = mean(value, na.rm = TRUE)) %>%
        dplyr::mutate(is_other = ifelse(!!sym(feature.level) == "Other", TRUE, FALSE)) %>%
        dplyr::arrange(is_other, overall_mean) %>%
        dplyr::mutate(!!feature.level := factor(!!sym(feature.level), levels = !!sym(feature.level)))

      # Update new_levels
      if (!is.na(other.abund.cutoff)) {
        new_levels <- c("Other", setdiff(levels(df_sorted[[feature.level]]), "Other"))
      } else {
        new_levels <- levels(df_sorted[[feature.level]])
      }

      # Apply new sorting
      sorted_merged_long_df <- sorted_merged_long_df %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels))

      # Modify the creation of df
      df <- mStat_prepare_stacked_positions(
        long.df = sorted_merged_long_df,
        feature.level = feature.level,
        id_var = "sample",
        ordered_levels = new_levels,
        terminal_ids = last_sample_ids$last_sample_id
      )

      if (!is.null(selected_features)) {
        df <- df %>% dplyr::filter(!!sym(feature.level) %in% selected_features)
      }

      # Update color palette
      color_pal <- setNames(pal[1:length(new_levels)], new_levels)

      bar_width <- 0.6
      bar_spacing <- bar_width / 2

      # The average barplot is plotted as follows
      if (is.factor(meta_tab[, time.var])) {
        last_time_ids <- dplyr::last(levels(meta_tab[, time.var]))
      } else if (is.numeric(meta_tab[, time.var])) {
        last_time_ids <- max(meta_tab[, time.var], na.rm = TRUE)
      } else {
        stop("The variable is neither factor nor numeric.")
      }

      has_group <- !is.null(group.var)
      placeholder_group <- mStat_ensure_group_placeholder(
        sorted_merged_long_df,
        group.var = group.var,
        value = "ALL",
        column_name = "ALL"
      )
      sorted_merged_long_df <- placeholder_group$df
      resolved_group_var <- placeholder_group$group.var
      if (!is.null(strata.var)) {
        sorted_merged_long_df <- sorted_merged_long_df %>%
          dplyr::mutate(!!sym(resolved_group_var) := interaction(!!sym(resolved_group_var), !!sym(strata.var), sep = .STRATA_SEP))
      }

      average_stack <- mStat_prepare_average_stack_data(
        long.df = sorted_merged_long_df,
        feature.level = feature.level,
        group.var = resolved_group_var,
        time.var = time.var,
        ordered_levels = new_levels,
        terminal_time_values = last_time_ids,
        bar_width = bar_width,
        bar_spacing = bar_spacing
      )
      df_average <- average_stack$df
      labels <- average_stack$labels

      # Separate group and strata variables and restore factor levels
      if(!is.null(strata.var)){
        group_var_for_restore <- if (has_group) group.var else resolved_group_var
        df_average <- df_average %>%
          tidyr::separate(!!sym(resolved_group_var), into = c(group_var_for_restore, strata.var), sep = .STRATA_SEP)
        df_average <- mStat_restore_factor_levels(df_average, fl$levels, group_var_for_restore, strata.var)
      }

      # Create the main stacked area plot
      stack_areaplot_average  <- 
        df_average %>%
        ggplot(aes(x = joint_factor_numeric, y = mean_value, fill = !!sym(feature.level))) +
        geom_area(stat = "identity", position = "fill") +
        scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
        scale_x_continuous(expand = c(0.01, 0.01), breaks = unique(df_average$joint_factor_numeric), labels = labels) +
        {
          # Add faceting if group or strata variables are present
          if (has_group && !is.null(strata.var)) {
            ggh4x::facet_nested(as.formula(paste(". ~", strata.var, "+", group.var)), drop = T, scales = "free", space = "free")
          } else if (has_group) {
            ggh4x::facet_nested(as.formula(paste(". ~", group.var)), drop = T, scales = "free", space = "free")
          } else if (!is.null(strata.var)) {
            ggh4x::facet_nested(as.formula(paste(". ~", strata.var)), drop = T, scales = "free", space = "free")
          }
        } +
        labs(fill = feature.level, y = "", x = "") +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        theme_to_use +
        theme(strip.background = element_rect(fill="white",color="black"),
              panel.spacing = unit(0,"lines"),
              strip.text.x = element_text(size= base.size,color="black"),
              axis.text.y=element_text(size= base.size,color="black"),
              axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size),
              axis.title.y = element_text(size= base.size,color="black"),
              legend.key=element_blank(),
              legend.title = element_text(size= base.size+2),
              legend.text = element_text(color="black",size= base.size),
              legend.spacing.x=unit(0.1,'cm'),
              legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.4,'cm'),
              legend.key.height=unit(0.4,'cm'),
              legend.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0("taxa_areaplot_pair",
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "feature_level_", feature.level,
                           "_",
                           "feature_number_", feature.number)
        pdf_name <- mStat_append_pdf_group_suffixes(
          pdf_name = pdf_name,
          group.var = if (has_group) group.var else NULL,
          strata.var = strata.var
        )
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }
        pdf_name <- paste0(pdf_name,"_average", ".pdf")
        ggsave(filename = pdf_name, plot = stack_areaplot_average, width = pdf.wid, height = pdf.hei)
      }

      return(stack_areaplot_average)
    })

    # Name the plots in the list
    names(plot_list_all) <- feature.level

    return(plot_list_all)
  }
