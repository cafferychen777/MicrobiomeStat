#' Generate Taxa Barplots for Longitudinal Data
#'
#' Creates stacked barplots showing taxonomic composition of samples over time.
#' Supports grouping and stratification for comparative visualization.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param features.plot Character vector of specific feature IDs to plot.
#'   If NULL, top features by mean abundance are displayed.
#' @param feature.number Integer specifying number of top features to display.
#'   Lower-ranked features are grouped into "Other". Default is 20.
#'
#' @return A list of ggplot objects for each feature level.
#'
#' @examples
#' \dontrun{
#' library(ggh4x)
#' library(vegan)
#' data(ecam.obj)
#' generate_taxa_barplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Family",
#'   feature.dat.type = "proportion",
#'   feature.number = 30,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_barplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Genus",
#'   feature.dat.type = "proportion",
#'   feature.number = 30,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' generate_taxa_barplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "delivery",
#'   strata.var = "diet",
#'   feature.level = "Genus",
#'   feature.dat.type = "proportion",
#'   features.plot = unique(ecam.obj$feature.ann[,"Genus"])[1:10],
#'   feature.number = 30,
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL
#' )
#'
#' data(subset_T2D.obj)
#' generate_taxa_barplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
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
#' generate_taxa_barplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "sample_body_site",
#'   strata.var = "subject_race",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   feature.number = 20,
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
#' generate_taxa_barplot_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   features.plot = unique(subset_T2D.obj$feature.ann[,"Genus"][1:6]),
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
#' @export
generate_taxa_barplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           feature.level = "original",
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

    # Validate input data and variables
    mStat_validate_data(data.obj)

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

    # Capture original factor levels for group.var and strata.var
    # so we can preserve ordering throughout the plotting pipeline
    if (!is.null(group.var)) {
      if (is.factor(data.obj$meta.dat[[group.var]])) {
        group_levels_original <- levels(data.obj$meta.dat[[group.var]])
      } else {
        group_levels_original <- unique(data.obj$meta.dat[[group.var]])
      }
      data.obj$meta.dat[[group.var]] <- factor(
        data.obj$meta.dat[[group.var]], levels = group_levels_original
      )
    }
    if (!is.null(strata.var)) {
      if (is.factor(data.obj$meta.dat[[strata.var]])) {
        strata_levels_original <- levels(data.obj$meta.dat[[strata.var]])
      } else {
        strata_levels_original <- unique(data.obj$meta.dat[[strata.var]])
      }
      data.obj$meta.dat[[strata.var]] <- factor(
        data.obj$meta.dat[[strata.var]], levels = strata_levels_original
      )
    }

    # Extract relevant metadata
    meta_tab <- data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

    # Get the appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Define color palette for plotting
    if (is.null(palette)){
      # If no palette is provided, use a default color scheme
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
      # Check if the provided palette is sufficient for the number of features
      if (feature.number > length(palette)) {
        stop("The number of unique features exceeds the length of the provided palette. Please provide a larger palette.")
      }

      pal = palette
    }

    # Note: For barplot with position="fill", normalization is not required
    # as ggplot2 automatically converts to proportions. We skip normalization
    # for all data types to preserve pre-computed feature.agg.list and improve
    # performance.
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw count format. For barplot visualization, ",
        "normalization is not required as ggplot2's position='fill' will ",
        "automatically convert to proportions."
      )
      # Skip normalization to preserve feature.agg.list
    } else if (feature.dat.type == "proportion") {
      message(
        "Your data is in proportion format. For barplot visualization, ",
        "additional normalization is not required as ggplot2's position='fill' ",
        "will automatically handle the proportions."
      )
      # Skip normalization to preserve feature.agg.list
    } else if (feature.dat.type == "other") {
      message("The 'other' type is suitable for situations where the user has ",
              "analyzed the data using a method not provided in ",
              "'mStat_normalize_data' method, and the 'barplot' is only ",
              "applicable to raw data that has not undergone any processing ",
              "or proportion data that adds up to 1. If you believe your data ",
              "falls into these two categories, please modify 'feature.dat.type'.")
    }

    # Generate plots for each feature level
    plot_list_all <- lapply(feature.level,function(feature.level){

      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Select the appropriate feature table
      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Subset features if specified (using %in% for robustness against NA or non-existent features)
      if (!is.null(features.plot)){
        otu_tax_agg <- otu_tax_agg[rownames(otu_tax_agg) %in% features.plot,]
      }

      # Prepare the feature table for analysis
      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        rownames_to_column(feature.level)

      # Transpose the data
      otu_tab_counts <- apply(t(otu_tax_agg %>% select(-all_of(feature.level))), 1, function(x) x)
      # Actually normalize to proportions for each sample (account for sequencing depth)
      otu_tab_norm <- sweep(otu_tab_counts, 2, colSums(otu_tab_counts), "/")
      rownames(otu_tab_norm) <- as.matrix(otu_tax_agg[, feature.level])

      meta_tab_sorted <- meta_tab[colnames(otu_tab_norm), ]

      # Calculate the average relative abundance of each taxon
      avg_abund <- rowMeans(otu_tab_norm)

      # Replace taxa with relative abundance lower than the threshold with "Other"
      otu_tab_other <- otu_tab_norm %>%
        as.data.frame() %>%
        rownames_to_column(feature.level)

      # Set the relative abundance threshold after feature.number
      other.abund.cutoff <- sort(avg_abund, decreasing=TRUE)[feature.number]

      if (!is.na(other.abund.cutoff)){
        otu_tab_other[, feature.level][avg_abund < other.abund.cutoff] <- "Other"
      }

      # Convert data frame to long format for plotting
      otu_tab_long <- otu_tab_other %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarize_all(sum) %>%
        tidyr::gather(key = "sample", value = "value", -feature.level)

      # Merge feature data with metadata
      merged_long_df <- otu_tab_long %>%
        dplyr::inner_join(meta_tab_sorted  %>% rownames_to_column("sample"), by = "sample")

      # Sort the merged data
      sorted_merged_long_df <- merged_long_df %>%
        dplyr::arrange(!!sym(subject.var), !!sym(time.var))

      # Identify the last sample for each subject
      last_sample_ids <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::summarize(last_sample_id = dplyr::last(sample))

      # Convert feature level to factor
      sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(feature.level) := as.factor(!!sym(feature.level)))

      # Calculate the average value of each feature and sort them
      df_sorted <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarise(overall_mean = mean(value, na.rm = TRUE)) %>%
        dplyr::mutate(is_other = ifelse(!!sym(feature.level) == "Other", TRUE, FALSE)) %>%
        dplyr::arrange(is_other, overall_mean) %>%
        dplyr::mutate(!!feature.level := factor(!!sym(feature.level), levels = !!sym(feature.level)))

      # Update feature levels
      if (!is.na(other.abund.cutoff)) {
        new_levels <- c("Other", setdiff(levels(df_sorted[[feature.level]]), "Other"))
      } else {
        new_levels <- levels(df_sorted[[feature.level]])
      }

      # Apply new sorting to the data
      sorted_merged_long_df <- sorted_merged_long_df %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels))

      # Create the main dataframe for plotting
      df <- sorted_merged_long_df %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels)) %>%
        dplyr::arrange(match(!!sym(feature.level), new_levels)) %>%
        dplyr::mutate(cumulative_value = (1-cumsum(value))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::mutate(next_cumulative_value = dplyr::if_else(sample %in% last_sample_ids$last_sample_id, NA_real_, dplyr::lead(cumulative_value))) %>%
        dplyr::ungroup()

      # Update color palette
      color_pal <- setNames(pal[1:length(new_levels)], new_levels)

      # Set bar width and spacing for the plot
      bar_width <- 0.6
      bar_spacing <- bar_width / 2

      # Identify the last time point for plotting
      if (is.factor(meta_tab[, time.var])) {
        last_time_ids <- dplyr::last(levels(meta_tab[, time.var]))
      } else if (is.numeric(meta_tab[, time.var])) {
        last_time_ids <- max(meta_tab[, time.var], na.rm = TRUE)
      } else {
        stop("The variable is neither factor nor numeric.")
      }

      # Handle grouping and stratification
      if (!is.null(strata.var)){
        if (!is.null(group.var)){
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var),!!sym(strata.var)))
        } else {
          group.var = ""
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := "")
        }
      } else {
        if (!is.null(group.var)){
        } else {
          group.var = ""
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := "")
        }
      }

      # Calculate average values for plotting
      df_average <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(feature.level),!!sym(group.var),!!sym(time.var)) %>%
        dplyr::summarise(mean_value  = mean(value)) %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = new_levels)) %>%
        dplyr::arrange(match(!!sym(feature.level), new_levels),!!sym(group.var),!!sym(time.var)) %>%
        dplyr::group_by(!!sym(group.var),!!sym(time.var)) %>%
        dplyr::mutate(cumulative_mean_value = (1-cumsum(mean_value))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::mutate(next_cumulative_mean_value = dplyr::if_else(!!sym(time.var) %in% last_time_ids, NA_real_, dplyr::lead(cumulative_mean_value))) %>%
        dplyr::ungroup()

      # Handle grouping for plotting
      if (group.var == ""){
        df_average <- df_average %>% dplyr::mutate(!!sym(group.var) := "")
        df_average <- df_average %>%
          dplyr::mutate(joint_factor = interaction(!!sym(time.var), ""))
      } else {
        df_average <- df_average %>%
          dplyr::mutate(joint_factor = interaction(!!sym(time.var), !!sym(group.var)))
      }

      # Calculate x-axis offsets for labels
      df_average <- df_average %>%
        dplyr::mutate(x_offset = ifelse(cumulative_mean_value == 0, (bar_width + bar_spacing) / 2, -(bar_width + bar_spacing) / 2))

      # Drop unused levels from joint_factor
      df_average$joint_factor <- droplevels(df_average$joint_factor)

      # Convert joint_factor to numeric for plotting
      df_average$joint_factor_numeric <- match(df_average$joint_factor, levels(df_average$joint_factor))

      # Extract labels for x-axis
      labels <- sub("\\..*", "", levels(df_average$joint_factor))

      # Separate group and strata variables if necessary
      if(!is.null(strata.var)){
        df_average <- df_average %>%
          tidyr::separate(!!sym(group.var), into = c(group.var, strata.var), sep = "\\.")
        # Restore factor levels after separate() which creates character columns
        df_average[[group.var]] <- factor(df_average[[group.var]], levels = group_levels_original)
        df_average[[strata.var]] <- factor(df_average[[strata.var]], levels = strata_levels_original)
      }

      # Create the main stacked barplot
      stack_barplot_average  <- 
        df_average %>%
        ggplot(aes(x = joint_factor_numeric, y = mean_value, fill = !!sym(feature.level))) +
        geom_bar(stat = "identity", position = "fill", width = bar_width) +
        geom_segment(aes(x = joint_factor_numeric + bar_spacing, xend = joint_factor_numeric + 1 - bar_spacing, y = cumulative_mean_value, yend = next_cumulative_mean_value, group = !!sym(feature.level), color = !!sym(feature.level)),linewidth = 0.6) +
        scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
        scale_x_continuous(expand = c(0.01, 0.01), breaks = unique(df_average$joint_factor_numeric), labels = labels) +
        {
          # Add faceting if group variable is present
          if (!is.null(group.var)){
            if (group.var == ""){
            } else {
              if (!is.null(strata.var)){
                ggh4x::facet_nested(as.formula(paste(". ~", group.var, "+", strata.var)), drop = T, scale = "free", space = "free")
              } else {
                ggh4x::facet_nested(as.formula(paste(". ~", group.var)), drop = T, scale = "free", space = "free")
              }
            }
          }
        } +
        labs(fill = feature.level, y = "", x = "") +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        theme_to_use +
        theme(strip.background = element_rect(fill="white"),
              panel.spacing = unit(0,"lines"),
              strip.text.x = element_text(size= base.size,color="black"),
              axis.text.y=element_text(size= base.size,color="black"),
              axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size),
              axis.title.y = element_text(size= base.size,color="black"),
              legend.key=element_blank(),
              legend.title = element_text(size= base.size), 
              legend.text = element_text(color="black",size= base.size),
              legend.spacing.x=unit(0.1,'cm'),
              legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.4,'cm'),
              legend.key.height=unit(0.4,'cm'),
              legend.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())

      # Save the stacked barplot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0("taxa_barplot_long",
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "feature_level_", feature.level,
                           "_",
                           "feature_number_", feature.number)
        if (!is.null(group.var)) {
          pdf_name <- paste0(pdf_name, "_", "group_", group.var)
        }
        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
        }
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }
        pdf_name <- paste0(pdf_name,"_average", ".pdf")
        ggsave(filename = pdf_name, plot = stack_barplot_average, width = pdf.wid, height = pdf.hei)
      }

      # Return the stacked barplot
      return(stack_barplot_average)
    })

    # Name the plots in the list according to feature levels
    names(plot_list_all) <- feature.level

    # Return the list of plots
    return(plot_list_all)
  }