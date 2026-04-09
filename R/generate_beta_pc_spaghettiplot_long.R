#' Generate PC Spaghetti Plots for Longitudinal Data
#'
#' Creates spaghetti plots showing individual PC trajectories over time.
#' Displays both individual subject lines and group mean trends.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param pc.ind Numeric vector specifying which PC axes to plot. Default c(1, 2).
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A list of ggplot objects for each distance measure and PC axis.
#' @examples
#' \dontrun{
#'
#' # Load required libraries and data
#' library(vegan)
#'
#' # Example with ecam.obj dataset
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, "BC")
#' pc.obj <- mStat_calculate_PC(dist.obj)
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:10],
#'   group.var = "diet",
#'   strata.var = "delivery",
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:10],
#'   group.var = "diet",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' # Example with peerj32.obj dataset
#' data(peerj32.obj)
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_pc_spaghettiplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = c('BC'),
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
generate_beta_pc_spaghettiplot_long <- function(data.obj = NULL,
                                                dist.obj = NULL,
                                                pc.obj = NULL,
                                                pc.ind = c(1, 2),
                                                subject.var,
                                                time.var,
                                                t0.level = NULL,
                                                ts.levels = NULL,
                                                group.var = NULL,
                                                strata.var = NULL,
                                                adj.vars = NULL,
                                                dist.name = c("BC"),
                                                base.size = 16,
                                                theme.choice = "bw",
                                                custom.theme = NULL,
                                                palette = NULL,
                                                pdf = TRUE,
                                                file.ann = NULL,
                                                pdf.wid = 11,
                                                pdf.hei = 8.5,
                                                ...) {

  # Check if distance metrics are provided
  if (is.null(dist.name)){
    return()
  }

  # If distance object is not provided, calculate it from the data object
  if (is.null(dist.obj)) {
    # Process time variable and extract relevant metadata
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        subject.var, time.var, group.var, strata.var
      )))
    
    # Calculate beta diversity
    dist.obj <-
      mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
    
    # If adjustment variables are provided, calculate adjusted distances
    if (!is.null(adj.vars)) {
      dist.obj <-
        mStat_calculate_adjusted_distance(
          data.obj = data.obj,
          dist.obj = dist.obj,
          adj.vars = adj.vars,
          dist.name = dist.name
        )
    }
  } else {
    prepared_context <- mStat_prepare_precomputed_beta_context(
      dist.obj = dist.obj,
      dist.name = dist.name,
      pc.obj = pc.obj,
      data.obj = data.obj,
      time.var = time.var,
      t0.level = t0.level,
      ts.levels = ts.levels,
      process_time = TRUE,
      required_pc_axes = max(pc.ind)
    )
    data.obj <- prepared_context$data.obj
    dist.obj <- prepared_context$dist.obj
    pc.obj <- prepared_context$pc.obj
    meta_tab <- mStat_extract_dist_metadata(
      dist.obj = dist.obj,
      dist.name = dist.name,
      vars = c(subject.var, time.var, group.var, strata.var),
      data.obj = data.obj
    )
  }

  # Get color palette
  col <- mStat_get_palette(palette)

  # Calculate new sizes based on base.size for consistent plot aesthetics
  title.size = base.size * 1.25
  axis.title.size = base.size * 0.75
  axis.text.size = base.size * 0.5
  legend.title.size = base.size * 1
  legend.text.size = base.size * 0.75

  # Get the appropriate theme
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Create a list to store plots for each distance metric
  plot_list <- lapply(dist.name, function(dist.name) {
    # If principal component object is not provided, calculate it using MDS
    if (is.null(pc.obj)) {
      message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
      message(
        "If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method."
      )
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj[dist.name],
          method = "mds",
          k = max(pc.ind),
          dist.name = dist.name
        )
    }

    # Extract principal component coordinates
    pc.mat <- pc.obj[[dist.name]]$points

    df <- mStat_prepare_pc_long_data(
      pc.points = pc.mat,
      pc.ind = pc.ind,
      meta.dat = meta_tab,
      vars = c(subject.var, time.var, group.var, strata.var),
      sample_col = "sample",
      join = "inner"
    )

    # If no group variable is provided, create a dummy group
    if (is.null(group.var)){
      df <- df %>% dplyr::mutate("ALL" = "ALL")
      group.var = "ALL"
    }

    # Calculate mean values for each group and time point
    if (is.null(strata.var)) {
      df.mean <- df %>%
        dplyr::group_by(!!sym(time.var),!!sym(group.var), PC) %>%
        dplyr::summarize(mean_value = mean(value, na.rm = TRUE))
      df <-
        dplyr::left_join(df, df.mean, by = c(time.var, group.var, "PC"))
    } else {
      df.mean <- df %>%
        dplyr::group_by(!!sym(time.var),!!sym(group.var),!!sym(strata.var), PC) %>%
        dplyr::summarize(mean_value = mean(value, na.rm = TRUE))
      df <-
        dplyr::left_join(df, df.mean, by = c(time.var, group.var, strata.var, "PC"))
    }

    # Create a list to store plots for each principal component
    sub_plot_list <- lapply(unique(df$PC), function(pc.index) {
      sub_df <- df %>% filter(PC == pc.index)

      # Create the spaghetti plot
      p <- ggplot() +
        geom_point(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "value",
            group = subject.var,
            color = group.var
          ),
          alpha = 0.3,
          size = 3
        ) +
        geom_line(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_value",
            group = group.var,
            color = group.var
          ),
          size = 2
        ) +
        geom_point(
          data = sub_df,
          aes_string(
            x = time.var,
            y = "mean_value",
            group = group.var,
            color = group.var
          ),
          size = 5
        ) +
        labs(x = time.var, y = paste("Distance:",
                                     dist.name,
                                     " - Axis",
                                     gsub("PC", "", pc.index)), color = group.var) +
        scale_color_manual(
          values = col
        ) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = axis.title.size*2),
          axis.title.y = element_text(size = axis.title.size*2),
          axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size * 2),
          axis.text.y = element_text(size = base.size*2),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16 * 2),
          legend.title = ggplot2::element_text(size = 16 * 2),
          legend.key.size = unit(10, "mm"),
          legend.key.spacing = unit(2, "mm")
        )

      # Remove legend if there's only one group
      if (group.var == "ALL"){
        p <- p + theme(legend.position = "none")
      }

      # Add faceting if strata variable is provided
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(
          cols = vars(!!sym(strata.var)),
          scale = "free",
          space = "free"
        )
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0(
          "beta_pc_spaghettiplot_long_",
          dist.name,
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var
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
          plot = p,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      }
      return(p)
    })

    names(sub_plot_list) <- unique(df$PC)
    return(sub_plot_list)
  })

  names(plot_list) <- dist.name

  return(plot_list)
}
