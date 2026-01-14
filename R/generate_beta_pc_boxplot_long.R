#' Generate PC Boxplots for Beta Diversity Over Time
#'
#' Creates boxplots of principal coordinate values across time points for
#' longitudinal microbiome data. Supports multiple dimension reduction methods.
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
#' # Example with ecam.obj dataset
#' data(ecam.obj)
#' generate_beta_pc_boxplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:4],
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
#'
#' # Example with peerj32.obj dataset
#' data(peerj32.obj)
#' generate_beta_pc_boxplot_long(
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
#'
#' }
#' @export
generate_beta_pc_boxplot_long <- function(data.obj = NULL,
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

  # Check if distance metrics are provided. If not, exit the function.
  if (is.null(dist.name)){
    return()
  }

  # This block handles the calculation or retrieval of distance matrices and metadata
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
    # Adjust distances if adjustment variables are provided
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
    # If distance object is provided, process metadata accordingly
    if (!is.null(data.obj) & !is.null(data.obj$meta.dat)) {
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <-
        data.obj$meta.dat %>% select(all_of(c(
          subject.var, time.var, group.var, strata.var
        )))
    } else {
      # Extract metadata from distance object if data object is not provided
      meta_tab <-
        attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(
          subject.var, time.var, group.var, strata.var
        )))
      data.obj <- list(meta.dat = meta_tab)
      data.obj <-
        mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat
      dist.obj <- mStat_subset_dist(dist.obj, colnames(meta_tab))
    }
  }

  # Determine the number of unique time points
  time.levels <-
    meta_tab %>% dplyr::select(all_of(c(time.var))) %>%
    pull() %>%
    as.factor() %>%
    levels() %>%
    length()

  # Get color palette for plotting
  col <- mStat_get_palette(palette)

  # Define aesthetic mappings for the plot
  # These determine how different variables will be represented visually
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

  # Define aesthetic mappings for the line plot
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

  # Get appropriate theme for plotting
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Generate plots for each distance metric
  plot_list <- lapply(dist.name, function(dist.name) {
    # Perform dimension reduction if not already done
    if (is.null(pc.obj)) {
      message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
      message(
        "If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method."
      )
      # Use Multidimensional Scaling (MDS) to reduce dimensions
      # This allows for visualization of high-dimensional distance data in lower-dimensional space
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj[dist.name],
          method = "mds",
          k = max(pc.ind),
          dist.name = dist.name
        )
    }

    # Extract principal coordinates
    pc.mat <- pc.obj[[dist.name]]$points

    # Rename columns to PC1, PC2, etc.
    colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

    # Convert to tibble for easier manipulation
    pc.mat <- pc.mat %>% as_tibble()

    # Combine principal coordinates with metadata
    df <-
      cbind(pc.mat[, paste0("PC", pc.ind)], meta_tab[, c(subject.var, time.var, group.var, strata.var)])

    # Reshape data from wide to long format
    df <-
      df %>%
      as_tibble() %>%
      tidyr::gather(
        key = "PC",
        value = "value",
        -all_of(subject.var, group.var, time.var, strata.var)
      )

    # Count unique subjects and time points
    n_subjects <- length(unique(df[[subject.var]]))
    n_times <- length(unique(df[[time.var]]))

    # Generate plots for each principal coordinate
    sub_plot_list <- lapply(unique(df$PC), function(pc.index) {
      sub_df <- df %>% filter(PC == pc.index)

      # Calculate average values if there are many subjects or time points
      # This helps to simplify the visualization when there's a lot of data
      average_sub_df <- NULL
      if (n_times > 10 || n_subjects > 25) {
        if (!is.null(strata.var) & !is.null(group.var)) {
          average_sub_df <- sub_df %>%
            dplyr::group_by(!!sym(strata.var),
                            !!sym(group.var),
                            !!sym(time.var)) %>%
            dplyr::summarise(dplyr::across(value, mean, na.rm = TRUE), .groups = "drop") %>%
            dplyr::ungroup() %>%
            dplyr::mutate(!!sym(subject.var) := "ALL")
        } else if (!is.null(group.var)) {
          average_sub_df <- sub_df %>%
            dplyr::group_by(!!sym(group.var),!!sym(time.var)) %>%
            dplyr::summarise(dplyr::across(value, mean, na.rm = TRUE), .groups = "drop") %>%
            dplyr::ungroup() %>%
            dplyr::mutate(!!sym(subject.var) := "ALL")
        } else {
          average_sub_df <- sub_df %>%
            dplyr::group_by(!!sym(time.var)) %>%
            dplyr::summarise(dplyr::across(value, mean, na.rm = TRUE), .groups = "drop") %>%
            dplyr::ungroup() %>%
            dplyr::mutate(!!sym(subject.var) := "ALL")
        }
      }

      # Create the boxplot
      boxplot <- ggplot(sub_df,
                        aes_function) +
        stat_boxplot(
          geom = "errorbar",
          position = position_dodge(width = 0.2),
          width = 0.3
        ) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.3,
        ) +
        geom_line(
          line_aes_function,
          alpha = 0.8,
          linewidth = 0.6,
          color = "black",
          linetype = "dashed",
          data = if (!is.null(average_sub_df))
            average_sub_df
          else
            sub_df
        ) +
        scale_fill_manual(values = col) +
        labs(x = time.var,
             y = paste("Distance:",
                       dist.name,
                       " - Axis",
                       gsub("PC", "", pc.index))) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = base.size),
          legend.title = ggplot2::element_text(size = base.size)
        )

      # Add faceting based on the number of time points and grouping variables
      if (time.levels > 2){
        if (!is.null(group.var)) {
          if (is.null(strata.var)) {
            boxplot <-
              boxplot + ggh4x::facet_nested(cols = vars(!!sym(group.var)),
                                            scales = "free",
                                            space = "free")
          } else {
            boxplot <-
              boxplot + ggh4x::facet_nested(
                cols = vars(!!sym(group.var)),
                rows = vars(!!sym(strata.var)),
                scales = "free",
                space = "free"
              )
          }
        }
      } else {
        if (!is.null(group.var)) {
          if (is.null(strata.var)) {
            boxplot <-
              boxplot + ggh4x::facet_nested(cols = vars(!!sym(group.var)),
                                            scales = "free",
                                            space = "free")
          } else {
            boxplot <-
              boxplot + ggh4x::facet_nested(
                cols = vars(!!sym(strata.var), !!sym(group.var)),
                scales = "free",
                space = "free"
              )
          }
        }
      }

      # Add jittered points if there are many subjects or time points
      # This helps to show the distribution of individual data points
      if (n_subjects > 10 || n_times > 10) {
        boxplot <- boxplot + geom_jitter(width = 0.1,
                                         alpha = 0.1,
                                         size = 1)
      }

      # Save the plot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "beta_pc_boxplot_long_",
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
          plot = boxplot,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      }
      return(boxplot)
    })

    # Assign names to the elements of sub_plot_list
    names(sub_plot_list) <- unique(df$PC)
    return(sub_plot_list)
  })

  # Assign names to the elements of plot_list
  names(plot_list) <- dist.name

  return(plot_list)
}