#' Generate PC Change Boxplots for Paired Samples
#'
#' Creates boxplots showing changes in PC coordinates between two time points.
#' Useful for visualizing individual-level shifts in ordination space.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param pc.ind Numeric vector specifying which PC axes to plot. Default c(1, 2).
#' @param change.base Character or numeric specifying the baseline time point.
#' @param change.func Function or "absolute change" for computing differences.
#' @param ... Additional arguments passed to the function.
#'
#' @return A named list of ggplot objects for each PC axis and distance metric.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(ggh4x)
#' data(peerj32.obj)
#' generate_beta_pc_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   change.base = "1",
#'   change.func = "absolute change",
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
#' data(peerj32.obj)
#' generate_beta_pc_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   change.func = "absolute change",
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
generate_beta_pc_change_boxplot_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           pc.ind = c(1, 2),
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           change.base = NULL,
           change.func = "absolute change",
           dist.name = c('BC', 'Jaccard'),
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

    # Ensure that at least one of data.obj or dist.obj is provided
    if (is.null(data.obj) & is.null(dist.obj)) {
      stop("Both data.obj and dist.obj cannot be NULL. Please provide at least one.")
    }

    # If distance object is not provided, calculate it from the data object
    if (is.null(dist.obj)) {
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      meta_tab <- data.obj$meta.dat
      if (is.null(meta_tab)) {
        stop("No metadata could be loaded from data.obj. Please ensure it contains the necessary metadata.")
      }
      # If adjustment variables are provided, calculate adjusted distances
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # Check if all requested distance metrics are available in the provided distance object
      if (!all(dist.name %in% names(dist.obj))) {
        stop(paste0("The requested dist.name(s) ", paste(dist.name[!dist.name %in% names(dist.obj)], collapse = ", "),
                    " are not available in the provided dist.obj. Please check again."))
      }
      # Extract metadata from the distance object
      meta_tab <- attr(dist.obj[[dist.name[1]]], "labels")
      if (is.null(meta_tab)) {
        message("No metadata found in dist.obj. Attempting to load metadata from data.obj.")
        if (is.null(data.obj)) {
          stop("No data.obj provided to load metadata from. Please ensure either dist.obj or data.obj contain the necessary metadata.")
        }
        meta_tab <- data.obj$meta.dat
        if (is.null(meta_tab)) {
          stop("No metadata could be loaded from data.obj. Please ensure it contains the necessary metadata.")
        }
      }
    }

    # If principal component object is not provided, calculate it using MDS
    if (is.null(pc.obj)) {
      message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
      message("If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method.")
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = max(pc.ind),
          dist.name = dist.name
        )
    }

    # Get color palette
    col <- mStat_get_palette(palette)

    # Get the appropriate theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Create a list to store plots for each distance metric
    plot_list <- lapply(dist.name, function(dist.name) {
      # Extract principal component coordinates
      pc.mat <- pc.obj[[dist.name]]$points

      colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

      pc.mat <- pc.mat %>% as_tibble()

      # Combine PC coordinates with metadata
      df <-
        cbind(pc.mat[, paste0("PC", pc.ind)], meta_tab[, c(subject.var, time.var, group.var, strata.var)])

      # Identify the time point to compare with the baseline
      change.after <-
        unique(df %>% select(all_of(c(time.var))))[unique(df %>% select(all_of(c(time.var)))) != change.base]

      # Reshape the data from wide to long format
      df <-
        df %>%
        as_tibble() %>%
        tidyr::gather(
          key = "PC",
          value = "value",-all_of(subject.var, group.var, time.var, strata.var)
        )

      # Split the data by time points
      split_data <-
        split(df, f = df %>%
                dplyr::group_by(!!sym(time.var)) %>% select(!!sym(time.var)))

      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # Combine data from two time points
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c("PC", subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # Calculate the change in PC coordinates
      combined_data <- combined_data %>%
        dplyr::mutate(value_diff = compute_taxa_change(
          value_after  = value_time_2,
          value_before = value_time_1,
          method       = change.func,
          verbose      = FALSE
        ))

      # Add metadata to the combined data
      combined_data <-
        combined_data %>% dplyr::left_join(meta_tab %>% select(all_of(
          c(subject.var, time.var, group.var, strata.var)
        )) %>% filter(!!sym(time.var) == change.after),
        by = subject.var)

      # If no group variable is provided, create a dummy group
      if (is.null(group.var)) {
        group.var = "ALL"
        combined_data$ALL <- "ALL"
      }

      # Create a list to store plots for each principal component
      sub.plot_list <- lapply(paste0("PC", pc.ind), function(pc.index) {
        # Create the boxplot
        boxplot <- ggplot(
          combined_data %>% filter(PC == pc.index),
          aes(
            x = !!sym(group.var),
            y = value_diff,
            fill = !!sym(group.var)
          )
        ) +
          #geom_violin(trim = F, alpha = 0.8) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.2),
            width = 0.3
          ) +
          geom_boxplot(
            position = position_dodge(width = 0.8),
            width = 0.3
            #fill = "white"
          ) +
          geom_jitter(width = 0.3,
                      alpha = 0.5,
                      size = 1.7) +
          scale_alpha_manual(values = c(0.5, 0.5)) +
          scale_fill_manual(values = col) +
          labs(x = group.var, y = paste("Change in ", "Axis ", gsub("PC", "", pc.index), " - ",
                                        if(is.function(change.func)){
                                          "custom function"
                                        } else {
                                          change.func
                                        })) +
          theme_to_use +
          theme(
            panel.spacing.x = unit(0, "cm"),
            panel.spacing.y = unit(0, "cm"),
            strip.text.x = element_text(size = 12, color = "black"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black", size = base.size),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = base.size),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            legend.text = ggplot2::element_text(size = base.size),
            legend.title = ggplot2::element_text(size = base.size)
          )

        # Add faceting if strata variable is provided
        if (is.null(strata.var)) {
          boxplot <- boxplot
        } else {
          boxplot <- boxplot +
            ggh4x::facet_nested(cols = vars(!!sym(strata.var)),
                                scales = "fixed",
                                space = "free")
        }

        # Adjust theme if there's only one group
        if (group.var == "ALL") {
          boxplot <- boxplot  + theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none",
            strip.text.x = element_blank()
          )
        }

        # Save the plot as a PDF if requested
        if (pdf) {
          pdf_name <- paste0(
            "beta_pc_change_boxplot_pair_",
            dist.name,
            "_",
            "subject_",
            subject.var,
            "_",
            "time_",
            time.var,
            "_",
            "change_base_",
            change.base
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
      names(sub.plot_list) <- paste0("PC", pc.ind)
      return(sub.plot_list)
    })
    names(plot_list) <- dist.name
    return(plot_list)
  }