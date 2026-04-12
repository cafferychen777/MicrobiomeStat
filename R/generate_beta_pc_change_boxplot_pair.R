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
    if (is.null(data.obj) && is.null(dist.obj)) {
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
      meta_tab <- mStat_extract_dist_metadata(
        dist.obj = dist.obj,
        dist.name = dist.name,
        vars = c(subject.var, time.var, group.var, strata.var),
        data.obj = data.obj
      )
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

      df <- mStat_prepare_pc_long_data(
        pc.points = pc.mat,
        pc.ind = pc.ind,
        meta.dat = meta_tab,
        vars = c(subject.var, time.var, group.var, strata.var),
        sample_col = "sample",
        join = "inner"
      )

      pair_change <- mStat_prepare_pc_pair_change_data(
        long.df = df,
        subject.var = subject.var,
        time.var = time.var,
        change.base = change.base,
        change.func = change.func,
        context = "beta PC change plotting"
      )
      combined_data <- pair_change$combined_data
      change.base <- pair_change$change.base
      change.after <- pair_change$change.after

      # Add metadata to the combined data
      combined_data <- mStat_attach_pair_metadata(
        df = combined_data,
        meta_tab = meta_tab,
        subject.var = subject.var,
        time.var = time.var,
        mode = "followup_time",
        change.after = change.after
      )

      placeholder_group <- mStat_ensure_group_placeholder(combined_data, group.var = group.var)
      combined_data <- placeholder_group$df
      resolved_group_var <- placeholder_group$group.var

      # Create a list to store plots for each principal component
      sub.plot_list <- lapply(paste0("PC", pc.ind), function(pc.index) {
        # Create the boxplot
        boxplot <- ggplot(
          combined_data %>% filter(PC == pc.index),
          aes(
            x = !!sym(resolved_group_var),
            y = value_diff,
            fill = !!sym(resolved_group_var)
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
          scale_fill_manual(values = col) +
          labs(
            x = if (is.null(group.var)) NULL else group.var,
            y = paste(
              "Change in ",
              "Axis ",
              gsub("PC", "", pc.index),
              " - ",
              if (is.function(change.func)) {
                "custom function"
              } else {
                change.func
              }
            )
          ) +
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
        if (!is.null(strata.var)) {
          boxplot <- boxplot +
            ggh4x::facet_nested(cols = vars(!!sym(strata.var)),
                                scales = "fixed",
                                space = "free")
        }

        # Adjust theme if there's only one group
        if (is.null(group.var)) {
          boxplot <- boxplot + theme(
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
