#' Generate Beta Diversity Ordination Plots for Longitudinal Data
#'
#' Creates ordination plots (PCoA/NMDS) showing sample trajectories over time,
#' with arrows connecting time points for each subject.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return A list of ggplot2 objects representing the beta ordination plots.
#' @seealso \code{\link{mStat_calculate_beta_diversity}}, \code{\link{mStat_calculate_PC}}
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' generate_beta_ordination_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = NULL,
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, "BC")
#' pc.obj <- mStat_calculate_PC(dist.obj)
#' generate_beta_ordination_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = "0",
#'   ts.levels = as.character(sort(as.numeric(unique(ecam.obj$meta.dat$month))))[2:10],
#'   group.var = "diet",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' generate_beta_ordination_long(
#'   data.obj = ecam.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = NULL,
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   dist.name = 'BC',
#'   base.size = 16,
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
generate_beta_ordination_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
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

    # Calculate beta diversity if not provided
    if (is.null(dist.obj)) {
      # Process time variable and extract relevant metadata
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      
      # Adjust distances if adjustment variables are provided
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # If distance object is provided, process metadata accordingly
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        data.obj <-
          mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      } else {
        # Extract metadata from distance object if data object is not provided
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
        data.obj <- list(meta.dat = meta_tab)
        data.obj <- mStat_process_time_variable(meta_tab, time.var, t0.level, ts.levels)
        meta_tab <- data.obj$meta.dat
      }
    }

    # Calculate principal coordinates if not provided
    if (is.null(pc.obj)) {
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = 2,
          dist.name = dist.name
        )
    }

    # Get color palette
    col <- mStat_get_palette(palette)

    # Define aesthetic mapping based on presence of group variable
    aes_function <- if (!is.null(group.var)) {
      aes(color = !!sym(group.var),
          alpha = !!sym(time.var))
    } else {
      aes(color = !!sym(time.var))
    }

    # Get appropriate theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each distance metric
    plot_list <- lapply(dist.name, function(dist.name) {
      # Extract principal coordinates
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]

      # Prepare data frame for plotting
      df <- as.data.frame(pc.mat) %>%
        setNames(c("PC1", "PC2")) %>%
        rownames_to_column("sample") %>%
        dplyr::inner_join(meta_tab %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var))) %>% rownames_to_column("sample"), by = "sample") %>%
        dplyr::mutate(x_start = PC1,
               y_start = PC2,
               x_end = NA,
               y_end = NA)

      # Calculate end points for arrows
      df <- df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::arrange(!!sym(time.var)) %>%
        dplyr::mutate(
          end_condition = if (is.factor(!!sym(time.var))) {
            levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
          } else {
            max(!!sym(time.var))
          },
          x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(PC1)),
          y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(PC2))
        ) %>%
        dplyr::ungroup()

      # Calculate mean positions for different grouping scenarios
      if (!is.null(strata.var) & !is.null(group.var)) {
        # Case: Both strata and group variables are present
        df_mean <- df %>%
          dplyr::group_by(!!sym(time.var), !!sym(group.var), !!sym(strata.var)) %>%
          dplyr::summarise(mean_PC1 = mean(PC1, na.rm = TRUE),
                           mean_PC2 = mean(PC2, na.rm = TRUE)) %>%
          dplyr::ungroup()

        df_mean <- df_mean %>%
          dplyr::group_by(!!sym(group.var), !!sym(strata.var)) %>%
          dplyr::arrange(!!sym(time.var)) %>%
          dplyr::mutate(
            end_condition = if (is.factor(!!sym(time.var))) {
              levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
            } else {
              max(!!sym(time.var))
            },
            x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC1)),
            y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC2))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(x_start = mean_PC1,
                        y_start = mean_PC2)

      } else if (!is.null(group.var)){
        # Case: Only group variable is present
        df_mean <- df %>%
          dplyr::group_by(!!sym(time.var), !!sym(group.var)) %>%
          dplyr::summarise(mean_PC1 = mean(PC1, na.rm = TRUE),
                           mean_PC2 = mean(PC2, na.rm = TRUE)) %>%
          dplyr::ungroup()

        df_mean <- df_mean %>%
          dplyr::group_by(!!sym(group.var)) %>%
          dplyr::arrange(!!sym(time.var)) %>%
          dplyr::mutate(
            end_condition = if (is.factor(!!sym(time.var))) {
              levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
            } else {
              max(!!sym(time.var))
            },
            x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC1)),
            y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC2))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(x_start = mean_PC1,
                        y_start = mean_PC2)

      } else {
        # Case: No grouping variables
        df_mean <- df %>%
          dplyr::group_by(!!sym(time.var)) %>%
          dplyr::summarise(mean_PC1 = mean(PC1, na.rm = TRUE),
                           mean_PC2 = mean(PC2, na.rm = TRUE)) %>%
          dplyr::ungroup()

        df_mean <- df_mean %>%
          dplyr::arrange(!!sym(time.var)) %>%
          dplyr::mutate(
            end_condition = if (is.factor(!!sym(time.var))) {
              levels(!!sym(time.var))[length(levels(!!sym(time.var)))]
            } else {
              max(!!sym(time.var))
            },
            x_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC1)),
            y_end = dplyr::if_else(!!sym(time.var) == end_condition, NA_real_, dplyr::lead(mean_PC2))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(x_start = mean_PC1,
                        y_start = mean_PC2)

      }

      # Create the plot
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(
          size = 5,
          aes_function,
          show.legend = T
        ) +
        ggplot2::geom_segment(
          aes(
            x = x_start,
            y = y_start,
            xend = x_end,
            yend = y_end
          ),
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
          size = 0.7,
          color = "gray70",
          alpha = 0.3
        ) +
        {
          # Add group-specific or time-specific arrows
          if (!is.null(group.var)){
            ggplot2::geom_segment(data = df_mean,
                                  aes(x = x_start,
                                      y = y_start,
                                      xend = x_end,
                                      yend = y_end,
                                      color = !!sym(group.var)),
                                  arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
                                  size = 1.5)
          } else {
            ggplot2::geom_segment(data = df_mean,
                                  aes(x = x_start,
                                      y = y_start,
                                      xend = x_end,
                                      yend = y_end,
                                      color = !!sym(time.var)),
                                  arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
                                  size = 1.5)
          }
        } +
        # Add axis labels with explained variance
        ggplot2::labs(
          x = ifelse(
            !is.null(pc.obj[[dist.name]]$eig),
            paste0("Axis 1 (", round(
              pc.obj[[dist.name]]$eig[1] / sum(pc.obj[[dist.name]]$eig) * 100, 2
            ), "%)"),
            "Axis 1"
          ),
          y = ifelse(
            !is.null(pc.obj[[dist.name]]$eig),
            paste0("Axis 2 (", round(
              pc.obj[[dist.name]]$eig[2] / sum(pc.obj[[dist.name]]$eig) * 100, 2
            ), "%)"),
            "Axis 2"
          )
        ) +
        # Set color scale based on grouping variables
        {
          if (!is.null(group.var) | !is.null(strata.var)){
            scale_color_manual(values = col)
          } else {
            ggplot2::scale_color_gradientn(colors = c("#92c5de", "#0571b0", "#f4a582", "#ca0020"))
          }
        } +
        # Add reference lines
        ggplot2::geom_vline(xintercept = 0,
                            linetype = "dashed",
                            color = "black") +
        ggplot2::geom_hline(yintercept = 0,
                            linetype = "dashed",
                            color = "black") +
        theme_to_use  +
        # Customize theme elements
        ggplot2::theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0, "cm"),
          axis.line.x = ggplot2::element_line(size = 1, colour = "black"),
          axis.line.y = ggplot2::element_line(size = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title = ggplot2::element_text(color = "black"),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color = "black",
                                            size = base.size),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      # Add faceting if strata variable is present
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested_wrap(as.formula(paste(".~", strata.var)), ncol = 3)
      }

      # Save the plot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "beta_ordination_long_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "dist.name_",
          dist.name
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

    # Assign names to the plot list
    names(plot_list) <- dist.name

    return(plot_list)
  }