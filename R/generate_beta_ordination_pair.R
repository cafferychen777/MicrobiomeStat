#' Generate Beta Diversity Ordination Plots for Paired Samples
#'
#' Creates PCoA ordination plots for paired samples showing change between
#' two time points with arrows connecting paired samples.
#'
#' @name generate_beta_ordination_pair
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param ... Additional arguments passed to the plotting function.
#'
#' @return A list of PCoA plots for each distance metric.
#' @seealso \code{\link{mStat_calculate_beta_diversity}}, \code{\link{mStat_calculate_PC}}
#'
#'
#' @examples
#' \dontrun{
#'
#' # Load necessary libraries and data
#' library(vegan)
#' data(peerj32.obj)
#'
#' # Perform beta ordination pair analysis using `generate_beta_ordination_pair`
#' generate_beta_ordination_pair(
#'   data.obj      = peerj32.obj,
#'   dist.obj      = NULL,
#'   pc.obj        = NULL,
#'   subject.var   = "subject",
#'   time.var      = "time",
#'   group.var     = "group",
#'   strata.var    = "sex",
#'   adj.vars      = "sex",
#'   dist.name     = c("BC"),
#'   base.size     = 16,
#'   theme.choice  = "bw",
#'   custom.theme  = NULL,
#'   palette       = NULL,
#'   pdf           = TRUE,
#'   file.ann      = NULL,
#'   pdf.wid       = 11,
#'   pdf.hei       = 8.5
#' )
#'
#' data(subset_pairs.obj)
#'
#' # Perform beta ordination pair analysis using `generate_beta_ordination_pair`
#' generate_beta_ordination_pair(
#'   data.obj      = subset_pairs.obj,
#'   dist.obj      = NULL,
#'   pc.obj        = NULL,
#'   subject.var   = "MouseID",
#'   time.var      = "Antibiotic",
#'   group.var     = "Sex",
#'   strata.var    = NULL,
#'   adj.vars      = NULL,
#'   dist.name     = c("BC"),
#'   base.size     = 16,
#'   theme.choice  = "bw",
#'   custom.theme  = NULL,
#'   palette       = NULL,
#'   pdf           = TRUE,
#'   file.ann      = NULL,
#'   pdf.wid       = 11,
#'   pdf.hei       = 8.5
#' )
#' }
#' @export
generate_beta_ordination_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           pc.obj = NULL,
           subject.var,
           time.var,
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

    # Check if distance metrics are provided. If not, exit the function.
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    if (is.null(dist.obj)) {
      # Compute beta diversity distances using specified metrics
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Extract relevant metadata
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      
      # Adjust distances if adjustment variables are provided
      # This step helps to account for potential confounding factors
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
      print(dist.obj)
    } else {
      # If distance object is provided, extract metadata accordingly
      if (is.null(data.obj)) {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      } else {
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var)))
      }
    }

    # Perform dimension reduction if not already done
    if (is.null(pc.obj)) {
      # Use Multidimensional Scaling (MDS) to reduce dimensions to 2
      # This allows for visualization of high-dimensional distance data in 2D space
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = 2,
          dist.name = dist.name
        )
    }

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Define aesthetic mapping based on presence of group variable
    # This determines how different variables will be represented visually
    aes_function <- if (!is.null(group.var)) {
      aes(color = !!sym(group.var),
          shape = !!sym(time.var))
    } else {
      aes(color = !!sym(time.var))
    }

    # Get appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each distance metric
    plot_list <- lapply(dist.name, function(dist.name) {
      # Extract the first two principal coordinates
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]
      
      # Prepare data frame for plotting
      df <- as.data.frame(pc.mat) %>%
        setNames(c("PC1", "PC2")) %>%
        dplyr::bind_cols(meta_tab[, c(subject.var, time.var, group.var, strata.var)]) %>%
        dplyr::mutate(
          x_start = PC1,
          y_start = PC2,
          x_end = NA,
          y_end = NA
        )

      # Get unique time points
      Time_choices <-
        df %>% dplyr::select(all_of(time.var)) %>% dplyr::pull() %>% unique()

      # Calculate end points for arrows
      # This creates the visual effect of change over time for each subject
      df <- df %>%
        dplyr::arrange(!!sym(subject.var),!!sym(time.var)) %>% 
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::mutate(x_end = dplyr::lead(PC1),
                      y_end = dplyr::lead(PC2)) %>%
        dplyr::ungroup()

      # Create a dataset with all points for each facet
      # This allows for comparison across strata while maintaining context
      if (!is.null(strata.var)) {
        all_strata <- unique(df[[strata.var]])
        df_all <- do.call(rbind, lapply(all_strata, function(s) {
          df_temp <- df
          df_temp$facet <- s
          df_temp$is_facet <- df_temp[[strata.var]] == s
          return(df_temp)
        }))
      } else {
        df_all <- df
        df_all$facet <- "All"
        df_all$is_facet <- TRUE
      }

      # Create base plot
      p <- ggplot2::ggplot(df_all, ggplot2::aes(PC1, PC2))

      # Add gray points and lines for all data
      # This provides context for the highlighted data in each facet
      p <- p +
        ggplot2::geom_point(data = subset(df_all, !is_facet),
                            aes_function,
                            size = 10, alpha = 0.5, color = "gray80") +
        ggplot2::geom_segment(
          data = subset(df_all, !is_facet),
          aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
          color = "gray80", size = 1, alpha = 0.5,
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open")
        )

      # Define aesthetic mapping for highlighted data
      aes_function <- if (!is.null(group.var)) {
        aes(shape = !!sym(time.var), color = !!sym(group.var))
      } else {
        aes(shape = !!sym(time.var), color = !!sym(time.var))
      }

      # Add highlighted points and arrows for the current facet
      p <- p +
        ggplot2::geom_point(data = subset(df_all, is_facet), aes_function, size = 10, show.legend = TRUE) +
        ggplot2::geom_segment(
          data = subset(df_all, is_facet),
          aes(
            x = x_start, y = y_start, xend = x_end, yend = y_end,
            color = if (!is.null(group.var)) !!sym(group.var) else !!sym(time.var)
          ),
          arrow = ggplot2::arrow(length = unit(0.25, "cm"), type = "open"),
          size = 1
        )

      # Add labels and customize plot appearance
      p <- p +
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
        scale_color_manual(values = col) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        theme_to_use +
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
          axis.text = ggplot2::element_text(color = "black", size = base.size),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      # Add faceting if strata variable is provided
      if (!is.null(strata.var)) {
        p <- p + ggh4x::facet_nested(as.formula(paste("~ facet")))
      }

      # Save the plot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0(
          "beta_ordination_pair_",
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
    names(plot_list) <- dist.name
    return(plot_list)
  }
