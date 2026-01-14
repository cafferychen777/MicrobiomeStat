#' Generate Beta Diversity Ordination Plot for Cross-Sectional Data
#'
#' Creates PCoA ordination plots for single time point or cross-sectional
#' microbiome data with optional group comparisons.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param t.level Character string specifying the time level to subset to.
#'   If NULL, uses all data.
#' @param point.size Numeric. Size of points in the scatter plot. Default is 4.
#' @param ... Additional arguments passed to underlying functions.
#'
#' @return A list of PCoA plots for each distance metric.
#' @seealso \code{\link{mStat_calculate_beta_diversity}}, \code{\link{mStat_calculate_PC}}
#'
#' @author Caffery Yang \email{cafferychen7850@@gmail.com}
#'
#' @examples
#' \dontrun{
#' library(aplot)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC','Jaccard'))
#' generate_beta_ordination_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   time.var = "time",
#'   t.level = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   point.size = 4,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_T2D.obj")
#' dist.obj <- mStat_calculate_beta_diversity(subset_T2D.obj, dist.name = c('BC', 'Jaccard'))
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC','Jaccard'))
#' generate_beta_ordination_single(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   time.var = "visit_number_num",
#'   t.level = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = "sample_body_site",
#'   dist.name = c("BC", 'Jaccard'),
#'   base.size = 20,
#'   point.size = 4,
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
generate_beta_ordination_single <-
  function(data.obj,
           time.var = NULL,
           t.level = NULL,
           group.var = NULL,
           adj.vars = NULL,
           strata.var = NULL,
           dist.obj = NULL,
           dist.name = c('BC', 'Jaccard'),
           pc.obj = NULL,
           base.size = 16,
           point.size = 4,
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
      # Process data based on time variable if it exists
      if (!is.null(time.var)){
        if (!is.null(t.level)){
          # Subset data to specific time point if t.level is provided
          condition <- paste(time.var, "== '", t.level, "'", sep = "")
          data.obj <- mStat_subset_data(data.obj, condition = condition)
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
          dist.obj <-
            mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        } else {
          # If no specific time point is provided, use all time points but warn if multiple exist
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
          if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
            message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
          }
          dist.obj <-
            mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        }
      } else {
        # If no time variable is provided, proceed with all data
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
        dist.obj <-
          mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      }
      # Adjust distances if adjustment variables are provided
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # If distance object is provided, extract metadata accordingly
      if (!is.null(data.obj)){
        if (!is.null(time.var)){
          if (!is.null(t.level)){
            condition <- paste(time.var, "== '", t.level, "'", sep = "")
            data.obj <- mStat_subset_data(data.obj, condition = condition)
            meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
          } else {
            meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
            if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
              message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
            }
          }
        } else {
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
        }
      }
      if (!is.null(attr(dist.obj[[dist.name[1]]], "labels"))){
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels")  %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
      }
    }

    # Perform dimension reduction if not already done
    # This step reduces the high-dimensional distance data to 2D for visualization
    if (is.null(pc.obj)) {
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

    # If no group variable is provided, create a dummy "ALL" group
    if (is.null(group.var)){
      group.var = "ALL"
      meta_tab$ALL <- "ALL"
    }

    # Define aesthetic mapping based on presence of strata variable
    aes_function <- if (!is.null(strata.var)) {
      aes(color = !!sym(group.var),
          shape = !!sym(strata.var))
    } else {
      aes(color = !!sym(group.var))
    }

    # Get appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each distance metric
    plot_list <- lapply(dist.name, function(dist.name) {

      # Extract the first two principal coordinates
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]

      # Prepare data frame for plotting
      df <-
        pc.mat %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        dplyr::left_join(meta_tab %>%
                           dplyr::select(all_of(c(time.var, group.var, strata.var))) %>%
                           rownames_to_column("sample"), by = "sample")

      # Filter out NA values in time variable if it exists
      if (!is.null(time.var)){
        df <- df %>% dplyr::filter(!is.na(!!sym(time.var)))
      }

      df <- df %>% column_to_rownames("sample")

      colnames(df)[1:2] <- c("PC1", "PC2")

      # Create the main scatter plot
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(size = point.size, aes_function, show.legend = T) +
        ggplot2::labs(
          x = ifelse(!is.null(pc.obj[[dist.name]]$eig),paste0("Axis 1 (", round(pc.obj[[dist.name]]$eig[1]/sum(pc.obj[[dist.name]]$eig)*100,2),"%)"),"Axis 1"),
          y = ifelse(!is.null(pc.obj[[dist.name]]$eig),paste0("Axis 2 (", round(pc.obj[[dist.name]]$eig[2]/sum(pc.obj[[dist.name]]$eig)*100,2),"%)"),"Axis 2")
        ) +
        # Add 95% confidence ellipses for each group
        ggplot2::stat_ellipse(ggplot2::aes(color = !!sym(group.var)),fill="white",geom = "polygon",
                              level=0.95,alpha = 0.01,show.legend = F) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed",
          color = "black"
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "black"
        ) +
        theme_to_use  +
        ggplot2::theme(
          axis.line.x = ggplot2::element_line(size = 1, colour = "black"),
          axis.line.y = ggplot2::element_line(size = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title = ggplot2::element_text(color = "black", size = 20),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(
            color = "black",
            size = base.size
          ),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = ggplot2::element_text(size = base.size),
          legend.title = ggplot2::element_text(size = base.size)
        )

      # Set color scale based on grouping
      if (group.var == "ALL") {
        p <- p + scale_color_manual(values = col, guide = "none")
      } else {
        p <- p + scale_color_manual(values = col)
      }

      # Create a boxplot for PC1 values
      Fig1a.taxa.pc1.boxplot <-
        ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x=!!sym(group.var), y=PC1, fill=!!sym(group.var)), color="black", alpha=0.5, show.legend = F) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand = c(0,0.001)) +
        ggplot2::labs(x=NULL, y=NULL) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.text.y =ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())+
        ggplot2::coord_flip()

      # Create a boxplot for PC2 values
      Fig1a.taxa.pc2.boxplot <-
        ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x=!!sym(group.var), y=PC2, fill=!!sym(group.var)), color="black", alpha=0.5, show.legend = F) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand = c(0,0.001)) +
        ggplot2::labs(x=NULL, y=NULL) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.text.y =ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())

      # Combine the main plot with the boxplots
      # This provides a comprehensive view of the data distribution along both principal coordinates
      p <- p %>%
        aplot::insert_top(Fig1a.taxa.pc1.boxplot, height = 0.2) %>%
        aplot::insert_right(Fig1a.taxa.pc2.boxplot, width=0.2) %>%
        as.ggplot()

      # Save the plot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0("beta_ordination_single_", "dist.name_", dist.name)
        if (!is.null(time.var)) {
          pdf_name <- paste0(pdf_name, "_", "time_", time.var)
        }
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

    # Assign names to the elements of the plot list
    names(plot_list) <- dist.name
    return(plot_list)
  }