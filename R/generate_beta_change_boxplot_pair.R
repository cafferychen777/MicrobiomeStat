#' Generate Beta Diversity Change Boxplot for Paired Samples
#'
#' Creates boxplots visualizing within-subject beta diversity changes between
#' two time points. Compares distributions across groups and strata.
#'
#' @name generate_beta_change_boxplot_pair
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param change.base Character or numeric specifying the baseline time point
#'   for computing beta diversity change.
#' @param ... Additional parameters passed on to ggsave()
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object.
#'
#' @return A named list of ggplot objects visualizing beta diversity change.
#' The list contains one plot for each distance metric specified in \code{dist.name}.
#' Each plot shows boxplots of the beta diversity changes, with samples faceted by
#' the \code{group_var} and \code{strata_var} variables if provided.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(vegan)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC")
#' generate_beta_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#' generate_beta_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   change.base = "Baseline",
#'   dist.name = c('BC'),
#'   base.size = 20,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_beta_change_boxplot_pair <-
  function(data.obj = NULL,
           dist.obj = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           change.base,
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

    # Exit the function if no distance metric is specified
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    # This step ensures we have the necessary distance matrices for the analysis
    if (is.null(dist.obj)) {
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Adjust distances for covariates if specified
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # Extract metadata if distance object is provided
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
        data.obj <- list(meta.dat = meta_tab)
        meta_tab <- data.obj$meta.dat
      }
    }

    # Add sample names to metadata
    meta_tab <- meta_tab %>% rownames_to_column("sample")

    # Set the baseline time point for change calculation if not provided
    if (is.null(change.base)){
      change.base <- unique(meta_tab %>% dplyr::select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'alpha_df' data frame: ", change.base)
    }

    # Identify the time point after the baseline for change calculation
    change.after <-
      unique(meta_tab %>% dplyr::select(all_of(c(time.var))))[unique(meta_tab %>% dplyr::select(all_of(c(time.var)))) != change.base]

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Set the theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each specified distance metric
    plot_list <- lapply(dist.name,function(dist.name){

      # Set y-axis label based on whether distances are adjusted for covariates
      if (is.null(adj.vars)) {
        y_label <- paste0("Distance from ", change.base, " to ", change.after)
      } else {
        y_label <- paste0("Distance from ", change.base, " to ", change.after, " (adjusted by: ", paste(adj.vars, collapse = ", "), ")")
      }

      # Convert distance matrix to long format for plotting
      dist.df <- as.matrix(dist.obj[[dist.name]]) %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      # Calculate pairwise distances between baseline and follow-up time points for each subject
      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      # Add group and strata information to the long-format data
      long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(-time.var) %>% dplyr::distinct(), by = subject.var)

      # Set up faceting formula based on presence of strata variable
      facet_formula <-
        if (!is.null(strata.var)) {
          paste(". ~", strata.var)
        } else {
          ". ~ 1"
        }

      # Create alternative x-axis variable for cases without grouping
      long.df$x_alternative <- "ALL"
      
      # Set up aesthetic mapping based on presence of grouping variable
      aes_function <- if (!is.null(group.var)){
        aes(
          x = !!sym(group.var),
          y = distance,
          fill = !!sym(group.var)
        )
      }else{
        aes(
          x = x_alternative,
          y = distance,
          fill = x_alternative
        )
      }

      # Create the boxplot
      p <-
        ggplot(long.df, aes_function) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.3) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.3,
        ) +
        geom_jitter(width = 0.3, alpha = 0.5, size = 1.7) +
        scale_fill_manual(values = col) +
        facet_wrap(as.formula(facet_formula), scales = "fixed") +
        xlab(group.var) +
        ylab(y_label) +
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

      # Adjust plot aesthetics based on presence of grouping and strata variables
      if(is.null(group.var) && is.null(strata.var)) {
        p <- p + theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          strip.text.x = element_blank()
        )
      }

      if (!is.null(group.var) & is.null(strata.var)){
        p <- p + theme(
          strip.text.x = element_blank()
        )
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0("beta_change_boxplot_pair_",
                           dist.name,
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "change_base_", change.base)

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
    
    # Name the plots in the list according to the distance metrics
    names(plot_list) <- dist.name
    return(plot_list)
  }