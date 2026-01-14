#' Generate Beta Diversity Change Boxplot for Longitudinal Data
#'
#' Creates boxplots comparing within-group and between-group beta diversity
#' across multiple time points. Useful for visualizing temporal changes in
#' community composition differences.
#'
#' @name generate_beta_change_boxplot_long
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param ... Additional parameters passed on to ggsave()
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object.
#'
#' @return A named list of ggplot objects visualizing within-group and between-group beta diversity across time points.
#' The list contains one plot for each distance metric specified in \code{dist.name}.
#' Each plot shows boxplots of the within-group and between-group distances at each time point,
#' with boxes dodge-positioned by the \code{group_var} variable if provided.
#' A line plot connecting the medians of boxes is overlaid to show the trajectory over time.
#'
#' @examples
#' \dontrun{
#' # Load required libraries and example data
#' library(vegan)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC")
#' generate_beta_change_boxplot_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   t0.level = "1",
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
#' generate_beta_change_boxplot_long(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   t0.level = "Baseline",
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
#' data("subset_T2D.obj")
#' generate_beta_change_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_gender",
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
#' generate_beta_change_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
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
#' data("ecam.obj")
#' generate_beta_change_boxplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "diet",
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
#' generate_beta_change_boxplot_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   group.var = "diet",
#'   strata.var = "delivery",
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
generate_beta_change_boxplot_long <-
  function(data.obj = NULL,
           dist.obj = NULL,
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

    # Exit the function if no distance metric is specified
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    # This step ensures we have the necessary distance matrices for the analysis
    if (is.null(dist.obj)) {
      # Process time variable to ensure proper ordering in longitudinal analysis
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(subject.var, time.var, group.var, strata.var, adj.vars)))
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Adjust distances for covariates if specified
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # Process time variable if distance object is provided
      data.obj <-
        mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
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

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # Set the theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each specified distance metric
    plot_list <- lapply(dist.name, function(dist.name) {

      # Convert distance object to a tibble for easier manipulation
      dist_tibble <- as_tibble(as.matrix(dist.obj[[dist.name]]), rownames = "Sample1")

      # Join metadata with distance information
      if (!is.null(strata.var)){
        dist_meta <- dist_tibble %>%
          tidyr::pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var), Strata = all_of(strata.var)), by = c("Sample1" = "sample")) %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var), Strata = all_of(strata.var)), by = c("Sample2" = "sample"))
      } else {
        dist_meta <- dist_tibble %>%
          tidyr::pivot_longer(cols = -Sample1, names_to = "Sample2", values_to = "Distance") %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var)), by = c("Sample1" = "sample")) %>%
          dplyr::left_join(meta_tab %>% select(sample, Group = all_of(group.var), Time = all_of(time.var)), by = c("Sample2" = "sample"))
      }

      # Handle different plotting scenarios based on the number of groups
      if (length(unique(dist_meta$Group.x)) <= 2) {
        # For two or fewer groups, calculate within- and between-group distances
        if (!is.null(strata.var)){
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = dplyr::if_else(Group.x == Group.y, "Within", "Between"),
              Group = dplyr::if_else(Type == "Within", Group.x, paste(Group.x, Group.y, sep = "-")),
              Time = Time.x
            ) %>%
            select(Group, Distance, Type, Time, Strata = Strata.x)
        } else {
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = dplyr::if_else(Group.x == Group.y, "Within", "Between"),
              Group = dplyr::if_else(Type == "Within", Group.x, paste(Group.x, Group.y, sep = "-")),
              Time = Time.x
            ) %>%
            select(Group, Distance, Type, Time)
        }

        # Create boxplot for within- and between-group distances
        p <- ggplot(within_between_dist, aes(x = Time, y = Distance, fill = Type)) +
          geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.8),
            width = 0.2
          ) +
          stat_summary(
            fun = mean,
            geom = "point",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 2
          ) +
          stat_summary(
            fun = mean,
            geom = "line",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 1
          ) +
          scale_fill_manual(values = col) +
          scale_color_manual(values = col) +
          labs(x = time.var, y = "Distance") +
          theme_to_use +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                plot.title = element_text(hjust = 0.5))

      } else {
        # For more than two groups, calculate pairwise distances between all groups
        if (!is.null(strata.var)){
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = paste(Group.x, Group.y, sep = "-"),
              Time = Time.x
            ) %>%
            select(Type, Distance, Time, Strata = Strata.x) %>%
            dplyr::mutate(Type = apply(select(., Type), 1, function(x) paste(sort(unlist(strsplit(x, "-"))), collapse = "-"))) %>%
            distinct(Type, Time, Distance, Strata, .keep_all = TRUE)
        } else {
          within_between_dist <- dist_meta %>%
            filter(Time.x == Time.y) %>%
            dplyr::mutate(
              Type = paste(Group.x, Group.y, sep = "-"),
              Time = Time.x
            ) %>%
            select(Type, Distance, Time) %>%
            dplyr::mutate(Type = apply(select(., Type), 1, function(x) paste(sort(unlist(strsplit(x, "-"))), collapse = "-"))) %>%
            distinct(Type, Time, Distance, .keep_all = TRUE)
        }

        # Create boxplot for pairwise distances between all groups
        p <- ggplot(within_between_dist, aes(x = Time, y = Distance, fill = Type)) +
          geom_boxplot(alpha = 0.7, position = position_dodge(0.8)) +
          stat_boxplot(
            geom = "errorbar",
            position = position_dodge(width = 0.8),
            width = 0.2
          ) +
          stat_summary(
            fun = mean,
            geom = "point",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 3
          ) +
          stat_summary(
            fun = mean,
            geom = "line",
            position = position_dodge(0.8),
            aes(group = Type, color = Type),
            size = 1
          ) +
          scale_fill_manual(values = col) +
          scale_color_manual(values = col) +
          labs(x = time.var, y = "Distance") +
          theme_to_use +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "right",
                plot.title = element_text(hjust = 0.5))
      }

      # Add faceting by strata if a strata variable is specified
      if (!is.null(strata.var)) {
        p <- p + facet_wrap2("Strata", scales = "free", nrow = length(unique(meta_tab[[strata.var]])))
      } else {
        p <- p
      }

      # Save the plot as a PDF if requested
      if (pdf) {
        pdf_name <- paste0("beta_change_boxplot_",
                           dist.name,
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "t0_level_", t0.level)

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