#' @title Generate Alpha Diversity Boxplot (Longitudinal)
#'
#' @description Generates boxplots of alpha diversity indices for longitudinal data,
#'   showing changes across time points with optional grouping and stratification.
#'
#' @name generate_alpha_boxplot_long
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @return A boxplot displaying the specified alpha diversity index dplyr::across different groupings and time points, stratified by the specified stratification variable (if provided). The boxplot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' generate_alpha_boxplot_long(
#'   data.obj = subset_T2D.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = c("sample_body_site"),
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#'
#' data("ecam.obj")
#' generate_alpha_boxplot_long(
#'   data.obj = ecam.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "subject.id",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "antiexposedall",
#'   strata.var = "diet",
#'   adj.vars = NULL,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#'
#' data(peerj32.obj)
#' generate_alpha_boxplot_long(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   t0.level = "1",
#'   ts.levels = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   base.size = 20,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5)
#'
#' data(subset_pairs.obj)
#' alpha.obj <- mStat_calculate_alpha_diversity(subset_pairs.obj$feature.tab,
#' c("shannon", "observed_species"))
#' generate_alpha_boxplot_long(
#'   data.obj = subset_pairs.obj,
#'   alpha.obj = alpha.obj,
#'   alpha.name = c("shannon", "observed_species"),
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   t0.level = "Baseline",
#'   ts.levels = NULL,
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 20,
#'   pdf.hei = 8.5)
#' }
#' @export
generate_alpha_boxplot_long <- function (data.obj,
                                         alpha.obj = NULL,
                                         alpha.name = c("shannon", "observed_species"),
                                         depth = NULL,
                                         subject.var,
                                         time.var,
                                         t0.level = NULL,
                                         ts.levels = NULL,
                                         group.var = NULL,
                                         strata.var = NULL,
                                         adj.vars = NULL,
                                         base.size = 12,
                                         theme.choice = "bw",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         ...) {
  # Check if alpha.name is provided
  if (is.null(alpha.name)){
    return()
  }

  # Check if alpha.obj is a list
  if (!is.null(alpha.obj) &&
      !is(alpha.obj, "list"))
    stop("`alpha.obj` should be a list or NULL.")
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

  # Calculate alpha diversity if not provided
  if (is.null(alpha.obj)) {
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
    if (!is.null(depth)) {
      message(
        "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <-
        mStat_rarefy_data(data.obj = data.obj, depth = depth)
    }
    otu_tab <- data.obj$feature.tab
    
    # Extract tree if faith_pd is requested
    tree <- NULL
    if ("faith_pd" %in% alpha.name) {
      tree <- data.obj$tree
    }
    
    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name, tree = tree)
  } else {
    # Verify that all alpha.name are present in alpha.obj
    if (!all(alpha.name %in% unlist(lapply(alpha.obj, function(x)
      colnames(x))))) {
      missing_alphas <- alpha.name[!alpha.name %in% names(alpha.obj)]
      stop(
        "The following alpha diversity indices are not available in alpha.obj: ",
        paste(missing_alphas, collapse = ", "),
        call. = FALSE
      )
    }
    data.obj <-
      mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)
  }

  # Prepare metadata and alpha diversity data
  meta_tab <- data.obj$meta.dat %>% 
    as.data.frame() %>% 
    dplyr::select(all_of(c(subject.var, group.var, time.var, strata.var, adj.vars)))

  time.levels <- meta_tab %>%
    dplyr::select(all_of(c(time.var))) %>%
    pull() %>%
    as.factor() %>%
    levels() %>%
    length()

  # Convert the alpha.obj list to a data frame
  alpha_df <-
    dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    dplyr::inner_join(meta_tab %>% rownames_to_column(var = "sample"),
                      by = c("sample"))

  # Set up theme and color palette
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Use mStat_get_palette to set the color palette
  col <- mStat_get_palette(palette)

  # Handle case when no grouping variable is provided
  if (is.null(group.var)) {
    alpha_df <- alpha_df %>% dplyr::mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  # Create a plot for each alpha diversity index
  plot_list <- lapply(alpha.name, function(index) {
    line_aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        group = !!sym(subject.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        group = !!sym(subject.var)
      )
    }

    aes_function <- if (!is.null(group.var)) {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        fill = !!sym(group.var)
      )
    } else {
      aes(
        x = !!sym(time.var),
        y = !!sym(index),
        fill = !!sym(time.var)
      )
    }

    # Adjust for covariates if specified
    if (!is.null(adj.vars)) {
      data_subset <- alpha_df %>%
        dplyr::select(all_of(adj.vars)) %>%
        dplyr::mutate(dplyr::across(where(is.character) &
                                      !is.factor, factor))

      M <-
        model.matrix(
          ~ 0 + .,
          data = data_subset,
          contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
        )

      # Center the covariates
      M_centered <- scale(M, scale = FALSE)

      # Fit regression model
      fit <- lm(alpha_df[[index]] ~ M_centered)

      # Compute the adjusted value
      adjusted_value <- fit$coefficients[1] + residuals(fit)

      # Update the alpha_df
      alpha_df[[index]] <- adjusted_value

      message(
        "Alpha diversity has been adjusted for the following covariates: ",
        paste(adj.vars, collapse = ", "),
        "."
      )
    }

    # Calculate average alpha diversity for large datasets
    average_alpha_df <- NULL
    if (length(unique(alpha_df[[time.var]])) > 10 ||
        length(unique(alpha_df[[subject.var]])) > 25) {
      if (!is.null(strata.var) & !is.null(group.var)) {
        average_alpha_df <- alpha_df %>%
          dplyr::group_by(!!sym(strata.var), !!sym(group.var), !!sym(time.var)) %>%
          dplyr::summarise(dplyr::across(!!sym(index), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL") %>%
          dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
      } else if (!is.null(group.var)) {
        average_alpha_df <- alpha_df %>%
          dplyr::group_by(!!sym(group.var), !!sym(time.var)) %>%
          dplyr::summarise(dplyr::across(!!sym(index), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL") %>%
          dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
      } else {
        average_alpha_df <- alpha_df %>%
          dplyr::group_by(!!sym(time.var)) %>%
          dplyr::summarise(dplyr::across(!!sym(index), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
          dplyr::ungroup() %>%
          dplyr::mutate(!!sym(subject.var) := "ALL") %>%
          dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))
      }
    }

    # Set up y-axis label
    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      y_label <- paste0(index, " index (adjusted by: ", covariates, ")")
    } else {
      y_label <- paste0(index, " index")
    }

    # Ensure time variable is a factor
    alpha_df <- alpha_df %>% dplyr::mutate(!!sym(time.var) := factor(!!sym(time.var)))

    boxplot <- ggplot(alpha_df,
                      aes_function) +
      #geom_violin(trim = FALSE, alpha = 0.8) +
      stat_boxplot(geom = "errorbar",
                   position = position_dodge(width = 0.2),
                   width = 0.3) +
      geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.3,
        #fill = "white"
      ) +
      geom_line(
        line_aes_function,
        alpha = 0.8,
        linewidth = 0.6,
        color = "black",
        linetype = "dashed",
        data = if (!is.null(average_alpha_df))
          average_alpha_df
        else
          alpha_df
      ) +
      scale_fill_manual(values = col) +
      {
        if (time.levels > 2) {
          if (!is.null(strata.var) & !is.null(group.var)) {
            ggh4x::facet_nested(
              cols = vars(!!sym(group.var)),
              rows = vars(!!sym(strata.var)),
              scale = "free",
              space = "free"
            )
          } else {
            if (group.var != "ALL") {
              ggh4x::facet_nested(
                cols = vars(!!sym(group.var)),
                scale = "free",
                space = "free"
              )
            }
          }
        } else {
          if (!is.null(strata.var) & !is.null(group.var)) {
            ggh4x::facet_nested(
              cols = vars(!!sym(strata.var), !!sym(group.var)),
              scale = "free",
              space = "free"
            )
          } else {
            if (group.var != "ALL") {
              ggh4x::facet_nested(
                cols = vars(!!sym(group.var)),
                scale = "free",
                space = "free"
              )
            }
          }
        }
      } +
      labs(x = time.var,
           y = y_label)  +
      theme_to_use +
      theme(
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0, "cm"),
        strip.text.x = element_text(size = base.size, color = "black"),
        strip.text.y = element_text(size = base.size, color = "black"),
        axis.text.x = element_text(
          angle = 90,
          color = "black",
          vjust = 0.5,
          size = base.size * 0.75
        ),
        axis.text.y = element_text(color = "black", size = base.size),
        axis.title.x = element_text(size = base.size),
        axis.title.y = element_text(size = base.size),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        legend.text = ggplot2::element_text(size = 16 * 2),
        legend.title = ggplot2::element_text(size = 16 * 2),
        legend.key.size = unit(10, "mm"),
        legend.key.spacing = unit(2, "mm")
      ) + {
        if (group.var == "ALL") {
          guides(fill = "none")
        }
      }

    # Add geom_jitter() if the number of unique time points or subjects is greater than 10
    if (length(unique(alpha_df[[time.var]])) > 10 ||
        length(unique(alpha_df[[subject.var]])) > 10) {
      boxplot <- boxplot + geom_jitter(width = 0.1,
                                       alpha = 0.1,
                                       size = 1)
    }

    return(boxplot)
  })

  # Save plots as PDF if requested
  if (pdf) {
    for (plot_index in seq_along(plot_list)) {
      plot <- plot_list[[plot_index]]
      current_alpha_name <- alpha.name[plot_index]

      pdf_name <- paste0(
        "alpha_boxplot_long_",
        current_alpha_name,
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

      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      print(plot)
      dev.off()
    }
  }

  # Name the plots in the list
  names(plot_list) <- alpha.name

  return(plot_list)
}
