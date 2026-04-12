#' @title Generate Alpha Diversity Change Boxplot (Paired)
#'
#' @description Generates boxplots comparing the change in alpha diversity indices
#'   between two time points (paired design), with optional grouping and stratification.
#'
#' @name generate_alpha_change_boxplot_pair
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param change.base The base time for calculating the change in alpha diversity.
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. Options include 'log fold change', 'absolute change',
#'   or a custom function taking two arguments (t, t0).
#' @param ... (Optional) Additional arguments to pass to the plotting function.
#'
#' @return A boxplot displaying the change in the specified alpha diversity index between two time points, stratified by the specified grouping and/or strata variables (if provided). The boxplot will be saved as a PDF if `pdf` is set to `TRUE`.
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#' # Example 1: Both group.var and strata.var are NULL
#' generate_alpha_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = "simpson",
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = NULL,
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   change.base = "1",
#'   alpha.change.func = "absolute change",
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = "no_groups",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' # Example 2: Only group.var is non-NULL
#' generate_alpha_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = "simpson",
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   change.base = "1",
#'   alpha.change.func = "log fold change",
#'   base.size = 16,
#'   theme.choice = "classic",
#'   palette = "npg",
#'   pdf = TRUE,
#'   file.ann = "group_only",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' # Example 3: Both group.var and strata.var are non-NULL
#' generate_alpha_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = "simpson",
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = NULL,
#'   change.base = "1",
#'   alpha.change.func = "log fold change",
#'   base.size = 16,
#'   theme.choice = "gray",
#'   palette = "aaas",
#'   pdf = TRUE,
#'   file.ann = "group_and_strata",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' # Example 4: Both group.var and adj.vars are non-NULL, strata.var is NULL
#' generate_alpha_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = "simpson",
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = NULL,
#'   adj.vars = c("sex"),
#'   change.base = "1",
#'   alpha.change.func = "absolute change",
#'   base.size = 16,
#'   theme.choice = "minimal",
#'   palette = "jama",
#'   pdf = TRUE,
#'   file.ann = "group_and_adj",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#' generate_alpha_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("simpson"),
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   alpha.change.func = "log fold change",
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_alpha_change_boxplot_pair <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = c("shannon", "observed_species"),
           depth = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           change.base = NULL,
           alpha.change.func = c("log fold change"),
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    if (is.null(alpha.name)){
      return()
    }

    prepared <- mStat_prepare_alpha_inputs(
      data.obj = data.obj,
      alpha.obj = alpha.obj,
      alpha.name = alpha.name,
      depth = depth
    )
    data.obj <- prepared$data.obj
    alpha.obj <- prepared$alpha.obj

    meta_tab <- mStat_prepare_alpha_meta_tab(
      data.obj = data.obj,
      vars = c(subject.var, group.var, time.var, strata.var, adj.vars)
    )

    alpha_df <- mStat_prepare_alpha_data(
      alpha.obj = alpha.obj,
      meta.dat = meta_tab,
      sample_col = "sample",
      join = "inner"
    )

    pair_change <- mStat_prepare_alpha_pair_change_data(
      alpha.df = alpha_df,
      alpha.name = alpha.name,
      subject.var = subject.var,
      time.var = time.var,
      change.base = change.base,
      change.func = alpha.change.func,
      context = "alpha change plotting",
      join_vars = subject.var
    )
    combined_alpha <- pair_change$combined_alpha
    change.base <- pair_change$change.base
    change.after <- pair_change$change.after

    combined_alpha <- mStat_attach_pair_metadata(
      df = combined_alpha,
      meta_tab = meta_tab,
      subject.var = subject.var,
      time.var = time.var,
      mode = "followup_time",
      change.after = change.after
    )

    placeholder_group <- mStat_ensure_group_placeholder(
      combined_alpha,
      group.var = group.var,
      value = "All",
      column_name = "group"
    )
    combined_alpha <- placeholder_group$df
    resolved_group_var <- placeholder_group$group.var

    col <- mStat_get_palette(palette)

    facet_formula <-
      if (!is.null(strata.var)) {
        paste(". ~", strata.var)
      } else {
        ". ~ 1"
      }

    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    ylab_label <- if (is.function(alpha.change.func)) {
      base_label <-
        paste0("Change from ", change.base, " (custom function)")
    } else {
      base_label <-
        paste0("Change from ", change.base, " (", alpha.change.func, ")")
    }

    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      ylab_label <-
        paste0(base_label, " (adjusted by: ", covariates, ")")
    } else {
      ylab_label <- base_label
    }

    plot_list <- lapply(alpha.name, function(index) {
      if (!is.null(adj.vars)) {
        data_subset <- combined_alpha %>%
          dplyr::select(all_of(adj.vars)) %>%
          dplyr::mutate(dplyr::across(where(~ is.character(.) & !is.factor(.)), factor))

        M <- model.matrix(
          ~ 0 + .,
          data = data_subset,
          contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
        )

        M_centered <- scale(M, scale = FALSE)

        fit <-
          lm(combined_alpha[[paste0(index, "_diff")]] ~ M_centered)

        adjusted_value <- fit$coefficients[1] + residuals(fit)

        combined_alpha[[paste0(index, "_diff")]] <- adjusted_value

        message(
          "Alpha diversity Change has been adjusted for the following covariates: ",
          paste(adj.vars, collapse = ", "),
          "."
        )
      }

      plot <-
        ggplot(combined_alpha, aes(
          x = !!sym(resolved_group_var),
          y = !!sym(paste0(index, "_diff")),
          fill = !!sym(resolved_group_var)
        )) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.3) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.3,
        ) +
        geom_jitter(width = 0.3,
                    alpha = 0.5,
                    size = 1.5) +
        scale_fill_manual(values = col) +
        ylab(ylab_label) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          strip.text.x = element_text(size = 15, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      if (is.null(group.var)) {
        plot <- plot +
          theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none"
          )
      }

      if (!is.null(strata.var)) {
        plot <- plot +
          facet_wrap(as.formula(facet_formula), scales = "fixed")
      }

      return(plot)
    })

    change_func_label <- if (is.function(alpha.change.func)) {
      "custom_function"
    } else {
      alpha.change.func
    }

    if (pdf) {
      lapply(seq_along(plot_list), function(i) {
        plot <- plot_list[[i]]
        alpha_index <- alpha.name[i]

        pdf_name <- paste0(
          "alpha_change_boxplot_pair_",
          alpha_index,
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "change_base_",
          change.base,
          "_",
          "change_func_",
          change_func_label
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
          plot = plot,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      })
    }

    names(plot_list) <- alpha.name

    return(plot_list)
  }
