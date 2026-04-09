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

    # This function generates boxplots to visualize changes in alpha diversity 
    # between two time points, with options for grouping and stratification.

    # Check if alpha diversity indices are specified. If not, exit the function.
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

    # Extract relevant metadata.
    # This step prepares the metadata for merging with the alpha diversity data.
    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
        subject.var, group.var, time.var, strata.var, adj.vars
      )))

    # Combine alpha diversity and metadata.
    # This creates a comprehensive dataset for our analysis.
    alpha_df <- mStat_prepare_alpha_data(
      alpha.obj = alpha.obj,
      meta.dat = meta_tab,
      sample_col = "sample",
      join = "inner"
    )

    pair_times <- mStat_resolve_pair_timepoints(
      values = alpha_df[[time.var]],
      time.var = time.var,
      change.base = change.base,
      context = "alpha change plotting"
    )
    change.base <- pair_times$change.base
    change.after <- pair_times$change.after

    # Split alpha diversity data by time points.
    # This separates the data into baseline and follow-up measurements.
    alpha_split <- split(alpha_df, f = as.character(alpha_df[[time.var]]))

    alpha_time_1 <- alpha_split[[change.base]]
    alpha_time_2 <- alpha_split[[change.after]]

    # Combine alpha diversity data from two time points.
    # This step pairs the baseline and follow-up measurements for each subject.
    combined_alpha <- alpha_time_1 %>%
      dplyr::inner_join(
        alpha_time_2,
        by = c(subject.var, group.var),
        suffix = c("_time_1", "_time_2")
      )

    # Calculate change in alpha diversity.
    # This is the core statistical operation of the function.
    diff_columns <- lapply(alpha.name, function(index) {
      diff_col_name <- paste0(index, "_diff")
      col_after  <- paste0(index, "_time_2")
      col_before <- paste0(index, "_time_1")

      combined_alpha %>%
        dplyr::mutate(!!diff_col_name := compute_alpha_change(
          value_after  = !!sym(col_after),
          value_before = !!sym(col_before),
          method       = alpha.change.func
        )) %>%
        dplyr::select(all_of(diff_col_name))
    })

    # Bind the calculated differences to the combined_alpha dataframe.
    # This step integrates the change calculations with the rest of our data.
    combined_alpha <- dplyr::bind_cols(combined_alpha, diff_columns)

    # Set up color palette for plotting.
    col <- mStat_get_palette(palette)

    # Set up faceting formula
    facet_formula <-
      if (!is.null(strata.var)) {
        paste(". ~", strata.var)
      } else {
        ". ~ 1"
      }

    # Handle grouping variable
    if (is.null(group.var)) {
      combined_alpha$group <- "All"
    } else{
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_df %>% dplyr::select(all_of(c(
          subject.var, group.var
        )))
        ,
        by = c(subject.var, group.var)) %>% dplyr::rename(group = group.var)
    }

    # Add strata variable if specified
    if (!is.null(strata.var)) {
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_time_1 %>% dplyr::select(all_of(c(
          subject.var, strata.var
        )))
        , by = c(subject.var))
    }

    # Add adjustment variables if specified
    if (!is.null(adj.vars) &&
        (is.null(strata.var) || strata.var != adj.vars)) {
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_time_1 %>% dplyr::select(all_of(c(
          subject.var, adj.vars
        )))
        , by = c(subject.var))
    }

    # Set up theme
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Set up y-axis label
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

    # Create plots for each alpha diversity index
    plot_list <- lapply(alpha.name, function(index) {
      # Adjust for covariates if specified
      if (!is.null(adj.vars)) {
        # Convert non-numerical covariates to factors
        data_subset <- combined_alpha %>%
          dplyr::select(all_of(adj.vars)) %>%
          dplyr::mutate(dplyr::across(where(~ is.character(.) & !is.factor(.)), factor))

        # Create a model matrix and set contrasts for non-numeric covariates
        M <- model.matrix(
          ~ 0 + .,
          data = data_subset,
          contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
        )

        # Center the covariates
        # Centering helps in interpreting the intercept as the expected value
        # when all covariates are at their mean
        M_centered <- scale(M, scale = FALSE)

        # Fit the regression model
        # This step performs a linear regression of the alpha diversity change
        # on the centered covariates
        fit <-
          lm(combined_alpha[[paste0(index, "_diff")]] ~ M_centered)

        # Calculate the adjusted alpha diversity change
        # The adjusted value is the sum of the intercept (expected value when
        # all covariates are at their mean) and the residuals (unexplained variation)
        adjusted_value <- fit$coefficients[1] + residuals(fit)

        # Update the alpha diversity change value in combined_alpha
        combined_alpha[[paste0(index, "_diff")]] <- adjusted_value

        # Display message indicating alpha diversity change has been adjusted for specific covariates
        message(
          "Alpha diversity Change has been adjusted for the following covariates: ",
          paste(adj.vars, collapse = ", "),
          "."
        )
      }

      # Create the plot
      # This step generates the boxplot visualization of alpha diversity changes
      plot <-
        ggplot(combined_alpha, aes(
          x = group,
          y = !!sym(paste0(index, "_diff")),
          fill = group
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

      # Adjust plot based on grouping
      # This step finalizes the plot appearance based on the data structure
      if (any(unique(combined_alpha$group) == "All")) {
        plot <- plot +
          theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none"
          )
      } else if (is.null(strata.var)) {
        # No additional modifications needed
      }
      else {
        plot <- plot +
          facet_wrap(as.formula(facet_formula), scales = "fixed")
      }

      return(plot)
    })

    # Set label for change function
    # This is used in file naming when saving plots
    change_func_label <- if (is.function(alpha.change.func)) {
      "custom_function"
    } else {
      alpha.change.func
    }

    # Save the plots as PDF files if requested
    if (pdf) {
      lapply(seq_along(plot_list), function(i) {
        plot <- plot_list[[i]]
        alpha_index <- alpha.name[i]

        # Construct PDF file name
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

        # Save plot as PDF
        ggsave(
          filename = pdf_name,
          plot = plot,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      })
    }

    # Name the plots in the list
    names(plot_list) <- alpha.name

    return(plot_list)
  }
