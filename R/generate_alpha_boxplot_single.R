#' @title Generate Alpha Diversity Boxplot (Single Time Point)
#'
#' @description Generates boxplots of alpha diversity indices for cross-sectional
#'   analysis or a single time point, with optional grouping and stratification.
#'
#' @name generate_alpha_boxplot_single
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#'
#' @param t.level Character string specifying the time level/value to subset data to,
#'   if a time variable is provided. Default NULL does not subset data.
#'
#' @return A list of boxplots displaying the specified alpha diversity indices at the specified time point dplyr::across different groupings, stratified by the specified stratification variable (if provided). Each boxplot in the list corresponds to one of the alpha diversity indices specified in `alpha.name`. The boxplots will be saved as PDF files if `pdf` is set to `TRUE`.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(ggh4x)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # First example with peerj32.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("simpson"),
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' alpha.obj <- mStat_calculate_alpha_diversity(peerj32.obj$feature.tab, "simpson")
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = alpha.obj,
#'   alpha.name   = c("simpson"),
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' # Load another dataset
#' data("subset_T2D.obj")
#'
#' # Second example with subset_T2D.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = subset_T2D.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("shannon"),
#'   time.var     = "visit_number",
#'   t.level      = "   3",
#'   group.var    = "subject_race",
#'   strata.var   = "subject_gender",
#'   adj.vars     = "sample_body_site",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = NULL,
#'   pdf          = TRUE,
#'   file.ann     = NULL,
#'   pdf.wid      = 20,
#'   pdf.hei      = 8.5
#' )
#'
#' }
#' library(vegan)
#' library(ggh4x)
#'
#' # Load data
#' data(peerj32.obj)
#'
#' # First example with peerj32.obj
#' generate_alpha_boxplot_single(
#'   data.obj     = peerj32.obj,
#'   alpha.obj    = NULL,
#'   alpha.name   = c("simpson"),
#'   time.var     = "time",
#'   t.level      = "2",
#'   group.var    = "group",
#'   strata.var   = "sex",
#'   adj.vars     = "sex",
#'   base.size    = 16,
#'   theme.choice = "bw",
#'   palette      = "lancet",
#'   pdf          = FALSE,
#'   file.ann     = NULL,
#'   pdf.wid      = 11,
#'   pdf.hei      = 8.5
#' )
#'
#' @export
generate_alpha_boxplot_single <- function (data.obj,
                                           alpha.obj = NULL,
                                           alpha.name = c("shannon",
                                                          "observed_species"),
                                           depth = NULL,
                                           time.var = NULL,
                                           t.level = NULL,
                                           group.var = NULL,
                                           strata.var = NULL,
                                           adj.vars = NULL,
                                           base.size = 16,
                                           theme.choice = "bw",
                                           custom.theme = NULL,
                                           palette = NULL,
                                           pdf = TRUE,
                                           file.ann = NULL,
                                           pdf.wid = 11,
                                           pdf.hei = 8.5,
                                           ...) {

  # Check if alpha diversity indices are specified
  if (is.null(alpha.name)){
    return()
  }

  # Calculate alpha diversity if not provided
  if (is.null(alpha.obj)) {
    # Perform rarefaction if depth is specified
    if (!is.null(depth)) {
      message(
        "Detected that the 'depth' parameter is not NULL. Proceeding with rarefaction. Call 'mStat_rarefy_data' to rarefy the data!"
      )
      data.obj <- mStat_rarefy_data(data.obj, depth = depth)
    }

    # Subset data to specific time point if specified
    if (!is.null(time.var) & !is.null(t.level)){
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    # Extract feature table and calculate alpha diversity
    otu_tab <- data.obj$feature.tab

    # Extract tree if faith_pd is requested
    tree <- NULL
    if ("faith_pd" %in% alpha.name) {
      tree <- data.obj$tree
    }

    alpha.obj <-
      mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name, tree = tree)
  } else {
    # Validate that all requested alpha.name are present in alpha.obj
    available_indices <- names(alpha.obj)
    missing_indices <- alpha.name[!alpha.name %in% available_indices]

    if (length(missing_indices) > 0) {
      stop(
        "The following alpha diversity indices are not available in alpha.obj: ",
        paste(missing_indices, collapse = ", "),
        ". Available indices: ",
        paste(available_indices, collapse = ", "),
        call. = FALSE
      )
    }

    # Subset data to specific time point if specified
    if (!is.null(time.var) & !is.null(t.level)){
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
  }

  # Extract metadata
  meta_tab <- data.obj$meta.dat

  # Combine alpha diversity and metadata
  alpha_df <-
    dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
    dplyr::inner_join(meta_tab %>% rownames_to_column(var = "sample"),
                      by = c("sample"))

  # Create a default group if not specified
  if (is.null(group.var)) {
    alpha_df <- alpha_df %>% dplyr::mutate("ALL" = "ALL")
    group.var <- "ALL"
  }

  # Set up the theme for plotting
  theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

  # Set up the color palette
  col <- mStat_get_palette(palette)

  # Create a plot for each alpha diversity index
  plot_list <- lapply(alpha.name, function(index) {
    # Define aesthetic mapping
    # Note: group.var is guaranteed to be non-NULL (set to "ALL" if not provided)
    # Note: index is guaranteed to exist in alpha_df (validated earlier)
    aes_function <- aes(
      x = !!sym(group.var),
      y = !!sym(index),
      fill = !!sym(group.var)
    )

    # Adjust for covariates if specified
    if (!is.null(adj.vars)) {
      # Prepare data for adjustment
      data_subset <- alpha_df %>%
        dplyr::select(all_of(adj.vars)) %>%
        dplyr::mutate(dplyr::across(where(is.character) & !is.factor, factor))

      # Create model matrix
      # This step converts categorical variables into dummy variables
      M <- model.matrix(
        ~ 0 + .,
        data = data_subset,
        contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
      )

      # Center the covariates
      # Centering helps to make the intercept interpretable as the expected value
      # when all covariates are at their mean
      M_centered <- scale(M, scale = FALSE)

      # Fit regression model
      # This step performs a linear regression of the alpha diversity index
      # on the centered covariates
      fit <- lm(alpha_df[[index]] ~ M_centered)

      # Calculate adjusted alpha diversity
      # The adjusted value is the sum of the intercept (expected value when
      # all covariates are at their mean) and the residuals (unexplained variation)
      adjusted_value <- fit$coefficients[1] + residuals(fit)

      # Update alpha diversity values
      alpha_df[[index]] <- adjusted_value

      # Inform user about adjustment
      message(
        "Alpha diversity has been adjusted for the following covariates: ",
        paste(adj.vars, collapse = ", "),
        "."
      )
    }

    # Set up y-axis label
    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      y_label <-
        paste0(index, " index (adjusted by: ", covariates, ")")
    } else {
      y_label <- paste0(index, " index")
    }

    # Create the boxplot
    boxplot <- ggplot(alpha_df,
                      aes_function) +
      geom_jitter(width = 0.3,
                  alpha = 0.5,
                  size = 2) +
      stat_boxplot(geom = "errorbar",
                   position = position_dodge(width = 0.2),
                   width = 0.3) +
      geom_boxplot(
        position = position_dodge(width = 0.8),
        width = 0.3,
      ) +
      scale_fill_manual(values = col) +
      # Add faceting if strata variable is provided
      {
        if (!is.null(strata.var) & !is.null(group.var)) {
          ggh4x::facet_nested(
            as.formula(paste(". ~", strata.var, "+", group.var)),
            drop = T,
            scale = "free",
            space = "free"
          )
        } else {
          if (group.var != "ALL") {
            ggh4x::facet_nested(
              as.formula(paste(". ~", group.var)),
              drop = T,
              scale = "free",
              space = "free"
            )
          }
        }
      } +
      # Add labels and title
      labs(y = y_label,
           title = dplyr::if_else(
             !is.null(time.var) &
               !is.null(t.level),
             paste0(time.var, " = ", t.level),
             ""
           ))  +
      theme_to_use +
      # Customize theme elements
      theme(
        panel.spacing.x = unit(0, "cm"),
        panel.spacing.y = unit(0, "cm"),
        strip.text = element_text(size = base.size),
        axis.text.y = element_text(color = "black", size = base.size),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = base.size),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
        legend.text = ggplot2::element_text(size = 16),
        legend.title = ggplot2::element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 20)
      ) + {
        if (group.var == "ALL") {
          guides(fill = "none")
        }
      }

    return(boxplot)
  })

  # Save the plots as PDF files if requested
  if (pdf) {
    plot_list <- lapply(seq_along(plot_list), function(plot_index) {
      plot <- plot_list[[plot_index]]
      current_alpha_name <- alpha.name[plot_index]

      # Construct PDF file name
      pdf_name <- paste0(
        "alpha_boxplot_single_",
        current_alpha_name,
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

      # Save plot as PDF
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      print(plot)
      dev.off()

      return(plot)
    })

  }

  # Name the plots in the list
  names(plot_list) <- alpha.name

  return(plot_list)
}