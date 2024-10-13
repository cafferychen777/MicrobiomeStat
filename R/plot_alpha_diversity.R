#' Plot Alpha Diversity
#'
#' This function, part of the MicrobiomeStat package, generates various types of plots
#' to visualize alpha diversity measures in microbiome data. It supports different
#' visualization schemes for single time point, paired time points, and longitudinal data.
#'
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param meta.dat A data frame containing metadata for the samples.
#' @param measure The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param group.var Character string specifying the grouping variable in meta.dat.
#' @param strata.var Character string specifying the stratification variable in meta.dat.
#' @param subject.var Character string specifying the subject variable in meta.dat.
#' @param time.var Character string specifying the time variable in meta.dat.
#' @param adj.vars Character vector specifying variables to adjust for. If specified, residuals will be used.
#' @param time.point.plot Character vector specifying time points to plot. The first point will be the reference.
#' @param is.plot.change Logical, whether to plot change from baseline. Only effective when there are multiple time points.
#' @param alpha.change.func Character string specifying the function to calculate alpha diversity change, default is "log fold change".
#' @param plot.type Character string specifying the plot type, one of c("boxplot", "spaghettiplot").
#' @param ... Additional arguments passed to the underlying plotting functions.
#'
#' @return A ggplot object representing the alpha diversity plot.
#'
#' @details
#' As a core component of the MicrobiomeStat package, this function leverages the package's
#' data structures and conventions to provide flexible and powerful visualizations of
#' alpha diversity in microbiome studies. It integrates seamlessly with other MicrobiomeStat
#' functions for comprehensive microbiome data analysis.
#'
#' This function provides a flexible framework for visualizing alpha diversity in microbiome data.
#' It can handle single time point, paired time points, and longitudinal data. The function offers
#' different plot types including box plots and spaghetti plots, and can visualize changes over time.
#'
#' This function provides a flexible framework for visualizing alpha diversity in microbiome data.
#' It can handle single time point, paired time points, and longitudinal data. The function offers
#' different plot types including box plots and spaghetti plots, and can visualize changes over time.
#'
#' The behavior of the function depends on the number of time points specified:
#' - For a single time point or when time variable is not specified, only boxplot is available.
#' - For two time points, both boxplot and spaghettiplot are available, with an option to plot change.
#' - For more than two time points, both boxplot and spaghettiplot are available for raw measures,
#'   but change plots are not supported.
#'
#' When `is.plot.change` is TRUE and multiple time points are specified, the function will calculate
#' changes from the baseline (first time point) using the method specified in `alpha.change.func`.
#'
#' @examples
#' data(peerj32.obj)
#' alpha.obj <- mStat_calculate_alpha_diversity(peerj32.obj$feature.tab,
#' alpha.name = c("shannon", "observed_species"))
#' plot_alpha_diversity(alpha.obj = alpha.obj,
#'                      meta.dat = peerj32.obj$meta.dat,
#'                      measure = c("shannon", "observed_species"),
#'                      group.var = "group",
#'                      strata.var = "sex",
#'                      time.var = "time",
#'                      time.point.plot = c("1"),
#'                      plot.type = "boxplot")
#'
#' plot_alpha_diversity(alpha.obj = alpha.obj,
#'                      meta.dat = peerj32.obj$meta.dat,
#'                      measure = c("shannon", "observed_species"),
#'                      group.var = "group",
#'                      subject.var = "subject",
#'                      strata.var = "sex",
#'                      time.var = "time",
#'                      time.point.plot = c("1", "2"),
#'                      is.plot.change = TRUE,
#'                      plot.type = "boxplot")
#'
#' plot_alpha_diversity(alpha.obj = alpha.obj,
#'                      meta.dat = peerj32.obj$meta.dat,
#'                      measure = c("shannon", "observed_species"),
#'                      group.var = "group",
#'                      subject.var = "subject",
#'                      strata.var = "sex",
#'                      time.var = "time",
#'                      time.point.plot = c("1", "2"),
#'                      is.plot.change = FALSE,
#'                      plot.type = "boxplot")
#'
#' plot_alpha_diversity(alpha.obj = alpha.obj,
#'                      meta.dat = peerj32.obj$meta.dat,
#'                      measure = c("shannon", "observed_species"),
#'                      group.var = "group",
#'                      subject.var = "subject",
#'                      strata.var = "sex",
#'                      time.var = "time",
#'                      time.point.plot = c("1", "2"),
#'                      is.plot.change = FALSE,
#'                      plot.type = "spaghettiplot")
#'
#' data(ecam.obj)
#' alpha.obj <- mStat_calculate_alpha_diversity(ecam.obj$feature.tab,
#' alpha.name = c("shannon", "observed_species"))
#' plot_alpha_diversity(alpha.obj = alpha.obj,
#'                      meta.dat = ecam.obj$meta.dat,
#'                      measure = c("shannon", "observed_species"),
#'                      group.var = "delivery",
#'                      subject.var = "antiexposedall",
#'                      strata.var = "diet",
#'                      time.var = "month_num",
#'                      time.point.plot = unique(ecam.obj$meta.dat$month_num),
#'                      is.plot.change = FALSE,
#'                      plot.type = "spaghettiplot")
#' plot_alpha_diversity(alpha.obj = alpha.obj,
#'                      meta.dat = ecam.obj$meta.dat,
#'                      measure = c("shannon", "observed_species"),
#'                      group.var = "delivery",
#'                      subject.var = "antiexposedall",
#'                      strata.var = "diet",
#'                      time.var = "month_num",
#'                      time.point.plot = unique(ecam.obj$meta.dat$month_num),
#'                      is.plot.change = FALSE,
#'                      plot.type = "boxplot")
#' @export
plot_alpha_diversity <- function (alpha.obj,
                                  meta.dat,
                                  measure,
                                  group.var = NULL,
                                  strata.var = NULL,
                                  subject.var = NULL,
                                  time.var = NULL,
                                  adj.vars = NULL,
                                  time.point.plot,
                                  is.plot.change = FALSE,
                                  alpha.change.func = "log fold change",
                                  plot.type = c("boxplot", "spaghettiplot"),
                                  ...) {

  # Validate and set the plot type
  # This step ensures that only supported plot types are used, enhancing function robustness
  plot.type <- match.arg(plot.type)

  # Create a standardized data object for analysis
  # This approach facilitates consistent data handling across different analysis functions
  data.obj <- list()
  data.obj$meta.dat <- meta.dat

  # The function adapts to various study designs based on time variables and number of time points
  # This flexibility allows for comprehensive analysis of different experimental setups

  if (is.null(time.var) | length(time.point.plot) == 1) {
    # Scenario 1: Cross-sectional data or single time point analysis
    # This scenario is applicable for studies comparing groups at a single time point or ignoring time effects

    if (plot.type == "boxplot"){
      # Generate a boxplot for single time point data
      # Boxplots provide a statistical summary of alpha diversity distribution:
      # - The box represents the interquartile range (IQR) containing the middle 50% of values
      # - The line in the box represents the median
      # - Whiskers typically extend to 1.5 times the IQR, representing the data's spread
      # - Points beyond the whiskers are potential outliers
      p <- generate_alpha_boxplot_single(
        data.obj = data.obj,
        alpha.obj = alpha.obj,
        alpha.name = measure,
        group.var = group.var,
        strata.var = strata.var,
        adj.vars = adj.vars,
        time.var = time.var,
        t.level = time.point.plot
      )
    } else {
      # Inform the user if the requested plot type is not supported for this scenario
      # This message helps users understand the function's capabilities and limitations
      message(paste(
        "Currently, we do not support",
        plot.type,
        "output for this scenario."
      ))
      return()
    }
  } else if (!is.null(time.var) & length(time.point.plot) == 2) {
    # Scenario 2: Paired data with two time points
    # This scenario is ideal for before-after studies or studies with baseline and follow-up measurements
    # It allows for the analysis of changes in alpha diversity over time

    # Ensure that the subject variable is specified for paired analysis
    # This is crucial for tracking individual changes over time
    if (is.null(subject.var)) {
      message("Subject variable not specified!")
      return()
    }

    if (is.plot.change) {
      # Generate a plot showing change in alpha diversity between two time points
      # This approach is useful for visualizing the magnitude and direction of changes over time
      if (plot.type == "boxplot"){
        # The change boxplot visualizes the distribution of alpha diversity changes
        # It allows for easy comparison of change magnitudes across different groups
        p <- generate_alpha_change_boxplot_pair(
          data.obj = data.obj,
          alpha.obj = alpha.obj,
          alpha.name = measure,
          subject.var = subject.var,
          time.var = time.var,
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars,
          change.base = time.point.plot[1],
          alpha.change.func = alpha.change.func,
        )
      } else {
        # Inform the user if the requested plot type is not supported for this scenario
        message(paste(
          "Currently, we do not support",
          plot.type,
          "output for this scenario."
        ))
        return()
      }
    } else {
      # Generate plots for two time points without explicitly showing change
      # This approach allows for comparison of alpha diversity at each time point
      if (plot.type == "boxplot") {
        # Boxplot for longitudinal data with two time points
        # This visualization allows for comparison of alpha diversity distributions between time points
        # It can reveal shifts in central tendency and spread of alpha diversity over time
        p <- generate_alpha_boxplot_long(
          data.obj = data.obj,
          alpha.obj = alpha.obj,
          alpha.name = measure,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars
        )
      } else if (plot.type == "spaghettiplot") {
        # Spaghetti plot for longitudinal data with two time points
        # This plot type visualizes individual trajectories of alpha diversity over time
        # It's particularly useful for identifying patterns of change at the individual level
        p <- generate_alpha_spaghettiplot_long(
          data.obj = data.obj,
          alpha.obj = alpha.obj,
          alpha.name = measure,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars
        )
      } else {
        # Inform the user if the requested plot type is not supported for this scenario
        message(paste(
          "Currently, we do not support",
          plot.type,
          "output for this scenario."
        ))
        return()
      }
    }
  } else if (!is.null(time.var) & length(time.point.plot) > 2) {
    # Scenario 3: Longitudinal data with more than two time points
    # This scenario is suitable for studies with multiple follow-up measurements
    # It allows for the analysis of alpha diversity trends over extended periods

    # Ensure that the subject variable is specified for longitudinal analysis
    # This is essential for tracking individual changes across multiple time points
    if (is.null(subject.var)) {
      message("Subject variable not specified!")
      return()
    }

    if (!is.plot.change) {
      # Generate plots for multiple time points
      # This approach allows for visualization of alpha diversity trends over time
      if (plot.type == "boxplot") {
        # Boxplot for longitudinal data with multiple time points
        # This visualization allows for comparison of alpha diversity distributions across all time points
        # It can reveal trends in central tendency and variability of alpha diversity over time
        p <- generate_alpha_boxplot_long(
          data.obj = data.obj,
          alpha.obj = alpha.obj,
          alpha.name = measure,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars
        )
      } else if (plot.type == "spaghettiplot") {
        # Spaghetti plot for longitudinal data with multiple time points
        # This plot type visualizes individual trajectories of alpha diversity over multiple time points
        # It's particularly useful for identifying patterns of change at the individual level and detecting potential subgroups with similar trajectories
        p <- generate_alpha_spaghettiplot_long(
          data.obj = data.obj,
          alpha.obj = alpha.obj,
          alpha.name = measure,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars
        )
      } else {
        # Inform the user if the requested plot type is not supported for this scenario
        message(paste(
          "Currently, we do not support",
          plot.type,
          "output for this scenario."
        ))
        return()
      }
    } else {
      # Inform the user that change plots are not supported for more than two time points
      # This limitation is due to the complexity of representing changes across multiple time points in a single plot
      message(paste(
        "Currently, we do not support",
        plot.type,
        "output for this scenario."
      ))
      return()
    }

  }

  # Return the generated plot
  # The returned plot object can be further customized or directly rendered
  return(p)

}