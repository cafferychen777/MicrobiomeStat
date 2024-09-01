#' Plot Beta Diversity
#'
#' @description
#' This function, part of the MicrobiomeStat package, generates various types of plots
#' to visualize beta diversity measures in microbiome data. It supports different
#' visualization schemes for single time point, paired time points, and longitudinal data.
#'
#' @param dist.obj A distance matrix or object containing beta diversity measures.
#' @param pc.obj A matrix or data frame containing principal coordinates.
#' @param meta.dat A data frame containing metadata for the samples.
#' @param measure Character string specifying the beta diversity measure to plot.
#' @param group.var Character string specifying the grouping variable in meta.dat.
#' @param strata.var Character string specifying the stratification variable in meta.dat.
#' @param subject.var Character string specifying the subject variable in meta.dat.
#' @param time.var Character string specifying the time variable in meta.dat.
#' @param adj.vars Character vector specifying variables to adjust for. If specified, adjusted PCoA will be performed.
#' @param time.point.plot Character vector specifying time points to plot. The first point will be the reference.
#' @param is.plot.change Logical, whether to plot change from baseline. Only effective when there are multiple time points.
#' @param plot.type Character string specifying the plot type, one of c("PCoA", "boxplot", "spaghettiplot").
#' @param ... Additional arguments passed to the underlying plotting functions.
#'
#' @return A ggplot object representing the beta diversity plot.
#'
#' @details
#' This function provides a flexible framework for visualizing beta diversity in microbiome data.
#' It can handle single time point, paired time points, and longitudinal data. The function offers
#' different plot types including PCoA plots, box plots, and spaghetti plots, and can visualize
#' changes over time.
#'
#' The behavior of the function depends on the number of time points specified:
#' - For a single time point or when time variable is not specified, only PCoA plot is available.
#' - For two time points, both PCoA and change boxplot are available.
#' - For more than two time points, PCoA, change boxplot, and change spaghettiplot are available.
#'
#' When `is.plot.change` is TRUE and multiple time points are specified, the function will calculate
#' changes from the baseline (first time point).
#'
#' @examples
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, c("BC", "Jaccard"))
#' plot_beta_diversity(dist.obj = dist.obj,
#'                     meta.dat = peerj32.obj$meta.dat,
#'                     measure = c("BC", "Jaccard"),
#'                     group.var = "group",
#'                     time.var = "time",
#'                     time.point.plot = c("1"),
#'                     plot.type = "PCoA")
#'
#' plot_beta_diversity(dist.obj = dist.obj,
#'                     meta.dat = peerj32.obj$meta.dat,
#'                     measure = c("BC", "Jaccard"),
#'                     group.var = "group",
#'                     subject.var = "subject",
#'                     time.var = "time",
#'                     time.point.plot = c("1","2"),
#'                     plot.type = "PCoA",
#'                     is.plot.change = FALSE)
#'
#' plot_beta_diversity(dist.obj = dist.obj,
#'                     meta.dat = peerj32.obj$meta.dat,
#'                     measure = c("BC", "Jaccard"),
#'                     group.var = "group",
#'                     subject.var = "subject",
#'                     time.var = "time",
#'                     time.point.plot = c("1","2"),
#'                     plot.type = "boxplot",
#'                     is.plot.change = TRUE)
#'
#' data(ecam.obj)
#' dist.obj <- mStat_calculate_beta_diversity(ecam.obj, c("BC", "Jaccard"))
#' plot_beta_diversity(dist.obj = dist.obj,
#'                     meta.dat = ecam.obj$meta.dat,
#'                     measure = c("BC", "Jaccard"),
#'                     group.var = "delivery",
#'                     subject.var = "subject.id",
#'                     time.var = "month_num",
#'                     time.point.plot = unique(ecam.obj$meta.dat$month_num),
#'                     plot.type = "boxplot",
#'                     is.plot.change = TRUE)
#'
#' plot_beta_diversity(dist.obj = dist.obj,
#'                     meta.dat = ecam.obj$meta.dat,
#'                     measure = c("BC", "Jaccard"),
#'                     group.var = "delivery",
#'                     subject.var = "subject.id",
#'                     time.var = "month_num",
#'                     time.point.plot = unique(ecam.obj$meta.dat$month_num),
#'                     plot.type = "spaghettiplot",
#'                     is.plot.change = TRUE)
#'
#' plot_beta_diversity(dist.obj = dist.obj,
#'                     meta.dat = ecam.obj$meta.dat,
#'                     measure = c("BC", "Jaccard"),
#'                     group.var = "delivery",
#'                     subject.var = "subject.id",
#'                     time.var = "month_num",
#'                     time.point.plot = unique(ecam.obj$meta.dat$month_num),
#'                     plot.type = "PCoA",
#'                     is.plot.change = FALSE)
#'
#' @export
plot_beta_diversity <- function (dist.obj,
                                 pc.obj = NULL,
                                 meta.dat,
                                 measure,
                                 group.var = NULL,
                                 strata.var = NULL,
                                 subject.var = NULL,
                                 time.var = NULL,
                                 adj.vars = NULL, # If specified, aPCoA will be performed.
                                 time.point.plot, # The first pt will be the reference.
                                 is.plot.change = TRUE,
                                 plot.type = c("PcoA", "boxplot", "spaghettiplot"), # dependent on settings
                                 ...
) {

  data.obj <- list()

  data.obj$meta.dat <- meta.dat

  if (is.null(time.var) | length(time.point.plot) == 1) {
    # Plot single time point or all samples if the time variable is not specified
    if (plot.type == "PCoA") {
      p <- generate_beta_ordination_single(
        data.obj = data.obj,
        subject.var = subject.var,
        time.var = time.var,
        t.level = time.point.plot,
        group.var = group.var,
        adj.vars = adj.vars,
        strata.var = strata.var,
        dist.obj = dist.obj,
        dist.name = measure,
        pc.obj = pc.obj
      )
    }
  } else if (!is.null(time.var) & length(time.point.plot) == 2) {
    # Plot two time points and sample pair
    if (is.null(subject.var)) {
      message("Subject variable not specified!")
      return()
    }

    if (!is.plot.change){
      if (plot.type == "PCoA"){
        p <- generate_beta_ordination_pair(
          data.obj = data.obj,
          dist.obj = dist.obj,
          pc.obj = pc.obj,
          subject.var = subject.var,
          time.var = time.var,
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars,
          dist.name = measure
        )
      } else {
        message(paste(
          "Currently, we do not support",
          plot.type,
          "output for this scenario."
        ))
        return()
      }
    } else {
      if (plot.type == "boxplot"){
        p <- generate_beta_change_boxplot_pair(
          data.obj = data.obj,
          dist.obj = dist.obj,
          subject.var = subject.var,
          time.var = time.var,
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars,
          change.base = time.point.plot[1],
          dist.name = measure
        )
      } else {
        message(paste(
          "Currently, we do not support",
          plot.type,
          "output for this scenario."
        ))
        return()
      }
    }
  } else if (!is.null(time.var) & length(time.point.plot) > 2) {
    # Plot more than two time points, which are truly longitudinal
    if (is.null(subject.var)) {
      message("Subject variable not specified!")
      return()
    }

    if (is.plot.change){
      if (plot.type == "boxplot"){
        p <- generate_beta_change_boxplot_long(
          data.obj = data.obj,
          dist.obj = dist.obj,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars,
          dist.name = measure
        )
      } else if (plot.type == "spaghettiplot"){
        p <- generate_beta_change_spaghettiplot_long(
          data.obj = data.obj,
          dist.obj = dist.obj,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars,
          dist.name = measure
        )
      } else {
        message(
          paste(
            "Currently, we do not support",
            plot.type,
            "output for this scenario."
          )
        )
        return()
      }
    } else {
      if (plot.type == "PCoA"){
        p <- generate_beta_ordination_long(
          data.obj = data.obj,
          dist.obj = dist.obj,
          pc.obj = pc.obj,
          subject.var = subject.var,
          time.var = time.var,
          t0.level = time.point.plot[1],
          ts.levels = time.point.plot[-1],
          group.var = group.var,
          strata.var = strata.var,
          adj.vars = adj.vars,
          dist.name = measure
        )
      } else {
        message(
          paste(
            "Currently, we do not support",
            plot.type,
            "output for this scenario."
          )
        )
        return()
      }
    }
  }


  return(p)
}


