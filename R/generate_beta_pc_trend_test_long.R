#' Trend Test on Principal Coordinates of Beta Diversity Metrics Over Time
#'
#' This function performs a trend test on the selected Principal Coordinates (PCs) of beta diversity distances over time for different groups.
#' Various methods of dimension reduction, such as Principal Coordinates Analysis (PCoA), non-metric multidimensional scaling (NMDS), t-SNE, or UMAP, can be used.
#' The default method for dimension reduction is PCoA.
#' The function relies on a numeric time variable and can make adjustments for distance calculations.
#'
#' @param data.obj A list object containing all data.
#' @param dist.obj List object containing distance matrices. If NULL, beta diversity is calculated from data.obj.
#' @param pc.obj List object containing Principal Coordinates matrices. If NULL, PCoA is performed by default.
#' @param pc.ind Numeric vector indicating which Principal Coordinates to include in the trend test.
#' @param subject.var Character string indicating the variable for subject identifiers.
#' @param time.var Character string indicating the variable for time points. Must be coded as numeric.
#' @param group.var Character string indicating the variable for group identifiers.
#' @param adj.vars Character string indicating the variables for adjustments in the distance calculations.
#' @param dist.name Vector of character strings indicating the distance measures to include in the trend test, such as "BC" for Bray-Curtis.
#' @param ... Additional arguments to be passed to the mixed-effects model function.
#'
#' @return A list of results for each distance measure and selected Principal Coordinate, including coefficients from the mixed-effects models.
#' @examples
#' \dontrun{
#' library(vegan)
#'
#'   data(ecam.obj)
#'   generate_beta_pc_trend_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   adj.vars = "delivery",
#'   dist.name = c('BC'),
#' )
#' }
#' @export
generate_beta_pc_trend_test_long <- function(data.obj = NULL,
                                             dist.obj = NULL,
                                             pc.obj = NULL,
                                             pc.ind = c(1, 2),
                                             subject.var,
                                             time.var,
                                             group.var = NULL,
                                             adj.vars = NULL,
                                             dist.name = c("BC"),
                                             ...) {
  if (is.null(dist.obj)) {
    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(
        subject.var, time.var, group.var, adj.vars
      )))
    dist.obj <-
      mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
    if (!is.null(adj.vars)){
      dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
    }
  } else {
    if (!is.null(data.obj) & !is.null(data.obj$meta.dat)) {
      meta_tab <-
        load_data_obj_metadata(data.obj) %>% select(all_of(c(
          subject.var, time.var, group.var, adj.vars
        )))
    } else {
      meta_tab <-
        attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(
          subject.var, time.var, group.var, adj.vars
        )))
      dist.obj <- mStat_subset_dist(dist.obj, colnames(meta_tab))
    }
  }

  message(
    "The trend test in 'generate_beta_trend_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
    "The time variable will be processed within the function if needed."
  )

  if (is.null(pc.obj)) {
    message("No pc.obj provided, using MDS (PCoA) for dimension reduction by default.")
    message(
      "If you prefer other methods such as NMDS, t-SNE or UMAP, you can use the mStat_calculate_PC function with a specified method."
    )
    pc.obj <-
      mStat_calculate_PC(
        dist.obj = dist.obj,
        method = "mds",
        k = max(pc.ind),
        dist.name = dist.name
      )
  }

  test.list <- lapply(dist.name, function(dist.name) {
    pc.mat <- pc.obj[[dist.name]]$points

    colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

    pc.mat <- pc.mat %>% as_tibble()

    df <-
      cbind(pc.mat[, paste0("PC", pc.ind)], meta_tab[, c(subject.var, time.var, group.var)])

    df <-
      df %>%
      as_tibble() %>%
      tidyr::gather(key = "PC",
                    value = "value",
                    -one_of(subject.var, group.var, time.var))

    sub_test.list <- lapply(unique(df$PC), function(pc.index) {
      sub_df <- df %>% filter(PC == pc.index)

      formula <-
        create_mixed_effects_formula(
          response.var = "value",
          time.var = time.var,
          group.var = group.var,
          subject.var = subject.var
        )

      sub_df <-
        sub_df %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

      model <- lmer(formula, data = sub_df, ...)

      coef.tab <- extract_coef(model)

      return(as_tibble(coef.tab))
    })

    names(sub_test.list) <- unique(df$PC)
    return(sub_test.list)
  })

  names(test.list) <- dist.name

  return(test.list)
}
