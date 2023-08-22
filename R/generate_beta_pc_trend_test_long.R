#' Trend Test on Principal Coordinates of Beta Diversity Metrics Over Time
#'
#' Performs linear trend tests on selected Principal Coordinates (PCs) of beta diversity
#' distance matrices over time, for different groups. Allows using PCoA, NMDS, t-SNE,
#' UMAP for dimension reduction, with PCoA as default.
#'
#' @param data.obj MicrobiomeStat data object containing OTU table, taxonomy, metadata etc.
#' Required if dist.obj is NULL.
#' @param dist.obj List of pre-calculated beta diversity distance matrices. Calculated
#' from data.obj if NULL. Includes distance matrix for each measure in dist.name.
#' @param pc.obj List of pre-calculated PCs matrices. Calculated by PCoA if NULL.
#' Includes PCs matrix for each distance matrix in dist.obj.
#' @param pc.ind Vector of indices (positive integers) indicating which PCs to include
#' in trend test. Should not exceed number of PCs calculated.
#' @param subject.var Name of subject ID variable in metadata (character).
#' @param time.var Name of time variable in metadata (numeric). Values should be
#' ordered sequentially.
#' @param group.var Name of grouping variable in metadata (character). Levels will be
#' used to fit separate models. Optional.
#' @param adj.vars Names of variables in metadata to adjust distance matrices for
#' (character). Optional.
#' @param dist.name Names of beta diversity distance measures to include, e.g.
#' "BC" for Bray-Curtis (character vector).
#' @param ... Additional named arguments to pass to lmer():
#' \itemize{
#' \item control: Control parameters for lmer.
#' \item weights: Prior weights for the residuals.
#' }
#'
#' @return A nested list by distance > PC. Each element contains:
#' \itemize{
#' \item coef: Data frame of coefficients from mixed effects model.
#' \item model: Fitted lmer model object.
#' }
#'
#' @details
#' This function allows performing linear trend tests on PCs of beta diversity
#' distances over time, across groups, with adjustments.
#'
#' It checks for pre-calculated distances and PCs, generating them from data if needed.
#' Sufficient PCs should be calculated to cover indices specified.
#'
#' For each distance, PCs are extracted for the selected indices. The metadata is
#' joined and formatted for mixed effects modeling with time as numeric.
#'
#' The model formula is created using the response, time, group, and subject variables.
#' Mixed effects models are fitted with lmer() and coefficients extracted.
#'
#' @seealso
#' \code{\link{mStat_calculate_beta_diversity}} to generate distance matrices.
#'
#' \code{\link{mStat_calculate_PC}} to generate PCs from distances.
#'
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
#'   group.var = "diet",
#'   adj.vars = "delivery",
#'   dist.name = c('BC')
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
