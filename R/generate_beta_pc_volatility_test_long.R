#' Generate Beta PC Volatility Test in Long Format
#'
#' This function computes the volatility test for principal coordinates (PC)
#' of beta diversity. It allows for a variety of distance measures and optional
#' adjustments. The function relies on a numeric time variable.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, containing components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). If dist.obj is provided, data.obj is not required.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param pc.obj A list containing the results of dimension reduction/Principal Component Analysis.
#' This should be the output from functions like \code{\link[MicrobiomeStat]{mStat_calculate_PC}},
#' containing the PC coordinates and other metadata. If NULL (default), dimension reduction
#' will be automatically performed using metric multidimensional scaling (MDS) via
#' \code{\link[MicrobiomeStat]{mStat_calculate_PC}}. The pc.obj list structure should contain:
#' \describe{
#'   \item{points}{A matrix with samples as rows and PCs as columns containing the coordinates.}
#'   \item{eig}{Eigenvalues for each PC dimension.}
#'   \item{vectors}{Loadings vectors for features onto each PC.}
#'   \item{Other metadata}{like method, dist.name, etc.}
#' }
#' See \code{\link[MicrobiomeStat]{mStat_calculate_PC}} function for details on output format.
#' @param pc.ind Numeric vector indicating which principal coordinate (PC) axes to use for
#'              volatility test, e.g. c(1,2) for PC1 and PC2. Defaults to c(1,2).
#' @param subject.var Character string specifying the column in metadata containing unique
#'                    subject IDs. This should uniquely identify each subject in the study.
#'                    Required to compute volatility within subjects.
#' @param time.var Character string specifying the column in metadata containing the time variable.
#'                This should be numeric with successive time points for each subject. Required to
#'                compute volatility over time.
#' @param group.var (Optional) Character string specifying the column in metadata containing a grouping
#'                  variable. This grouping variable will be used to test for differences in volatility
#'                  between groups. Can be left NULL.
#' @param adj.vars (Optional) Character vector specifying columns in metadata containing variables to
#'                 adjust for in the volatility models. Can be left NULL.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param ... (Optional) Additional arguments to pass to internal functions.
#'
#' @return A list of results for each distance measure and selected Principal Coordinate, including coefficients from the mixed-effects models.
#' @examples
#' \dontrun{
#' library(vegan)
#' data(ecam.obj)
#' generate_beta_pc_volatility_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "studyid",
#'   time.var = "month",
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   dist.name = c('BC')
#' )
#' }
#' @export
generate_beta_pc_volatility_test_long <- function(data.obj,
                                             dist.obj = NULL,
                                             pc.obj = NULL,
                                             pc.ind = c(1, 2),
                                             subject.var,
                                             time.var,
                                             group.var = NULL,
                                             adj.vars = NULL,
                                             dist.name = c("BC"),
                                             ...) {

  if (is.null(dist.name)){
    return()
  }

  if (is.null(dist.obj)) {
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
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
        data.obj$meta.dat %>% select(all_of(c(
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
    "The volatility test in 'generate_beta_pc_volatility_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
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

    df <- df %>%
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

      volatility_df <- sub_df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::arrange(!!sym(time.var)) %>%
        dplyr::mutate(
          value_diff = abs(value - dplyr::lag(value)),
          time_diff = !!sym(time.var) - dplyr::lag(!!sym(time.var))
        ) %>%
        dplyr::filter(!is.na(value_diff), !is.na(time_diff)) %>%
        dplyr::filter(time_diff != 0) %>%
        dplyr::summarise(
          volatility = mean(value_diff / time_diff, na.rm = TRUE)
        )

      test_df <- volatility_df %>%
        dplyr::left_join(meta_tab %>%
                    select(all_of(c(subject.var, group.var))) %>%
                    dplyr::distinct(),
                  by = subject.var)

      # Test the association between the volatility and the grp.var
      formula <- as.formula(paste("volatility ~", group.var))
      test_result <- lm(formula, data = test_df)

      coef.tab <- extract_coef(test_result)

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(unique(test_df[[group.var]])) > 1) {
        anova <- anova(test_result)
        anova.tab <- anova %>% as.data.frame() %>%
          rownames_to_column("Term") %>%
          dplyr::select(
                 Term,
                 Statistic = `F value`,
                 P.Value = `Pr(>F)`) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          as_tibble() %>%
          select(
            Term,
            Estimate,
            Std.Error,
            Statistic,
            P.Value
          )

        coef.tab <-
          rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
      }

      return(as_tibble(coef.tab))
    })

    names(sub_test.list) <- unique(df$PC)
    return(sub_test.list)
  })

  names(test.list) <- dist.name

  return(test.list)
}
