#' Volatility Test on Principal Coordinates for Longitudinal Data
#'
#' Tests association between PC coordinate volatility and group variable
#' for longitudinal microbiome data.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param pc.ind Numeric vector specifying which PC axes to test. Default c(1, 2).
#' @param group.var Required. Character string specifying the grouping variable in metadata.
#' @param ... Additional arguments passed to internal functions.
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
                                             group.var,
                                             adj.vars = NULL,
                                             dist.name = c("BC"),
                                             ...) {

  # Check if distance metrics are provided
  if (is.null(dist.name)){
    return()
  }

  # If distance object is not provided, calculate it from the data object
  if (is.null(dist.obj)) {
    data.obj <- mStat_process_time_variable(data.obj, time.var)

    # Extract relevant metadata
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        subject.var, time.var, group.var, adj.vars
      )))
    # Calculate beta diversity
    dist.obj <-
      mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
    # If adjustment variables are provided, calculate adjusted distances
    if (!is.null(adj.vars)){
      dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
    }
  } else {
    prepared_context <- mStat_prepare_precomputed_beta_context(
      dist.obj = dist.obj,
      dist.name = dist.name,
      pc.obj = pc.obj,
      data.obj = data.obj,
      time.var = time.var,
      process_time = TRUE,
      required_pc_axes = max(pc.ind)
    )
    data.obj <- prepared_context$data.obj
    dist.obj <- prepared_context$dist.obj
    pc.obj <- prepared_context$pc.obj
    meta_tab <- mStat_extract_dist_metadata(
      dist.obj = dist.obj,
      dist.name = dist.name,
      vars = c(subject.var, time.var, group.var, adj.vars),
      data.obj = data.obj
    )
  }

  mStat_validate_group_var_contract(
    meta.dat = meta_tab,
    group.var = group.var,
    subject.var = subject.var,
    context = "beta PC volatility testing"
  )

  # Inform the user about the importance of numeric time variable
  message(
    "The volatility test in 'generate_beta_pc_volatility_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
    "The time variable will be processed within the function if needed."
  )

  # If principal component object is not provided, calculate it using MDS
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

  # Perform volatility test for each distance metric
  test.list <- lapply(dist.name, function(dist.name) {
    # Extract principal component coordinates
    pc.mat <- pc.obj[[dist.name]]$points

    df <- mStat_prepare_pc_long_data(
      pc.points = pc.mat,
      pc.ind = pc.ind,
      meta.dat = meta_tab,
      vars = c(subject.var, time.var, group.var, adj.vars),
      sample_col = "sample",
      join = "inner"
    )

    # Perform volatility test for each principal component
    sub_test.list <- lapply(unique(df$PC), function(pc.index) {
      sub_df <- df %>% filter(PC == pc.index)


      # Ensure time variable is numeric
      sub_df[[time.var]] <- mStat_coerce_time_to_numeric(
        sub_df[[time.var]],
        time.var = time.var,
        context = "beta PC volatility analysis"
      )

      # Calculate volatility for each subject
      # Volatility is defined as the mean of absolute differences in PC values divided by time differences
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

      # Join volatility data with group information
      test_df <- volatility_df %>%
        dplyr::left_join(meta_tab %>%
                    select(all_of(c(subject.var, group.var))) %>%
                    dplyr::distinct(),
                  by = subject.var,
                  relationship = "many-to-one")

      model_terms <- c()
      if (!is.null(group.var)) {
        model_terms <- c(model_terms, group.var)
      }
      if (!is.null(adj.vars)) {
        model_terms <- c(model_terms, adj.vars)
      }
      model_rhs <- if (length(model_terms) > 0) {
        paste(model_terms, collapse = " + ")
      } else {
        "1"
      }

      test_result <- lm(as.formula(paste("volatility ~", model_rhs)), data = test_df)

      coef.tab <- extract_coef(test_result)

      if (!is.null(group.var) && length(unique(test_df[[group.var]])) > 1) {
        anova_result <- anova(test_result)
        group_row <- as.data.frame(anova_result)
        terms <- rownames(group_row)
        group_index <- which(terms == group.var)
        if (length(group_index) > 0) {
          selected <- group_row[group_index[[1]], , drop = FALSE]
          coef.tab <- rbind(
            coef.tab,
            data.frame(
              Term = group.var,
              Estimate = NA_real_,
              Std.Error = NA_real_,
              Statistic = selected$`F value`,
              P.Value = selected$`Pr(>F)`,
              check.names = FALSE
            )
          )
        }
      }

      return(tibble::as_tibble(coef.tab))
    })

    names(sub_test.list) <- unique(df$PC)
    return(sub_test.list)
  })

  names(test.list) <- dist.name

  return(test.list)
}
