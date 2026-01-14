#' Trend Test on Principal Coordinates for Longitudinal Data
#'
#' Performs linear trend tests on PC axes of beta diversity over time using
#' mixed effects models.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param pc.obj A list containing dimension reduction results from
#'   \code{\link{mStat_calculate_PC}}. If NULL, PCoA is performed automatically.
#' @param pc.ind Numeric vector specifying which PC axes to test. Default c(1, 2).
#' @param ... Additional arguments passed to lmer().
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
#' data("ecam.obj")
#' generate_beta_pc_trend_test_long(
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
#'
#' data("subset_T2D.obj")
#' generate_beta_pc_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   pc.ind = c(1, 2),
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   group.var = "subject_race",
#'   adj.vars = "subject_gender",
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

  # Check if distance metrics are provided
  if (is.null(dist.name)){
    return()
  }

  # If distance object is not provided, calculate it from the data object
  if (is.null(dist.obj)) {
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
    # If data object is provided with metadata, extract relevant information
    if (!is.null(data.obj) & !is.null(data.obj$meta.dat)) {
      meta_tab <-
        data.obj$meta.dat %>% select(all_of(c(
          subject.var, time.var, group.var, adj.vars
        )))
    } else {
      # If no data object, extract metadata from distance object
      meta_tab <-
        attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(
          subject.var, time.var, group.var, adj.vars
        )))
      dist.obj <- mStat_subset_dist(dist.obj, colnames(meta_tab))
    }
  }

  # Inform the user about the importance of numeric time variable
  message(
    "The trend test in 'generate_beta_trend_test_long' relies on a numeric time variable.\n",
    "Please ensure that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the trend test.\n",
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

  # Perform trend test for each distance metric
  test.list <- lapply(dist.name, function(dist.name) {
    # Extract principal component coordinates
    pc.mat <- pc.obj[[dist.name]]$points

    colnames(pc.mat) <- paste0("PC", 1:ncol(pc.mat))

    pc.mat <- pc.mat %>% as_tibble()

    # Combine PC coordinates with metadata
    df <-
      cbind(pc.mat[, paste0("PC", pc.ind)], meta_tab[, c(subject.var, time.var, group.var)])

    # Reshape data from wide to long format
    df <-
      df %>%
      as_tibble() %>%
      tidyr::gather(key = "PC",
                    value = "value",
                    -all_of(subject.var, group.var, time.var))

    # Perform trend test for each principal component
    sub_test.list <- lapply(unique(df$PC), function(pc.index) {
      sub_df <- df %>% filter(PC == pc.index)

      # Create formula for mixed effects model
      formula <-
        create_mixed_effects_formula(
          response.var = "value",
          time.var = time.var,
          group.var = group.var,
          subject.var = subject.var
        )

      # Ensure time variable is numeric
      sub_df <-
        sub_df %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

      # Fit linear mixed effects model
      model <- lmer(formula, data = sub_df, ...)

      # Check if group variable has more than two categories
      if (length(unique(sub_df[[group.var]])) > 2) {
        # Perform Type III ANOVA for multi-category group variable
        anova_result <- anova(model, type = "III")

        # Extract coefficients from the model
        coef.tab <- extract_coef(model)
        
        # Extract the last row of ANOVA results (overall effect of group variable)
        last_row <- utils::tail(anova_result, 1)
        var_name <- rownames(last_row)[1]

        # Adjust the last row to match the format of coefficient table
        adjusted_last_row <- data.frame(
          Term = var_name,
          Estimate = NA,
          Std.Error = NA,
          Statistic = last_row$`F value`,
          P.Value = last_row$`Pr(>F)`
        )

        # Combine coefficient table with adjusted last row
        coef.tab <- rbind(coef.tab, adjusted_last_row)

      } else {
        # For binary group variable, extract coefficients directly
        coef.tab <- extract_coef(model)
      }

      return(as_tibble(coef.tab))
    })

    names(sub_test.list) <- unique(df$PC)
    return(sub_test.list)
  })

  names(test.list) <- dist.name

  return(test.list)
}