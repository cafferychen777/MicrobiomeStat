#' Beta Diversity Volatility Test for Longitudinal Data
#'
#' Tests association between beta diversity volatility and group variable
#' for longitudinal microbiome data.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing the result of the volatility test for each specified beta diversity index.
#'   Each element in the list is a tibble with the coefficients extracted from the linear model fitted
#'   for each distance, and an ANOVA table if the group variable is multi-categorical.
#'
#' @details
#' The function starts by validating the input data, processing the time variable, and calculating the beta diversity if necessary.
#' Adjustments are made based on the provided adjusting variables. The volatility of the beta diversity is
#' computed for each subject, and linear models are fitted to test the association between volatility and
#' the specified group variable. The coefficients and ANOVA results are extracted and returned for each
#' beta diversity index specified.
#'
#' @note
#' A warning message will be displayed to ensure that the time variable is coded as numeric.
#' Non-numeric coding may lead to issues in the volatility test computation.
#'
#' @seealso
#' mStat_calculate_beta_diversity, mStat_calculate_adjusted_distance
#'
#' @examples
#' \dontrun{
#' data(ecam.obj)
#' generate_beta_volatility_test_long(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   dist.name = c("BC", "Jaccard")
#' )
#'
#' data(subset_T2D.obj)
#' generate_beta_volatility_test_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   dist.name = c("BC", "Jaccard")
#' )
#' }
#' @export
generate_beta_volatility_test_long <-
  function(data.obj,
           dist.obj = NULL,
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

    # Validate the input data object
    mStat_validate_data(data.obj)

    # Inform the user about the importance of numeric time variable
    message(
      "The volatility test in 'generate_beta_volatility_test_long' relies on a numeric time variable.\n",
      "Please ensure that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
      "The time variable will be processed within the function if needed."
    )

    # If distance object is not provided, calculate it from the data object
    if (is.null(dist.obj)) {
      # Extract relevant metadata
      meta_tab <- data.obj$meta.dat %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      # Calculate beta diversity
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # If adjustment variables are provided, calculate adjusted distances
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # If distance object is provided, extract metadata from the appropriate source
      if (!is.null(data.obj) & !is.null(data.obj$meta.dat)){
        meta_tab <- data.obj$meta.dat %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      } else {
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels") %>% select(all_of(c(subject.var, time.var, group.var, adj.vars)))
      }
    }

    # Ensure the distance object and metadata have matching dimensions
    if (nrow(as.matrix(dist.obj[[dist.name[1]]])) > nrow(meta_tab)){
      samIDs <- rownames(meta_tab)
      dist.obj <- mStat_subset_dist(dist.obj = dist.obj, samIDs = samIDs)
    }

    # Perform volatility test for each distance metric
    test.list <- lapply(dist.name,function(dist.name){

      # Convert distance matrix to long format
      dist.df <- as.matrix(dist.obj[[dist.name]])
      dist.df <- dist.df %>%
        as.data.frame() %>%
        rownames_to_column("sample")
      meta_tab <- meta_tab %>% rownames_to_column("sample")

      # Prepare data for volatility analysis
      # This step calculates the distance between consecutive time points for each subject
      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        dplyr::mutate(min_time_level = min(!!sym(paste0(time.var, ".subject")))[1]) %>%
        dplyr::arrange(!!sym(paste0(time.var, ".sample"))) %>%
        dplyr::mutate(prev_time_level = dplyr::lag(!!sym(paste0(time.var, ".subject")))) %>%
        filter(!!sym(paste0(time.var, ".sample")) == !!sym("prev_time_level")) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        filter(!!sym(paste0(time.var, ".subject")) != !!sym("min_time_level")) %>%
        dplyr::ungroup() %>%
        select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), !!sym(paste0(time.var, ".sample")) ,distance) %>%
        dplyr::rename(
          !!sym(subject.var) := !!sym(paste0(subject.var, ".subject")),
          !!sym(time.var) := !!sym(paste0(time.var, ".subject")),
          !!sym(paste0(time.var, ".before")) := !!sym(paste0(time.var, ".sample"))
        )

      # Ensure time variables are numeric
      long.df <- long.df %>%
        dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)),
               !!sym(paste0(time.var, ".before")) := as.numeric(!!sym(paste0(time.var, ".before"))))

      # Calculate volatility for each subject
      # Volatility is defined as the mean of distances divided by time differences
      volatility_df <- long.df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::arrange(!!sym(time.var)) %>%
        dplyr::mutate(
          time_diff = !!sym(time.var) - !!sym(paste0(time.var, ".before"))
        ) %>%
        dplyr::filter(!is.na(distance), !is.na(time_diff)) %>%
        dplyr::filter(time_diff != 0) %>%
        dplyr::summarise(
          volatility = mean(distance / time_diff, na.rm = TRUE)
        )

      # Join volatility data with group information
      test_df <- volatility_df %>%
        dplyr::left_join(meta_tab %>%
                           select(all_of(c(subject.var, group.var))) %>%
                           dplyr::distinct(), by = subject.var, relationship = "many-to-many")

      # Test the association between volatility and group variable using linear regression
      formula <- as.formula(paste("volatility ~", group.var))
      test_result <- lm(formula, data = test_df)

      # Extract coefficients from the linear model
      coef.tab <- extract_coef(test_result)

      # If group variable has more than one level, perform ANOVA
      if (length(unique(test_df[[group.var]])) > 1) {
        anova <- anova(test_result)
        # Format ANOVA results to match coefficient table structure
        anova.tab <- anova %>% as.data.frame() %>%
          rownames_to_column("Term") %>%
          select(
            Term,
            Statistic = `F value`,
            P.Value = `Pr(>F)`
          ) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          as_tibble() %>%
          select(
            Term,
            Estimate,
            Std.Error,
            Statistic,
            P.Value
          )

        # Combine coefficient table with ANOVA results
        coef.tab <-
          rbind(coef.tab, anova.tab)
      }

      return(as_tibble(coef.tab))
    })

    # Assign names to the elements of test.list
    names(test.list) <- dist.name

    return(test.list)
  }