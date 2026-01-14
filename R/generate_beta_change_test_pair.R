#' Test Beta Diversity Change Between Time Points
#'
#' Tests within-subject beta diversity changes between two time points using
#' linear models. Supports group comparisons and covariate adjustment.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#'   If NULL, uses the first time point.
#'
#' @examples
#' \dontrun{
#'
#' # Load packages
#' library(vegan)
#'
#' # Load data
#' data(peerj32.obj)
#' generate_beta_change_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = NULL,
#'   change.base = "1",
#'   dist.name = c('BC', 'Jaccard')
#' )
#' generate_beta_change_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#' data("subset_pairs.obj")
#' generate_beta_change_test_pair(
#'   data.obj = subset_pairs.obj,
#'   dist.obj = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   dist.name = c('BC', 'Jaccard')
#' )
#' }
#' @return A named list containing linear modeling results for each beta diversity metric.

#' Each list element corresponds to one of the distance metrics specified in \code{dist.name}.

#' It contains a coefficient table from fitting a linear model with the beta diversity change as
#' response and the \code{group_var} and \code{adj_vars} as predictors.

#' If \code{group_var} has multiple levels, ANOVA results are also included after the coefficients.

#' Column names include:
#' Term, Estimate, Std.Error, Statistic, P.Value
#' @export
#' @name generate_beta_change_test_pair
generate_beta_change_test_pair <-
  function(data.obj,
           dist.obj = NULL,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars = NULL,
           change.base = NULL,
           dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

    # Check if dist.name is provided, if not, return early
    if (is.null(dist.name)){
      return()
    }

    # Calculate beta diversity if not provided
    if (is.null(dist.obj) & !is.null(data.obj)) {
      # Calculate beta diversity using specified distance metrics
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      # Extract relevant metadata
      meta_vars <- c(subject.var, group.var, time.var)
      if (!is.null(adj.vars)) {
        meta_vars <- c(meta_vars, adj.vars)
      }
      meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(meta_vars)) %>% rownames_to_column("sample")
    } else {
      # Extract metadata from data.obj or dist.obj
      if (!is.null(data.obj)){
        meta_vars <- c(subject.var, group.var, time.var)
        if (!is.null(adj.vars)) {
          meta_vars <- c(meta_vars, adj.vars)
        }
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(meta_vars)) %>% rownames_to_column("sample")
      } else {
        meta_vars <- c(subject.var, group.var, time.var)
        if (!is.null(adj.vars)) {
          meta_vars <- c(meta_vars, adj.vars)
        }
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels")  %>% dplyr::select(all_of(meta_vars)) %>% rownames_to_column("sample")
      }
    }

    # Initialize variable to store time-varying information
    time_varying_info <- NULL

    # Handle time-varying covariates
    if (!is.null(adj.vars)){
      # Identify time-varying and non-time-varying variables
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # Adjust distances for non-time-varying variables
      if (length(time_varying_info$non_time_varying_vars) > 0){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj, dist.obj, time_varying_info$non_time_varying_vars, dist.name)
      }
    }

    # Set default change.base if not provided
    if (is.null(change.base)){
      change.base <- unique(meta_tab %>% dplyr::select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'meta.dat' data frame: ", change.base)
    }

    # Identify time points after the baseline
    change.after <-
      unique(meta_tab %>% dplyr::select(all_of(c(time.var))))[unique(meta_tab %>% dplyr::select(all_of(c(time.var)))) != change.base]

    # Perform analysis for each distance metric
    test.list <- lapply(dist.name, function(dist.name){

      # Convert distance matrix to long format
      dist.df <- as.matrix(dist.obj[[dist.name]]) %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      # Prepare data for analysis
      long.df <- dist.df %>%
        tidyr::gather(key = "sample2", value = "distance", -sample) %>%
        dplyr::left_join(meta_tab, by = "sample") %>%
        dplyr::left_join(meta_tab, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        # Filter for within-subject comparisons
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        dplyr::group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        # Filter for baseline and follow-up time points
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        dplyr::ungroup() %>%
        dplyr::select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      # Join with metadata
      long.df <- long.df %>% dplyr::left_join(meta_tab %>% dplyr::select(-any_of(c(time.var, "sample"))) %>% dplyr::distinct(), by = subject.var)

      # CRITICAL FIX: Check group levels AFTER data processing and BEFORE fitting model
      predictors <- c(time_varying_info$time_varying_vars, group.var)
      predictors <- predictors[!is.null(predictors)]
      
      # Check if group variable still has multiple levels after filtering
      if (!is.null(group.var) && group.var %in% names(long.df)) {
        remaining_group_levels <- unique(long.df[[group.var]])
        remaining_group_levels <- remaining_group_levels[!is.na(remaining_group_levels)]
        
        if (length(remaining_group_levels) < 2) {
          warning("After data filtering, group variable '", group.var, "' has only ", 
                  length(remaining_group_levels), " level(s): ", 
                  paste(remaining_group_levels, collapse = ", "), 
                  ". Removing from model and proceeding with coefficient estimates only.")
          # Remove group variable from predictors
          predictors <- predictors[predictors != group.var]
        }
      }
      
      # Create formula for linear model
      if (length(predictors) > 0) {
        formula <- stats::as.formula(paste0("distance", "~", paste(predictors, collapse = "+")))
      } else {
        formula <- stats::as.formula("distance ~ 1")  # Intercept-only model
      }

      # Fit linear model
      lm.model <- lm(formula, data = long.df)

      # Extract model summary
      summary <- summary(lm.model)

      # Create coefficient table
      coef.tab <- summary$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("Term") %>%
        dplyr::select(
                Term,
                Estimate,
                Std.Error = `Std. Error`,
                Statistic = `t value`,
                P.Value = `Pr(>|t|)`) %>%
        as_tibble()

      # CRITICAL FIX: Only perform ANOVA if group variable is in the final model
      # Check if group variable is actually included in the fitted model
      if (!is.null(group.var) && group.var %in% attr(lm.model$terms, "term.labels")) {
        # Group variable is in the model, safe to perform ANOVA
        anova <- anova(lm.model)
        # Create ANOVA table
        anova.tab <- anova %>%
          as.data.frame() %>%
          rownames_to_column("Term") %>%
          dplyr::select(Term,
                        Statistic = `F value`,
                        P.Value = `Pr(>F)`) %>%
          dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
          as_tibble() %>%
          dplyr::select(
            Term,
            Estimate,
            Std.Error,
            Statistic,
            P.Value
          )

        # Combine coefficient and ANOVA tables
        coef.tab <-
          rbind(coef.tab, anova.tab)
      }

      return(coef.tab)
    })

    # Assign names to the elements of test.list
    names(test.list) <- dist.name

    return(test.list)

  }
