#' @title Taxa Change Test for Paired Data
#'
#' @description Computes taxa abundance changes between two time points and tests for
#' differential change between groups using linear models (lm) or ANOVA.
#'
#' @inheritParams mStat_data_obj_doc
#' @param change.base Character or numeric specifying the baseline time point.
#'   If NULL, uses the first unique value from time.var.
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param ref.level Character specifying the reference level for group comparisons.
#'   If NULL, uses the first level alphabetically.
#' @param winsor.qt Numeric (0-1) specifying the quantile for winsorization. Default 0.97.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' generate_taxa_change_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   feature.dat.type = "count"
#' )
#'
#' data(subset_pairs.obj)
#' generate_taxa_change_test_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   feature.dat.type = "count"
#' )
#' }
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested comparisons between groups
#'         \itemize{
#'           \item Elements named as "Level vs Reference (Reference)"
#'           \item If \code{group.var} has >2 levels, includes ANOVA results
#'         }
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Effect size of the change between time points
#'                 (interpretation depends on \code{feature.change.func}:
#'                 for "log fold change", represents difference in log2 abundances;
#'                 for "relative change", represents relative difference;
#'                 for "absolute change", represents absolute difference)
#'           \item \code{SE}: Standard error of the coefficient from the linear model
#'           \item \code{P.Value}: Raw p-value from standard linear model (lm)
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value using Benjamini-Hochberg method
#'           \item \code{Mean.Abundance}: Mean abundance of the feature across all samples
#'           \item \code{Prevalence}: Proportion of samples where the feature is present (non-zero)
#'         }
#' }
#'
#' This function analyzes CHANGE SCORES (differences between two time points) rather than
#' raw abundances, using standard linear models rather than LinDA mixed-effects models.
#' @export
#' @name generate_taxa_change_test_pair
generate_taxa_change_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var,
           ref.level = NULL,
           adj.vars = NULL,
           change.base = NULL,
           feature.change.func = "relative change",
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion", "other"),
           winsor.qt = 0.97) {
    # Validate the input data object
    data.obj <- mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Extract relevant columns from the metadata
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    mStat_validate_group_var_contract(
      meta.dat = meta_tab,
      group.var = group.var,
      subject.var = subject.var,
      context = "taxa change testing"
    )

    # Set reference level for group variable
    if (!is.null(group.var)) {
      # Convert to factor (this function always treats group.var as categorical)
      meta_tab[[group.var]] <- as.factor(meta_tab[[group.var]])

      # Get available levels
      available_levels <- levels(meta_tab[[group.var]])

      # Validate and set reference level
      if (!is.null(ref.level)) {
        if (!(ref.level %in% available_levels)) {
          stop(
            "ref.level '", ref.level, "' not found in group.var '", group.var, "'. ",
            "Available levels: ", paste(available_levels, collapse = ", ")
          )
        }
        meta_tab[[group.var]] <- relevel(meta_tab[[group.var]], ref = ref.level)
        message("Reference level for '", group.var, "': ", ref.level)
      } else {
        message(
          "Reference level for '", group.var, "': ", available_levels[1],
          " (alphabetically first)"
        )
      }
    }

    # Check for time-varying covariates
    if (!is.null(adj.vars)){
      # Identify time-varying variables in the metadata
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # If time-varying variables are found, stop the analysis
      if (length(time_varying_info$time_varying_vars) > 0) {
        stop("Feature-level analysis does not yet support adjustment for time-varying variables. Found time-varying variables: ",
             paste(time_varying_info$time_varying_vars, collapse = ", "),
             ". Future versions will support this feature.")
      }
    }

    # Extract levels of the grouping variable
    group_level <-
      meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

    # Set the reference level for the grouping variable
    reference_level <- group_level[1]

    formula <- mStat_build_formula(
      response = "value",
      terms = c(group.var, adj.vars)
    )

    # Adjust filters for 'other' data type
    if (feature.dat.type == "other") {
      prev.filter <- 0
      abund.filter <- 0
    }

    # Normalize count data if necessary.
    data.obj <- mStat_normalize_count_data_if_needed(data.obj, feature.dat.type)

    # Perform analysis for each feature level
    test.list <- lapply(feature.level, function(feature.level) {
      # Aggregate data by taxonomy if necessary
      otu_tax_agg_filter <- get_taxa_data(data.obj, feature.level, prev.filter, abund.filter)

      # Perform data imputation and winsorization for count or proportion data
      if (feature.dat.type %in% c("count", "proportion")) {
        # Calculate half of the minimum non-zero value for each taxon across all samples and time points
        # This ensures the same pseudocount is used for each taxon at both time points,
        # eliminating time-point specific bias in log fold change calculations
        half_nonzero_min <- apply(otu_tax_agg_filter[, -1], 1, function(x) {
          nonzero_values <- x[x > 0]
          if (length(nonzero_values) > 0) {
            min(nonzero_values) / 2
          } else {
            1e-10  # Fallback for taxa with all zeros
          }
        })

        # Create a logical matrix identifying zero values
        zero_matrix <- otu_tax_agg_filter[, -1] == 0

        # Create a matrix of imputation values (per-taxon pseudocount applied across all samples)
        half_nonzero_min_matrix <- matrix(half_nonzero_min, nrow = nrow(zero_matrix),
                                          ncol = ncol(zero_matrix), byrow = FALSE)

        # Impute zero values with half of the minimum non-zero value
        otu_tax_agg_filter[, -1][zero_matrix] <- half_nonzero_min_matrix[zero_matrix]

        message("Zero-handling: Per-taxon half-minimum pseudocount calculated across ALL samples and time points combined. ",
                "This ensures unbiased change calculations by using the same pseudocount at both time points.")

        # Apply winsorization to limit extreme values
        otu_tax_agg_filter[, -1] <- apply(otu_tax_agg_filter[, -1], 2, function(x) {
          qt <- quantile(x, probs = c((1 - winsor.qt) / 2, 1 - (1 - winsor.qt) / 2))
          x[x < qt[1]] <- qt[1]
          x[x > qt[2]] <- qt[2]
          return(x)
        })

      }
      # "other": no preprocessing — respect user's pre-processed data

      merged_data <- mStat_prepare_taxa_long_data(
        feature.dat = otu_tax_agg_filter,
        feature.level = feature.level,
        value_col = "value",
        meta.dat = meta_tab,
        feature_in_column = TRUE,
        join = "inner"
      )

      pair_change <- mStat_prepare_taxa_pair_change_data(
        long.df = merged_data,
        feature.level = feature.level,
        subject.var = subject.var,
        time.var = time.var,
        change.base = change.base,
        feature.change.func = feature.change.func,
        context = "taxa change testing"
      )
      combined_data <- pair_change$combined_data

      message("Note: For repeated measurements of the same subject at the same time point, the average will be taken.")

      # Create a matrix of value differences
      value_diff_matrix <- combined_data %>%
        select(all_of(feature.level), !!sym(subject.var), value_diff) %>%
        dplyr::group_by(!!sym(feature.level), !!sym(subject.var)) %>%
        dplyr::summarise(value_diff = mean(value_diff, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = !!sym(subject.var), values_from = value_diff) %>%
        tibble::column_to_rownames(var = feature.level) %>%
        as.matrix()

      # Prepare metadata for analysis
      aligned_subjects <- mStat_align_subject_metadata_to_matrix(
        value_matrix = value_diff_matrix[, colnames(na.omit(value_diff_matrix)), drop = FALSE],
        meta_tab = meta_tab,
        subject.var = subject.var,
        keep_vars = c(group.var, adj.vars)
      )
      value_diff_matrix <- aligned_subjects$value_matrix
      subject_meta <- aligned_subjects$sorted_meta

      # Calculate average abundance and prevalence for each feature
      prop_prev_data <- mStat_summarize_taxa_features(
        feature.dat = otu_tax_agg_filter,
        feature.level = feature.level,
        feature_in_column = TRUE
      )

      # Convert the value difference matrix to long format
      value_diff_long <- value_diff_matrix %>%
        as.data.frame() %>%
        tibble::rownames_to_column(feature.level) %>%
        tidyr::pivot_longer(
          cols = -all_of(feature.level),
          names_to = subject.var,
          values_to = "value"
        )

      # Perform statistical tests for each feature
      sub_test.list <-
        lapply(value_diff_long %>% select(all_of(feature.level)) %>% pull() %>% unique(), function(taxon) {
          # Prepare data for the current feature
          test_df <- mStat_attach_subject_level_metadata(
            df = value_diff_long %>% dplyr::filter(!!sym(feature.level) == taxon),
            meta.dat = subject_meta,
            subject.var = subject.var,
            vars = c(group.var, adj.vars)
          )

          # Fit the linear model
          test_result <- lm(formula, data = test_df)

          # Extract coefficients from the linear model
          coef.tab <- extract_coef(test_result)

          # Perform ANOVA if the grouping variable has more than two levels
          if (length(unique(stats::na.omit(test_df[[group.var]]))) > 2) {
            anova <- anova(test_result)
            anova.tab <- anova %>%
              as.data.frame() %>%
              tibble::rownames_to_column("Term") %>%
              select(
                Term,
                Statistic = `F value`,
                P.Value = `Pr(>F)`
              ) %>%
              tibble::as_tibble() %>%
              dplyr::mutate(Estimate = NA, Std.Error = NA)

            # Reorder the columns to match coef.tab
            anova.tab <- anova.tab %>%
              select(
                Term,
                Estimate,
                Std.Error,
                Statistic,
                P.Value
              )

            # Combine the coefficient table and ANOVA table
            coef.tab <-
              rbind(coef.tab, anova.tab)
          }
          return(tibble::as_tibble(coef.tab))
        })

      # Assign names to the elements of sub_test.list
      names(sub_test.list) <- value_diff_long %>%
        select(all_of(feature.level)) %>%
        pull() %>%
        unique()

      # Extract unique terms related to the grouping variable
      group_prefix_pattern <- paste0("^", mStat_escape_regex(group.var))
      unique_terms <-
        grep(paste0(group_prefix_pattern, "$|", group_prefix_pattern, ".*"),
             unique(unlist(
               lapply(sub_test.list, function(df)
                 unique(df$Term))
             )),
             value = TRUE)

      # Compile results for each term
      result_list <- lapply(unique_terms, function(term) {
        do.call(rbind, lapply(sub_test.list, function(df) {
          df %>% dplyr::filter(Term == term)
        })) %>%
          dplyr::mutate(!!sym(feature.level) := names(sub_test.list)) %>%
          dplyr::left_join(prop_prev_data, by = feature.level) %>%
          dplyr::mutate(Adjusted.P.Value = p.adjust(P.Value, method = "fdr")) %>%
          dplyr::select(all_of(
            c(
              feature.level,
              "Estimate",
              "Std.Error",
              "P.Value",
              "Adjusted.P.Value",
              "avg_abundance",
              "prevalence"
            )
          )) %>%
          dplyr::rename(
            Coefficient = Estimate,
            SE = Std.Error,
            Variable = feature.level,
            Mean.Abundance = avg_abundance,
            Prevalence = prevalence
          )
      })

      # Assign names to the result list
      names(result_list) <- unique_terms

      # Modify result names to include reference level information
      new_names <- sapply(names(result_list), function(name) {

        if (grepl(group_prefix_pattern, name) &&
            !grepl(paste0(group_prefix_pattern, "$"), name)) {
          sub_name <- sub(group_prefix_pattern, "", name)
          return(paste(sub_name, "vs", reference_level, "(Reference)"))
        }
        return(name)
      })

      names(result_list) <- new_names

      return(result_list)

    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)

  }
