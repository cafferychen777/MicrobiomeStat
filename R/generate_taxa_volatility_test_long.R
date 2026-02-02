#' Longitudinal Taxa Abundance Volatility Test
#'
#' Calculates and tests taxa abundance volatility (variability over time) between groups.
#' Volatility is the mean absolute difference between consecutive time points,
#' normalized by time difference.
#'
#' @inheritParams mStat_data_obj_doc
#'
#' @param transform Character; transformation method before volatility calculation.
#'   "CLR" applies CLR transform (default). Set to other value to skip transformation.
#' @param ... Additional arguments passed to other methods.
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Second level: Named by tested comparisons between groups
#'         (e.g., "Level vs Reference (Reference)")
#'   \item Each element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Effect size for volatility differences between groups
#'           \item \code{SE}: Standard error of the coefficient from the linear model
#'           \item \code{P.Value}: Raw p-value from standard linear model (lm)
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value (Benjamini-Hochberg)
#'           \item \code{Mean.Abundance}: Mean abundance across all samples
#'           \item \code{Prevalence}: Proportion of samples where feature is present (non-zero)
#'         }
#' }
#'
#' This function analyzes VOLATILITY (variability over time) using standard linear models,
#' NOT LinDA. Volatility is calculated as mean absolute differences between consecutive time points.
#'
#' @examples
#' \dontrun{
#' data("subset_T2D.obj")
#' test.list <- generate_taxa_volatility_test_long(
#' data.obj = subset_T2D.obj,
#' time.var = "visit_number",
#' subject.var = "subject_id",
#' group.var = "subject_race",
#' adj.vars = "sample_body_site",
#' prev.filter = 0.1,
#' abund.filter = 0.0001,
#' feature.level = c("Genus"),
#' feature.dat.type = "count",
#' transform = "CLR"
#' )
#' plot.list <- generate_taxa_volatility_volcano_long(data.obj = subset_T2D.obj,
#'                                                    group.var = "subject_race",
#'                                                    test.list = test.list,
#'                                                    feature.sig.level = 0.1,
#'                                                    feature.mt.method = "none")
#'
#' data("ecam.obj")
#' test.list <- generate_taxa_volatility_test_long(
#'   data.obj = ecam.obj,
#'   time.var = "month_num",
#'   subject.var = "subject.id",
#'   group.var = "antiexposedall",
#'   adj.vars = "delivery",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   feature.level = c("Order", "Family", "Genus"),
#'   feature.dat.type = "proportion",
#'   transform = "CLR"
#' )
#' plot.list <- generate_taxa_volatility_volcano_long(
#'   data.obj = ecam.obj,
#'   group.var = "antiexposedall",
#'   test.list = test.list,
#'   feature.sig.level = 0.2,
#'   feature.mt.method = "none"
#' )
#' }
#' @export
generate_taxa_volatility_test_long <- function(data.obj,
                                               time.var,
                                               subject.var,
                                               group.var,
                                               adj.vars = NULL,
                                               prev.filter = 0,
                                               abund.filter = 0,
                                               feature.level,
                                               feature.dat.type = c("count", "proportion", "other"),
                                               transform = "CLR",
                                               ...) {
  # Validate the input data object
  mStat_validate_data(data.obj)

  # Match the feature data type argument
  feature.dat.type <- match.arg(feature.dat.type)

  # Inform the user about the importance of numeric time variable
  message(
    "The volatility calculation relies on a numeric time variable.\n",
    "Please check that your time variable is coded as numeric.\n",
    "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
    "You can ensure the time variable is numeric by mutating it in the metadata."
  )

  # Convert time variable to numeric
  data.obj$meta.dat <-
    data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

  # Extract relevant variables from metadata
  meta_tab <- data.obj$meta.dat %>%
    dplyr::select(all_of(c(
      time.var, subject.var, group.var, adj.vars
    ))) %>% rownames_to_column("sample")

  # Extract group levels and set reference level
  group_level <- meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels
  reference_level <- group_level[1]

  # If CLR transformation is used, set abundance filter to 0
  if (transform == "CLR"){
    abund.filter <- 0
  }

  # Create a formula for the linear model
  formula_str <- paste("volatility ~", group.var)
  if (!is.null(adj.vars)) {
    formula_str <- paste(formula_str, "+", paste(adj.vars, collapse = " + "))
  }
  formula <- as.formula(formula_str)

  # Perform analysis for each taxonomic level
  test.list <- lapply(feature.level, function(feature.level) {

    # Normalize count data if necessary
    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    # Aggregate data to the specified taxonomic level if necessary
    otu_tax_agg <- get_taxa_data(data.obj, feature.level, feature.col = FALSE)

    # Calculate average abundance and prevalence for each feature
    prop_prev_data <-
      otu_tax_agg %>%
      as.matrix() %>%
      as.table() %>%
      as.data.frame() %>%
      dplyr::group_by(Var1) %>%
      dplyr::summarise(
        avg_abundance = mean(Freq),
        prevalence = sum(Freq > 0) / dplyr::n()
      ) %>% column_to_rownames("Var1") %>%
      rownames_to_column(feature.level)

    # Function for performing CLR transformation on data
    clr_transform <- function(x) {
      gm = exp(mean(log(x)))
      return(log(x / gm))
    }

    # Function to impute zeros with half the minimum non-zero value
    impute_zeros_rowwise <- function(row) {
      min_nonzero <- min(row[row > 0])
      row[row == 0] <- min_nonzero / 2
      return(row)
    }

    # Impute zeros and apply CLR transformation
    otu_tax_agg_imputed_temp <- apply(otu_tax_agg, 1, impute_zeros_rowwise)
    otu_tax_agg_imputed <- t(otu_tax_agg_imputed_temp)

    # Apply CLR transformation per SAMPLE (column-wise)
    # CLR must be applied within each sample, not across samples for each taxon
    # Each sample's CLR values should sum to 0 (compositional data property)
    otu_tax_agg_clr <- apply(otu_tax_agg_imputed, 2, clr_transform)

    # Convert to long format and apply filters
    otu_tax_agg_clr_long <- otu_tax_agg_clr %>%
      as.data.frame() %>%
      mStat_filter(prev.filter = prev.filter,
                   abund.filter = abund.filter) %>%
      rownames_to_column(feature.level) %>%
      tidyr::gather(key = "sample", value = "value",-feature.level)

    taxa.levels <- otu_tax_agg_clr_long %>% select(all_of(feature.level)) %>% pull() %>% unique()

    # Perform volatility analysis for each taxon
    sub_test.list <-
      lapply(taxa.levels, function(taxon) {

        taxa_df <- otu_tax_agg_clr_long %>%
          dplyr::filter(!!sym(feature.level) == taxon) %>%
          dplyr::left_join(meta_tab, by = "sample")

        # Calculate volatility for each subject
        volatility_df <- taxa_df %>%
          dplyr::arrange(!!sym(subject.var),!!sym(time.var)) %>%
          dplyr::group_by(!!sym(subject.var)) %>%
          dplyr::mutate(
            diff_value = abs(value - dplyr::lag(value)),
            diff_time = !!sym(time.var) - dplyr::lag(!!sym(time.var))
          ) %>%
          dplyr::filter(!is.na(diff_value),!is.na(diff_time)) %>%
          dplyr::filter(diff_time != 0) %>%
          dplyr::summarize(volatility = mean(diff_value / diff_time),
                    .groups = 'drop')

        test_df <- volatility_df %>%
          dplyr::left_join(meta_tab %>%
                             select(all_of(c(subject.var, group.var, adj.vars))) %>%
                             dplyr::distinct(),
                           by = subject.var)

        # Fit linear model
        test_result <- lm(formula, data = test_df)

        coef.tab <- extract_coef(test_result)

        # Perform ANOVA if group variable has more than two levels
        if (length(unique(taxa_df[[group.var]])) > 2) {
          anova <- anova(test_result)
          anova.tab <- anova %>% as.data.frame() %>%
            rownames_to_column("Term") %>%
            dplyr::select(
                   Term,
                   Statistic = `F value`,
                   P.Value = `Pr(>F)`) %>%
            as_tibble() %>%
            dplyr::mutate(Estimate = NA, Std.Error = NA) %>%
            select(
              Term,
              Estimate,
              Std.Error,
              Statistic,
              P.Value
            )

          coef.tab <- rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
        }
        return(as_tibble(coef.tab))
      })

    # Assign names to the elements of test.list
    names(sub_test.list) <- otu_tax_agg_clr_long %>% select(all_of(feature.level)) %>% pull() %>% unique()

    # Find all unique terms related to the group variable
    unique_terms <- grep(paste0("^", group.var, "$|^", group.var, ".*"), unique(unlist(lapply(sub_test.list, function(df) unique(df$Term)))), value = TRUE)

    # Extract and process results for each term
    result_list <- lapply(unique_terms, function(term) {
      do.call(rbind, lapply(sub_test.list, function(df) {
        df %>% dplyr::filter(Term == term)
      })) %>%
        dplyr::mutate(!!sym(feature.level) := names(sub_test.list)) %>%
        dplyr::left_join(prop_prev_data, by = feature.level) %>%
        dplyr::mutate(Adjusted.P.Value = p.adjust(P.Value, method = "fdr")) %>%
        dplyr::select(all_of(c(feature.level, "Estimate", "Std.Error", "P.Value", "Adjusted.P.Value", "avg_abundance", "prevalence"))) %>%
        dplyr::rename(Coefficient = Estimate,
                      SE = Std.Error,
                      Variable = feature.level,
                      Mean.Abundance = avg_abundance,
                      Prevalence = prevalence)
    })

    # Name the result list
    names(result_list) <- unique_terms

    # Rename the result list elements to include reference level information
    new_names <- sapply(names(result_list), function(name) {
      if (grepl(paste0("^", group.var), name) && !grepl(paste0("^", group.var, "$"), name)) {
        sub_name <- sub(paste0(group.var), "", name)
        return(paste(sub_name, "vs", reference_level, "(Reference)"))
      }
      return(name)
    })

    names(result_list) <- new_names

    return(result_list)
  })

  names(test.list) <- feature.level

  return(test.list)
}
