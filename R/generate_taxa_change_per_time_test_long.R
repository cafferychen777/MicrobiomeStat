#' Longitudinal Taxa Test in Microbiome Data
#'
#' This function performs a comprehensive analysis of microbiome data, focusing on the longitudinal trends of various taxa. It is specifically designed to work with data where the primary interest is to understand how the abundance of microbial taxa changes over time, across different groups, or under different conditions.
#'
#' Utilizing a standard linear model approach, the function is adept at identifying significant temporal variations in taxa abundance. It provides a robust framework for comparing microbial communities at different time points, thereby offering valuable insights into the dynamics of these communities over extended periods.
#'
#' The function's ability to handle both count and proportion data types, along with its features for adjusting covariates, makes it a versatile tool for microbiome research. It is particularly beneficial for studies investigating the effects of treatments, environmental changes, or other interventions on the microbiome's composition and behavior over time.
#'
#' @param data.obj A MicrobiomeStat data object containing microbiome data and metadata.
#' @param subject.var A string specifying the column name in meta.dat that uniquely identifies each subject.
#' @param time.var Optional; a string representing the time variable in the meta.dat. If provided, enables longitudinal analysis.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param group.var Optional; a string specifying the group variable in meta.dat for
#'                  between-group comparisons of change scores. This variable is ALWAYS
#'                  treated as categorical (factor) by the underlying `generate_taxa_change_test_pair`
#'                  function, regardless of input type. The first level (alphabetically)
#'                  is used as the reference group.
#' @param adj.vars Optional; a vector of strings representing covariates in meta.dat for adjustment in the analysis.
#' @param feature.level A string or vector of strings indicating the taxonomic level(s) for analysis (e.g., "Phylum", "Class").
#' @param prev.filter Numeric; a minimum prevalence threshold for taxa inclusion in the analysis.
#' @param abund.filter Numeric; a minimum abundance threshold for taxa inclusion in the analysis.
#' @param feature.dat.type Character; the type of feature data, passed to
#'                         \code{generate_taxa_change_test_pair}. Should be one of:
#' \itemize{
#'   \item "count": Raw count data. Will be TSS-normalized, zero-imputed using global
#'                  half-minimum pseudocount, then analyzed with standard lm()
#'   \item "proportion": Pre-normalized proportional data. Will be zero-imputed using
#'                       global half-minimum pseudocount, then analyzed with standard lm()
#'   \item "other": Pre-transformed data. Analyzed with standard lm() without normalization
#' }
#' Default is "count". All processing is done by \code{generate_taxa_change_test_pair}.
#' @param feature.change.func A function or character string specifying how to calculate
#' the change from baseline value. This allows flexible options:
#' - If a function is provided, it will be applied to each row to calculate change.
#'   The function should take 2 arguments: value at timepoint t and value at baseline t0.
#' - If a character string is provided, following options are supported:
#'   - 'relative change': (value_t - value_t0) / (value_t + value_t0)
#'   - 'absolute change': value_t - value_t0
#'   - 'log fold change': log2(value_t + 1e-5) - log2(value_t0 + 1e-5)
#' - Default is 'relative change'.
#'
#' If none of the above options are matched, an error will be thrown indicating
#' the acceptable options or prompting the user to provide a custom function.
#' @param ... Additional arguments passed to other methods.
#' @details
#' The function integrates various data manipulations, normalization procedures, and statistical tests to assess the significance of taxa changes over time or between groups. It allows for the adjustment of covariates and is capable of handling both count and proportion data types.
#'
#' The function uses a standard linear model (lm) to analyze the data. It handles fixed effects to account for the influence of different variables on the taxa. Filtering is performed based on prevalence and abundance thresholds, and normalization and aggregation procedures are applied as necessary.
#'
#' A key feature of the function is its ability to conduct differential abundance analysis separately for each time point in the longitudinal data. This method is particularly effective for identifying significant changes in taxa at specific time points, offering insights into the temporal dynamics of the microbiome.
#'
#' @return A nested list structure where:
#' \itemize{
#'   \item First level: Named by follow-up time points (\code{ts.levels})
#'   \item Second level: Named by \code{feature.level} (e.g., "Phylum", "Genus")
#'   \item Third level: Named by tested comparisons between groups
#'         (e.g., "Level vs Reference (Reference)")
#'   \item Each final element is a data.frame with the following columns:
#'         \itemize{
#'           \item \code{Variable}: Feature/taxon name
#'           \item \code{Coefficient}: Effect size of the change from baseline
#'                 (interpretation depends on \code{feature.change.func})
#'           \item \code{SE}: Standard error of the coefficient from the linear model
#'           \item \code{P.Value}: Raw p-value from standard linear model (lm)
#'           \item \code{Adjusted.P.Value}: FDR-adjusted p-value (Benjamini-Hochberg)
#'           \item \code{Mean.Abundance}: Mean abundance across all samples
#'           \item \code{Prevalence}: Proportion of samples where feature is present
#'         }
#' }
#'
#' This function is a wrapper that calls \code{generate_taxa_change_test_pair} separately
#' for each follow-up time point, analyzing changes from baseline at each time point.
#'
#' This function is especially useful for longitudinal microbiome studies, facilitating the exploration of temporal patterns in microbial communities. By analyzing different time points against a baseline, it helps to uncover significant temporal shifts in the abundance of various taxa.
#'
#' The function is tailored for investigations that aim to monitor changes in microbial communities over time, such as in response to treatments or environmental changes. The structured output assists in interpreting temporal trends and identifying key taxa that contribute to these changes.
#' @examples
#' \dontrun{
#' # Example1: Analyzing the Type 2 Diabetes dataset
#' data("subset_T2D.obj")
#' # Longitudinal analysis of microbial changes in different racial groups
#' result <- generate_taxa_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   adj.vars = "subject_gender",
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus", "Family"),
#'   feature.dat.type = "count"
#' )
#' # Visualizing the results for the Type 2 Diabetes dataset
#' dotplot_T2D <- generate_taxa_per_time_dotplot_long(
#'   data.obj = subset_T2D.obj,
#'   test.list = result,
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   time.var = "visit_number",
#'   feature.level = c("Genus", "Family")
#' )
#' result2 <- generate_taxa_change_per_time_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   t0.level = unique(subset_T2D.obj$meta.dat$visit_number)[1],
#'   ts.levels = unique(subset_T2D.obj$meta.dat$visit_number)[-1],
#'   group.var = "subject_race",
#'   adj.vars = NULL,
#'   prev.filter = 0.1,
#'   abund.filter = 0.001,
#'   feature.level = c("Genus", "Family"),
#'   feature.dat.type = "count"
#' )
#' }
#' @export
generate_taxa_change_per_time_test_long <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           adj.vars = NULL,
           feature.level,
           feature.change.func = "relative change",
           feature.dat.type = c("count", "proportion", "other"),
           prev.filter = 0.001,
           abund.filter = 0.001,
           ...) {

    # Match the feature data type argument to ensure it's one of the allowed options
    feature.dat.type <- match.arg(feature.dat.type)

    # Validate the input data object to ensure it meets the required format
    mStat_validate_data(data.obj)

    # Check if the input variables are of the correct type
    if (!is.character(subject.var))
      stop("`subject.var` should be a character string.")
    if (!is.character(time.var))
      stop("`time.var` should be a character string.")
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")

    # Process the time variable to ensure it's in the correct format for analysis
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

    # Extract relevant metadata for the analysis
    meta_tab <- data.obj$meta.dat %>%
      as.data.frame() %>%
      select(all_of(c(subject.var, group.var, time.var, adj.vars)))

    # Check for time-varying covariates, which are not currently supported
    if (!is.null(adj.vars)){
      # Identify any time-varying variables among the adjustment variables
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # If time-varying variables are found, stop the analysis and inform the user
      if (length(time_varying_info$time_varying_vars) > 0) {
        stop("Feature-level analysis does not yet support adjustment for time-varying variables. Found time-varying variables: ",
             paste(time_varying_info$time_varying_vars, collapse = ", "),
             ". Future versions will support this feature.")
      }
    }

    # If baseline time point is not specified, determine it from the data
    if (is.null(t0.level)) {
      if (is.numeric(meta_tab[, time.var])) {
        t0.level <- sort(unique(meta_tab[, time.var]))[1]
      } else {
        t0.level <- levels(meta_tab[, time.var])[1]
      }
    }

    # If follow-up time points are not specified, determine them from the data
    if (is.null(ts.levels)) {
      if (is.numeric(meta_tab[, time.var])) {
        ts.levels <- sort(unique(meta_tab[, time.var]))[-1]
      } else {
        ts.levels <- levels(meta_tab[, time.var])[-1]
      }
    }

    # Perform the longitudinal analysis for each follow-up time point
    test_list <- lapply(ts.levels, function(ts.level){
        # Subset the data to include only the baseline and current follow-up time point
        subset.ids <- get_sample_ids(data.obj, time.var, c(t0.level, ts.level))
        subset_data.obj <- mStat_subset_data(data.obj, samIDs = subset.ids)

        # Perform the pairwise analysis between baseline and current follow-up time point
        subset.test.list <- generate_taxa_change_test_pair(data.obj = subset_data.obj,
                                       subject.var = subject.var,
                                       time.var = time.var,
                                       group.var = group.var,
                                       adj.vars = adj.vars,
                                       change.base = t0.level,
                                       feature.change.func = feature.change.func,
                                       feature.level = feature.level,
                                       prev.filter = prev.filter,
                                       abund.filter = abund.filter,
                                       feature.dat.type = feature.dat.type)

        return(subset.test.list)
      })

    # Name the results list with the follow-up time points for easy reference
    names(test_list) <- ts.levels

    # Return the complete list of test results for all follow-up time points
    return(test_list)
  }