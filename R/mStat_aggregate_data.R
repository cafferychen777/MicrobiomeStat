#' @title Aggregate Data Function
#' @description This function aggregates data by a specified subject and optional strata.
#' It provides flexible strategies for handling metadata conflicts in replicates.
#' @name mStat_aggregate_data
#' @param data.obj A list object in a format specific to MicrobiomeStat...
#' @param subject.var A character specifying the name of the subject variable.
#' @param strata.var A character specifying the name of the strata variable. If NULL, no stratification is performed.
#' @param meta.handle.conflict A character string specifying how to handle conflicting metadata in replicates.
#'   Options are:
#'   - `"first"` (default): Use the metadata from the first record encountered for each group. Issues a warning about conflicts.
#'   - `"stop"`: Stop execution with an error if any metadata conflict is found.
#'   - `"summarise"`: For numeric variables, compute the mean. For non-numeric (character/factor) variables, if they are inconsistent, stop with an error. Issues a warning for summarized numeric variables.
#'
#' @return A new list object with aggregated data.
#'
#'
#' @examples
#' \dontrun{
#' # Prepare data for the function
#' data(peerj32.obj)
#'
#' # Call the function with the default subject variable "subject"
#' aggregated_data <- mStat_aggregate_data(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   strata.var = NULL
#' )
#' 
#' # Example with a different subject variable name
#' # Let's pretend the subject ID column is called "participant"
#' # peerj32.obj$meta.dat$participant <- peerj32.obj$meta.dat$subject
#' # aggregated_data_2 <- mStat_aggregate_data(
#' #   data.obj = peerj32.obj,
#' #   subject.var = "participant",
#' #   strata.var = "group"
#' # )
#' }
#' @export
mStat_aggregate_data <- function(data.obj,
                                  subject.var,
                                  strata.var = NULL,
                                  meta.handle.conflict = c("first", "stop", "summarise")) {

  # --- 0. Argument and Data Validation ---
  meta.handle.conflict <- match.arg(meta.handle.conflict)
  message(paste("Aggregating data. Metadata conflict handling strategy:", meta.handle.conflict))

  grouping_vars <- subject.var
  if (!is.null(strata.var)) {
    grouping_vars <- c(subject.var, strata.var)
  }

  if (!all(grouping_vars %in% names(data.obj$meta.dat))) {
    stop("One or more specified 'subject.var' or 'strata.var' not found in data.obj$meta.dat.")
  }
  if (!"feature.tab" %in% names(data.obj)) {
    stop("Input data object must contain 'feature.tab'.")
  }
  if (!identical(rownames(data.obj$meta.dat), colnames(data.obj$feature.tab))) {
    stop("Data integrity issue: rownames of meta.dat do not match colnames of feature.tab.")
  }

  # --- 1. Detect Metadata Conflicts ---
  # This check is performed regardless of the strategy.
  meta_consistency_check <- data.obj$meta.dat %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~dplyr::n_distinct(.)), .groups = "drop")

  # Find all groups that have at least one inconsistent variable
  inconsistent_groups <- meta_consistency_check %>%
    dplyr::filter(dplyr::if_any(-dplyr::all_of(grouping_vars), ~ . > 1))

  has_conflicts <- nrow(inconsistent_groups) > 0

  # --- 2. Handle Conflicts Based on Selected Strategy ---
  if (has_conflicts) {
    # Build a detailed message about conflicts
    conflict_message <- "Inconsistent metadata found in replicates. Please check the following:\n"
    for (i in 1:nrow(inconsistent_groups)) {
      group_info <- paste(inconsistent_groups[i, grouping_vars], collapse = "_")
      inconsistent_vars <- names(inconsistent_groups)[which(inconsistent_groups[i, ] > 1)]
      inconsistent_vars <- setdiff(inconsistent_vars, grouping_vars)
      if (length(inconsistent_vars) > 0) {
        conflict_message <- paste0(conflict_message, "- Group '", group_info, "': inconsistent variable(s) '", paste(inconsistent_vars, collapse = ", "), "'\n")
      }
    }

    if (meta.handle.conflict == "stop") {
      stop(conflict_message, call. = FALSE)
    } else {
      # For both "first" and "summarise", we issue a warning.
      warning(paste("Strategy:", meta.handle.conflict, "\n", conflict_message), call. = FALSE)
    }
  }

  # --- 3. Aggregate Feature Table (Main and List) ---
  # This logic is now safe because we've handled the conflict reporting.
  data.obj.new <- list()

  # Function to aggregate a single feature table
  aggregate_features <- function(abund_matrix) {
    abund_df <- as.data.frame(t(abund_matrix))
    original_feature_names <- colnames(abund_df)

    meta_subset <- data.obj$meta.dat[, grouping_vars, drop = FALSE]
    agg_data <- cbind(meta_subset, abund_df)

    agg_feature_df <- agg_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::summarise(dplyr::across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

    if (is.null(strata.var)) {
      new_colnames <- agg_feature_df[[subject.var]]
    } else {
      new_colnames <- paste(agg_feature_df[[subject.var]], agg_feature_df[[strata.var]], sep = "_")
    }

    feature.tab <- t(as.matrix(agg_feature_df[, -which(names(agg_feature_df) %in% grouping_vars)]))
    colnames(feature.tab) <- new_colnames
    rownames(feature.tab) <- original_feature_names
    return(feature.tab)
  }

  data.obj.new$feature.tab <- aggregate_features(data.obj$feature.tab)

  if (!is.null(data.obj$feature.agg.list)) {
    data.obj.new$feature.agg.list <- lapply(data.obj$feature.agg.list, aggregate_features)
  }

  # --- 4. Aggregate Metadata Based on Strategy ---
  if (meta.handle.conflict == "summarise") {
    # Strategy: Summarise numerics, check for conflicts in non-numerics
    meta.dat <- data.obj$meta.dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::summarise(
        # For numeric columns, calculate the mean.
        dplyr::across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
        # For non-numeric columns, check for consistency. If not consistent, stop. Otherwise, take the first.
        dplyr::across(where(function(x) !is.numeric(x)), ~{
          if (dplyr::n_distinct(.) > 1) {
            stop(paste("Cannot summarise non-numeric variable '", dplyr::cur_column(), "' for group '", paste(dplyr::cur_group(), collapse="_"), "' as it has conflicting values: ", paste(unique(.), collapse=", ")), call. = FALSE)
          }
          .[1]
        }),
        .groups = "drop"
      )
  } else {
    # Strategy: "first" (or "stop", which would have already exited)
    meta.dat <- data.obj$meta.dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), ~ .[1]), .groups = "drop")
  }

  # --- 5. Finalize the New Data Object ---
  if (is.null(strata.var)) {
    new_rownames <- meta.dat[[subject.var]]
  } else {
    new_rownames <- paste(meta.dat[[subject.var]], meta.dat[[strata.var]], sep = "_")
  }

  meta.dat <- as.data.frame(meta.dat)
  rownames(meta.dat) <- new_rownames
  data.obj.new$meta.dat <- meta.dat

  data.obj.new$tree <- data.obj$tree
  data.obj.new$feature.ann <- data.obj$feature.ann

  message("Finished!")
  return(data.obj.new)
}