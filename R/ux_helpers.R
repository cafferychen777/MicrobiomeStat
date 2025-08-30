# User Experience Enhancement Functions for NewInterface.R
# User experience optimization helper functions

# Enhanced error message helper
mStat_user_friendly_error <- function(error_msg, function_name = "", param_name = "", suggested_fix = "") {
  # Create user-friendly error messages with suggestions
  
  # Common error patterns and solutions
  error_patterns <- list(
    "object not found" = list(
      pattern = "object.*not found|could not find",
      message = "Variable or data object not found",
      suggestion = "Please check if the variable name is spelled correctly and confirm the data is loaded properly"
    ),
    "column_missing" = list(
      pattern = "Column.*doesn't exist|subscript out of bounds",
      message = "Specified column/variable does not exist in the data",
      suggestion = "Please use names(data.obj$meta.dat) to view available variables"
    ),
    "time_level_missing" = list(
      pattern = "t0.level|ts.levels.*not found",
      message = "Specified time point does not exist in the data", 
      suggestion = "Please use unique(data.obj$meta.dat$time.var) to view available time points"
    ),
    "feature_missing" = list(
      pattern = "features.plot.*not found|No features found",
      message = "Specified features do not exist in the data",
      suggestion = "Please use rownames(data.obj$feature.tab) to view available features"
    ),
    "dimension_error" = list(
      pattern = "array of at least two dimensions|Matrices must have same dimensions",
      message = "Data dimension mismatch or data structure issue",
      suggestion = "Please check data integrity and ensure sufficient sample size for analysis"
    ),
    "alpha_calculation" = list(
      pattern = "alpha diversity|colSums.*array",
      message = "Alpha diversity calculation failed",
      suggestion = "This may be an underlying function issue, consider providing alpha diversity data manually"
    )
  )
  
  # Match error pattern
  matched_pattern <- NULL
  for (pattern_name in names(error_patterns)) {
    if (grepl(error_patterns[[pattern_name]]$pattern, error_msg, ignore.case = TRUE)) {
      matched_pattern <- error_patterns[[pattern_name]]
      break
    }
  }
  
  # Build friendly error message
  if (!is.null(matched_pattern)) {
    friendly_msg <- sprintf(
      "\n[ERROR] %sError: %s\n[SUGGESTION] Suggested solution: %s\n[FUNCTION] Function: %s\n[DETAILS] Original error: %s",
      ifelse(function_name != "", paste0("[", function_name, "] "), ""),
      matched_pattern$message,
      matched_pattern$suggestion,
      function_name,
      substr(error_msg, 1, 100)
    )
  } else {
    # Generic friendly message
    friendly_msg <- sprintf(
      "\n[ERROR] %sExecution error\n[SUGGESTION] Suggestion: Please check if input parameters are correct\n[FUNCTION] Function: %s\n[DETAILS] Detailed error: %s",
      ifelse(function_name != "", paste0("[", function_name, "] "), ""),
      function_name,
      substr(error_msg, 1, 150)
    )
  }
  
  return(friendly_msg)
}

# Parameter validation helper
mStat_validate_parameters <- function(data.obj, subject.var = NULL, time.var = NULL, group.var = NULL, 
                                      strata.var = NULL, feature.level = NULL, function_name = "") {
  
  warnings <- c()
  errors <- c()
  
  # Check data object structure
  if (is.null(data.obj)) {
    errors <- c(errors, "Data object (data.obj) cannot be null")
  } else {
    if (!is.list(data.obj)) {
      errors <- c(errors, "data.obj must be a list object")
    } else {
      # Check required components
      required_components <- c("feature.tab", "meta.dat")
      missing_components <- setdiff(required_components, names(data.obj))
      if (length(missing_components) > 0) {
        errors <- c(errors, paste("Data object missing required components:", paste(missing_components, collapse = ", ")))
      }
    }
  }
  
  # Check metadata variables
  if (!is.null(data.obj) && !is.null(data.obj$meta.dat)) {
    meta_vars <- names(data.obj$meta.dat)
    
    # Check subject.var
    if (!is.null(subject.var) && !subject.var %in% meta_vars) {
      errors <- c(errors, sprintf("subject.var '%s' does not exist in metadata. Available variables: %s", 
                                  subject.var, paste(meta_vars, collapse = ", ")))
    }
    
    # Check time.var
    if (!is.null(time.var) && !time.var %in% meta_vars) {
      errors <- c(errors, sprintf("time.var '%s' does not exist in metadata. Available variables: %s", 
                                  time.var, paste(meta_vars, collapse = ", ")))
    }
    
    # Check group.var
    if (!is.null(group.var) && !group.var %in% meta_vars) {
      errors <- c(errors, sprintf("group.var '%s' does not exist in metadata. Available variables: %s", 
                                  group.var, paste(meta_vars, collapse = ", ")))
    }
    
    # Check strata.var
    if (!is.null(strata.var) && !strata.var %in% meta_vars) {
      errors <- c(errors, sprintf("strata.var '%s' does not exist in metadata. Available variables: %s", 
                                  strata.var, paste(meta_vars, collapse = ", ")))
    }
    
    # Check time variable values
    if (!is.null(time.var) && time.var %in% meta_vars) {
      time_values <- unique(data.obj$meta.dat[[time.var]])
      if (length(time_values) < 2) {
        warnings <- c(warnings, sprintf("Time variable '%s' has only %d time points, some analyses may not be possible", 
                                        time.var, length(time_values)))
      }
    }
    
    # Check group variable balance
    if (!is.null(group.var) && group.var %in% meta_vars) {
      group_counts <- table(data.obj$meta.dat[[group.var]])
      if (min(group_counts) < 3) {
        warnings <- c(warnings, sprintf("Group variable '%s' has some groups with small sample sizes (minimum: %d), may affect statistical analysis", 
                                        group.var, min(group_counts)))
      }
    }
  }
  
  # Check feature level
  if (!is.null(feature.level) && !is.null(data.obj)) {
    if (feature.level != "original" && !is.null(data.obj$feature.agg.list)) {
      available_levels <- names(data.obj$feature.agg.list)
      if (!feature.level %in% available_levels) {
        errors <- c(errors, sprintf("feature.level '%s' is not available. Available levels: %s", 
                                    feature.level, paste(available_levels, collapse = ", ")))
      }
    }
  }
  
  # Return validation results
  result <- list(
    is_valid = length(errors) == 0,
    errors = errors,
    warnings = warnings
  )
  
  if (!result$is_valid) {
    result$error_message <- mStat_user_friendly_error(
      paste(errors, collapse = "; "), 
      function_name = function_name,
      suggested_fix = "Please check the above parameter settings"
    )
  }
  
  if (length(warnings) > 0) {
    result$warning_message <- paste("[WARNING] Warning:", paste(warnings, collapse = "\n[WARNING] Warning: "))
  }
  
  return(result)
}

# Progress indicator helper
mStat_progress_indicator <- function(current_step, total_steps, step_name = "", show_time = TRUE) {
  
  if (total_steps <= 1) return(invisible(NULL))  # Don't show for simple operations
  
  # Calculate progress percentage
  progress_pct <- round((current_step / total_steps) * 100)
  
  # Create progress bar
  bar_length <- 20
  filled_length <- round((progress_pct / 100) * bar_length)
  bar <- paste0(
    "[", 
    paste(rep("=", filled_length), collapse = ""),
    paste(rep("-", bar_length - filled_length), collapse = ""),
    "]"
  )
  
  # Time estimate (simple)
  time_info <- ""
  if (show_time && current_step > 1) {
    # This is a simplified time estimate - in practice you'd track actual time
    estimated_remaining <- ifelse(current_step < total_steps, 
                                  paste("~", round((total_steps - current_step) * 0.1, 1), "seconds"), 
                                  "Complete")
    time_info <- paste(" |", estimated_remaining)
  }
  
  # Create status message
  status_msg <- sprintf("Progress: %s %d/%d (%d%%) %s%s", 
                        bar, current_step, total_steps, progress_pct, step_name, time_info)
  
  # Print progress (use cat to overwrite previous line)
  if (current_step == 1) {
    cat("\n")  # New line for first progress update
  }
  cat("\r", status_msg)
  
  # Add newline when complete
  if (current_step == total_steps) {
    cat("\n[COMPLETE] Processing complete!\n")
  }
  
  flush.console()  # Ensure immediate output
}

# Enhanced wrapper for commonly failing operations
mStat_safe_execute <- function(expr, function_name = "", operation_name = "", 
                               show_progress = FALSE, progress_info = NULL) {
  
  if (show_progress && !is.null(progress_info)) {
    mStat_progress_indicator(progress_info$current, progress_info$total, operation_name)
  }
  
  tryCatch({
    result <- expr
    return(result)
  }, error = function(e) {
    friendly_error <- mStat_user_friendly_error(
      as.character(e), 
      function_name = function_name,
      suggested_fix = paste("Operation failed:", operation_name)
    )
    stop(friendly_error, call. = FALSE)
  }, warning = function(w) {
    cat("[WARNING] Warning [", function_name, "]:", as.character(w), "\n")
    invokeRestart("muffleWarning")
  })
}

# Data quality checker
mStat_check_data_quality <- function(data.obj, function_name = "") {
  
  issues <- c()
  
  if (!is.null(data.obj$feature.tab)) {
    # Check for empty features
    feature_sums <- rowSums(data.obj$feature.tab, na.rm = TRUE)
    empty_features <- sum(feature_sums == 0)
    if (empty_features > 0) {
      issues <- c(issues, sprintf("Found %d empty features (all-zero rows)", empty_features))
    }
    
    # Check for empty samples
    sample_sums <- colSums(data.obj$feature.tab, na.rm = TRUE)
    empty_samples <- sum(sample_sums == 0)
    if (empty_samples > 0) {
      issues <- c(issues, sprintf("Found %d empty samples (all-zero columns)", empty_samples))
    }
    
    # Check for excessive zeros
    zero_ratio <- sum(data.obj$feature.tab == 0, na.rm = TRUE) / length(data.obj$feature.tab)
    if (zero_ratio > 0.9) {
      issues <- c(issues, sprintf("Data is too sparse (%.1f%% zeros)", zero_ratio * 100))
    }
  }
  
  if (length(issues) > 0) {
    warning_msg <- paste("[DATA QUALITY] Data quality notice [", function_name, "]:\n",
                         paste("  * ", issues, collapse = "\n"), "\n")
    cat(warning_msg)
  }
  
  return(length(issues) == 0)
}