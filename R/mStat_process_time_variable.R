#' Handle Time Variable in Meta Data Table
#'
#' This function handles a time variable in a meta data table. It sets the levels of the time variable based on provided parameters and subsets the data if necessary.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' table (`feature.tab`), feature names, a full feature name list, and a feature aggregation list
#' (`feature.agg.list`).
#' @param time.var A character string that specifies the name of the time variable in the meta data table.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#'
#' @return A data frame that is the modified version of the input meta data table.
#'
#' @examples
#' # Assuming 'meta_data' is your meta data table and 'time' is your time variable
#' # meta_data_handled <- mStat_process_time_variable(meta_data, 'time')
#' data(ecam.obj)
#' mStat_process_time_variable(data.obj = ecam.obj,time.var = "month",
#' t0.level = "0",ts.levels = c("1","2"))
#'
#' @details
#' The function first gets the unique values of the time variable. Then, it sets the levels of the time variable based on the provided t0.level and ts.levels. If these parameters are not provided, the levels are set to the unique values of the time variable. The function then modifies the time variable in the meta data table based on its type (factor, numeric, or character). If the time variable is a character, it is converted to numeric. If there are any NA values generated during this conversion, a message is displayed. Finally, if the levels set do not include all the unique values of the time variable, the function subsets the data to exclude the missing levels and displays a message.
#'
#' @author Your Name
#' @seealso \code{\link[dplyr]{mutate}}, \code{\link[rlang]{sym}}
#' @export
mStat_process_time_variable <-
  function(data.obj,
           time.var,
           t0.level = NULL,
           ts.levels = NULL) {
    # Extract all unique values of the time variable and sort them.
    # This step is crucial for understanding the range and distribution of time points in the dataset.
    unique_vals <- sort(unique(data.obj$meta.dat[[time.var]]))

    # Determine the levels for the time variable based on provided parameters.
    # This allows for flexible handling of time points, accommodating different study designs.
    levels <- if (!is.null(t0.level) && !is.null(ts.levels)) {
      # If both baseline and follow-up time points are specified, use them as levels.
      c(t0.level, ts.levels)
    } else if (!is.null(t0.level) && is.null(ts.levels)) {
      # If only baseline is specified, use it as the first level and include all other time points.
      unique_vals <- as.character(unique_vals)
      c(t0.level, unique_vals[!unique_vals %in% t0.level])
    } else {
      # If no specific levels are provided, use all unique values as levels.
      unique_vals
    }

    # Process the time variable in the metadata of the data object.
    # This step ensures that the time variable is properly formatted for downstream analyses.
    data.obj$meta.dat <- data.obj$meta.dat %>%
      dplyr::mutate(!!sym(time.var) := {
        # Retrieve the current values of the time variable.
        time.var_val <- data.obj$meta.dat[[time.var]]

        # Handle the time variable based on its data type.
        # This ensures appropriate treatment of different time formats (factor, numeric, character).
        if (is.factor(time.var_val)) {
          # For factor variables, reorder levels if necessary.
          factor(time.var_val, levels = levels)
        } else if (is.numeric(time.var_val)) {
          # For numeric variables, convert to factor if specific levels are provided.
          if (!is.null(t0.level) && !is.null(ts.levels)) {
            factor(time.var_val, levels = levels)
          } else {
            time.var_val
          }
        } else if (is.character(time.var_val)) {
          # For character variables, attempt to convert to numeric and handle accordingly.
          num_time.var <- suppressWarnings(as.numeric(time.var_val))
          if (any(is.na(num_time.var))) {
            # Warn about potential issues with character to numeric conversion.
            message(
              "NA values were generated when converting a character variable to numeric. Please note that no levels, t0.level, and ts.levels were set for the current character variable. As a result, visualizations may exhibit disorderly behavior."
            )
          }
          if (!is.null(t0.level) && !is.null(ts.levels)) {
            factor(time.var_val, levels = levels)
          } else {
            # Create an ordered factor based on the numeric values.
            factor(time.var_val, levels = as.character(sort(
              unique(num_time.var), decreasing = FALSE
            )))
          }
        } else {
          # Stop execution if the time variable is not of an expected type.
          stop("time.var must be factor, numeric or character.")
        }
      })

    # Handle cases where some levels are not present in the specified t0.level or ts.levels.
    # This ensures that the dataset only includes the time points of interest.
    if (any(is.na(data.obj$meta.dat[[time.var]]))) {
      # Identify missing levels and inform the user.
      missing_levels <- setdiff(unique_vals, levels)
      message(
        paste(
          "The data is subsetted to exclude the following levels as they are not in t0.level or ts.levels:",
          paste(missing_levels, collapse = ", ")
        )
      )
      # Create a condition to subset the data, excluding NA values in the time variable.
      condition <- paste("!is.na(", time.var, ")")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    # Return the processed data object.
    # The time variable in this object is now properly formatted for longitudinal analyses.
    return(data.obj)
  }
