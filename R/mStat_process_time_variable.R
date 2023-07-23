#' Handle Time Variable in Meta Data Table
#'
#' This function handles a time variable in a meta data table. It sets the levels of the time variable based on provided parameters and subsets the data if necessary.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' table (`feature.tab`), feature names, a full feature name list, and a feature aggregation list
#' (`feature.agg.list`).
#' @param time.var A character string that specifies the name of the time variable in the meta data table.
#' @param t0.level A single value that specifies the level to be used as the reference level (optional).
#' @param ts.levels A vector of values that specifies the levels to be used as the time series levels (optional).
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
    # 获取 time.var 的所有唯一值
    unique_vals <- sort(unique(data.obj$meta.dat[[time.var]]))

    # 根据提供的 t0.level 和 ts.levels，创建 levels
    levels <- if (!is.null(t0.level) && !is.null(ts.levels)) {
      c(t0.level, ts.levels)
    } else if (!is.null(t0.level) && is.null(ts.levels)) {
      # 如果 t0.level 不为 NULL，而 ts.levels 为 NULL
      c(t0.level, unique_vals[!unique_vals %in% t0.level])
    } else {
      unique_vals
    }

    # 对 data.obj 中的 meta.dat 的 time.var 进行处理
    data.obj$meta.dat <- data.obj$meta.dat %>%
      dplyr::mutate(!!sym(time.var) := {
        # 先获取 time.var 的值
        time.var_val <- data.obj$meta.dat[[time.var]]

        # 判断 time.var 的类型
        if (is.factor(time.var_val)) {
          factor(time.var_val, levels = levels)
        } else if (is.numeric(time.var_val)) {
          if (!is.null(t0.level) && !is.null(ts.levels)) {
            factor(time.var_val, levels = levels)
          } else {
            time.var_val
          }
        } else if (is.character(time.var_val)) {
          num_time.var <- suppressWarnings(as.numeric(time.var_val))
          if (any(is.na(num_time.var))) {
            message(
              "NA values were generated when converting a character variable to numeric. Please note that no levels, t0.level, and ts.levels were set for the current character variable. As a result, visualizations may exhibit disorderly behavior."
            )
          }
          if (!is.null(t0.level) && !is.null(ts.levels)) {
            factor(time.var_val, levels = levels)
          } else {
            factor(time.var_val, levels = as.character(sort(
              unique(num_time.var), decreasing = FALSE
            )))
          }
        } else {
          stop("time.var must be factor, numeric or character.")
        }
      })

    # 如果 levels 并不包含所有的 levels，则对数据做子集处理
    if (any(is.na(data.obj$meta.dat[[time.var]]))) {
      missing_levels <- setdiff(unique_vals, levels)
      message(
        paste(
          "The data is subsetted to exclude the following levels as they are not in t0.level or ts.levels:",
          paste(missing_levels, collapse = ", ")
        )
      )
      condition <- paste("!is.na(", time.var, ")")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }

    return(data.obj)
  }
