#' Report Generation Utility Functions
#'
#' Internal utility functions for report generation.
#' These functions reduce code duplication across report generators.
#'
#' @name report_utils
#' @keywords internal
NULL

#' Ensure Output File Has Correct Extension
#'
#' Adds the appropriate file extension if not present.
#'
#' @param output.file The output file path
#' @param output.format The output format ("pdf" or "html")
#' @return The file path with correct extension
#' @noRd
ensure_file_extension <- function(output.file, output.format) {
  if (output.format == "pdf" && !grepl("\\.pdf$", output.file, ignore.case = TRUE)) {
    output.file <- paste0(output.file, ".pdf")
  } else if (output.format == "html" && !grepl("\\.html$", output.file, ignore.case = TRUE)) {
    output.file <- paste0(output.file, ".html")
  }
  output.file
}

#' Ensure Output Directory Exists
#'
#' Creates the output directory if it doesn't exist.
#'
#' @param output.file The output file path
#' @return The directory path (invisibly)
#' @noRd
ensure_output_dir <- function(output.file) {
  output_dir <- dirname(output.file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  invisible(output_dir)
}

#' Get Result Output Mode
#'
#' Determines the knitr result output mode based on pdf setting.
#'
#' @param pdf Logical indicating if PDF output
#' @return "asis" for PDF, "markup" otherwise
#' @noRd
get_result_output_mode <- function(pdf) {
  if (pdf) "asis" else "markup"
}

#' Create Parameter Status String
#'
#' Converts a parameter value to a display string for reports.
#'
#' @param param The parameter value to check
#' @param show_value If TRUE and param is not NULL, show the actual value
#' @return A string representation of the parameter status
#' @noRd
param_status <- function(param, show_value = FALSE) {
  if (is.null(param)) {
    "NULL"
  } else if (show_value && (is.character(param) || is.numeric(param))) {
    toString(param)
  } else {
    "Not NULL"
  }
}

#' Filter Results by P-value Method
#'
#' Filters a data frame based on the multiple testing correction method.
#'
#' @param df Data frame with P.Value and Adjusted.P.Value columns
#' @param method Either "fdr" or "none"
#' @param threshold Significance threshold
#' @return A list with filtered results and p-value description string
#' @noRd
filter_by_p_method <- function(df, method, threshold) {
  if (method == "fdr") {
    list(
      results = dplyr::filter(df, Adjusted.P.Value < threshold),
      p_str = "adjusted p-value"
    )
  } else {
    list(
      results = dplyr::filter(df, P.Value < threshold),
      p_str = "p-value"
    )
  }
}

#' Generate YAML Output Format String
#'
#' Creates the YAML output format section based on output format.
#'
#' @param output.format The output format ("pdf" or "html")
#' @return YAML string for the output format
#' @noRd
generate_yaml_output <- function(output.format) {
  if (output.format == "pdf") {
    paste0(
      "output:\n",
      "  pdf_document:\n",
      "    latex_engine: pdflatex\n",
      "    keep_tex: false"
    )
  } else {
    paste0(
      "output:\n",
      "  html_document:\n",
      "    toc: true\n",
      "    toc_float: true"
    )
  }
}

#' Save Plots as PDF
#'
#' Saves a list of plots to a PDF file.
#'
#' @param plots List of plot objects
#' @param file_path Full path to the PDF file
#' @param width PDF width in inches
#' @param height PDF height in inches
#' @return The file path (invisibly)
#' @noRd
save_plots_pdf <- function(plots, file_path, width = 11, height = 8.5) {
  ensure_output_dir(file_path)
  pdf(file_path, width = width, height = height)
  lapply(plots, print)
  dev.off()
  invisible(file_path)
}

#' Select Metadata Variables Safely
#'
#' Selects variables from metadata, filtering out NULL variable names.
#' This is useful when optional parameters like time.var may be NULL.
#'
#' @param meta.dat The metadata data frame
#' @param ... Variable names (can include NULL values which will be filtered out)
#' @return Data frame with selected columns
#' @noRd
select_meta_vars <- function(meta.dat, ...) {
 vars <- c(...)
 vars <- vars[!sapply(vars, is.null)]
 meta.dat %>% as.data.frame() %>% dplyr::select(dplyr::all_of(vars))
}

#' Calculate Sample Count for Plot Width
#'
#' Calculates the sample count for determining plot width.
#' Handles NULL time.var for single time point analysis.
#'
#' @param meta.dat The metadata data frame
#' @param time.var Time variable name (can be NULL)
#' @param group.var Group variable name
#' @return Integer count of unique combinations
#' @noRd
calc_sample_count <- function(meta.dat, time.var, group.var) {
 group_count <- length(unique(meta.dat[[group.var]]))
 if (!is.null(time.var) && time.var %in% names(meta.dat)) {
   time_count <- length(unique(meta.dat[[time.var]]))
   return(time_count * group_count)
 }
 group_count
}

#' Build PDF Filename for Report Plots
#'
#' Constructs a standardized PDF filename based on analysis parameters.
#'
#' @param prefix Base name prefix (e.g., "taxa_boxplot")
#' @param params Named list of parameters to include in filename
#' @param file.ann Optional file annotation
#' @return The constructed filename
#' @noRd
build_pdf_filename <- function(prefix, params, file.ann = NULL) {
  parts <- prefix

  for (name in names(params)) {
    value <- params[[name]]
    if (!is.null(value)) {
      parts <- paste0(parts, "_", name, "_", paste(value, collapse = "-"))
    }
  }

  if (!is.null(file.ann)) {
    parts <- paste0(parts, "_", file.ann)
  }

  paste0(parts, ".pdf")
}
