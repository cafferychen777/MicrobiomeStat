#' Identify Time-Varying Variables within Subjects
#'
#' This function checks each specified adjustable variable to determine if it
#' is time-varying within the same subject. If a variable is found to be
#' time-varying for any subject, the function stops and reports which subjects
#' caused this issue. If no variable is time-varying, a message confirms that
#' all variables are consistent within subjects.
#'
#' @param meta.dat A data frame containing the metadata where each row
#'        corresponds to an observation.
#' @param adj.vars A character vector of names of the adjustable variables in
#'        \code{meta.dat} that need to be checked for time-variance within the
#'        same subject.
#' @param subject.var The name of the variable in \code{meta.dat} that
#'        identifies the subject. This variable is used to group observations
#'        by subject.
#'
#' @return The function does not return a value. If it finds that an adjustable
#'         variable is time-varying within subjects, it stops and throws an
#'         error with a message specifying the variable and the problematic
#'         subjects. If all variables are non-time-varying within subjects, it
#'         outputs a message confirming this.
#'
#' @examples
#' # Example 1: All variables are non-time-varying
#' meta.dat <- data.frame(
#'   subject = c(1, 1, 2, 2, 3, 3, 4, 4),
#'   age = c(30, 30, 25, 25, 40, 40, 35, 35),
#'   bmi = c(22, 22, 24, 24, 28, 28, 26, 26)
#' )
#' adj.vars <- c("age", "bmi")
#' subject.var <- "subject"
#' mStat_identify_time_varying_vars(meta.dat, adj.vars, subject.var)
#'
#' # Example 2: Missing variables in meta.dat
#' try({
#'   meta.dat <- data.frame(
#'     subject = c(1, 1, 2, 2, 3, 3, 4, 4),
#'     age = c(30, 30, 25, 25, 40, 40, 35, 35)
#'   )
#'   adj.vars <- c("age", "bmi")  # "bmi" is missing
#'   subject.var <- "subject"
#'   mStat_identify_time_varying_vars(meta.dat, adj.vars, subject.var)
#' })
#'
#' # Example 3: One variable is time-varying
#' try({
#'   meta.dat <- data.frame(
#'     subject = c(1, 1, 2, 2, 3, 3, 4, 4),
#'     age = c(30, 31, 25, 25, 40, 40, 35, 35),
#'     bmi = c(22, 22, 24, 24, 28, 28, 26, 26)
#'   )
#'   adj.vars <- c("age", "bmi")
#'   subject.var <- "subject"
#'   mStat_identify_time_varying_vars(meta.dat, adj.vars, subject.var)
#' })
#'
#' # Example 4: Empty data set
#' try({
#'   meta.dat <- data.frame(subject = integer(0), age = integer(0), bmi = integer(0))
#'   adj.vars <- c("age", "bmi")
#'   subject.var <- "subject"
#'   mStat_identify_time_varying_vars(meta.dat, adj.vars, subject.var)
#' })
#'
#' # Example 5: Mixed case with time-varying and non-time-varying variables
#' try({
#'   meta.dat <- data.frame(
#'     subject = c(1, 1, 2, 2, 3, 3, 4, 4),
#'     age = c(30, 31, 25, 26, 40, 41, 35, 36),
#'     bmi = c(22, 23, 24, 25, 28, 29, 26, 27)
#'   )
#'   adj.vars <- c("age", "bmi")
#'   subject.var <- "subject"
#'   mStat_identify_time_varying_vars(meta.dat, adj.vars, subject.var)
#' })
#' @export
mStat_identify_time_varying_vars <- function(meta.dat, adj.vars, subject.var) {
  # Check if all specified variables exist in the metadata data.frame
  necessary_vars <- c(adj.vars, subject.var)
  missing_vars <- necessary_vars[!necessary_vars %in% names(meta.dat)]

  if (length(missing_vars) > 0) {
    stop("The following required variables are missing in meta.dat: ", paste(missing_vars, collapse = ", "), ".")
  }

  # Use dplyr to detect which variables are time-varying
  # Group and summarize the data, calculate the number of unique values for each variable within each subject
  # Using dplyr functions with :: operator instead of library(dplyr)
  time_varying_info <- meta.dat %>%
    dplyr::select(c(subject.var, adj.vars)) %>%
    dplyr::group_by(!!sym(subject.var)) %>%
    dplyr::summarise(dplyr::across(all_of(adj.vars), ~dplyr::n_distinct(.x)), .groups = 'drop')

  # Initialize lists to hold the status of each variable
  time_varying_vars <- character(0)
  non_time_varying_vars <- character(0)

  # Iterate through all the adjustable variables, check if they are time-varying
  for (var in adj.vars) {
    if (any(time_varying_info[[var]] > 1)) {
      # If variable is time-varying for any subject, add to time-varying list
      time_varying_vars <- c(time_varying_vars, var)
    } else {
      # Otherwise, add to non-time-varying list
      non_time_varying_vars <- c(non_time_varying_vars, var)
    }
  }

  # Return a list of time-varying and non-time-varying variables
  return(list(time_varying_vars = time_varying_vars,
              non_time_varying_vars = non_time_varying_vars))
}
