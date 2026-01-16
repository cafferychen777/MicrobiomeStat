#!/usr/bin/env Rscript
# Batch update parameter documentation to use @inheritParams
# This script replaces verbose repeated parameter docs with @inheritParams

library(stringr)

# Define the long data.obj pattern to replace
long_data_obj_pattern <- "#' @param data\\.obj A list object in a format specific to MicrobiomeStat.*?(?=#'\\s*@(?:param|return|examples|export|name|title|description|details|seealso|author|inheritParams))"

# Parameters that should be inherited from mStat_data_obj_doc
data_obj_params <- c(
  "data.obj", "subject.var", "time.var", "group.var", "strata.var",
 "adj.vars", "feature.level", "feature.dat.type", "prev.filter",
 "abund.filter", "t0.level", "ts.levels"
)

# Parameters that should be inherited from mStat_plot_params_doc
plot_params <- c(
  "base.size", "theme.choice", "custom.theme", "palette",
  "pdf", "file.ann", "pdf.wid", "pdf.hei"
)

# Parameters that should be inherited from mStat_test_params_doc
test_params <- c(
  "alpha.obj", "alpha.name", "dist.obj", "dist.name", "depth"
)

# Function to update a single file
update_file <- function(filepath) {
  content <- readLines(filepath, warn = FALSE)
  full_content <- paste(content, collapse = "\n")

  # Check if already updated
  if (grepl("@inheritParams mStat_data_obj_doc", full_content)) {
    message(sprintf("  [SKIP] %s - already updated", basename(filepath)))
    return(FALSE)
  }

  # Check if contains the long pattern
  if (!grepl("@param data\\.obj A list object in a format specific to MicrobiomeStat", full_content)) {
    message(sprintf("  [SKIP] %s - no matching pattern", basename(filepath)))
    return(FALSE)
  }

  # Find which inherit params are needed based on what params exist in file
  needs_data_obj_doc <- any(sapply(data_obj_params, function(p) {
    grepl(paste0("@param ", p, "\\b"), full_content)
  }))

  needs_plot_doc <- any(sapply(plot_params, function(p) {
    grepl(paste0("@param ", p, "\\b"), full_content)
  }))

  needs_test_doc <- any(sapply(test_params, function(p) {
    grepl(paste0("@param ", p, "\\b"), full_content)
  }))

  # Build inheritParams block
  inherit_lines <- c()
  if (needs_data_obj_doc) inherit_lines <- c(inherit_lines, "#' @inheritParams mStat_data_obj_doc")
  if (needs_test_doc) inherit_lines <- c(inherit_lines, "#' @inheritParams mStat_test_params_doc")
  if (needs_plot_doc) inherit_lines <- c(inherit_lines, "#' @inheritParams mStat_plot_params_doc")

  inherit_block <- paste(inherit_lines, collapse = "\n")

  # Patterns to remove (will be inherited)
  all_inherit_params <- c(data_obj_params, plot_params, test_params)

  # Process line by line
  new_lines <- c()
  i <- 1
  inherit_added <- FALSE
  in_param_to_remove <- FALSE

  while (i <= length(content)) {
    line <- content[i]

    # Check if this is a @param line for an inherited parameter
    is_inherit_param <- FALSE
    for (p in all_inherit_params) {
      if (grepl(paste0("^#'\\s*@param\\s+", p, "\\b"), line)) {
        is_inherit_param <- TRUE
        in_param_to_remove <- TRUE
        break
      }
    }

    # Check if we're starting a new @param or other tag (end of multi-line param)
    if (in_param_to_remove && grepl("^#'\\s*@(param|return|examples|export|name|details|seealso|author)", line) && !is_inherit_param) {
      in_param_to_remove <- FALSE
    }

    # Check if we've left roxygen block
    if (in_param_to_remove && !grepl("^#'", line)) {
      in_param_to_remove <- FALSE
    }

    if (is_inherit_param) {
      # Add inheritParams block before first removed param
      if (!inherit_added) {
        new_lines <- c(new_lines, "#'")
        new_lines <- c(new_lines, inherit_block)
        new_lines <- c(new_lines, "#'")
        inherit_added <- TRUE
      }
      # Skip this line (inherited param)
    } else if (in_param_to_remove) {
      # Skip continuation of inherited param
    } else {
      new_lines <- c(new_lines, line)
    }

    i <- i + 1
  }

  # Write updated content
  writeLines(new_lines, filepath)
  message(sprintf("  [UPDATED] %s", basename(filepath)))
  return(TRUE)
}

# Get all files to update
files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
files <- files[!grepl("mStat_param_docs\\.R$", files)]  # Exclude the template file

message("Starting batch update of parameter documentation...")
message(sprintf("Found %d R files to process\n", length(files)))

updated_count <- 0
for (f in files) {
  if (update_file(f)) {
    updated_count <- updated_count + 1
  }
}

message(sprintf("\nCompleted! Updated %d files.", updated_count))
