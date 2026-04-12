#!/usr/bin/env Rscript
# Batch update roxygen parameter documentation to use @inheritParams tags.

get_current_file <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      if (!is.null(frame$ofile)) {
        return(frame$ofile)
      }
      NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]

  if (length(frame_files) > 0) {
    return(normalizePath(frame_files[length(frame_files)], winslash = "/", mustWork = TRUE))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_args <- grep("^--file=", args, value = TRUE)
  if (length(file_args) > 0) {
    script_path <- sub("^--file=", "", file_args[length(file_args)])
    return(normalizePath(script_path, winslash = "/", mustWork = TRUE))
  }

  stop("Unable to determine current script path.")
}

script_path <- get_current_file()
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
r_dir <- file.path(repo_root, "R")

data_obj_params <- c(
  "data.obj", "subject.var", "time.var", "group.var", "strata.var",
  "adj.vars", "feature.level", "feature.dat.type", "prev.filter",
  "abund.filter", "t0.level", "ts.levels"
)

plot_params <- c(
  "base.size", "theme.choice", "custom.theme", "palette",
  "pdf", "file.ann", "pdf.wid", "pdf.hei"
)

test_params <- c(
  "alpha.obj", "alpha.name", "dist.obj", "dist.name", "depth"
)

update_file <- function(filepath) {
  content <- readLines(filepath, warn = FALSE)
  full_content <- paste(content, collapse = "\n")

  if (grepl("@inheritParams mStat_data_obj_doc", full_content, fixed = TRUE)) {
    message(sprintf("  [SKIP] %s - already updated", basename(filepath)))
    return(FALSE)
  }

  if (!grepl("@param data\\.obj A list object in a format specific to MicrobiomeStat", full_content)) {
    message(sprintf("  [SKIP] %s - no matching pattern", basename(filepath)))
    return(FALSE)
  }

  needs_data_obj_doc <- any(vapply(data_obj_params, function(param_name) {
    grepl(paste0("@param ", param_name, "\\b"), full_content)
  }, logical(1)))

  needs_plot_doc <- any(vapply(plot_params, function(param_name) {
    grepl(paste0("@param ", param_name, "\\b"), full_content)
  }, logical(1)))

  needs_test_doc <- any(vapply(test_params, function(param_name) {
    grepl(paste0("@param ", param_name, "\\b"), full_content)
  }, logical(1)))

  inherit_lines <- character(0)
  if (needs_data_obj_doc) {
    inherit_lines <- c(inherit_lines, "#' @inheritParams mStat_data_obj_doc")
  }
  if (needs_test_doc) {
    inherit_lines <- c(inherit_lines, "#' @inheritParams mStat_test_params_doc")
  }
  if (needs_plot_doc) {
    inherit_lines <- c(inherit_lines, "#' @inheritParams mStat_plot_params_doc")
  }

  inherit_block <- paste(inherit_lines, collapse = "\n")
  inherited_params <- unique(c(data_obj_params, plot_params, test_params))

  new_lines <- character(0)
  inherit_added <- FALSE
  in_param_to_remove <- FALSE

  for (line in content) {
    is_inherited_param <- any(vapply(inherited_params, function(param_name) {
      grepl(paste0("^#'\\s*@param\\s+", param_name, "\\b"), line)
    }, logical(1)))

    if (in_param_to_remove &&
        grepl("^#'\\s*@(param|return|examples|export|name|details|seealso|author)", line) &&
        !is_inherited_param) {
      in_param_to_remove <- FALSE
    }

    if (in_param_to_remove && !grepl("^#'", line)) {
      in_param_to_remove <- FALSE
    }

    if (is_inherited_param) {
      in_param_to_remove <- TRUE

      if (!inherit_added) {
        new_lines <- c(new_lines, "#'", inherit_block, "#'")
        inherit_added <- TRUE
      }

      next
    }

    if (!in_param_to_remove) {
      new_lines <- c(new_lines, line)
    }
  }

  writeLines(new_lines, filepath)
  message(sprintf("  [UPDATED] %s", basename(filepath)))
  TRUE
}

files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
files <- files[!grepl("mStat_param_docs\\.R$", files)]

message("Starting batch update of parameter documentation...")
message(sprintf("Processing %d R files under %s\n", length(files), r_dir))

updated_count <- 0L
for (filepath in files) {
  if (update_file(filepath)) {
    updated_count <- updated_count + 1L
  }
}

message(sprintf("\nCompleted! Updated %d files.", updated_count))
