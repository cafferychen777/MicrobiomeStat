test_that("difference tests fail fast when group.var is missing", {
  data(peerj32.obj)

  expect_error(
    test_alpha(peerj32.obj, test.type = "difference", alpha.name = "shannon"),
    "group.var is required"
  )

  expect_error(
    test_beta(peerj32.obj, test.type = "difference", dist.name = "BC"),
    "group.var is required"
  )

  expect_error(
    test_taxa(peerj32.obj, test.type = "difference", feature.level = "Genus"),
    "group.var is required"
  )
})

test_that("test_taxa propagates multiple-testing settings into formatted output", {
  data(peerj32.obj)

  default_result <- suppressWarnings(
    suppressMessages(
      test_taxa(
        peerj32.obj,
        test.type = "difference",
        group.var = "group",
        feature.level = "Genus"
      )
    )
  )

  none_result <- suppressWarnings(
    suppressMessages(
      test_taxa(
        peerj32.obj,
        test.type = "difference",
        group.var = "group",
        feature.level = "Genus",
        feature.mt.method = "none",
        feature.sig.level = 1
      )
    )
  )

  expect_false(identical(default_result, none_result))
  expect_equal(
    none_result$Genus[[1]]$Adjusted.P.Value,
    none_result$Genus[[1]]$P.Value
  )
  expect_true(all(none_result$Genus[[1]]$Significant))
})


test_that("unified interface scopes pair time.points before routing", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(5, 5, 5, 5, 5, 5),
      nrow = 1,
      dimnames = list("f1", sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2"), each = 3),
      time = rep(c("T1", "T2", "T10"), times = 2),
      group = rep(c("A", "B"), each = 3),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )
  alpha.obj <- list(
    shannon = data.frame(
      shannon = c(1, 2, 3, 4, 5, 6),
      row.names = sample_ids
    )
  )

  result <- suppressWarnings(
    test_alpha(
      data.obj,
      test.type = "difference",
      alpha.obj = alpha.obj,
      alpha.name = "shannon",
      subject.var = "subject",
      time.var = "time",
      time.points = c("T1", "T2"),
      group.var = "group"
    )
  )

  expect_true("shannon" %in% names(result))
  expect_true(nrow(result$shannon) > 0)
})


test_that("unified interface rejects unsupported design instead of falling back", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        5, 4, 3, 2, 1, 0,
        0, 1, 2, 3, 4, 5
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2"), each = 3),
      time = rep(c("T1", "T2", "T3"), times = 2),
      group = rep(c("A", "B"), each = 3),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  expect_error(
    test_taxa(
      data.obj,
      test.type = "difference",
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      feature.level = "original"
    ),
    "Supported designs"
  )
})

test_that("unified interface validates strata.var before dispatch", {
  data(peerj32.obj)

  expect_error(
    plot_taxa(
      peerj32.obj,
      plot.type = "barplot",
      feature.level = "Phylum",
      strata.var = "not_a_column"
    ),
    "strata.var 'not_a_column' not found in meta.dat"
  )
})

test_that("data summary uses sample count for average reads per sample", {
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 20,
        30, 40,
        50, 60
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("tax", 1:3), c("s1", "s2"))
    )
  )

  summary_tbl <- suppressMessages(mStat_summarize_data_obj(data.obj))
  avg_reads <- summary_tbl$Value[
    summary_tbl$Variable == "Average reads per sample"
  ]

  expect_equal(as.numeric(avg_reads), 105)
})

test_that("data summary works without feature annotations", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("tax1", "tax2"), c("s1", "s2"))
    )
  )

  expect_no_error(suppressMessages(mStat_summarize_data_obj(data.obj)))
})

test_that("data summary handles missing group.var and Date time variables", {
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 20,
        30, 40
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("tax1", "tax2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      visit = as.Date(c("2024-01-01", "2024-01-10")),
      row.names = c("s1", "s2")
    )
  )

  summary_tbl <- suppressMessages(
    suppressWarnings(mStat_summarize_data_obj(data.obj, time.var = "visit"))
  )

  expect_true(all(c(
    "Earliest sample time-point",
    "Latest sample time-point"
  ) %in% summary_tbl$Variable))
})

test_that("resolve_params preserves explicit pdf values for plot targets", {
  design_info <- list(design = "single", t.level = NULL)

  resolved_default <- resolve_params(
    args = list(theme = "bw"),
    target_func = "generate_taxa_barplot_single",
    design_info = design_info
  )
  resolved_explicit <- resolve_params(
    args = list(theme = "bw", pdf = TRUE),
    target_func = "generate_taxa_barplot_single",
    design_info = design_info
  )

  expect_false(resolved_default$pdf)
  expect_true(resolved_explicit$pdf)
})

test_that("alpha preparation keeps sample ids when selecting covariates", {
  data(peerj32.obj)

  meta.dat <- peerj32.obj$meta.dat %>%
    tibble::rownames_to_column("sample")

  alpha.obj <- suppressWarnings(
    mStat_calculate_alpha_diversity(
      x = peerj32.obj$feature.tab,
      alpha.name = "shannon"
    )
  )

  prepared <- mStat_prepare_alpha_data(
    alpha.obj = alpha.obj,
    meta.dat = meta.dat,
    vars = "group",
    sample_col = "sample",
    join = "inner"
  )

  expect_equal(nrow(prepared), nrow(meta.dat))
  expect_true(all(c("sample", "shannon", "group") %in% names(prepared)))
  expect_setequal(prepared$sample, meta.dat$sample)
})

test_that("change metadata casts join keys to match long-form distances", {
  change.df <- tibble::tibble(
    subject = c("s1", "s1", "s2"),
    month = c("1", "2", "1"),
    distance = c(0.1, 0.2, 0.3)
  )

  meta.dat <- data.frame(
    subject = c("s1", "s1", "s2"),
    month = c(1, 2, 1),
    diet = c("A", "B", "A"),
    row.names = c("sample-1", "sample-2", "sample-3"),
    stringsAsFactors = FALSE
  )

  attached <- mStat_attach_change_metadata(
    change.df = change.df,
    meta.dat = meta.dat,
    by = c("subject", "month"),
    vars = "diet"
  )

  expect_equal(attached$diet, c("A", "B", "A"))
})

test_that("plot_taxa fallback clears stale change.type before dispatch", {
  ns <- asNamespace("MicrobiomeStat")
  original_validate_inputs <- get("validate_inputs", envir = ns)
  original_filter_target_params <- get("filter_target_params", envir = ns)
  original_generate_taxa_barplot_single <- get("generate_taxa_barplot_single", envir = ns)
  call_count <- 0L
  captured <- NULL

  unlockBinding("validate_inputs", ns)
  assign(
    "validate_inputs",
    function(...) {
      call_count <<- call_count + 1L
      if (call_count == 1L) {
        return(list(is_valid = FALSE, errors = "change unsupported", warnings = character()))
      }
      list(
        is_valid = TRUE,
        errors = character(),
        warnings = character(),
        target_func = "generate_taxa_barplot_single"
      )
    },
    envir = ns
  )
  lockBinding("validate_inputs", ns)

  unlockBinding("filter_target_params", ns)
  assign(
    "filter_target_params",
    function(resolved, target_func = NULL) resolved,
    envir = ns
  )
  lockBinding("filter_target_params", ns)

  unlockBinding("generate_taxa_barplot_single", ns)
  assign(
    "generate_taxa_barplot_single",
    function(...) {
      captured <<- list(...)
      list(mock_plot = ggplot2::ggplot())
    },
    envir = ns
  )
  lockBinding("generate_taxa_barplot_single", ns)

  on.exit({
    unlockBinding("validate_inputs", ns)
    assign("validate_inputs", original_validate_inputs, envir = ns)
    lockBinding("validate_inputs", ns)

    unlockBinding("filter_target_params", ns)
    assign("filter_target_params", original_filter_target_params, envir = ns)
    lockBinding("filter_target_params", ns)

    unlockBinding("generate_taxa_barplot_single", ns)
    assign("generate_taxa_barplot_single", original_generate_taxa_barplot_single, envir = ns)
    lockBinding("generate_taxa_barplot_single", ns)
  }, add = TRUE)

  data(peerj32.obj)

  result <- suppressMessages(
    plot_taxa(
      peerj32.obj,
      plot.type = "barplot",
      feature.level = "Phylum",
      change.type = "relative"
    )
  )

  expect_type(result, "list")
  expect_true(is.null(captured$change.type))
  expect_false("feature.change.func" %in% names(captured))
})

test_that("test_taxa exposes per_time through the unified interface", {
  ns <- asNamespace("MicrobiomeStat")
  original_generate_taxa_per_time_test_long <- get("generate_taxa_per_time_test_long", envir = ns)
  captured <- NULL

  unlockBinding("generate_taxa_per_time_test_long", ns)
  assign(
    "generate_taxa_per_time_test_long",
    function(...) {
      captured <<- list(...)
      list(mock = data.frame(ok = TRUE))
    },
    envir = ns
  )
  lockBinding("generate_taxa_per_time_test_long", ns)

  on.exit({
    unlockBinding("generate_taxa_per_time_test_long", ns)
    assign("generate_taxa_per_time_test_long", original_generate_taxa_per_time_test_long, envir = ns)
    lockBinding("generate_taxa_per_time_test_long", ns)
  }, add = TRUE)

  data.obj <- list(
    feature.tab = matrix(
      c(
        5, 6, 7, 8, 9, 10,
        2, 3, 4, 5, 6, 7
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), paste0("s", 1:6))
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2"), each = 3),
      time = rep(c("T1", "T2", "T3"), times = 2),
      group = rep(c("A", "B"), each = 3),
      row.names = paste0("s", 1:6),
      stringsAsFactors = FALSE
    )
  )

  result <- suppressMessages(
    suppressWarnings(
      test_taxa(
        data.obj,
        test.type = "per_time",
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        feature.level = "original"
      )
    )
  )

  expect_type(result, "list")
  expect_false(is.null(captured))
  expect_identical(attr(result, "test.type"), "per_time")
})

test_that("single report overview uses sample-id scoping for time level filtering", {
  section_text <- generate_single_report_section_overview()

  expect_match(section_text, "mStat_subset_data\\(data.obj, samIDs = subset_ids\\)")
  expect_no_match(section_text, "condition <- paste\\(time.var, '==', t.level")
})

test_that("single helper files use shared sample-id time scoping", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  helper_paths <- c(
    "R/generate_taxa_barplot_single.R",
    "R/generate_taxa_dotplot_single.R",
    "R/generate_taxa_heatmap_single.R",
    "R/generate_taxa_indiv_boxplot_single.R",
    "R/generate_taxa_boxplot_single.R",
    "R/generate_taxa_test_single.R",
    "R/generate_beta_test_single.R",
    "R/generate_beta_ordination_single.R"
  )

  for (path in helper_paths) {
    contents <- paste(readLines(file.path(package_root, path), warn = FALSE), collapse = "\n")
    expect_true(
      grepl("mStat_subset_by_meta_values", contents, fixed = TRUE) ||
        grepl("mStat_prepare_taxa_single_context", contents, fixed = TRUE)
    )
    expect_no_match(contents, "condition <- paste")
  }
})

test_that("time variable processing subsets by sample ids instead of parsed condition strings", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  contents <- paste(
    readLines(file.path(package_root, "R/mStat_process_time_variable.R"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(contents, "samIDs = rownames\\(data.obj\\$meta.dat\\)\\[keep_rows\\]")
  expect_no_match(contents, "condition <- paste\\(\"!is.na")
})

test_that("single and long taxa helpers reuse shared context preparation", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  helper_paths <- c(
    "R/generate_taxa_dotplot_single.R",
    "R/generate_taxa_boxplot_single.R",
    "R/generate_taxa_indiv_boxplot_single.R",
    "R/generate_taxa_test_single.R",
    "R/generate_taxa_barplot_single.R",
    "R/generate_taxa_spaghettiplot_long.R",
    "R/generate_taxa_barplot_long.R",
    "R/generate_taxa_boxplot_long.R"
  )

  for (path in helper_paths) {
    contents <- paste(readLines(file.path(package_root, path), warn = FALSE), collapse = "\n")
    expect_true(
      grepl("mStat_prepare_taxa_single_context", contents, fixed = TRUE) ||
        grepl("mStat_prepare_taxa_long_context", contents, fixed = TRUE)
    )
  }
})

test_that("pair helpers reuse shared metadata attachment and feature selection helpers", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  helper_paths <- c(
    "R/generate_taxa_change_boxplot_pair.R",
    "R/generate_taxa_indiv_change_boxplot_pair.R",
    "R/generate_taxa_change_scatterplot_pair.R",
    "R/generate_taxa_indiv_change_scatterplot_pair.R",
    "R/generate_taxa_change_heatmap_pair.R",
    "R/generate_taxa_change_test_pair.R"
  )

  for (path in helper_paths) {
    contents <- paste(readLines(file.path(package_root, path), warn = FALSE), collapse = "\n")
    expect_match(contents, "mStat_prepare_taxa_pair_change_data")
  }

  expect_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_boxplot_pair.R"), warn = FALSE), collapse = "\n"),
    "mStat_attach_pair_metadata"
  )
  expect_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_test_pair.R"), warn = FALSE), collapse = "\n"),
    "mStat_attach_subject_level_metadata"
  )
  expect_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_heatmap_pair.R"), warn = FALSE), collapse = "\n"),
    "mStat_attach_subject_level_metadata"
  )
  expect_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_scatterplot_pair.R"), warn = FALSE), collapse = "\n"),
    "mStat_resolve_selected_features"
  )
})

test_that("long helper files avoid empty-string group placeholders and feature leakage patterns", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)

  barplot_long <- paste(readLines(file.path(package_root, "R/generate_taxa_barplot_long.R"), warn = FALSE), collapse = "\n")
  areaplot_long <- paste(readLines(file.path(package_root, "R/generate_taxa_areaplot_long.R"), warn = FALSE), collapse = "\n")
  dotplot_single <- paste(readLines(file.path(package_root, "R/generate_taxa_dotplot_single.R"), warn = FALSE), collapse = "\n")

  expect_no_match(barplot_long, 'group.var = ""')
  expect_no_match(areaplot_long, 'group.var = ""')
  expect_match(barplot_long, "mStat_ensure_group_placeholder")
  expect_match(areaplot_long, "mStat_ensure_group_placeholder")
  expect_match(dotplot_single, "current_features_plot")
  expect_no_match(dotplot_single, 'features.plot <- names\\(sort')
})

test_that("remaining taxa helpers use local feature selection variables and grouped summary helper", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)

  boxplot_single <- paste(readLines(file.path(package_root, "R/generate_taxa_boxplot_single.R"), warn = FALSE), collapse = "\n")
  heatmap_single <- paste(readLines(file.path(package_root, "R/generate_taxa_heatmap_single.R"), warn = FALSE), collapse = "\n")
  change_dotplot_pair <- paste(readLines(file.path(package_root, "R/generate_taxa_change_dotplot_pair.R"), warn = FALSE), collapse = "\n")
  heatmap_pair <- paste(readLines(file.path(package_root, "R/generate_taxa_heatmap_pair.R"), warn = FALSE), collapse = "\n")
  indiv_boxplot_long <- paste(readLines(file.path(package_root, "R/generate_taxa_indiv_boxplot_long.R"), warn = FALSE), collapse = "\n")
  indiv_spaghetti_long <- paste(readLines(file.path(package_root, "R/generate_taxa_indiv_spaghettiplot_long.R"), warn = FALSE), collapse = "\n")
  heatmap_long <- paste(readLines(file.path(package_root, "R/generate_taxa_heatmap_long.R"), warn = FALSE), collapse = "\n")
  spaghetti_long <- paste(readLines(file.path(package_root, "R/generate_taxa_spaghettiplot_long.R"), warn = FALSE), collapse = "\n")
  change_heatmap_long <- paste(readLines(file.path(package_root, "R/generate_taxa_change_heatmap_long.R"), warn = FALSE), collapse = "\n")
  ma_single <- paste(readLines(file.path(package_root, "R/generate_taxa_ma_plot_single.R"), warn = FALSE), collapse = "\n")
  volcano_single <- paste(readLines(file.path(package_root, "R/generate_taxa_volcano_single.R"), warn = FALSE), collapse = "\n")
  trend_volcano_long <- paste(readLines(file.path(package_root, "R/generate_taxa_trend_volcano_long.R"), warn = FALSE), collapse = "\n")
  volatility_volcano_long <- paste(readLines(file.path(package_root, "R/generate_taxa_volatility_volcano_long.R"), warn = FALSE), collapse = "\n")
  per_time_dotplot_long <- paste(readLines(file.path(package_root, "R/generate_taxa_per_time_dotplot_long.R"), warn = FALSE), collapse = "\n")

  expect_match(boxplot_single, "current_features_plot")
  expect_match(heatmap_single, "current_features_plot")
  expect_match(heatmap_single, "heatmap_title <- if")
  expect_match(change_dotplot_pair, "current_features_plot")
  expect_match(heatmap_pair, "current_features_plot")
  expect_match(indiv_boxplot_long, "current_features_plot")
  expect_match(indiv_spaghetti_long, "current_features_plot")
  expect_match(heatmap_long, "current_features_plot")
  expect_match(change_heatmap_long, "selected_features")
  expect_match(heatmap_long, "mStat_prepare_taxa_long_data")
  expect_match(heatmap_long, "mStat_summarize_mean_by_groups")
  expect_match(change_heatmap_long, "mStat_prepare_taxa_long_data")
  expect_match(change_heatmap_long, "mStat_summarize_mean_by_groups")
  expect_match(change_heatmap_long, "mStat_prepare_change_heatmap_long_matrix")
  expect_no_match(change_heatmap_long, "tidyr::spread")
  expect_match(heatmap_pair, "mStat_prepare_taxa_long_data")
  expect_match(heatmap_pair, "mStat_summarize_mean_by_groups")
  expect_match(spaghetti_long, "mStat_summarize_mean_by_groups")
  expect_match(heatmap_single, "heatmap_title <- if")

  expect_no_match(heatmap_long, "tidyr::gather")
  expect_no_match(heatmap_long, "tidyr::spread")
  expect_match(ma_single, "mStat_filter_test_result_features")
  expect_match(volcano_single, "mStat_filter_test_result_features")
  expect_match(trend_volcano_long, "mStat_filter_test_result_features")
  expect_match(volatility_volcano_long, "mStat_filter_test_result_features")
  expect_match(per_time_dotplot_long, "mStat_filter_test_result_features")

  expect_no_match(boxplot_single, 'features.plot <- names\\(sort')
  expect_no_match(heatmap_single, 'features.plot <- names\\(sort')
  expect_no_match(change_dotplot_pair, 'features.plot <- names\\(sort')
  expect_no_match(change_dotplot_pair, 'tidyr::spread')
  expect_no_match(change_dotplot_pair, 'tidyr::gather')
  expect_no_match(heatmap_pair, 'features.plot <- names\\(sort')
  expect_no_match(indiv_boxplot_long, 'features.plot <- names\\(sort')
  expect_no_match(indiv_spaghetti_long, 'features.plot <- names\\(sort')
})
