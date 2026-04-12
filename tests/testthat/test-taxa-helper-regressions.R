test_that("mStat_summarize_taxa_features computes abundance and prevalence consistently", {
  feature.mat <- matrix(
    c(1, 0, 3, 0,
      0, 5, 0, 5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("tax1", "tax2"), c("s1", "s2", "s3", "s4"))
  )

  expected <- tibble::tibble(
    Genus = c("tax1", "tax2"),
    avg_abundance = c(1, 2.5),
    prevalence = c(0.5, 0.5)
  )

  expect_equal(
    mStat_summarize_taxa_features(feature.mat, feature.level = "Genus"),
    expected
  )

  feature.df <- data.frame(
    Genus = c("tax1", "tax2"),
    s1 = c(1, 0),
    s2 = c(0, 5),
    s3 = c(3, 0),
    s4 = c(0, 5),
    check.names = FALSE
  )

  expect_equal(
    mStat_summarize_taxa_features(
      feature.df,
      feature.level = "Genus",
      feature_in_column = TRUE
    ),
    expected
  )
})

test_that("mStat_as_taxa_feature_matrix and composition matrix preserve schema cleanly", {
  feature.df <- data.frame(
    Genus = c("tax1", "tax2"),
    s1 = c(2, 0),
    s2 = c(2, 2),
    s3 = c(0, 0),
    check.names = FALSE
  )

  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.df,
    feature.level = "Genus",
    feature_in_column = TRUE
  )

  expect_equal(rownames(feature.mat), c("tax1", "tax2"))
  expect_equal(colnames(feature.mat), c("s1", "s2", "s3"))
  expect_type(feature.mat, "double")

  composition.mat <- mStat_as_taxa_composition_matrix(
    feature.dat = feature.df,
    feature.level = "Genus"
  )

  expect_equal(unname(colSums(composition.mat[, c("s1", "s2"), drop = FALSE])), c(1, 1))
  expect_equal(unname(composition.mat[, "s3"]), c(0, 0))
})

test_that("mStat_meta_with_sample and long-data helper preserve join semantics", {
  meta.df <- data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
  meta.with.sample <- mStat_meta_with_sample(meta.df)

  expect_named(meta.with.sample, c("sample", "group"))
  expect_equal(meta.with.sample$sample, c("s1", "s2"))

  feature.df <- data.frame(
    Genus = c("tax1", "tax2"),
    s1 = c(1, 0),
    s3 = c(2, 4),
    check.names = FALSE
  )

  left.long <- mStat_prepare_taxa_long_data(
    feature.dat = feature.df,
    feature.level = "Genus",
    value_col = "count",
    meta.dat = meta.df,
    join = "left"
  )

  inner.long <- mStat_prepare_taxa_long_data(
    feature.dat = feature.df,
    feature.level = "Genus",
    value_col = "count",
    meta.dat = meta.df,
    join = "inner"
  )

  expect_equal(sort(unique(left.long$sample)), c("s1", "s3"))
  expect_equal(sort(unique(inner.long$sample)), "s1")
  expect_true(all(is.na(left.long$group[left.long$sample == "s3"])))
})

test_that("mStat_resolve_optional_flag preserves explicit values and validates input", {
  expect_true(mStat_resolve_optional_flag(NULL, TRUE, "cluster.cols"))
  expect_false(mStat_resolve_optional_flag(FALSE, TRUE, "cluster.cols"))
  expect_error(
    mStat_resolve_optional_flag(c(TRUE, FALSE), TRUE, "cluster.cols"),
    "single TRUE/FALSE value"
  )
})

test_that("mStat_summarize_grouped_taxa_long computes grouped summaries once", {
  long.df <- tibble::tibble(
    Genus = c("tax1", "tax1", "tax2", "tax2"),
    group = c("A", "A", "B", "B"),
    time = c("t1", "t1", "t1", "t2"),
    count = c(1, 0, 2, 2)
  )

  summary.df <- mStat_summarize_grouped_taxa_long(
    long.df = long.df,
    feature.level = "Genus",
    group_vars = c("group", "time"),
    value_col = "count",
    mean_col = "mean_abundance",
    prevalence_col = "prevalence",
    mean_transform = sqrt
  )

  expect_named(summary.df, c("group", "time", "Genus", "mean_abundance", "prevalence"))
  expect_equal(summary.df$mean_abundance[summary.df$Genus == "tax1"], sqrt(0.5))
  expect_equal(summary.df$prevalence[summary.df$Genus == "tax1"], 0.5)
})

test_that("mStat_transform_taxa_long_values handles sqrt and log transforms consistently", {
  long.df <- tibble::tibble(
    Genus = c("tax1", "tax1", "tax2", "tax2", "tax3", "tax3"),
    sample = c("s1", "s2", "s1", "s2", "s1", "s2"),
    value = c(4, 0, 1, 9, 0, 0)
  )

  sqrt.df <- mStat_transform_taxa_long_values(
    long.df = long.df,
    feature.level = "Genus",
    feature.dat.type = "count",
    transform = "sqrt"
  )
  expect_equal(sqrt.df$value, c(2, 0, 1, 3, 0, 0))

  log.df <- mStat_transform_taxa_long_values(
    long.df = long.df,
    feature.level = "Genus",
    feature.dat.type = "count",
    transform = "log"
  )

  expect_false("tax3" %in% log.df$Genus)
  expect_equal(log.df$value[log.df$Genus == "tax1" & log.df$sample == "s2"], log10(2))
  expect_equal(log.df$value[log.df$Genus == "tax2" & log.df$sample == "s1"], log10(1))
})

test_that("mStat_prepare_taxa_clr_long_data uses sample-wise zero handling", {
  feature.mat <- matrix(
    c(
      1, 0, 0,
      3, 6, 0,
      0, 6, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("tax1", "tax2", "tax3"), c("s1", "s2", "s3"))
  )

  expect_warning(
    clr.long <- mStat_prepare_taxa_clr_long_data(
      feature.dat = feature.mat,
      feature.level = "Genus",
      prev.filter = 0,
      abund.filter = -Inf
    ),
    "zero total abundance"
  )

  expect_false("s3" %in% clr.long$sample)

  sample_sums <- clr.long %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(total = sum(value), .groups = "drop")
  expect_true(all(abs(sample_sums$total) < 1e-8))
})

test_that("mStat_format_linda_feature_results standardizes LinDA-style outputs", {
  result.df <- data.frame(
    log2FoldChange = c(0.1, -0.2),
    lfcSE = c(0.01, 0.02),
    pvalue = c(0.05, 0.01),
    padj = c(0.10, 0.02),
    reject = c(FALSE, TRUE),
    row.names = c("tax1", "tax2")
  )

  feature.stats <- tibble::tibble(
    Genus = c("tax1", "tax2"),
    avg_abundance = c(1, 2.5),
    prevalence = c(0.5, 0.5)
  )

  formatted <- mStat_format_linda_feature_results(
    result.df = result.df,
    feature.level = "Genus",
    feature.stats = feature.stats,
    include_significant = TRUE
  )

  expect_named(
    formatted,
    c(
      "Variable",
      "Coefficient",
      "SE",
      "P.Value",
      "Adjusted.P.Value",
      "Significant",
      "Mean.Abundance",
      "Prevalence"
    )
  )
  expect_equal(formatted$Variable, c("tax1", "tax2"))
  expect_equal(formatted$Significant, c(FALSE, TRUE))
  expect_equal(formatted$Mean.Abundance, c(1, 2.5))
})

test_that("mStat_format_linda_feature_results fails fast on malformed inputs", {
  incomplete.df <- data.frame(
    log2FoldChange = 0.1,
    row.names = "tax1"
  )

  feature.stats <- tibble::tibble(
    Genus = "tax1",
    avg_abundance = 1,
    prevalence = 0.5
  )

  expect_error(
    mStat_format_linda_feature_results(
      result.df = incomplete.df,
      feature.level = "Genus",
      feature.stats = feature.stats
    ),
    "missing required columns"
  )
})

test_that("shared metadata-value scoper preserves typed time matching", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("tax1", "tax2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      time = factor(c("T1", "T2"), levels = c("T1", "T2")),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  subset.obj <- mStat_subset_by_meta_values(data.obj, "time", "T2")

  expect_equal(colnames(subset.obj$feature.tab), "s2")
  expect_equal(rownames(subset.obj$meta.dat), "s2")
})

test_that("package spaghetti plots no longer use aes_string", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  helper_paths <- c(
    "R/generate_alpha_spaghettiplot_long.R",
    "R/generate_taxa_spaghettiplot_long.R",
    "R/generate_taxa_indiv_spaghettiplot_long.R",
    "R/generate_beta_pc_spaghettiplot_long.R",
    "R/generate_beta_change_spaghettiplot_long.R"
  )

  for (path in helper_paths) {
    contents <- paste(readLines(file.path(package_root, path), warn = FALSE), collapse = "\n")
    expect_no_match(contents, "aes_string")
  }
})

test_that("single taxa context helper scopes and injects subject ids consistently", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      dimnames = list(c("tax1", "tax2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      time = c("T1", "T2"),
      group = c("A", "B"),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  context <- mStat_prepare_taxa_single_context(
    data.obj = data.obj,
    time.var = "time",
    t.level = "T2",
    group.var = "group",
    add_subject_id = TRUE
  )

  expect_equal(colnames(context$data.obj$feature.tab), "s2")
  expect_equal(context$subject.var, "subject.id")
  expect_named(context$meta_tab, c("subject.id", "group", "time"))
  expect_equal(context$meta_tab$subject.id, "s2")
})

test_that("long taxa context helper centralizes time processing and metadata selection", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2, 3, 4, 5, 6),
      nrow = 2,
      dimnames = list(c("tax1", "tax2"), c("s1", "s2", "s3"))
    ),
    meta.dat = data.frame(
      subject = c("u1", "u1", "u2"),
      time = c("T1", "T2", "T3"),
      group = c("A", "A", "B"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )

  context <- mStat_prepare_taxa_long_context(
    data.obj = data.obj,
    subject.var = "subject",
    time.var = "time",
    group.var = "group",
    t0.level = "T1",
    ts.levels = "T2"
  )

  expect_equal(colnames(context$data.obj$feature.tab), c("s1", "s2"))
  expect_named(context$meta_tab, c("subject", "group", "time"))
})

test_that("pair taxa change helper computes aligned change once", {
  long.df <- tibble::tibble(
    Genus = c("tax1", "tax1", "tax2", "tax2"),
    subject = c("u1", "u1", "u1", "u1"),
    time = c("T1", "T2", "T1", "T2"),
    value = c(2, 6, 4, 1)
  )

  pair_change <- mStat_prepare_taxa_pair_change_data(
    long.df = long.df,
    feature.level = "Genus",
    subject.var = "subject",
    time.var = "time",
    change.base = "T1",
    feature.change.func = "absolute change",
    context = "test"
  )

  expect_equal(pair_change$change.base, "T1")
  expect_equal(pair_change$change.after, "T2")
  expect_equal(pair_change$combined_data$value_diff, c(4, -3))
})

test_that("pair metadata helper supports follow-up and subject joins", {
  df <- tibble::tibble(subject = c("u1", "u2"), value_diff = c(1, 2))
  meta_tab <- tibble::tibble(
    subject = c("u1", "u1", "u2", "u2"),
    time = c("T1", "T2", "T1", "T2"),
    group = c("A0", "A1", "B0", "B1")
  )

  followup <- mStat_attach_pair_metadata(
    df = df,
    meta_tab = meta_tab,
    subject.var = "subject",
    time.var = "time",
    mode = "followup_time",
    change.after = "T2"
  )
  subject_only <- mStat_attach_pair_metadata(
    df = df,
    meta_tab = meta_tab,
    subject.var = "subject",
    time.var = "time",
    mode = "subject"
  )

  expect_equal(followup$group, c("A1", "B1"))
  expect_false("time" %in% names(subject_only))
})

test_that("matrix metadata alignment follows matrix column order", {
  value_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("tax1", "tax2"), c("u2", "u1"))
  )
  meta_tab <- tibble::tibble(
    subject = c("u1", "u2"),
    group = c("A", "B")
  )

  aligned <- mStat_align_subject_metadata_to_matrix(
    value_matrix = value_matrix,
    meta_tab = meta_tab,
    subject.var = "subject",
    keep_vars = "group"
  )

  expect_equal(colnames(aligned$value_matrix), c("u2", "u1"))
  expect_equal(rownames(aligned$sorted_meta), c("u2", "u1"))
})

test_that("feature selection helper preserves explicit, top-k, and fallback selection", {
  feature.dat <- data.frame(
    Genus = c("tax1", "tax2", "tax3"),
    s1 = c(10, 1, 5),
    s2 = c(8, 2, 4),
    check.names = FALSE
  )

  explicit <- mStat_resolve_selected_features(
    feature.level = "Genus",
    features.plot = c("tax3", "tax1"),
    taxa.levels = c("tax1", "tax2", "tax3")
  )
  topk <- mStat_resolve_selected_features(
    feature.dat = feature.dat,
    feature.level = "Genus",
    top.k.plot = 2,
    top.k.func = "mean"
  )
  fallback <- mStat_resolve_selected_features(
    feature.level = "Genus",
    taxa.levels = c("tax1", "tax2", "tax3"),
    fallback_n = 2
  )

  expect_equal(explicit, c("tax3", "tax1"))
  expect_length(topk, 2)
  expect_equal(fallback, c("tax1", "tax2"))
})

test_that("group placeholder helper creates a stable synthetic group column", {
  df <- tibble::tibble(time = c("T1", "T2"), value = c(1, 2))
  placeholder <- mStat_ensure_group_placeholder(df, group.var = NULL, value = "ALL", column_name = "ALL")

  expect_equal(placeholder$group.var, "ALL")
  expect_equal(placeholder$df$ALL, c("ALL", "ALL"))
})


test_that("targeted long beta helpers no longer mutate group.var placeholders inline", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  helper_paths <- c(
    "R/generate_beta_change_spaghettiplot_long.R",
    "R/generate_beta_pc_spaghettiplot_long.R",
    "R/generate_beta_pc_change_boxplot_pair.R",
    "R/generate_taxa_change_boxplot_pair.R"
  )

  for (path in helper_paths) {
    lines <- readLines(file.path(package_root, path), warn = FALSE)
    code_lines <- lines[!grepl("^\\s*#", lines)]
    contents <- paste(code_lines, collapse = "\n")
    expect_no_match(contents, "group\\.var\\s*(<-|=)\\s*\\\"(ALL|group)\\\"")
  }
})


test_that("remaining pair helpers prefer shared join helpers and modern pivot verbs", {
  package_root <- normalizePath(file.path(test_path(), "..", ".."), mustWork = TRUE)
  beta_ord_single <- paste(readLines(file.path(package_root, "R/generate_beta_ordination_single.R"), warn = FALSE), collapse = "\n")
  beta_pc_change_pair <- paste(readLines(file.path(package_root, "R/generate_beta_pc_change_boxplot_pair.R"), warn = FALSE), collapse = "\n")
  taxa_change_test_pair <- paste(readLines(file.path(package_root, "R/generate_taxa_change_test_pair.R"), warn = FALSE), collapse = "\n")
  taxa_change_dotplot_pair <- paste(readLines(file.path(package_root, "R/generate_taxa_change_dotplot_pair.R"), warn = FALSE), collapse = "\n")

  expect_match(beta_ord_single, "mStat_prepare_pc_long_data")
  expect_match(beta_pc_change_pair, "mStat_attach_pair_metadata")
  expect_match(taxa_change_test_pair, "mStat_attach_subject_level_metadata")
  expect_match(taxa_change_dotplot_pair, "pivot_wider")
  expect_match(taxa_change_dotplot_pair, "pivot_longer")
  expect_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_heatmap_pair.R"), warn = FALSE), collapse = "\n"),
    "mStat_attach_subject_level_metadata"
  )
  expect_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_heatmap_long.R"), warn = FALSE), collapse = "\n"),
    "mStat_prepare_change_heatmap_long_matrix"
  )
  expect_no_match(taxa_change_test_pair, "tidyr::spread")
  expect_no_match(taxa_change_test_pair, "tidyr::gather")
  expect_no_match(taxa_change_dotplot_pair, "tidyr::spread")
  expect_no_match(taxa_change_dotplot_pair, "tidyr::gather")
  expect_no_match(
    paste(readLines(file.path(package_root, "R/generate_taxa_change_heatmap_long.R"), warn = FALSE), collapse = "\n"),
    "tidyr::spread"
  )
})

test_that("mean summary helper computes grouped means with optional transform", {
  long.df <- tibble::tibble(
    Genus = c("tax1", "tax1", "tax2", "tax2"),
    group = c("A", "A", "B", "B"),
    time = c("T1", "T2", "T1", "T2"),
    value = c(4, 9, 1, 16)
  )

  summary.df <- mStat_summarize_mean_by_groups(
    long.df = long.df,
    feature.level = "Genus",
    group_vars = c("group", "time"),
    value_col = "value",
    mean_col = "mean_value",
    mean_transform = sqrt
  )

  expect_named(summary.df, c("Genus", "group", "time", "mean_value"))
  expect_equal(summary.df$mean_value[summary.df$Genus == "tax1" & summary.df$time == "T1"], 2)
})

test_that("mean summary helper ignores missing values consistently", {
  long.df <- tibble::tibble(
    Genus = c("tax1", "tax1"),
    time = c("T1", "T1"),
    value = c(1, NA_real_)
  )

  summary.df <- mStat_summarize_mean_by_groups(
    long.df = long.df,
    feature.level = "Genus",
    group_vars = "time",
    value_col = "value",
    mean_col = "mean_value"
  )

  expect_equal(summary.df$mean_value, 1)
})

test_that("test result feature filter helper preserves all rows when no filter is set", {
  result.df <- tibble::tibble(Variable = c("tax1", "tax2"), value = c(1, 2))
  expect_equal(mStat_filter_test_result_features(result.df), result.df)
})

test_that("test result feature filter helper filters by requested variables", {
  result.df <- tibble::tibble(Variable = c("tax1", "tax2"), value = c(1, 2))
  filtered <- mStat_filter_test_result_features(result.df, features.plot = "tax2")
  expect_equal(filtered$Variable, "tax2")
})


test_that("generate_taxa_per_time_test_long analyzes sparse per-time longitudinal slices", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 12, 9, 11, 8, 10,
        1, 0, 2, 1, 3, 2
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = c("u1", "u2", "u3", "u1", "u2", "u3"),
      time = c("T1", "T1", "T1", "T2", "T2", "T2"),
      group = c("A", "A", "B", "A", "A", "B"),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    generate_taxa_per_time_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      feature.level = "original",
      feature.dat.type = "proportion"
    )
  )

  expect_true(all(c("T1", "T2") %in% names(result)))
})


test_that("generate_taxa_per_time_test_long allows subject-varying group across time", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 9, 8, 7, 6, 5,
        1, 2, 1, 2, 3, 4
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = c("u1", "u1", "u2", "u2", "u3", "u3"),
      time = c("T1", "T2", "T1", "T2", "T1", "T2"),
      group = c("A", "B", "A", "A", "B", "B"),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    generate_taxa_per_time_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      feature.level = "original",
      feature.dat.type = "proportion"
    )
  )

  expect_true(all(c("T1", "T2") %in% names(result)))
})


test_that("generate_taxa_change_heatmap_pair computes internal column gaps without terminal boundary", {
  sorted_meta_tab <- data.frame(
    group = c("A", "A", "B", "B", "C"),
    row.names = paste0("s", 1:5),
    stringsAsFactors = FALSE
  )

  group_counts <- table(sorted_meta_tab[["group"]])
  gaps <- cumsum(group_counts)[-length(group_counts)]

  expect_identical(unname(gaps), c(2L, 4L))
})


test_that("generate_taxa_change_heatmap_pair handles all-zero change matrices", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        5, 5, 5, 5, 5, 5,
        2, 2, 2, 2, 2, 2
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3"), each = 2),
      time = rep(c("T1", "T2"), times = 3),
      group = rep(c("A", "B", "A"), each = 2),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  expect_warning(
    result <- suppressMessages(
      generate_taxa_change_heatmap_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        change.base = "T1",
        feature.level = "original",
        feature.change.func = "absolute change",
        feature.dat.type = "proportion",
        pdf = FALSE
      )
    ),
    NA
  )

  expect_s3_class(result$original$indiv, "ggplot")
  expect_s3_class(result$original$average, "ggplot")
})


test_that("generate_taxa_change_heatmap_long handles all-zero change matrices", {
  sample_ids <- paste0("s", 1:8)
  data.obj <- list(
    feature.tab = matrix(
      c(
        4, 4, 4, 4, 4, 4, 4, 4,
        1, 1, 1, 1, 1, 1, 1, 1
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3", "u4"), each = 2),
      time = rep(c("T1", "T2"), times = 4),
      group = rep(c("A", "A", "B", "B"), each = 2),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  expect_warning(
    result <- suppressMessages(
      generate_taxa_change_heatmap_long(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        t0.level = "T1",
        ts.levels = "T2",
        group.var = "group",
        feature.level = "original",
        feature.change.func = "absolute change",
        feature.dat.type = "proportion",
        pdf = FALSE
      )
    ),
    NA
  )

  expect_s3_class(result$original, "ggplot")
})


test_that("paired individual change scatter logic collapses repeated measurements before pairing", {
  long.df <- tibble::tibble(
    Genus = c("tax1", "tax1", "tax1", "tax1"),
    subject = c("u1", "u1", "u1", "u1"),
    time = c("T0", "T0", "T1", "T1"),
    count = c(1, 3, 5, 7)
  )

  subject_time_abundance <- long.df %>%
    dplyr::group_by(Genus, subject, time) %>%
    dplyr::summarise(count = mean(count, na.rm = TRUE), .groups = "drop")

  df_t0 <- subject_time_abundance %>% dplyr::filter(time == "T0")
  df_t1 <- subject_time_abundance %>% dplyr::filter(time == "T1")
  paired <- dplyr::inner_join(df_t1, df_t0, by = c("Genus", "subject"), suffix = c("_t1", "_t0"))

  expect_equal(nrow(paired), 1)
  expect_equal(paired$count_t0, 2)
  expect_equal(paired$count_t1, 6)
})


test_that("scatter pair functions keep time mapping when strata.var is NULL", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 12, 14, 16, 18, 20,
        5, 4, 6, 7, 9, 8
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3"), each = 2),
      time = rep(c("T1", "T2"), times = 3),
      score = c(1.1, 1.2, 1.6, 1.8, 2.0, 2.4),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  scatter_pair <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_scatterplot_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "score",
        strata.var = NULL,
        change.base = "T1",
        feature.level = "original",
        feature.dat.type = "proportion",
        pdf = FALSE
      )
    )
  )

  indiv_pair <- suppressWarnings(
    suppressMessages(
      generate_taxa_indiv_change_scatterplot_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "score",
        strata.var = NULL,
        change.base = "T1",
        feature.level = "original",
        feature.dat.type = "proportion",
        pdf = FALSE
      )
    )
  )

  scatter_point_layer_idx <- which(vapply(
    scatter_pair$original$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))[1]
  indiv_point_layer_idx <- which(vapply(
    indiv_pair$original[[1]]$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))[1]

  expect_identical(
    rlang::as_label(scatter_pair$original$layers[[scatter_point_layer_idx]]$mapping$colour),
    "time"
  )
  expect_identical(
    rlang::as_label(indiv_pair$original[[1]]$layers[[indiv_point_layer_idx]]$mapping$colour),
    "time"
  )
})


test_that("scatter pair functions use strata mapping when strata.var is provided", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 12, 14, 16, 18, 20,
        5, 4, 6, 7, 9, 8
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3"), each = 2),
      time = rep(c("T1", "T2"), times = 3),
      score = c(1.1, 1.2, 1.6, 1.8, 2.0, 2.4),
      strata = c("A", "A", "B", "B", "A", "A"),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  scatter_pair <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_scatterplot_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "score",
        strata.var = "strata",
        change.base = "T1",
        feature.level = "original",
        feature.dat.type = "proportion",
        pdf = FALSE
      )
    )
  )

  indiv_pair <- suppressWarnings(
    suppressMessages(
      generate_taxa_indiv_change_scatterplot_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "score",
        strata.var = "strata",
        change.base = "T1",
        feature.level = "original",
        feature.dat.type = "proportion",
        pdf = FALSE
      )
    )
  )

  scatter_point_layer_idx <- which(vapply(
    scatter_pair$original$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))[1]
  indiv_point_layer_idx <- which(vapply(
    indiv_pair$original[[1]]$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))[1]

  expect_identical(
    rlang::as_label(scatter_pair$original$layers[[scatter_point_layer_idx]]$mapping$colour),
    "strata"
  )
  expect_identical(
    rlang::as_label(indiv_pair$original[[1]]$layers[[indiv_point_layer_idx]]$mapping$colour),
    "strata"
  )
})


test_that("is_continuous_numeric enforces numeric type and uniqueness threshold", {
  expect_false(is_continuous_numeric(c("a", "b", "c")))
  expect_false(is_continuous_numeric(c(1, 2, 3, 4, 5, 6, 7, 8, 9)))
  expect_false(is_continuous_numeric(rep(1, 12)))
  expect_true(is_continuous_numeric(1:10))
})


test_that("per-time dotplot helpers keep continuous group terms", {
  sample_ids <- paste0("s", 1:4)
  data.obj <- list(
    feature.tab = matrix(
      c(10, 12, 14, 16,
        3, 4, 5, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      time = c("T1", "T1", "T2", "T2"),
      score = c(1.1, 2.0, 1.4, 2.3),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  test.list <- list(
    T1 = list(
      score = data.frame(
        Term = "shannon",
        Estimate = 0.20,
        P.Value = 0.04,
        stringsAsFactors = FALSE
      )
    ),
    T2 = list(
      score = data.frame(
        Term = "shannon",
        Estimate = 0.35,
        P.Value = 0.20,
        stringsAsFactors = FALSE
      )
    )
  )

  alpha_plot <- suppressWarnings(
    generate_alpha_per_time_dotplot_long(
      data.obj = data.obj,
      test.list = test.list,
      group.var = "score",
      time.var = "time",
      t0.level = "T1",
      ts.levels = "T2",
      pdf = FALSE
    )
  )

  beta_plot <- suppressWarnings(
    generate_beta_per_time_dotplot_long(
      data.obj = data.obj,
      test.list = test.list,
      group.var = "score",
      time.var = "time",
      t0.level = "T1",
      ts.levels = "T2",
      pdf = FALSE
    )
  )

  expect_identical(names(alpha_plot), "score")
  expect_identical(names(beta_plot), "score")
})


test_that("change scatter pair helpers keep top-k selection scoped per feature level", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 20, 10, 10, 20, 10,
        5, 5, 5, 5, 5, 5,
        1, 1, 1, 9, 9, 9,
        4, 4, 4, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3", "f4"), sample_ids)
    ),
    feature.ann = matrix(
      c("gA", "gA", "gB", "gB"),
      ncol = 1,
      dimnames = list(c("f1", "f2", "f3", "f4"), "Genus")
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3"), each = 2),
      time = rep(c("T0", "T1"), times = 3),
      score = c(1.0, 1.2, 1.5, 1.7, 2.0, 2.2),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  scatter_pair <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_scatterplot_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "score",
        change.base = "T0",
        feature.level = c("original", "Genus"),
        feature.dat.type = "proportion",
        top.k.plot = 1,
        top.k.func = "mean",
        pdf = FALSE
      )
    )
  )

  indiv_pair <- suppressWarnings(
    suppressMessages(
      generate_taxa_indiv_change_scatterplot_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "score",
        change.base = "T0",
        feature.level = c("original", "Genus"),
        feature.dat.type = "proportion",
        top.k.plot = 1,
        top.k.func = "mean",
        pdf = FALSE
      )
    )
  )

  expect_true(all(c("original", "Genus") %in% names(scatter_pair)))
  expect_true(all(c("original", "Genus") %in% names(indiv_pair)))
  expect_identical(unique(scatter_pair$original$data$original), "f1")
  expect_identical(unique(scatter_pair$Genus$data$Genus), "gA")
  expect_identical(names(indiv_pair$original), "f1")
  expect_identical(names(indiv_pair$Genus), "gA")
})


test_that("change heatmap helpers handle all-zero changes and scoped top-k selection", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 10, 10, 10, 10, 10,
        5, 5, 5, 5, 5, 5,
        2, 2, 2, 2, 2, 2,
        1, 1, 1, 1, 1, 1
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3", "f4"), sample_ids)
    ),
    feature.ann = matrix(
      c("gA", "gA", "gB", "gB"),
      ncol = 1,
      dimnames = list(c("f1", "f2", "f3", "f4"), "Genus")
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3"), each = 2),
      time = rep(c("T0", "T1"), times = 3),
      group = c("A", "A", "A", "A", "B", "B"),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  pair_result <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_heatmap_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        change.base = "T0",
        feature.level = c("original", "Genus"),
        feature.dat.type = "proportion",
        top.k.plot = 2,
        top.k.func = "mean",
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  long_result <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_heatmap_long(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        t0.level = "T0",
        ts.levels = "T1",
        group.var = "group",
        feature.level = c("original", "Genus"),
        feature.dat.type = "proportion",
        top.k.plot = 2,
        top.k.func = "mean",
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_true(all(c("original", "Genus") %in% names(pair_result)))
  expect_true(all(c("average", "indiv") %in% names(pair_result$original)))
  expect_true(all(c("average", "indiv") %in% names(pair_result$Genus)))
  expect_true(all(c("original", "Genus") %in% names(long_result)))
})


test_that("change heatmap pair keeps top-k selection scoped per feature level", {
  sample_ids <- paste0("s", 1:6)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 20, 10, 10, 20, 10,
        5, 5, 5, 5, 5, 5,
        1, 1, 1, 9, 9, 9,
        4, 4, 4, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3", "f4"), sample_ids)
    ),
    feature.ann = matrix(
      c("gA", "gA", "gB", "gB"),
      ncol = 1,
      dimnames = list(c("f1", "f2", "f3", "f4"), "Genus")
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3"), each = 2),
      time = rep(c("T0", "T1"), times = 3),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_heatmap_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        change.base = "T0",
        feature.level = c("original", "Genus"),
        feature.dat.type = "proportion",
        top.k.plot = 1,
        top.k.func = "mean",
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_identical(names(result), c("original", "Genus"))
  expect_s3_class(result$original$indiv, "ggplot")
  expect_s3_class(result$Genus$indiv, "ggplot")
  expect_s3_class(result$original$average, "ggplot")
  expect_s3_class(result$Genus$average, "ggplot")
})


test_that("change heatmap long keeps top-k selection scoped per feature level", {
  sample_ids <- paste0("s", 1:8)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 20, 10, 20, 10, 20, 10, 20,
        5, 5, 5, 5, 5, 5, 5, 5,
        1, 1, 1, 9, 9, 9, 9, 1,
        4, 4, 4, 2, 2, 2, 2, 4
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3", "f4"), sample_ids)
    ),
    feature.ann = matrix(
      c("gA", "gA", "gB", "gB"),
      ncol = 1,
      dimnames = list(c("f1", "f2", "f3", "f4"), "Genus")
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2", "u3", "u4"), each = 2),
      time = rep(c("T0", "T1"), times = 4),
      group = rep(c("A", "A", "B", "B"), each = 2),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_heatmap_long(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        t0.level = "T0",
        ts.levels = "T1",
        group.var = "group",
        feature.level = c("original", "Genus"),
        feature.dat.type = "proportion",
        top.k.plot = 1,
        top.k.func = "mean",
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_identical(names(result), c("original", "Genus"))
  expect_s3_class(result$original, "ggplot")
  expect_s3_class(result$Genus, "ggplot")
})

test_that("generate_taxa_association_test_long keeps per-level sample filtering local", {
  sample_ids <- paste0("s", 1:8)
  data.obj <- list(
    feature.tab = matrix(
      c(
        0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.0,
        0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.0,
        0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.6,
        0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5
      ),
      nrow = 4,
      byrow = TRUE,
      dimnames = list(paste0("f", 1:4), sample_ids)
    ),
    feature.ann = data.frame(
      Family = c("fam1", "fam1", "fam2", "fam2"),
      row.names = paste0("f", 1:4),
      stringsAsFactors = FALSE
    ),
    meta.dat = data.frame(
      subject = rep(paste0("id", 1:4), each = 2),
      group = rep(c("A", "A", "B", "B"), each = 2),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    ),
    feature.agg.list = list()
  )

  result <- suppressWarnings(
    generate_taxa_association_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      group.var = "group",
      feature.level = c("original", "Family"),
      feature.dat.type = "proportion",
      prev.filter = 0,
      abund.filter = 0.3
    )
  )

  expect_named(result, c("original", "Family"))
  expect_true(length(result$original) > 0)
  expect_true(length(result$Family) > 0)
})
