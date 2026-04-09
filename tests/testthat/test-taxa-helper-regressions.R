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
