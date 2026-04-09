test_that("mStat_prepare_precomputed_beta_context preserves metadata-only time processing", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 1, 3, 2,
        0, 5, 1, 4,
        2, 2, 2, 2),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3", "s4"))
    ),
    meta.dat = data.frame(
      time = c("T0", "T1", "T2", "T2"),
      group = c("A", "A", "B", "B"),
      row.names = c("s1", "s2", "s3", "s4"),
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(data.obj, dist.name = "BC"))

  prepared <- mStat_prepare_precomputed_beta_context(
    dist.obj = dist.obj,
    dist.name = "BC",
    data.obj = NULL,
    time.var = "time",
    t0.level = "T0",
    ts.levels = "T1",
    process_time = TRUE
  )

  expect_identical(rownames(prepared$data.obj$meta.dat), c("s1", "s2"))
  expect_identical(mStat_get_dist_labels(prepared$dist.obj$BC), c("s1", "s2"))

  extracted <- mStat_extract_dist_metadata(prepared$dist.obj, "BC", vars = c("time", "group"))
  expect_identical(rownames(extracted), c("s1", "s2"))
  expect_identical(as.character(extracted$time), c("T0", "T1"))
})

test_that("mStat_prepare_precomputed_beta_context reorders pc.obj to distance labels", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 1, 3,
        2, 4, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    ),
    meta.dat = data.frame(
      group = c("A", "A", "B"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(data.obj, dist.name = "BC"))
  pc.obj <- list(
    BC = list(
      points = matrix(
        c(3, 30,
          1, 10,
          2, 20),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("s3", "s1", "s2"), c("PC1", "PC2"))
      )
    )
  )

  prepared <- mStat_prepare_precomputed_beta_context(
    dist.obj = dist.obj,
    dist.name = "BC",
    pc.obj = pc.obj,
    data.obj = NULL,
    required_pc_axes = 2
  )

  expect_identical(
    rownames(prepared$pc.obj$BC$points),
    mStat_get_dist_labels(prepared$dist.obj$BC)
  )
  expect_equal(unname(prepared$pc.obj$BC$points[, "PC1"]), c(1, 2, 3))
})

test_that("mStat_prepare_precomputed_beta_context fails fast on malformed pc.obj", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 1, 3,
        2, 4, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    ),
    meta.dat = data.frame(
      group = c("A", "A", "B"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(data.obj, dist.name = "BC"))

  wrong_labels_pc <- list(
    BC = list(
      points = matrix(
        c(1, 10,
          2, 20,
          4, 40),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("s1", "s2", "s4"), c("PC1", "PC2"))
      )
    )
  )

  expect_error(
    mStat_prepare_precomputed_beta_context(
      dist.obj = dist.obj,
      dist.name = "BC",
      pc.obj = wrong_labels_pc,
      data.obj = NULL,
      required_pc_axes = 2
    ),
    "sample labels must match"
  )

  too_short_pc <- list(
    BC = list(
      points = matrix(
        c(1, 2, 3),
        ncol = 1,
        dimnames = list(c("s1", "s2", "s3"), "PC1")
      )
    )
  )

  expect_error(
    mStat_prepare_precomputed_beta_context(
      dist.obj = dist.obj,
      dist.name = "BC",
      pc.obj = too_short_pc,
      data.obj = NULL,
      required_pc_axes = 2
    ),
    "has only 1 axis/axes"
  )
})

test_that("mStat_extract_dist_metadata fails fast when metadata is unavailable", {
  raw_dist <- list(BC = stats::dist(matrix(c(1, 2, 3, 4), nrow = 2)))

  expect_error(
    mStat_extract_dist_metadata(raw_dist, "BC", vars = "group"),
    "Metadata is required"
  )
})

test_that("mStat_prepare_alpha_inputs keeps alpha.obj aligned after time subsetting", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 5, 5,
        1, 2, 3),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
    ),
    meta.dat = data.frame(
      time = c("T1", "T1", "T2"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )

  alpha.obj <- list(
    shannon = data.frame(
      shannon = c(1, 2, 99),
      row.names = c("s1", "s2", "s3")
    )
  )

  prepared <- mStat_prepare_alpha_inputs(
    data.obj = data.obj,
    alpha.obj = alpha.obj,
    alpha.name = "shannon",
    time.var = "time",
    t.level = "T1"
  )

  expect_identical(rownames(prepared$data.obj$meta.dat), c("s1", "s2"))
  expect_identical(rownames(prepared$alpha.obj$shannon), c("s1", "s2"))
  expect_equal(prepared$alpha.obj$shannon$shannon, c(1, 2))
})

test_that("mStat_prepare_alpha_inputs validates precomputed alpha.obj sample labels", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 5,
        1, 2),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      time = c("T1", "T1"),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  alpha.obj <- list(
    shannon = data.frame(
      shannon = 1,
      row.names = "s1"
    )
  )

  expect_error(
    mStat_prepare_alpha_inputs(
      data.obj = data.obj,
      alpha.obj = alpha.obj,
      alpha.name = "shannon"
    ),
    "missing the following samples: s2"
  )
})

test_that("winsorization uses one shared implementation", {
  Y <- matrix(
    c(1, 100,
      2, 50),
    nrow = 2,
    byrow = TRUE
  )

  shared <- winsor_feature_table(Y, 0.5, "proportion")
  alias <- winsor.fun(Y, 0.5, "proportion")

  expect_identical(shared, alias)
  expect_true(all(shared <= matrix(c(50.5, 50.5, 26, 26), nrow = 2, byrow = TRUE)))
})

test_that("generate_beta_ordination_single works with precomputed dist metadata only", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 8, 7, 0, 1, 2,
        0, 1, 2, 9, 8, 7,
        3, 3, 3, 3, 3, 3),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3", "s4", "s5", "s6"))
    ),
    meta.dat = data.frame(
      time = c("T1", "T1", "T1", "T2", "T2", "T2"),
      row.names = c("s1", "s2", "s3", "s4", "s5", "s6"),
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(data.obj, dist.name = "BC"))

  result <- suppressWarnings(
    generate_beta_ordination_single(
      data.obj = NULL,
      dist.obj = dist.obj,
      time.var = "time",
      t.level = "T1",
      group.var = NULL,
      dist.name = "BC",
      pdf = FALSE
    )
  )

  expect_true("BC" %in% names(result))
})

test_that("generate_beta_change_test_pair works with precomputed dist metadata only", {
  samples <- paste0("s", 1:8)
  data.obj <- list(
    feature.tab = matrix(
      c(10, 9, 8, 7, 1, 2, 3, 4,
        1, 2, 3, 4, 10, 9, 8, 7,
        5, 5, 5, 5, 5, 5, 5, 5),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), samples)
    ),
    meta.dat = data.frame(
      subject = rep(paste0("id", 1:4), each = 2),
      time = rep(c("T0", "T1"), times = 4),
      group = rep(c("A", "A", "B", "B"), each = 2),
      row.names = samples,
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(data.obj, dist.name = "BC"))

  result <- suppressWarnings(
    generate_beta_change_test_pair(
      data.obj = NULL,
      dist.obj = dist.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      change.base = "T0",
      dist.name = "BC"
    )
  )

  expect_true("BC" %in% names(result))
  expect_true(any(result$BC$Term == "groupB"))
})

test_that("test_alpha routes longitudinal difference requests to change-per-time implementation", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 8, 6, 4, 2, 1, 3, 5, 7,
        1, 2, 3, 4, 5, 6, 7, 8, 9),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"))
    ),
    meta.dat = data.frame(
      subject = rep(c("id1", "id2", "id3"), each = 3),
      time = rep(c("T0", "T1", "T2"), times = 3),
      group = c("A", "A", "A", "B", "B", "B", "A", "A", "A"),
      row.names = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9"),
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    test_alpha(
      data.obj = data.obj,
      test.type = "difference",
      alpha.name = "shannon",
      subject.var = "subject",
      time.var = "time",
      group.var = "group"
    )
  )

  expect_identical(attr(result, "design"), "long")
  expect_identical(names(result), c("T1", "T2"))
})

test_that("test_beta routes longitudinal difference requests to change-per-time implementation", {
  samples <- paste0("s", 1:9)
  data.obj <- list(
    feature.tab = matrix(
      c(10, 9, 8, 1, 2, 3, 7, 6, 5,
        1, 2, 3, 10, 9, 8, 4, 5, 6,
        5, 5, 5, 5, 5, 5, 5, 5, 5),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), samples)
    ),
    meta.dat = data.frame(
      subject = rep(c("id1", "id2", "id3"), each = 3),
      time = rep(c("T0", "T1", "T2"), times = 3),
      group = c("A", "A", "A", "B", "B", "B", "A", "A", "A"),
      row.names = samples,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    test_beta(
      data.obj = data.obj,
      test.type = "difference",
      dist.name = "BC",
      subject.var = "subject",
      time.var = "time",
      group.var = "group"
    )
  )

  expect_identical(attr(result, "design"), "long")
  expect_identical(names(result), c("T1", "T2"))
})

test_that("longitudinal change and volatility helpers reject subject-varying group assignments", {
  alpha_data.obj <- list(
    feature.tab = matrix(
      c(10, 9, 8, 7,
        1, 2, 3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4"))
    ),
    meta.dat = data.frame(
      subject = c("id1", "id1", "id2", "id2"),
      time = c("T0", "T1", "T0", "T1"),
      group = c("A", "B", "A", "A"),
      row.names = c("s1", "s2", "s3", "s4"),
      stringsAsFactors = FALSE
    )
  )

  beta_data.obj <- list(
    feature.tab = matrix(
      c(10, 9, 8, 7,
        1, 2, 3, 4,
        5, 5, 5, 5),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3", "s4"))
    ),
    meta.dat = alpha_data.obj$meta.dat
  )

  expect_error(
    suppressWarnings(
      generate_alpha_change_test_pair(
        data.obj = alpha_data.obj,
        alpha.name = "shannon",
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        change.base = "T0"
      )
    ),
    "must be constant within each subject"
  )

  expect_error(
    suppressWarnings(
      generate_alpha_volatility_test_long(
        data.obj = alpha_data.obj,
        alpha.name = "shannon",
        time.var = "time",
        subject.var = "subject",
        group.var = "group"
      )
    ),
    "must be constant within each subject"
  )

  expect_error(
    suppressWarnings(
      generate_beta_volatility_test_long(
        data.obj = beta_data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        dist.name = "BC"
      )
    ),
    "must be constant within each subject"
  )

  expect_error(
    suppressWarnings(
      generate_beta_pc_volatility_test_long(
        data.obj = beta_data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        dist.name = "BC"
      )
    ),
    "must be constant within each subject"
  )
})

test_that("per-time helpers respect ordered character time values", {
  data.obj <- list(
    feature.tab = matrix(
      c(9, 8, 7, 6, 5, 4,
        1, 2, 3, 4, 5, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2", "s3", "s4", "s5", "s6"))
    ),
    meta.dat = data.frame(
      subject = rep(c("id1", "id2"), each = 3),
      time = rep(c("T10", "T2", "T1"), times = 2),
      group = rep(c("A", "B"), each = 3),
      row.names = c("s1", "s2", "s3", "s4", "s5", "s6"),
      stringsAsFactors = FALSE
    )
  )

  alpha_result <- suppressWarnings(
    generate_alpha_change_per_time_test_long(
      data.obj = data.obj,
      alpha.name = "shannon",
      time.var = "time",
      t0.level = NULL,
      ts.levels = NULL,
      subject.var = "subject",
      group.var = "group"
    )
  )

  beta_result <- suppressWarnings(
    generate_beta_change_per_time_test_long(
      data.obj = data.obj,
      time.var = "time",
      t0.level = NULL,
      ts.levels = NULL,
      subject.var = "subject",
      group.var = "group",
      dist.name = "BC"
    )
  )

  expect_identical(names(alpha_result), c("T2", "T10"))
  expect_identical(names(beta_result), c("T2", "T10"))
})

test_that("mStat_import_qiime2_as_data_obj handles optional taxa and refseq cleanly", {
  ns <- asNamespace("MicrobiomeStat")
  original_read_qza <- get("read_qza", envir = ns)

  unlockBinding("read_qza", ns)
  assign(
    "read_qza",
    function(file, temp = tempdir()) {
      if (identical(file, "otu")) {
        return(matrix(1, nrow = 1, ncol = 1, dimnames = list("f1", "s1")))
      }
      if (identical(file, "ref")) {
        return("mock-refseq")
      }
      stop("unexpected qza input")
    },
    envir = ns
  )
  lockBinding("read_qza", ns)

  on.exit({
    unlockBinding("read_qza", ns)
    assign("read_qza", original_read_qza, envir = ns)
    lockBinding("read_qza", ns)
  }, add = TRUE)

  imported <- suppressMessages(
    mStat_import_qiime2_as_data_obj(
      otu_qza = "otu",
      taxa_qza = NULL,
      refseq_qza = "ref"
    )
  )

  expect_null(imported$feature.ann)
  expect_identical(imported$refseq, "mock-refseq")
})


test_that("mStat_validate_dist_object rejects cross-distance label mismatch", {
  dist.obj <- list(
    BC = stats::as.dist(matrix(
      c(0, 1, 2,
        1, 0, 3,
        2, 3, 0),
      nrow = 3,
      dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
    )),
    Jaccard = stats::as.dist(matrix(
      c(0, 4, 5,
        4, 0, 6,
        5, 6, 0),
      nrow = 3,
      dimnames = list(c("s2", "s1", "s3"), c("s2", "s1", "s3"))
    ))
  )

  expect_error(
    mStat_validate_dist_object(dist.obj, c("BC", "Jaccard")),
    "share identical sample labels and order"
  )
})


test_that("mStat_meta_to_tibble handles existing sample column", {
  meta.dat <- data.frame(
    sample = c("legacy-1", "legacy-2"),
    group = c("A", "B"),
    row.names = c("s1", "s2"),
    stringsAsFactors = FALSE
  )

  meta.tbl <- mStat_meta_to_tibble(meta.dat, sample_col = "sample")

  expect_identical(meta.tbl$sample, c("s1", "s2"))
  expect_identical(meta.tbl$group, c("A", "B"))
  expect_false(anyDuplicated(colnames(meta.tbl)) > 0)
})


test_that("generate_beta_pc_trend_test_long supports group.var = NULL", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 8, 6, 7, 5, 3,
        1, 2, 3, 2, 3, 4,
        3, 3, 3, 3, 3, 3),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3", "s4", "s5", "s6"))
    ),
    meta.dat = data.frame(
      subject = rep(c("id1", "id2"), each = 3),
      time = rep(c("0", "1", "2"), times = 2),
      row.names = c("s1", "s2", "s3", "s4", "s5", "s6"),
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    generate_beta_pc_trend_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = NULL,
      dist.name = "BC"
    )
  )

  expect_true("BC" %in% names(result))
  expect_true(all(c("PC1", "PC2") %in% names(result$BC)))
})


test_that("generate_beta_pc_volatility_test_long fails clearly without group.var", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 8, 6, 7, 5, 3,
        1, 2, 3, 2, 3, 4,
        3, 3, 3, 3, 3, 3),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3", "s4", "s5", "s6"))
    ),
    meta.dat = data.frame(
      subject = rep(c("id1", "id2"), each = 3),
      time = rep(c("0", "1", "2"), times = 2),
      row.names = c("s1", "s2", "s3", "s4", "s5", "s6"),
      stringsAsFactors = FALSE
    )
  )

  expect_error(
    suppressWarnings(
      generate_beta_pc_volatility_test_long(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = NULL,
        dist.name = "BC"
      )
    ),
    "`group.var` is required"
  )
})


test_that("generate_alpha_volatility_test_long ignores invariant covariates and keeps group test", {
  sample_ids <- paste0("s", 1:8)
  data.obj <- list(
    feature.tab = matrix(
      c(10, 11, 12, 13, 8, 7, 6, 5,
        1, 2, 1, 2, 3, 4, 3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("id1", "id2", "id3", "id4"), each = 2),
      time = rep(c("0", "1"), times = 4),
      group = rep(c("A", "A", "B", "B"), each = 2),
      invariant = 1,
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    generate_alpha_volatility_test_long(
      data.obj = data.obj,
      alpha.name = "shannon",
      time.var = "time",
      subject.var = "subject",
      group.var = "group",
      adj.vars = "invariant"
    )
  )

  expect_true("shannon" %in% names(result))
  expect_true(any(grepl("group", result$shannon$Term, fixed = TRUE)))
})


test_that("group and time contracts reject missing or degenerate inputs", {
  meta <- data.frame(
    subject = c("id1", "id1", "id2", "id2"),
    group = c("A", "A", "B", "B"),
    time = c("T1", "T2", "T1", "T2"),
    stringsAsFactors = FALSE
  )

  expect_error(
    mStat_validate_group_var_contract(meta, group.var = NULL, context = "contract test"),
    "`group.var` is required"
  )

  expect_error(
    mStat_validate_group_var_contract(meta, group.var = "missing", context = "contract test"),
    "was not found in metadata"
  )

  expect_error(
    mStat_validate_group_var_contract(
      transform(meta, group = "A"),
      group.var = "group",
      context = "contract test"
    ),
    "must have at least two observed levels/values"
  )

  expect_error(
    mStat_validate_time_var_contract(meta, time.var = NULL, context = "contract test"),
    "`time.var` is required"
  )

  expect_error(
    mStat_validate_time_var_contract(meta, time.var = "missing", context = "contract test"),
    "was not found in metadata"
  )

  expect_error(
    mStat_validate_time_var_contract(
      transform(meta, time = NA_character_),
      time.var = "time",
      context = "contract test"
    ),
    "has no observed values"
  )

  expect_error(
    mStat_validate_time_var_contract(
      transform(meta, time = "T1"),
      time.var = "time",
      context = "contract test",
      require_variation = TRUE
    ),
    "must have at least two observed levels/values"
  )
})


test_that("group-required public helpers keep non-NULL signature contract", {
  expect_identical(
    formals(generate_alpha_per_time_test_long)$group.var,
    quote(expr = )
  )
  expect_identical(
    formals(generate_beta_change_per_time_test_long)$group.var,
    quote(expr = )
  )
  expect_identical(
    formals(generate_beta_volatility_test_long)$group.var,
    quote(expr = )
  )
  expect_identical(
    formals(generate_beta_pc_volatility_test_long)$group.var,
    quote(expr = )
  )
})
