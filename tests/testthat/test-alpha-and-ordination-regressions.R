test_that("mStat_calculate_adjusted_alpha_diversity aligns by sample labels", {
  alpha.obj <- list(
    shannon = data.frame(
      shannon = c(11, 12, 13),
      row.names = c("s1", "s2", "s3")
    )
  )
  meta.dat <- data.frame(
    age = c(30, 20, 10),
    row.names = c("s3", "s2", "s1")
  )

  adjusted <- mStat_calculate_adjusted_alpha_diversity(
    alpha.obj = alpha.obj,
    meta.dat = meta.dat,
    adj.vars = "age"
  )

  fit <- stats::lm(adjusted$shannon$shannon ~ meta.dat[rownames(adjusted$shannon), "age"])

  expect_identical(rownames(adjusted$shannon), c("s1", "s2", "s3"))
  expect_lt(abs(unname(stats::coef(fit)[2])), 1e-10)
})

test_that("mStat_calculate_adjusted_alpha_diversity errors on missing covariates", {
  alpha.obj <- list(
    shannon = data.frame(shannon = c(1, 2), row.names = c("s1", "s2"))
  )
  meta.dat <- data.frame(group = c("A", "B"), row.names = c("s1", "s2"))

  expect_error(
    mStat_calculate_adjusted_alpha_diversity(alpha.obj, meta.dat, "age"),
    "missing"
  )
})

test_that("mStat_calculate_PC accepts legacy method vectors and defaults dist.name", {
  dist.obj <- list(
    BC = stats::dist(
      matrix(
        c(1, 2,
          3, 4,
          5, 6),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("s1", "s2", "s3"), c("x", "y"))
      )
    )
  )

  expect_warning(
    pc.obj <- mStat_calculate_PC(
      dist.obj = dist.obj,
      method = c("mds", "nmds"),
      k = 2,
      dist.name = NULL
    ),
    "supports one method per call"
  )

  expect_identical(names(pc.obj), "BC")
  expect_identical(rownames(pc.obj$BC$points), c("s1", "s2", "s3"))
})

test_that("mStat_calculate_beta_diversity maps UniFrac variants correctly", {
  skip_if_not_installed("GUniFrac")

  otu <- matrix(
    c(10, 0, 5,
      0, 10, 5,
      5, 5, 5),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("t1", "t2", "t3"), c("s1", "s2", "s3"))
  )

  tree <- ape::read.tree(text = "((t1:1,t2:1):1,t3:1);")
  data.obj <- list(
    feature.tab = otu,
    meta.dat = data.frame(
      group = c("A", "B", "C"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    ),
    tree = tree
  )

  dist.obj <- suppressWarnings(
    mStat_calculate_beta_diversity(data.obj, dist.name = c("UniFrac", "WUniFrac"))
  )
  expected <- GUniFrac::GUniFrac(t(otu), tree, alpha = 1)$unifracs

  expect_equal(as.matrix(dist.obj$UniFrac), expected[, , "d_UW"], tolerance = 1e-8)
  expect_equal(as.matrix(dist.obj$WUniFrac), expected[, , "d_1"], tolerance = 1e-8)
})

test_that("mStat_calculate_beta_diversity uses binary Jaccard semantics", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 10,
        1, 1),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(data.obj, dist.name = "Jaccard"))

  expect_equal(as.matrix(dist.obj$Jaccard)[1, 2], 0)
  expect_identical(rownames(attr(dist.obj$Jaccard, "metadata")), c("s1", "s2"))
})

test_that("mStat_calculate_adjusted_distance preserves metadata and sample order", {
  data.obj <- list(
    meta.dat = data.frame(
      group = c("B", "A", "B"),
      row.names = c("s3", "s1", "s2"),
      stringsAsFactors = FALSE
    )
  )
  distance_matrix <- matrix(
    c(0, 1, 2,
      1, 0, 3,
      2, 3, 0),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
  )
  dist.obj <- list(
    BC = mStat_attach_dist_metadata(stats::as.dist(distance_matrix), data.obj$meta.dat)
  )

  adjusted <- mStat_calculate_adjusted_distance(
    data.obj = data.obj,
    dist.obj = dist.obj,
    adj.vars = "group",
    dist.name = "BC"
  )

  expect_identical(attr(adjusted$BC, "Labels"), c("s1", "s2", "s3"))
  expect_identical(rownames(attr(adjusted$BC, "metadata")), c("s1", "s2", "s3"))
})

test_that("mStat_calculate_adjusted_alpha_diversity rejects missing covariates clearly", {
  alpha.obj <- list(
    shannon = data.frame(shannon = c(1, 2, 3), row.names = c("s1", "s2", "s3"))
  )
  meta.dat <- data.frame(
    age = c(10, NA, 30),
    row.names = c("s1", "s2", "s3")
  )

  expect_error(
    mStat_calculate_adjusted_alpha_diversity(alpha.obj, meta.dat, "age"),
    "missing covariates"
  )
})

test_that("mStat_calculate_adjusted_distance rejects missing covariates clearly", {
  data.obj <- list(
    meta.dat = data.frame(
      group = c("A", NA, "B"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )
  distance_matrix <- matrix(
    c(0, 1, 2,
      1, 0, 3,
      2, 3, 0),
    nrow = 3,
    dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
  )
  dist.obj <- list(
    BC = mStat_attach_dist_metadata(stats::as.dist(distance_matrix), data.obj$meta.dat)
  )

  expect_error(
    mStat_calculate_adjusted_distance(
      data.obj = data.obj,
      dist.obj = dist.obj,
      adj.vars = "group",
      dist.name = "BC"
    ),
    "missing covariates"
  )
})

test_that("generate_alpha_change_test_pair rejects more than one follow-up level", {
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

  expect_error(
    generate_alpha_change_test_pair(
      data.obj = data.obj,
      alpha.obj = alpha.obj,
      alpha.name = "shannon",
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      change.base = "T1"
    ),
    "exactly one follow-up level"
  )
})


test_that("paired change tests fail fast when group.var changes within subject", {
  sample_ids <- paste0("s", 1:4)
  data.obj <- list(
    feature.tab = matrix(
      c(
        5, 6, 7, 8,
        8, 7, 6, 5
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), sample_ids)
    ),
    meta.dat = data.frame(
      subject = rep(c("u1", "u2"), each = 2),
      time = rep(c("T1", "T2"), times = 2),
      group = c("A", "B", "C", "C"),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )

  expect_error(
    suppressWarnings(
      generate_taxa_change_test_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        feature.level = "original",
        feature.dat.type = "proportion"
      )
    ),
    "must be constant within each subject"
  )

  expect_error(
    suppressWarnings(
      generate_taxa_test_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        feature.level = "original",
        feature.dat.type = "proportion"
      )
    ),
    "must be constant within each subject"
  )

  expect_error(
    suppressWarnings(
      generate_beta_change_test_pair(
        data.obj = data.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        dist.name = "BC"
      )
    ),
    "must be constant within each subject"
  )
})


test_that("generate_alpha_trend_test_long falls back from random slopes", {
  data(peerj32.obj)

  expect_message(
    result <- suppressWarnings(
      generate_alpha_trend_test_long(
        data.obj = peerj32.obj,
        alpha.name = "shannon",
        time.var = "time",
        subject.var = "subject",
        group.var = "group"
      )
    ),
    "Simplifying the random-effects structure"
  )

  expect_true("shannon" %in% names(result))
  expect_true("Term" %in% colnames(result$shannon))
})


test_that("create_mixed_effects_formula includes adjustment covariates", {
  formula <- create_mixed_effects_formula(
    response.var = "distance",
    time.var = "time",
    group.var = "group",
    subject.var = "subject",
    adj.vars = "batch",
    random_slopes = TRUE
  )

  formula_str <- paste(deparse(formula), collapse = "")
  expect_match(formula_str, "batch")
})


test_that("ordination axis labels are robust when eigenvalues include negatives", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2, 3,
        2, 1, 3,
        3, 3, 1),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3"))
    ),
    meta.dat = data.frame(
      group = c("A", "A", "B"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )

  dist.obj <- list(
    BC = stats::as.dist(matrix(
      c(0, 1, 2,
        1, 0, 1,
        2, 1, 0),
      nrow = 3,
      dimnames = list(c("s1", "s2", "s3"), c("s1", "s2", "s3"))
    ))
  )
  attr(dist.obj$BC, "metadata") <- data.obj$meta.dat

  pc.obj <- list(
    BC = list(
      points = matrix(
        c(1, 0,
          0, 1,
          -1, -1),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(c("s1", "s2", "s3"), c("V1", "V2"))
      ),
      eig = c(0.4, -0.3, -0.1)
    )
  )

  plot_single <- suppressWarnings(
    generate_beta_ordination_single(
      data.obj = data.obj,
      dist.obj = dist.obj,
      pc.obj = pc.obj,
      group.var = "group",
      dist.name = "BC",
      pdf = FALSE
    )
  )

  pair_data <- data.obj
  pair_data$meta.dat$subject <- c("id1", "id1", "id2")
  pair_data$meta.dat$time <- c("T0", "T1", "T0")

  plot_pair <- suppressWarnings(
    generate_beta_ordination_pair(
      data.obj = pair_data,
      dist.obj = dist.obj,
      pc.obj = pc.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      dist.name = "BC",
      pdf = FALSE
    )
  )

  axis_labels <- mStat_build_axis_labels_from_eig(pc.obj$BC$eig)

  expect_match(axis_labels$x, "Axis 1")
  expect_match(axis_labels$y, "Axis 2")
  expect_false(grepl("NaN", axis_labels$x, fixed = TRUE))
  expect_false(grepl("NaN", axis_labels$y, fixed = TRUE))

  expect_s3_class(plot_single$BC, "ggplot")
  expect_s3_class(plot_pair$BC, "ggplot")
})

test_that("generate_beta_ordination_pair orders character time labels by time semantics", {
  sample_ids <- c("s1", "s2", "s3", "s4")
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2, 3, 4),
      nrow = 1,
      dimnames = list("f1", sample_ids)
    ),
    meta.dat = data.frame(
      subject = c("u1", "u1", "u2", "u2"),
      time = c("10", "2", "10", "2"),
      group = c("A", "A", "B", "B"),
      row.names = sample_ids,
      stringsAsFactors = FALSE
    )
  )
  dist.obj <- list(
    BC = stats::dist(
      matrix(
        c(0, 0,
          1, 0,
          0, 1,
          1, 1),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(sample_ids, c("x", "y"))
      )
    )
  )
  pc.obj <- list(
    BC = list(
      points = matrix(
        c(10, 0,
          2, 0,
          11, 0,
          3, 0),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(sample_ids, c("PC1", "PC2"))
      ),
      eig = c(1, 1)
    )
  )

  plot_pair <- suppressWarnings(
    generate_beta_ordination_pair(
      data.obj = data.obj,
      dist.obj = dist.obj,
      pc.obj = pc.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      dist.name = "BC",
      pdf = FALSE
    )
  )

  built <- ggplot_build(plot_pair$BC)
  segment_layers <- built$data[vapply(built$data, function(layer) {
    all(c("x", "xend", "y", "yend") %in% names(layer))
  }, logical(1))]

  subject_segments <- unique(segment_layers[[1]][, c("x", "xend")])
  expect_equal(subject_segments$x, c(2, 10, 3, 11))
  expect_equal(subject_segments$xend, c(10, NA, 11, NA))
})
