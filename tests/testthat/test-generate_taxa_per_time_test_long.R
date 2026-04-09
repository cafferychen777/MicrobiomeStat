test_that("generate_taxa_per_time_test_long drops time points with no group contrasts", {
  feature.tab <- matrix(
    c(10, 12, 14,
      3, 4, 5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
  )

  meta.dat <- data.frame(
    subject = c("u1", "u2", "u3"),
    time = c("T1", "T2", "T2"),
    group = c("A", "A", "B"),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  feature.ann <- matrix(
    c("g1", "g2"),
    ncol = 1,
    dimnames = list(c("f1", "f2"), "Genus")
  )

  data.obj <- list(
    feature.tab = feature.tab,
    meta.dat = meta.dat,
    feature.ann = feature.ann
  )

  result <- suppressWarnings(
    generate_taxa_per_time_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      feature.level = "Genus",
      feature.dat.type = "proportion"
    )
  )

  expect_identical(names(result), "T2")
  expect_true("Genus" %in% names(result$T2))
  expect_gt(length(result$T2$Genus), 0)
})


test_that("generate_taxa_per_time_test_long keeps continuous group labels", {
  set.seed(123)

  feature.tab <- matrix(
    rpois(20, lambda = 8),
    nrow = 2,
    dimnames = list(c("f1", "f2"), paste0("s", 1:10))
  )

  meta.dat <- data.frame(
    subject = paste0("u", 1:10),
    time = rep(c("T1", "T2"), each = 5),
    score = c(1.2, 1.6, 2.0, 2.5, 3.0, 1.1, 1.7, 2.1, 2.4, 3.2),
    row.names = colnames(feature.tab),
    stringsAsFactors = FALSE
  )

  feature.ann <- matrix(
    c("g1", "g2"),
    ncol = 1,
    dimnames = list(c("f1", "f2"), "Genus")
  )

  data.obj <- list(
    feature.tab = feature.tab,
    meta.dat = meta.dat,
    feature.ann = feature.ann
  )

  result <- suppressWarnings(
    generate_taxa_per_time_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "score",
      feature.level = "Genus",
      feature.dat.type = "proportion"
    )
  )

  expect_true(all(c("T1", "T2") %in% names(result)))
  expect_true("score" %in% names(result$T1$Genus))
  expect_true("score" %in% names(result$T2$Genus))
})


test_that("generate_taxa_per_time_test_long fails fast on invalid group/time inputs", {
  feature.tab <- matrix(
    c(10, 12, 14,
      3, 4, 5),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), c("s1", "s2", "s3"))
  )

  feature.ann <- matrix(
    c("g1", "g2"),
    ncol = 1,
    dimnames = list(c("f1", "f2"), "Genus")
  )

  base_meta <- data.frame(
    subject = c("u1", "u2", "u3"),
    time = c("T1", "T2", "T2"),
    group = c("A", "A", "B"),
    row.names = c("s1", "s2", "s3"),
    stringsAsFactors = FALSE
  )

  base_data <- list(feature.tab = feature.tab, meta.dat = base_meta, feature.ann = feature.ann)

  expect_error(
    generate_taxa_per_time_test_long(
      data.obj = base_data,
      subject.var = "subject",
      time.var = NULL,
      group.var = "group",
      feature.level = "Genus",
      feature.dat.type = "proportion"
    ),
    "`time.var` is required"
  )

  expect_error(
    generate_taxa_per_time_test_long(
      data.obj = base_data,
      subject.var = "subject",
      time.var = "missing_time",
      group.var = "group",
      feature.level = "Genus",
      feature.dat.type = "proportion"
    ),
    "was not found in metadata"
  )

  one_group_data <- list(
    feature.tab = feature.tab,
    meta.dat = transform(base_meta, group = "A"),
    feature.ann = feature.ann
  )

  expect_error(
    suppressWarnings(
      generate_taxa_per_time_test_long(
        data.obj = one_group_data,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        feature.level = "Genus",
        feature.dat.type = "proportion"
      )
    ),
    "No time points could be analyzed"
  )
})
