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
