test_that("mStat_normalize_data TMM uses effective library sizes", {
  skip_if_not_installed("edgeR")

  data.obj <- list(
    feature.tab = matrix(
      c(10, 20,
        100, 200),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
  )

  norm_result <- mStat_normalize_data(data.obj, method = "TMM")

  expect_equal(
    unname(norm_result$data.obj.norm$feature.tab[, "s1"]),
    unname(norm_result$data.obj.norm$feature.tab[, "s2"]),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(colSums(norm_result$data.obj.norm$feature.tab)),
    c(1, 1),
    tolerance = 1e-8
  )
})

test_that("mStat_normalize_data DESeq uses median-ratio size factors", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 20,
        100, 200),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
  )

  norm_result <- mStat_normalize_data(data.obj, method = "DESeq")

  expect_equal(
    unname(norm_result$data.obj.norm$feature.tab[, "s1"]),
    unname(norm_result$data.obj.norm$feature.tab[, "s2"]),
    tolerance = 1e-8
  )
  expect_equal(
    unname(norm_result$scale_factor["s2"] / norm_result$scale_factor["s1"]),
    10,
    tolerance = 1e-8
  )
})

test_that("mStat_normalize_data DESeq falls back cleanly for sparse counts", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 0,
        0, 10),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
  )

  expect_warning(
    norm_result <- mStat_normalize_data(data.obj, method = "DESeq"),
    "Falling back to a positive-count geometric mean variant"
  )

  expect_true(all(is.finite(norm_result$scale_factor)))
  expect_true(all(norm_result$scale_factor > 0))
  expect_true(all(is.finite(norm_result$data.obj.norm$feature.tab)))
})

test_that("mStat_normalize_data validates rarefaction depth before calling vegan", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 4,
        3, 2),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
  )

  expect_error(
    mStat_normalize_data(data.obj, method = "Rarefy", depth = Inf),
    "single non-negative finite number"
  )
  expect_error(
    mStat_normalize_data(data.obj, method = "Rarefy", depth = c(1, 2)),
    "single non-negative finite number"
  )
})

test_that("mStat_rarefy_data validates rarefaction depth before calling vegan", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 4,
        3, 2),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    )
  )

  expect_error(
    mStat_rarefy_data(data.obj, depth = Inf),
    "single non-negative finite number"
  )
  expect_error(
    mStat_rarefy_data(data.obj, depth = c(1, 2)),
    "single non-negative finite number"
  )
})
