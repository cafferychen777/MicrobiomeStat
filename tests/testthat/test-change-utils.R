test_that("compute_taxa_change: absolute change works", {
  after  <- c(5, 10, 0.3)
  before <- c(2,  8, 0.1)
  result <- compute_taxa_change(after, before, "absolute change")
  expect_equal(result, c(3, 2, 0.2))
})

test_that("compute_taxa_change: relative change works", {
  after  <- c(6, 0, 0)
  before <- c(4, 0, 0)
  result <- compute_taxa_change(after, before, "relative change")
  # (6-4)/(6+4) = 0.2; 0/0 -> 0
  expect_equal(result, c(0.2, 0, 0))
})

test_that("compute_taxa_change: log fold change with no zeros", {
  after  <- c(8, 4)
  before <- c(2, 4)
  result <- compute_taxa_change(after, before, "log fold change",
                                feature_id = c("A", "A"), verbose = FALSE)
  expect_equal(result, c(log2(8) - log2(2), log2(4) - log2(4)))
  expect_equal(result, c(2, 0))
})

test_that("compute_taxa_change: log fold change applies per-feature pseudocount for zeros", {
  # Feature A: values are 0.1, 0.2, 0 -> half-min of non-zero = 0.1/2 = 0.05

  # Feature B: values are 0.4, 0.6 -> no zeros, no pseudocount needed
  after  <- c(0,   0.6)
  before <- c(0.1, 0.4)
  fid    <- c("A", "B")
  result <- compute_taxa_change(after, before, "log fold change",
                                feature_id = fid, verbose = FALSE)

  # Feature A: zero in after -> replaced with 0.05, so log2(0.05) - log2(0.1)
  expect_equal(result[1], log2(0.05) - log2(0.1))
  # Feature B: no zeros, straightforward

  expect_equal(result[2], log2(0.6) - log2(0.4))
})

test_that("compute_taxa_change: log fold change all-zero feature uses 1e-10 fallback", {
  after  <- c(0, 0)
  before <- c(0, 0)
  fid    <- c("A", "A")
  result <- compute_taxa_change(after, before, "log fold change",
                                feature_id = fid, verbose = FALSE)
  # Both replaced with 1e-10, so log2(1e-10) - log2(1e-10) = 0
  expect_equal(result, c(0, 0))
})

test_that("compute_taxa_change: log fold change without feature_id uses global pseudocount", {
  after  <- c(0,   0.4)
  before <- c(0.2, 0.6)
  result <- compute_taxa_change(after, before, "log fold change",
                                verbose = FALSE)
  # Global half-min = min(0.2, 0.4, 0.6) / 2 = 0.1
  expect_equal(result[1], log2(0.1) - log2(0.2))
  expect_equal(result[2], log2(0.4) - log2(0.6))
})

test_that("compute_taxa_change: custom function passthrough", {
  my_func <- function(a, b) a * b
  result <- compute_taxa_change(c(3, 4), c(2, 5), my_func)
  expect_equal(result, c(6, 20))
})

test_that("compute_taxa_change: unknown method defaults to absolute change", {
  result <- compute_taxa_change(c(10, 20), c(3, 8), "unknown method")
  expect_equal(result, c(7, 12))
})

test_that("compute_alpha_change: log fold change works", {
  after  <- c(3.0, 2.0)
  before <- c(1.5, 2.0)
  result <- compute_alpha_change(after, before, "log fold change")
  expect_equal(result, c(log2(3.0 / 1.5), log2(2.0 / 2.0)))
  expect_equal(result, c(1, 0))
})

test_that("compute_alpha_change: absolute change works", {
  after  <- c(3.0, 2.0)
  before <- c(1.5, 2.0)
  result <- compute_alpha_change(after, before, "absolute change")
  expect_equal(result, c(1.5, 0))
})

test_that("compute_alpha_change: custom function passthrough", {
  result <- compute_alpha_change(c(4, 6), c(2, 3), function(a, b) a / b)
  expect_equal(result, c(2, 2))
})

test_that("compute_alpha_change: unknown method defaults to absolute change with message", {
  expect_message(
    result <- compute_alpha_change(c(5, 10), c(2, 7), "bogus"),
    "Defaulting to"
  )
  expect_equal(result, c(3, 3))
})

test_that("describe_change_method: taxa context", {
  desc <- describe_change_method("log fold change", context = "taxa")
  expect_true(grepl("log2 fold changes", desc))
  expect_true(grepl("The changes were", desc))

  desc <- describe_change_method("relative change", context = "taxa")
  expect_true(grepl("relative changes", desc))

  desc <- describe_change_method("absolute change", context = "taxa")
  expect_true(grepl("absolute changes", desc))
})

test_that("describe_change_method: alpha context", {
  desc <- describe_change_method("log fold change", context = "alpha")
  expect_true(grepl("t0.level", desc))
})

test_that("describe_change_method: custom function", {
  desc <- describe_change_method(function(a, b) a - b)
  expect_true(grepl("custom function", desc))
})

test_that("constants are defined correctly", {
  expect_equal(.CHANGE_LOG_FOLD, "log fold change")
  expect_equal(.CHANGE_RELATIVE, "relative change")
  expect_equal(.CHANGE_ABSOLUTE, "absolute change")
  expect_equal(.CHANGE_METHODS, c("log fold change", "relative change", "absolute change"))
})
