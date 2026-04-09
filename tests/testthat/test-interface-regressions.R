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
