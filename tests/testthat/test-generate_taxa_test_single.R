test_that("generate_taxa_test_single handles empty samples correctly", {
  # Load test data
  data(peerj32.obj)

  # Create a case with empty samples
  peerj32.obj$feature.tab[, 1] <- 0
  peerj32.obj$feature.tab[14, 1] <- 29438 # sum of column 1

  # Test with only group.var
  test.list <- generate_taxa_test_single(
    data.obj = peerj32.obj,
    group.var = "group",
    feature.dat.type = "count",
    feature.level = c("Genus"),
    prev.filter = 0.1,
    abund.filter = 0.01
  )

  # Expectations
  expect_type(test.list, "list")
  expect_named(test.list, "Genus")
  expect_false(any(is.null(test.list$Genus)))

  # Verify the structure of the results
  expect_true(all(c("Variable", "Coefficient", "SE", "P.Value",
                    "Adjusted.P.Value", "Mean.Abundance", "Prevalence") %in%
                    names(test.list$Genus[[1]])))
})

test_that("generate_taxa_test_single handles single column metadata correctly", {
  # Load test data
  data(peerj32.obj)

  # Test with only group.var (single column metadata)
  test.list <- generate_taxa_test_single(
    data.obj = peerj32.obj,
    group.var = "group",
    feature.dat.type = "count",
    feature.level = c("Genus"),
    prev.filter = 0.1,
    abund.filter = 0.01
  )

  # Expectations
  expect_type(test.list, "list")
  expect_named(test.list, "Genus")
  expect_false(any(is.null(test.list$Genus)))
})

test_that("generate_taxa_test_single handles multiple metadata variables correctly", {
  # Load test data
  data(peerj32.obj)

  # Test with group.var and adj.vars
  test.list <- generate_taxa_test_single(
    data.obj = peerj32.obj,
    group.var = "group",
    adj.vars = "sex",
    feature.dat.type = "count",
    feature.level = c("Genus"),
    prev.filter = 0.1,
    abund.filter = 0.01
  )

  # Expectations
  expect_type(test.list, "list")
  expect_named(test.list, "Genus")
  expect_false(any(is.null(test.list$Genus)))
})

test_that("generate_taxa_test_single handles filtering edge cases", {
  # Load test data
  data(peerj32.obj)

  # Test with less extreme filtering parameters
  test.list <- generate_taxa_test_single(
    data.obj = peerj32.obj,
    group.var = "group",
    feature.dat.type = "count",
    feature.level = c("Genus"),
    prev.filter = 0.5,  # Lower threshold
    abund.filter = 0.05  # Lower threshold
  )

  # Expectations for normal case
  expect_type(test.list, "list")
  expect_named(test.list, "Genus")

  # Test extreme filtering case
  expect_warning(
    test.list_extreme <- generate_taxa_test_single(
      data.obj = peerj32.obj,
      group.var = "group",
      feature.dat.type = "count",
      feature.level = c("Genus"),
      prev.filter = 0.9,  # Extreme threshold
      abund.filter = 0.1   # Extreme threshold
    ),
    "No features remain after filtering"
  )

  # Check structure of empty result
  expect_type(test.list_extreme, "list")
  expect_named(test.list_extreme, "Genus")
  expect_true(length(test.list_extreme$Genus) == 0)
})

test_that("generate_taxa_test_single handles time variable correctly", {
  # Load test data
  data(peerj32.obj)

  # Test with time variable
  test.list <- generate_taxa_test_single(
    data.obj = peerj32.obj,
    time.var = "time",
    t.level = "1",
    group.var = "group",
    feature.dat.type = "count",
    feature.level = c("Genus"),
    prev.filter = 0.1,
    abund.filter = 0.01
  )

  # Expectations
  expect_type(test.list, "list")
  expect_named(test.list, "Genus")
  expect_false(any(is.null(test.list$Genus)))
})
