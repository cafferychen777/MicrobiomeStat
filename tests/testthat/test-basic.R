test_that("basic validation returns a usable data object", {
  data.obj <- list(
    feature.tab = matrix(
      c(3, 1,
        2, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      group = c("A", "B"),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  validated <- suppressMessages(mStat_validate_data(data.obj))

  expect_type(validated, "list")
  expect_identical(colnames(validated$feature.tab), rownames(validated$meta.dat))
})
