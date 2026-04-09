make_test_data_obj <- function(sample_order = c("s1", "s2")) {
  feature.tab <- matrix(
    c(5, 1,
      2, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), sample_order)
  )

  meta.dat <- data.frame(
    group = c("A", "B"),
    row.names = c("s1", "s2"),
    stringsAsFactors = FALSE
  )

  feature.ann <- matrix(
    c("Phylum1", "Phylum2"),
    ncol = 1,
    dimnames = list(c("f1", "f2"), "Phylum")
  )

  list(
    feature.tab = feature.tab,
    meta.dat = meta.dat,
    feature.ann = feature.ann
  )
}

test_that("mStat_validate_data returns an aligned object and keeps feature.ann optional", {
  data.obj <- list(
    feature.tab = matrix(
      c(5, 1,
        2, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s2", "s1"))
    ),
    meta.dat = data.frame(
      group = c("A", "B"),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  validated <- suppressMessages(mStat_validate_data(data.obj))

  expect_type(validated, "list")
  expect_identical(rownames(validated$meta.dat), c("s2", "s1"))
  expect_null(validated$feature.ann)
})

test_that("mStat_subset_data is a no-op when no subset is requested", {
  data.obj <- make_test_data_obj()

  subsetted <- suppressMessages(mStat_subset_data(data.obj))

  expect_identical(rownames(subsetted$meta.dat), c("s1", "s2"))
  expect_identical(colnames(subsetted$feature.tab), c("s1", "s2"))
  expect_identical(rownames(subsetted$feature.ann), c("f1", "f2"))
})

test_that("mStat_subset_data preserves requested sample order across metadata and feature table", {
  data.obj <- make_test_data_obj()

  subsetted <- suppressMessages(mStat_subset_data(data.obj, samIDs = c("s2", "s1")))

  expect_identical(rownames(subsetted$meta.dat), c("s2", "s1"))
  expect_identical(colnames(subsetted$feature.tab), c("s2", "s1"))
})

test_that("mStat_subset_alpha resolves numeric indices against row names", {
  alpha.obj <- list(
    shannon = data.frame(
      shannon = c(1.1, 2.2),
      row.names = c("s1", "s2")
    )
  )

  subsetted <- mStat_subset_alpha(alpha.obj, c(2, 1))

  expect_identical(rownames(subsetted$shannon), c("s2", "s1"))
  expect_equal(subsetted$shannon$shannon, c(2.2, 1.1))
})

test_that("mStat_process_time_variable keeps non-numeric character labels intact", {
  data.obj <- list(
    feature.tab = matrix(
      c(10, 20,
        30, 40),
      nrow = 2,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      tp = c("baseline", "followup"),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    )
  )

  processed <- suppressMessages(mStat_process_time_variable(data.obj, "tp"))

  expect_identical(rownames(processed$meta.dat), c("s1", "s2"))
  expect_identical(levels(processed$meta.dat$tp), c("baseline", "followup"))
})

test_that("detect_study_design orders numeric-like and natural character times consistently", {
  numeric_like <- list(
    meta.dat = data.frame(
      tp = c("1", "2", "10"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )
  natural_like <- list(
    meta.dat = data.frame(
      tp = c("T1", "T2", "T10"),
      row.names = c("s1", "s2", "s3"),
      stringsAsFactors = FALSE
    )
  )

  numeric_design <- detect_study_design(numeric_like, time.var = "tp")
  natural_design <- detect_study_design(natural_like, time.var = "tp")

  expect_identical(numeric_design$t0.level, "1")
  expect_identical(numeric_design$ts.levels, c("2", "10"))
  expect_identical(natural_design$t0.level, "T1")
  expect_identical(natural_design$ts.levels, c("T2", "T10"))
})

test_that("mStat_coerce_time_to_numeric uses labels instead of factor codes", {
  expect_equal(
    mStat_coerce_time_to_numeric(factor(c("0", "2", "4")), "time"),
    c(0, 2, 4)
  )

  expect_error(
    mStat_coerce_time_to_numeric(c("baseline", "followup"), "time"),
    "must be numeric or coercible to numeric"
  )
})

test_that("mStat_subset_PC resolves numeric indices from points and preserves eig", {
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
  pc.obj <- mStat_calculate_PC(dist.obj, method = "mds", k = 2, dist.name = "BC")

  subsetted <- mStat_subset_PC(pc.obj, c(2, 1))

  expect_identical(rownames(subsetted$BC$points), c("s2", "s1"))
  expect_identical(subsetted$BC$eig, pc.obj$BC$eig)
})

test_that("mStat_prepare_alpha_data aligns alpha tables by sample labels", {
  alpha.obj <- list(
    shannon = data.frame(
      shannon = c(1.1, 2.2),
      row.names = c("s1", "s2")
    ),
    simpson = data.frame(
      simpson = c(20, 10),
      row.names = c("s2", "s1")
    )
  )

  alpha_df <- mStat_prepare_alpha_data(alpha.obj, sample_col = "sample")

  expect_identical(alpha_df$sample, c("s1", "s2"))
  expect_equal(alpha_df$shannon, c(1.1, 2.2))
  expect_equal(alpha_df$simpson, c(10, 20))
})

test_that("mStat_prepare_alpha_data rejects inconsistent alpha sample sets", {
  alpha.obj <- list(
    shannon = data.frame(
      shannon = c(1.1, 2.2),
      row.names = c("s1", "s2")
    ),
    simpson = data.frame(
      simpson = c(3.3, 4.4),
      row.names = c("s1", "s3")
    )
  )

  expect_error(
    mStat_prepare_alpha_data(alpha.obj, sample_col = "sample"),
    "same sample set"
  )
})

test_that("mStat_prepare_pc_long_data joins metadata by sample labels", {
  pc.points <- matrix(
    c(20, 200,
      10, 100),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(c("s2", "s1"), c("Axis1", "Axis2"))
  )
  meta.dat <- data.frame(
    group = c("A", "B"),
    row.names = c("s1", "s2"),
    stringsAsFactors = FALSE
  )

  pc_long <- mStat_prepare_pc_long_data(
    pc.points = pc.points,
    pc.ind = c(1, 2),
    meta.dat = meta.dat,
    vars = "group",
    sample_col = "sample",
    join = "inner"
  )

  s1_pc1 <- pc_long[pc_long$sample == "s1" & pc_long$PC == "PC1", ]
  s2_pc2 <- pc_long[pc_long$sample == "s2" & pc_long$PC == "PC2", ]

  expect_equal(s1_pc1$value, 10)
  expect_identical(s1_pc1$group, "A")
  expect_equal(s2_pc2$value, 200)
  expect_identical(s2_pc2$group, "B")
})

test_that("mStat_resolve_pair_timepoints uses ordered labels and rejects multiple follow-ups", {
  resolved <- suppressMessages(
    mStat_resolve_pair_timepoints(
      values = c("T10", "T2"),
      time.var = "time",
      change.base = NULL,
      context = "test"
    )
  )

  expect_identical(resolved$change.base, "T2")
  expect_identical(resolved$change.after, "T10")

  expect_error(
    mStat_resolve_pair_timepoints(
      values = c("T1", "T2", "T10"),
      time.var = "time",
      change.base = "T1",
      context = "test"
    ),
    "exactly one follow-up level"
  )
})

test_that("mStat_filter ignores NA values and preserves matrix shape", {
  x <- matrix(
    c(
      1, NA, 0,
      0, 0, 2,
      NA, NA, NA
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("f1", "f2", "f3"), c("s1", "s2", "s3"))
  )

  filtered <- mStat_filter(x, prev.filter = 0.5, abund.filter = -Inf)

  expect_true(is.matrix(filtered))
  expect_identical(rownames(filtered), "f1")
  expect_identical(colnames(filtered), c("s1", "s2", "s3"))
  expect_equal(dim(filtered), c(1L, 3L))
})

test_that("mStat_remove_feature returns the input object when level is missing or unmatched", {
  data.obj <- make_test_data_obj()

  unchanged_missing <- suppressMessages(
    mStat_remove_feature(data.obj, featureIDs = "Phylum1", feature.level = NULL)
  )
  unchanged_unmatched <- suppressMessages(
    mStat_remove_feature(data.obj, featureIDs = "Phylum1", feature.level = "Genus")
  )

  expect_identical(unchanged_missing, data.obj)
  expect_identical(unchanged_unmatched, data.obj)
})

test_that("mStat_remove_feature preserves matrix dimensions when one feature remains", {
  data.obj <- make_test_data_obj()

  removed <- suppressMessages(
    mStat_remove_feature(data.obj, featureIDs = "f1", feature.level = "original")
  )

  expect_true(is.matrix(removed$feature.tab))
  expect_true(is.matrix(removed$feature.ann))
  expect_equal(dim(removed$feature.tab), c(1L, 2L))
  expect_equal(dim(removed$feature.ann), c(1L, 1L))
  expect_identical(rownames(removed$feature.tab), "f2")
})

test_that("mStat_combine_data merges feature/meta/taxonomy without suffix corruption", {
  data.obj1 <- list(
    feature.tab = matrix(
      c(1, 2,
        3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    meta.dat = data.frame(
      group = c("A", "B"),
      row.names = c("s1", "s2"),
      stringsAsFactors = FALSE
    ),
    feature.ann = matrix(
      c("p1", "p2"),
      ncol = 1,
      dimnames = list(c("f1", "f2"), "Phylum")
    )
  )

  data.obj2 <- list(
    feature.tab = matrix(
      c(1, 5,
        3, 6),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f3"), c("s1", "s3"))
    ),
    meta.dat = data.frame(
      group = c("A", "C"),
      row.names = c("s1", "s3"),
      stringsAsFactors = FALSE
    ),
    feature.ann = matrix(
      c("p1", "p3"),
      ncol = 1,
      dimnames = list(c("f1", "f3"), "Phylum")
    )
  )

  combined <- suppressMessages(mStat_combine_data(data.obj1, data.obj2))

  expect_identical(colnames(combined$feature.tab), c("s1", "s2", "s3"))
  expect_identical(rownames(combined$feature.tab), c("f1", "f2", "f3"))
  expect_equal(combined$feature.tab["f2", "s3"], 0)
  expect_equal(combined$feature.tab["f3", "s2"], 0)

  expect_true(is.data.frame(combined$meta.dat))
  expect_identical(rownames(combined$meta.dat), c("s1", "s2", "s3"))
  expect_identical(colnames(combined$meta.dat), "group")
  expect_identical(combined$meta.dat$group, c("A", "B", "C"))

  expect_true(is.matrix(combined$feature.ann))
  expect_identical(rownames(combined$feature.ann), c("f1", "f2", "f3"))
  expect_identical(colnames(combined$feature.ann), "Phylum")
  expect_identical(as.character(combined$feature.ann[, "Phylum"]), c("p1", "p2", "p3"))
})

test_that("mStat_aggregate_by_taxonomy fails early on disjoint feature IDs", {
  data.obj <- list(
    feature.tab = matrix(
      c(1, 2,
        3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("f1", "f2"), c("s1", "s2"))
    ),
    feature.ann = matrix(
      c("p1", "p2"),
      ncol = 1,
      dimnames = list(c("g1", "g2"), "Phylum")
    )
  )

  expect_error(
    suppressMessages(mStat_aggregate_by_taxonomy(data.obj, feature.level = "Phylum")),
    "No overlapping feature IDs"
  )
})

test_that("mStat_import_dada2_as_data_obj does not require Biostrings attachment", {
  seq_tab <- matrix(
    c(10, 20,
      30, 40),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("s1", "s2"), c("ACGT", "TGCA"))
  )

  tax_tab <- matrix(
    c("k__Bacteria", "k__Bacteria"),
    ncol = 1,
    dimnames = list(c("ACGT", "TGCA"), "Kingdom")
  )

  data.obj <- mStat_import_dada2_as_data_obj(seq_tab = seq_tab, tax_tab = tax_tab)

  expect_identical(rownames(data.obj$feature.tab), c("ASV1", "ASV2"))
  expect_identical(colnames(data.obj$feature.tab), c("s1", "s2"))
  expect_identical(rownames(data.obj$feature.ann), c("ASV1", "ASV2"))
  expect_identical(
    unname(data.obj$feature.ann[, "Kingdom"]),
    c("k__Bacteria", "k__Bacteria")
  )
})

test_that("mStat_import_mothur_as_data_obj creates metadata rownames from feature.tab", {
  shared_file <- tempfile(fileext = ".shared")
  writeLines(
    c(
      "label\tGroup\tnumOtus\tOtu1\tOtu2",
      "0.03\ts1\t2\t1\t0",
      "0.03\ts2\t2\t0\t1"
    ),
    con = shared_file
  )

  data.obj <- suppressMessages(
    mStat_import_mothur_as_data_obj(mothur_shared_file = shared_file)
  )

  expect_true(is.data.frame(data.obj$meta.dat))
  expect_identical(rownames(data.obj$meta.dat), colnames(data.obj$feature.tab))
  expect_equal(ncol(data.obj$meta.dat), 0)
})
