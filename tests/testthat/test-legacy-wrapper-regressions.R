test_that("plot_taxa cladogram drops unsupported unified params before dispatch", {
  skip_if_not_installed("ggtree")
  skip_if_not_installed("ggtreeExtra")

  data(peerj32.obj)

  result <- suppressWarnings(
    suppressMessages(
      plot_taxa(
        peerj32.obj,
        plot.type = "cladogram",
        group.var = "group",
        strata.var = "sex",
        feature.level = c("Phylum", "Family", "Genus"),
        feature.select = 5,
        theme = list(base.size = 11, choice = "classic")
      )
    )
  )

  expect_type(result, "list")
  expect_true(length(result) > 0)
  expect_true(all(vapply(result, inherits, logical(1), "ggplot")))
})

test_that("legacy alpha wrapper requires precomputed alpha input and stays in-memory", {
  data(peerj32.obj)

  expect_error(
    plot_alpha_diversity(
      meta.dat = peerj32.obj$meta.dat,
      measure = "shannon",
      group.var = "group",
      time.point.plot = "1"
    ),
    "`alpha.obj` is required"
  )

  alpha.obj <- mStat_calculate_alpha_diversity(
    peerj32.obj$feature.tab,
    alpha.name = "shannon"
  )

  tmpdir <- tempfile("plot-alpha-")
  dir.create(tmpdir)
  old_wd <- setwd(tmpdir)
  on.exit(setwd(old_wd), add = TRUE)

  result <- suppressWarnings(
    suppressMessages(
      plot_alpha_diversity(
        alpha.obj = alpha.obj,
        meta.dat = peerj32.obj$meta.dat,
        measure = "shannon",
        group.var = "group",
        time.var = "time",
        time.point.plot = "1"
      )
    )
  )

  expect_type(result, "list")
  expect_length(list.files(pattern = "\\.pdf$"), 0)
})

test_that("legacy beta wrapper defaults to PCoA and stays in-memory", {
  data(peerj32.obj)
  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC"))

  tmpdir <- tempfile("plot-beta-")
  dir.create(tmpdir)
  old_wd <- setwd(tmpdir)
  on.exit(setwd(old_wd), add = TRUE)

  result <- suppressWarnings(
    suppressMessages(
      plot_beta_diversity(
        dist.obj = dist.obj,
        meta.dat = peerj32.obj$meta.dat,
        measure = "BC",
        group.var = "group",
        time.var = "time",
        time.point.plot = "1"
      )
    )
  )

  expect_type(result, "list")
  expect_true("BC" %in% names(result))
  expect_length(list.files(pattern = "\\.pdf$"), 0)
})

test_that("generate_beta_ordination_single can save composed plots to PDF", {
  data(peerj32.obj)
  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC"))

  tmpdir <- tempfile("beta-ordination-save-")
  dir.create(tmpdir)
  old_wd <- setwd(tmpdir)
  on.exit(setwd(old_wd), add = TRUE)

  result <- suppressWarnings(
    generate_beta_ordination_single(
      data.obj = peerj32.obj,
      dist.obj = dist.obj,
      dist.name = "BC",
      time.var = "time",
      t.level = "1",
      group.var = "group",
      pdf = TRUE
    )
  )

  expect_true("BC" %in% names(result))
  expect_true(file.exists("beta_ordination_single_dist.name_BC_time_time_group_group.pdf"))
})

test_that("tree smoothing preserves the caller RNG state", {
  skip_if_not_installed("ape")

  tree <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
  tax.names <- tree$tip.label

  set.seed(123)
  seed_before <- .Random.seed
  smoothing_info <- get_tree_smoothing_info(tree, tax.names, lambda = 0.1, k.neighbors = 2)

  expect_identical(.Random.seed, seed_before)
  expect_equal(smoothing_info$n.matched, length(tax.names))
})

test_that("alpha mixed-model helper keeps every adjustment covariate in fallback formulas", {
  formula <- construct_formula(
    index = "shannon",
    group.var = "group",
    time.var = "time",
    subject.var = "subject",
    adj.vars = c("sex", "age"),
    random_slopes = FALSE
  )

  expect_identical(
    paste(deparse(formula), collapse = " "),
    "shannon ~ group * time + (1 | subject) + sex + age"
  )
})

test_that("beta trend warning names the beta helper", {
  data(peerj32.obj)
  dist.obj <- suppressWarnings(mStat_calculate_beta_diversity(peerj32.obj, dist.name = "BC"))

  expect_message(
    suppressWarnings(
      try(
        generate_beta_trend_test_long(
          data.obj = peerj32.obj,
          dist.obj = dist.obj,
          subject.var = "subject",
          time.var = "time",
          group.var = "group",
          dist.name = "BC"
        ),
        silent = TRUE
      )
    ),
    "generate_beta_trend_test_long"
  )
})

test_that("taxa pair test accepts NULL time without interaction parsing errors", {
  data(peerj32.obj)

  result <- suppressWarnings(
    generate_taxa_test_pair(
      data.obj = peerj32.obj,
      subject.var = "subject",
      time.var = NULL,
      group.var = "group",
      feature.level = "original",
      feature.dat.type = "count"
    )
  )

  expect_true("original" %in% names(result))
  expect_true(length(result$original) > 0)
  expect_true(all(grepl("\\[Main Effect\\]$", names(result$original))))
})

test_that("taxa trend test returns time, group main effects, and interactions", {
  samples <- paste0("s", 1:12)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        5, 4, 6, 5, 7, 6, 8, 7, 9, 8, 10, 9,
        2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6, 5
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), samples)
    ),
    meta.dat = data.frame(
      subject = rep(paste0("id", 1:4), each = 3),
      time = rep(c(0, 1, 2), times = 4),
      group = rep(c("A", "A", "B", "B"), each = 3),
      row.names = samples,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    generate_taxa_trend_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group",
      feature.level = "original",
      feature.dat.type = "count"
    )
  )

  level_names <- names(result$original)
  expect_true("time" %in% level_names)
  expect_true(any(grepl("\\[Main Effect\\]$", level_names)))
  expect_true(any(grepl("\\[Interaction\\]$", level_names)))
})

test_that("taxa trend test keeps continuous group labels explicit", {
  samples <- paste0("s", 1:12)
  data.obj <- list(
    feature.tab = matrix(
      c(
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
        5, 4, 6, 5, 7, 6, 8, 7, 9, 8, 10, 9,
        2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6, 5
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), samples)
    ),
    meta.dat = data.frame(
      subject = rep(paste0("id", 1:4), each = 3),
      time = rep(c(0, 1, 2), times = 4),
      group_num = rep(c(1, 1, 2, 2), each = 3),
      row.names = samples,
      stringsAsFactors = FALSE
    )
  )

  result <- suppressWarnings(
    generate_taxa_trend_test_long(
      data.obj = data.obj,
      subject.var = "subject",
      time.var = "time",
      group.var = "group_num",
      feature.level = "original",
      feature.dat.type = "count"
    )
  )

  level_names <- names(result$original)
  expect_true("time" %in% level_names)
  expect_true("group_num" %in% level_names)
  expect_true("group_num:time" %in% level_names)
  expect_false(any(grepl(" vs ", level_names, fixed = TRUE)))
})

test_that("taxa trend helper returns main and interaction labels for synthetic output", {
  linda_output <- list(
    time = data.frame(log2FoldChange = 0.1, row.names = "tax1"),
    groupB = data.frame(log2FoldChange = 0.2, row.names = "tax1"),
    "groupB:time" = data.frame(log2FoldChange = 0.3, row.names = "tax1")
  )

  extracted <- .mStat_extract_trend_linda_outputs(
    linda_output = linda_output,
    group_var = "group",
    time_var = "time",
    reference_level = "A"
  )

  expect_true(all(c("time", "B vs A (Reference) [Main Effect]", "B vs A (Reference) [Interaction]") %in% names(extracted)))
})

test_that("taxa trend helper keeps continuous group labels explicit", {
  linda_output <- list(
    time = data.frame(log2FoldChange = 0.1, row.names = "tax1"),
    group_num = data.frame(log2FoldChange = 0.2, row.names = "tax1"),
    "group_num:time" = data.frame(log2FoldChange = 0.3, row.names = "tax1")
  )

  extracted <- .mStat_extract_trend_linda_outputs(
    linda_output = linda_output,
    group_var = "group_num",
    time_var = "time",
    reference_level = NULL
  )

  expect_true(all(c("time", "group_num", "group_num:time") %in% names(extracted)))
  expect_false(any(grepl(" vs ", names(extracted), fixed = TRUE)))
})
