load_all_helpers <- function() {
  pkgload::load_all(test_path("..", ".."), helpers = FALSE, quiet = TRUE)
}

expect_plot_renders <- function(plot) {
  file <- tempfile(fileext = ".pdf")

  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplotGrob(plot))

  grDevices::pdf(file = file)
  on.exit(grDevices::dev.off(), add = TRUE)
  on.exit(unlink(file), add = TRUE)
  expect_no_error(print(plot))
}

expect_all_plots_render <- function(x) {
  if (inherits(x, "ggplot")) {
    expect_plot_renders(x)
    return(invisible(NULL))
  }

  if (is.list(x)) {
    for (item in x) {
      expect_all_plots_render(item)
    }
    return(invisible(NULL))
  }

  stop("Expected a ggplot object or nested list of ggplot objects.")
}

test_that("single taxa heatmaps build and render", {
  load_all_helpers()
  data("peerj32.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_heatmap_single(
        data.obj = peerj32.obj,
        time.var = "time",
        t.level = "2",
        group.var = "group",
        strata.var = "sex",
        feature.level = "Family",
        feature.dat.type = "count",
        top.k.plot = 12,
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "Family")
  expect_all_plots_render(plots)
})

test_that("paired taxa heatmaps build and render", {
  load_all_helpers()
  data("subset_pairs.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_heatmap_pair(
        data.obj = subset_pairs.obj,
        subject.var = "MouseID",
        time.var = "Antibiotic",
        group.var = "Sex",
        feature.level = "Genus",
        feature.dat.type = "count",
        top.k.plot = 12,
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "Genus")
  expect_named(plots$Genus, c("indiv", "average"))
  expect_all_plots_render(plots)
})

test_that("longitudinal taxa heatmaps build and render", {
  load_all_helpers()
  data("ecam.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_heatmap_long(
        data.obj = ecam.obj,
        subject.var = "subject.id",
        time.var = "month_num",
        t0.level = "0",
        ts.levels = c("1", "2", "3", "4", "5", "6"),
        group.var = "delivery",
        strata.var = "diet",
        feature.level = "Family",
        feature.dat.type = "count",
        top.k.plot = 12,
        top.k.func = mean,
        abund.filter = 0.001,
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "Family")
  expect_named(plots$Family, c("indiv", "average"))
  expect_all_plots_render(plots)
})

test_that("paired change heatmaps build and render", {
  load_all_helpers()
  data("subset_pairs.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_heatmap_pair(
        data.obj = subset_pairs.obj,
        subject.var = "MouseID",
        time.var = "Antibiotic",
        group.var = "Sex",
        change.base = "Baseline",
        feature.level = "Genus",
        feature.dat.type = "count",
        top.k.plot = 12,
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "Genus")
  expect_named(plots$Genus, c("average", "indiv"))
  expect_all_plots_render(plots)
})

test_that("longitudinal change heatmaps build and render", {
  load_all_helpers()
  data("ecam.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_heatmap_long(
        data.obj = ecam.obj,
        subject.var = "subject.id",
        time.var = "month_num",
        t0.level = "0",
        ts.levels = c("1", "2", "3", "4", "5", "6"),
        group.var = "delivery",
        strata.var = "diet",
        change.base = "0",
        feature.level = "Family",
        feature.dat.type = "count",
        top.k.plot = 12,
        abund.filter = 0.001,
        cluster.rows = FALSE,
        cluster.cols = FALSE,
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "Family")
  expect_all_plots_render(plots)
})

test_that("longitudinal beta change boxplots build with strata faceting", {
  load_all_helpers()
  data("peerj32.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_beta_change_boxplot_long(
        data.obj = peerj32.obj,
        dist.obj = NULL,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        strata.var = "sex",
        adj.vars = "sex",
        t0.level = "1",
        dist.name = "BC",
        base.size = 20,
        theme.choice = "bw",
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "BC")
  expect_all_plots_render(plots)
})

test_that("longitudinal beta ordination renders without grouping", {
  load_all_helpers()
  data("subset_T2D.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_beta_ordination_long(
        data.obj = subset_T2D.obj,
        dist.obj = NULL,
        pc.obj = NULL,
        subject.var = "subject_id",
        time.var = "visit_number_num",
        t0.level = NULL,
        ts.levels = NULL,
        group.var = NULL,
        strata.var = NULL,
        adj.vars = NULL,
        dist.name = "BC",
        base.size = 12,
        theme.choice = "bw",
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "BC")
  expect_all_plots_render(plots)
})

test_that("paired taxa change dotplots render with strata faceting", {
  load_all_helpers()
  data("peerj32.obj", package = "MicrobiomeStat")

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_change_dotplot_pair(
        data.obj = peerj32.obj,
        subject.var = "subject",
        time.var = "time",
        group.var = "group",
        strata.var = "sex",
        change.base = "1",
        feature.change.func = "log fold change",
        feature.level = "Family",
        feature.dat.type = "count",
        top.k.plot = 8,
        top.k.func = "mean",
        prev.filter = 0.01,
        abund.filter = 1e-04,
        base.size = 12,
        theme.choice = "bw",
        pdf = FALSE
      )
    )
  )

  expect_named(plots, "Family")
  expect_all_plots_render(plots)
})

test_that("taxa cladogram single renders with ggplot tile fruit layer", {
  load_all_helpers()
  data("peerj32.obj", package = "MicrobiomeStat")

  test.list <- suppressWarnings(
    suppressMessages(
      generate_taxa_test_single(
        data.obj = peerj32.obj,
        time.var = "time",
        group.var = "group",
        adj.vars = "sex",
        feature.level = c("Phylum", "Family", "Genus"),
        feature.dat.type = "count",
        prev.filter = 0.1,
        abund.filter = 1e-04
      )
    )
  )

  plots <- suppressWarnings(
    suppressMessages(
      generate_taxa_cladogram_single(
        data.obj = peerj32.obj,
        test.list = test.list,
        group.var = "group",
        feature.level = c("Phylum", "Family", "Genus"),
        cutoff = 0.3,
        color.group.level = "Family",
        pdf = FALSE
      )
    )
  )

  expect_true(length(plots) > 0)
  expect_all_plots_render(plots)
})
