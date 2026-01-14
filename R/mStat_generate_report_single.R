#' Generate Single Time Point Microbiome Analysis Report
#'
#' Generates a comprehensive PDF/HTML report for cross-sectional microbiome analysis
#' including alpha diversity, beta diversity, and taxonomic composition.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param test.adj.vars Character vector of covariate names for statistical adjustment.
#' @param vis.adj.vars Character vector of covariate names to visualize in plots.
#' @param t.level Time point to subset data to (if time.var provided). Default NULL uses all data.
#' @param pc.obj Pre-calculated PCoA/PCA results from mStat_calculate_PC. If NULL, computed automatically.
#' @param bar.area.feature.no Number of top features to show in bar/area plots (default 40).
#' @param heatmap.feature.no Number of top features to show in heatmaps (default 40).
#' @param dotplot.feature.no Number of top features to show in dotplots (default 40).
#' @param vis.feature.level Taxonomic level(s) for visualization.
#' @param test.feature.level Taxonomic level(s) for statistical testing.
#' @param feature.analysis.rarafy Logical, whether to rarefy data for feature analysis (default TRUE).
#' @param feature.mt.method Multiple testing correction: "fdr" or "none" (default "fdr").
#' @param feature.sig.level Significance threshold for highlighting features (default 0.1).
#' @param feature.box.axis.transform Y-axis transformation for boxplots: "identity", "sqrt", or "log".
#' @param output.file Output report filename (required).
#' @param output.format Output format: "pdf" or "html".
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A report file containing the microbial ecology analysis results.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' library(aplot)
#' data(peerj32.obj)
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   vis.adj.vars = NULL,
#'   test.adj.vars = NULL,
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = NULL,
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = c("Phylum", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/peerj32.obj_report.pdf",
#'   output.format = c("pdf")
#' )
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   vis.adj.vars = NULL,
#'   test.adj.vars = NULL,
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = NULL,
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = c("Phylum", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/peerj32.obj_report.html",
#'   output.format = c("html")
#' )
#' data(subset_T2D.obj)
#' mStat_generate_report_single(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "subject_race",
#'   vis.adj.vars = "sample_body_site",
#'   test.adj.vars = "sample_body_site",
#'   time.var = "visit_number_num",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = 2000,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = 1,
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "subject_gender",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/subset_T2D.obj_report.pdf"
#' )
#' mStat_generate_report_single(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "subject_race",
#'   vis.adj.vars = "sample_body_site",
#'   test.adj.vars = "sample_body_site",
#'   time.var = "visit_number_num",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = 2000,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = 1,
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "subject_gender",
#'   vis.feature.level = c("Order", "Family", "Genus"),
#'   test.feature.level = c("Order", "Family", "Genus"),
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/subset_T2D.obj_report.html",
#'   output.format = c("html")
#' )
#' mStat_generate_report_single(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "subject_race",
#'   vis.adj.vars = "sample_body_site",
#'   test.adj.vars = "sample_body_site",
#'   time.var = "visit_number_num",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = 2000,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = 1,
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "subject_gender",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/report.pdf"
#' )
#' mStat_generate_report_single(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "subject_race",
#'   vis.adj.vars = "sample_body_site",
#'   test.adj.vars = "sample_body_site",
#'   time.var = "visit_number_num",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = 2000,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = 1,
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "subject_gender",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/report.pdf"
#' )
#' data(ecam.obj)
#' mStat_generate_report_single(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "delivery",
#'   vis.adj.vars = "diet",
#'   test.adj.vars = "diet",
#'   time.var = "month_num",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = 1,
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "antiexposedall",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "proportion",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/ecam.obj_report.pdf"
#' )
#' mStat_generate_report_single(
#'   data.obj = ecam.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "delivery",
#'   vis.adj.vars = "diet",
#'   test.adj.vars = "diet",
#'   time.var = "month_num",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = 1,
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "antiexposedall",
#'   vis.feature.level = c("Order", "Family", "Genus"),
#'   test.feature.level = c("Order", "Family", "Genus"),
#'   feature.dat.type = "proportion",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/ecam.obj_report.html",
#'   output.format = c("html")
#' )
#' }
#' @export
mStat_generate_report_single <- function(data.obj,
                                         group.var,
                                         vis.adj.vars = NULL,
                                         test.adj.vars = NULL,
                                         strata.var = NULL,
                                         time.var = NULL,
                                         t.level = NULL,
                                         alpha.obj = NULL,
                                         alpha.name = c("shannon", "observed_species"),
                                         depth = NULL,
                                         dist.obj = NULL,
                                         dist.name = c('BC', 'Jaccard'),
                                         pc.obj = NULL,
                                         prev.filter = 0.1,
                                         abund.filter = 0.0001,
                                         bar.area.feature.no = 40,
                                         heatmap.feature.no = 40,
                                         dotplot.feature.no = 40,
                                         vis.feature.level = NULL,
                                         test.feature.level = NULL,
                                         feature.dat.type = c("count", "proportion", "other"),
                                         feature.analysis.rarafy = TRUE,
                                         feature.mt.method = c("fdr", "none"),
                                         feature.sig.level = 0.1,
                                         feature.box.axis.transform = c("sqrt"),
                                         base.size = 16,
                                         theme.choice = "bw",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         output.file,
                                         output.format = c("pdf", "html"),
                                         ...) {

  output.format <- match.arg(output.format)
  feature.dat.type <- match.arg(feature.dat.type)
  feature.mt.method <- match.arg(feature.mt.method)

  # Calculate sample count for determining plot width
  sample_count <- calc_sample_count(data.obj$meta.dat, time.var, group.var)
  has_many_samples <- sample_count > 8

  if (has_many_samples) {
    width <- 7
  } else {
    width <- 3
  }

  # Ensure output.format is either "pdf" or "html"
  output.format <- match.arg(output.format)

  # Set pdf to TRUE if output.format is "pdf", FALSE otherwise
  pdf <- output.format == "pdf"

  # Ensure output.file has the correct extension
  output.file <- ensure_file_extension(output.file, output.format)

  # Set result output mode based on format
  result.output <- get_result_output_mode(pdf)

  # Build YAML output format string
  if (output.format == "pdf") {
    yaml_output <- "
output:
  pdf_document:
    toc: true
    toc_depth: 3
    latex_engine: lualatex
"
  } else {
    yaml_output <- "
output:
  html_document:
    toc: true
    toc_depth: 3
"
  }

  # Build template from section generators
  template <- paste0(
    generate_single_report_yaml_header(yaml_output),
    generate_single_report_section_overview(),
    generate_single_report_section_alpha(),
    generate_single_report_section_beta(),
    generate_single_report_section_taxa()
  )

  rmd_code <- knitr::knit_expand(
                          text = template,
                          data.obj = data.obj,
                          group.var = group.var,
                          vis.adj.vars = vis.adj.vars,
                          test.adj.vars = test.adj.vars,
                          strata.var = strata.var,
                          time.var = time.var,
                          t.level = t.level,
                          alpha.obj = alpha.obj,
                          alpha.name = alpha.name,
                          depth = depth,
                          dist.obj = dist.obj,
                          dist.name = dist.name,
                          pc.obj = pc.obj,
                          prev.filter = prev.filter,
                          abund.filter = abund.filter,
                          bar.area.feature.no = bar.area.feature.no,
                          heatmap.feature.no = heatmap.feature.no,
                          dotplot.feature.no = dotplot.feature.no,
                          vis.feature.level = vis.feature.level,
                          test.feature.level = test.feature.level,
                          feature.dat.type = feature.dat.type,
                          feature.analysis.rarafy = feature.analysis.rarafy,
                          feature.mt.method = feature.mt.method,
                          feature.sig.level = feature.sig.level,
                          feature.box.axis.transform = feature.box.axis.transform,
                          base.size = base.size,
                          theme.choice = theme.choice,
                          custom.theme = custom.theme,
                          palette = palette,
                          pdf = pdf,
                          file.ann = file.ann,
                          pdf.wid = pdf.wid,
                          pdf.hei = pdf.hei)

  rmd_file <- tempfile(fileext = ".Rmd")
  writeLines(rmd_code, con = rmd_file)

  report_file <-
    rmarkdown::render(input = rmd_file,
                      output_file = output.file,
                      quiet = FALSE,
                      clean = FALSE)

  return(report_file)
}
