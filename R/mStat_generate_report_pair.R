#' Generate Paired Microbiome Analysis Report
#'
#' Generates a comprehensive PDF/HTML report for paired (two time point) microbiome analysis
#' including changes in alpha diversity, beta diversity, and taxonomic composition.
#'
#' @importFrom pander pander
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param test.adj.vars Character vector of covariate names for statistical adjustment.
#' @param vis.adj.vars Character vector of covariate names to visualize in plots.
#' @param change.base The baseline time point for calculating changes.
#' @param alpha.change.func Method for calculating alpha diversity change: "log fold change" or "absolute change".
#' @param pc.obj Pre-calculated PCoA/PCA results from mStat_calculate_PC. If NULL, computed automatically.
#' @param feature.change.func Method for calculating change: "relative change", "absolute change", or "log fold change".
#' @param bar.area.feature.no Number of top features to show in bar/area plots (default 30).
#' @param heatmap.feature.no Number of top features to show in heatmaps (default 30).
#' @param dotplot.feature.no Number of top features to show in dotplots (default 30).
#' @param vis.feature.level Taxonomic level(s) for visualization.
#' @param test.feature.level Taxonomic level(s) for statistical testing.
#' @param feature.analysis.rarafy Logical, whether to rarefy data for feature analysis (default TRUE).
#' @param feature.mt.method Multiple testing correction: "fdr" or "none" (default "fdr").
#' @param feature.sig.level Significance threshold for highlighting features (default 0.1).
#' @param feature.box.axis.transform Y-axis transformation for boxplots: "identity", "sqrt", or "log".
#' @param output.file Output report filename (required).
#' @param output.format Output format: "pdf" or "html".
#'
#' @return A report file containing the microbial ecology analysis results for paired data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' mStat_generate_report_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   dist.name = c("BC",'Jaccard'),
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   strata.var = NULL,
#'   vis.feature.level = c("Phylum","Family","Genus"),
#'   test.feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
#'   theme.choice = "bw",
#'   base.size = 18,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/pair_peerj32.obj_report.pdf",
#'   output.format = "pdf"
#' )
#' mStat_generate_report_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   dist.name = c("BC",'Jaccard'),
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   strata.var = NULL,
#'   vis.feature.level = c("Phylum","Family","Genus"),
#'   test.feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
#'   theme.choice = "bw",
#'   base.size = 18,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/pair_peerj32.obj_report.html",
#'   output.format = "html"
#' )
#' mStat_generate_report_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   dist.name = c("BC",'Jaccard'),
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   strata.var = NULL,
#'   vis.feature.level = c("Phylum","Family","Genus"),
#'   test.feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
#'   theme.choice = "bw",
#'   base.size = 18,
#'   output.file = "/Users/apple/MicrobiomeStat/report.pdf"
#' )
#' data(peerj32.obj)
#' mStat_generate_report_pair(
#'   data.obj = peerj32.obj,
#'   group.var = "group",
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   strata.var = "sex",
#'   subject.var = "subject",
#'   time.var = "time",
#'   change.base = "1",
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   dist.obj = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   feature.change.func = "relative change",
#'   vis.feature.level = c("Genus"),
#'   test.feature.level = c("Genus"),
#'   bar.area.feature.no = 30,
#'   heatmap.feature.no = 30,
#'   dotplot.feature.no = 20,
#'   feature.dat.type = "count",
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
#'   theme.choice = "bw",
#'   base.size = 18,
#'   output.file = "/Users/apple/MicrobiomeStat/report.pdf")
#'
#' data(subset_pairs.obj)
#' mStat_generate_report_pair(
#'   data.obj = subset_pairs.obj,
#'   group.var = "Sex",
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   strata.var = NULL,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   change.base = "Baseline",
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon", "observed_species"),
#'   dist.obj = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   feature.change.func = "relative change",
#'   vis.feature.level = c("Genus", "Family"),
#'   test.feature.level = c("Genus", "Family"),
#'   bar.area.feature.no = 30,
#'   heatmap.feature.no = 30,
#'   dotplot.feature.no = 20,
#'   feature.dat.type = "count",
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
#'   theme.choice = "bw",
#'   base.size = 18,
#'   output.file = "/Users/apple/MicrobiomeStat/report.pdf"
#' )
#'
#' }
#' @export
mStat_generate_report_pair <- function(data.obj,
                                       group.var = NULL,
                                       strata.var = NULL,
                                       test.adj.vars = NULL,
                                       vis.adj.vars = NULL,
                                       subject.var,
                                       time.var,
                                       change.base,
                                       alpha.obj = NULL,
                                       alpha.name = c("shannon", "observed_species"),
                                       alpha.change.func = "log fold change",
                                       depth = NULL,
                                       dist.obj = NULL,
                                       dist.name = c('BC', 'Jaccard'),
                                       pc.obj = NULL,
                                       prev.filter = 0.1,
                                       abund.filter = 0.0001,
                                       bar.area.feature.no = 30,
                                       heatmap.feature.no = 30,
                                       dotplot.feature.no = 30,
                                       vis.feature.level,
                                       test.feature.level,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       feature.analysis.rarafy = TRUE,
                                       feature.change.func = "relative change",
                                       feature.mt.method = c("fdr"),
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
                                       output.format = c("pdf", "html")) {

  # Use pander package to ensure it is actually loaded
  pander::pander(list("MicrobiomeStat Report Generation", Sys.time(), "Parameters", paste("Output Format:", output.format)))
  
  # Ensure output.format is either "pdf" or "html"
  output.format <- match.arg(output.format)
  feature.dat.type <- match.arg(feature.dat.type)
  feature.mt.method <- match.arg(feature.mt.method)

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
    generate_pair_report_yaml_header(yaml_output),
    generate_pair_report_section_overview(),
    generate_pair_report_section_alpha(),
    generate_pair_report_section_beta(),
    generate_pair_report_section_taxa()
  )

rmd_code <- knitr::knit_expand(
  text = template,
  data.obj = data.obj,
  group.var = group.var,
  vis.adj.vars = vis.adj.vars,
  test.adj.vars = test.adj.vars,
  strata.var = strata.var,
  subject.var = subject.var,
  time.var = time.var,
  change.base = change.base,
  alpha.obj = alpha.obj,
  alpha.name = alpha.name,
  depth = depth,
  alpha.change.func = alpha.change.func,
  dist.obj = dist.obj,
  dist.name = dist.name,
  pc.obj = pc.obj,
  prev.filter = prev.filter,
  abund.filter = abund.filter,
  feature.dat.type = feature.dat.type,
  feature.analysis.rarafy = feature.analysis.rarafy,
  feature.change.func = feature.change.func,
  bar.area.feature.no = bar.area.feature.no,
  heatmap.feature.no = heatmap.feature.no,
  dotplot.feature.no = dotplot.feature.no,
  vis.feature.level = vis.feature.level,
  test.feature.level = test.feature.level,
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
  pdf.hei = pdf.hei
)

rmd_file <- tempfile(fileext = ".Rmd")
writeLines(rmd_code, con = rmd_file)

report_file <-
  rmarkdown::render(input = rmd_file,
                    output_file = output.file,
                    quiet = FALSE,
                    clean = FALSE)

return(report_file)
}
