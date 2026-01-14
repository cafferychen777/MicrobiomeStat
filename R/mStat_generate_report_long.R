#' Generate Longitudinal Microbiome Analysis Report
#'
#' Generates a comprehensive PDF/HTML report for longitudinal microbiome analysis
#' including alpha diversity, beta diversity, and taxonomic composition over time.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_test_params_doc
#' @inheritParams mStat_plot_params_doc
#' @param test.adj.vars Character vector of covariate names for statistical adjustment.
#' @param vis.adj.vars Character vector of covariate names to visualize in plots.
#' @param pc.obj Pre-calculated PCoA/PCA results from mStat_calculate_PC. If NULL, computed automatically.
#' @param feature.change.func Method for calculating change: "relative change", "absolute change", or "log fold change".
#' @param bar.area.feature.no Number of top features to show in bar/area plots (default 20).
#' @param heatmap.feature.no Number of top features to show in heatmaps (default 20).
#' @param vis.feature.level Taxonomic level(s) for visualization.
#' @param test.feature.level Taxonomic level(s) for statistical testing.
#' @param feature.analysis.rarafy Logical, whether to rarefy data for feature analysis (default TRUE).
#' @param feature.mt.method Multiple testing correction: "fdr" or "none" (default "fdr").
#' @param feature.sig.level Significance threshold for highlighting features (default 0.1).
#' @param feature.box.axis.transform Y-axis transformation for boxplots: "identity", "sqrt", or "log".
#' @param output.file Output report filename (required).
#' @param output.format Output format: "pdf" or "html".
#'
#' @return A PDF report containing:
#'
#' - Summary statistics table: Sample size, covariates, time points
#'
#' - Alpha diversity boxplots: Boxplots colored by time and groups
#'
#' - Alpha diversity spaghetti plots: Line plots of trajectories
#'
#' - Alpha diversity trends: Tables with LME model results
#'
#' - Alpha volatility: Tables with LM model results
#'
#' - Beta diversity PCoA plots: Ordination plots colored by time and groups
#'
#' - Beta distance boxplots: Boxplots colored by time and groups
#'
#' - Beta spaghetti plots: Individual line plots
#'
#' - Beta diversity trends: Tables with LME results on distances
#'
#' - Beta volatility: Tables with LM results on distances
#'
#' - Feature area plots: Stacked area plots showing composition
#'
#' - Feature heatmaps: Heatmaps colored by relative abundance
#'
#' - Feature volcano plots: With trend/volatility significance
#'
#' - Feature boxplots: Distribution by time and groups
#'
#' - Feature spaghetti plots: Individual line plots
#'
#' @details This function generates a comprehensive longitudinal microbiome report with:
#'
#' 1. Summary Statistics
#' - Table with sample size, number of timepoints, covariates
#'
#' 2. Alpha Diversity Analysis
#' - Boxplots of alpha diversity vs time, colored by groups
#' - Spaghetti plots of alpha trajectories for each subject
#' - Linear mixed effects model results for alpha diversity trend
#' - Linear model results for alpha diversity volatility
#'
#' 3. Beta Diversity Analysis
#' - PCoA plots colored by time points and groups
#' - Boxplots of PCoA coordinate 1 vs time, colored by groups
#' - Spaghetti plots of distance from baseline vs time
#' - Linear mixed effects models for beta diversity distance trend
#' - Linear models for beta diversity volatility
#'
#' 4. Taxonomic Composition Analysis
#' - Stacked area plots of average composition by time and groups
#' - Heatmaps of relative abundance colored from low (blue) to high (red)
#' - Volcano plots highlighting significant taxa in trend and volatility
#' - Boxplots of significant taxa by time and groups
#' - Spaghetti plots for significant taxa vs time
#'
#'
#' @author Chen Yang \email{cafferychen7850@email.com}
#'
#' @seealso \code{\link{mStat_calculate_alpha_diversity}},
#'          \code{\link{mStat_calculate_beta_diversity}},
#'          \code{\link{mStat_calculate_PC}}
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#' mStat_generate_report_long(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   strata.var = NULL,
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","observed_species"),
#'   dist.obj = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   pc.obj = NULL,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.3,
#'   vis.feature.level = c("Family","Genus"),
#'   test.feature.level = c("Family"),
#'   feature.change.func = "relative change",
#'   feature.dat.type = "count",
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   bar.area.feature.no = 40,
#'   heatmap.feature.no = 40,
#'   feature.box.axis.transform = "sqrt",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/per_time_test_report.pdf",
#'   output.format = "pdf"
#' )
#'
#' mStat_generate_report_long(
#'   data.obj = subset_T2D.obj,
#'   group.var = "subject_race",
#'   strata.var = NULL,
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","observed_species"),
#'   dist.obj = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   pc.obj = NULL,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.3,
#'   vis.feature.level = c("Family","Genus"),
#'   test.feature.level = c("Family"),
#'   feature.change.func = "relative change",
#'   feature.dat.type = "count",
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   bar.area.feature.no = 40,
#'   heatmap.feature.no = 40,
#'   feature.box.axis.transform = "sqrt",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/per_time_test_report.html",
#'   output.format = "html"
#' )
#' data(ecam.obj)
#' mStat_generate_report_long(
#'   data.obj = ecam.obj,
#'   group.var = "antiexposedall",
#'   strata.var = NULL,
#'   test.adj.vars = "delivery",
#'   vis.adj.vars = "delivery",
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","observed_species"),
#'   dist.obj = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   pc.obj = NULL,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.3,
#'   vis.feature.level = c("Family","Genus"),
#'   test.feature.level = c("Family"),
#'   feature.change.func = "relative change",
#'   feature.dat.type = "proportion",
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   bar.area.feature.no = 40,
#'   heatmap.feature.no = 40,
#'   feature.box.axis.transform = "sqrt",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/ecam_obj_report_long.pdf"
#' )
#'
#' mStat_generate_report_long(
#'   data.obj = ecam.obj,
#'   group.var = "antiexposedall",
#'   strata.var = NULL,
#'   test.adj.vars = "delivery",
#'   vis.adj.vars = "delivery",
#'   subject.var = "subject.id",
#'   time.var = "month_num",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   alpha.obj = NULL,
#'   alpha.name = c("shannon","observed_species"),
#'   dist.obj = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   pc.obj = NULL,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.3,
#'   vis.feature.level = c("Family","Genus"),
#'   test.feature.level = c("Family"),
#'   feature.change.func = "relative change",
#'   feature.dat.type = "proportion",
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   bar.area.feature.no = 40,
#'   heatmap.feature.no = 40,
#'   feature.box.axis.transform = "sqrt",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   output.file = "/Users/apple/Research/MicrobiomeStat/result/ecam_obj_report_long.html",
#'   output.format = "html"
#' )
#' }
#' @export
mStat_generate_report_long <- function(data.obj,
                                       group.var,
                                       vis.adj.vars = NULL,
                                       test.adj.vars = NULL,
                                       strata.var = NULL,
                                       subject.var,
                                       time.var,
                                       t0.level = NULL,
                                       ts.levels = NULL,
                                       depth = NULL,
                                       alpha.obj = NULL,
                                       alpha.name = c("shannon", "observed_species"),
                                       dist.obj = NULL,
                                       dist.name = c('BC', 'Jaccard'),
                                       pc.obj = NULL,
                                       prev.filter = 0.1,
                                       abund.filter = 0.0001,
                                       bar.area.feature.no = 30,
                                       heatmap.feature.no = 30,
                                       vis.feature.level = NULL,
                                       test.feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       feature.analysis.rarafy = TRUE,
                                       feature.change.func = "relative change",
                                       feature.mt.method = "none",
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

  # Ensure output.format is either "pdf" or "html"
  output.format <- match.arg(output.format)

  # Set pdf to TRUE if output.format is "pdf", FALSE otherwise
  pdf <- output.format == "pdf"

  # Ensure output.file has the correct extension
  output.file <- ensure_file_extension(output.file, output.format)

  # Set result output mode based on format
  result.output <- get_result_output_mode(pdf)

  # Adjust the YAML front matter based on output.format
  if (grepl("\\.pdf$", output.file)) {
    # Check if tinytex is installed and can compile PDFs
    if (!tinytex::is_tinytex()) {
      message("TinyTeX is not installed. PDF output may not work correctly.")
      message("Consider installing TinyTeX with tinytex::install_tinytex()")
    }
    yaml_output <- paste0("output: 
    pdf_document:
    toc: true
    toc_depth: 3
")
  } else {
    yaml_output <- paste0("output: 
    html_document:
    toc: true
    toc_depth: 3
")
  }

  # Build template by combining section generators
  template <- paste0(
    generate_report_yaml_header(yaml_output),
    generate_report_section_overview(),
    generate_report_section_alpha(),
    generate_report_section_beta(),
    generate_report_section_taxa()
  )

  # Template content defined in mStat_generate_report_long_sections.R

rmd_code <- knitr::knit_expand(
  text = template,
  data.obj = data.obj,
  group.var = group.var,
  strata.var = strata.var,
  test.adj.vars = test.adj.vars,
  vis.adj.vars = vis.adj.vars,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  alpha.obj = alpha.obj,
  alpha.name = alpha.name,
  depth = depth,
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
  pdf.hei = pdf.hei,
  output.file = output.file
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
