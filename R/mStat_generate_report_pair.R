#' Generate a report for microbial ecology analysis of paired data
#'
#' This function generates a comprehensive report for microbial ecology analysis,
#' including changes in alpha diversity, beta diversity, and taxonomic features between paired data.
#' The function is specifically designed for analysis of paired data.
#' 
#' @importFrom pander pander
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var Variable name used for grouping samples.
#' @param strata.var Character, column name in metadata containing stratification variable, e.g "sex". Optional.
#' @param test.adj.vars Character vector, names of columns in the metadata containing covariates to be adjusted for in statistical tests and models, such as linear mixed effects models for longitudinal data analysis. This allows the user to account for the effects of additional variables in assessing the effects of primary variables of interest such as time and groups. Default is NULL, which indicates no covariates are adjusted for in statistical testing.
#' @param vis.adj.vars Character vector, names of columns in the metadata containing covariates to visualize in plots, in addition to the primary variables of interest such as groups. For example, if sex is provided in vis.adj.vars, plots will display facets or colors for different sex groups. This allows visualization of effects across multiple covariates. Default is NULL, which indicates only the primary variables of interest will be visualized without additional covariates.
#' @param subject.var Variable name used for subject identification.
#' @param time.var Character, column name in metadata containing time variable, e.g. "week". Required.
#' @param change.base The base level for calculating changes in paired data.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name Names of alpha diversity indices to include in the analysis. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", "pielou", and "faith_pd".
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. This allows flexible options to quantify change:
#'
#'   - If a function is provided: The function will be applied to compare alpha diversity
#'     at timepoint t vs baseline t0. The function should take two arguments
#'     representing the alpha diversity values at t and t0. For instance, a custom function to
#'     calculate the percentage change might look like:
#'     \preformatted{
#'       percentage_change <- function(t, t0) {
#'         return ((t - t0) / t0) * 100
#'       }
#'     }
#'     You can then pass this function as the value for `alpha.change.func`.
#'
#'   - If a string is provided, the following options are supported:
#'     - 'log fold change': Calculates the log2 fold change of alpha diversity at t compared to t0.
#'     - 'absolute change': Calculates the absolute difference in alpha diversity at t compared to t0.
#'     - Any other value: A warning will be given that the provided method is not recognized,
#'       and the default method ('absolute change') will be used.
#'
#'   - Default behavior (if no recognized string or function is provided) is to compute the absolute difference between t and t0.
#' @param depth An integer specifying the sequencing depth for the "Rarefy" and "Rarefy-TSS" methods.
#' If NULL, no rarefaction is performed.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param pc.obj A list containing the results of dimension reduction/Principal Component Analysis.
#' This should be the output from functions like \code{\link[MicrobiomeStat]{mStat_calculate_PC}},
#' containing the PC coordinates and other metadata. If NULL (default), dimension reduction
#' will be automatically performed using metric multidimensional scaling (MDS) via
#' \code{\link[MicrobiomeStat]{mStat_calculate_PC}}. The pc.obj list structure should contain:
#' \describe{
#'   \item{points}{A matrix with samples as rows and PCs as columns containing the coordinates.}
#'   \item{eig}{Eigenvalues for each PC dimension.}
#'   \item{vectors}{Loadings vectors for features onto each PC.}
#'   \item{Other metadata}{like method, dist.name, etc.}
#' }
#' See \code{\link[MicrobiomeStat]{mStat_calculate_PC}} function for details on output format.
#' @param feature.change.func A function or character string specifying how to calculate
#' the change from baseline value. This allows flexible options:
#' - If a function is provided, it will be applied to each row to calculate change.
#'   The function should take 2 arguments: value at timepoint t and value at baseline t0.
#' - If a character string is provided, following options are supported:
#'   - 'relative change': (value_t - value_t0) / (value_t + value_t0)
#'   - 'absolute change': value_t - value_t0
#'   - 'log fold change': log2(value_t + 1e-5) - log2(value_t0 + 1e-5)
#' - Default is 'relative change'.
#'
#' If none of the above options are matched, an error will be thrown indicating
#' the acceptable options or prompting the user to provide a custom function.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param bar.area.feature.no A numeric value indicating the number of top abundant features to retain in both barplot and areaplot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 40.
#' @param heatmap.feature.no A numeric value indicating the number of top abundant features to retain in the heatmap. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 40.
#' @param dotplot.feature.no A numeric value indicating the number of top abundant features to retain in the dotplot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 40.
#' @param vis.feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for visualization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param test.feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for testing or analytical purposes. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be analyzed separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param feature.analysis.rarafy Logical, indicating whether to rarefy the data at the feature-level for analysis.
#' If TRUE, the feature data will be rarefied before analysis. Default is TRUE.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param feature.mt.method Character, multiple testing method for features, "fdr" or "none", default is "fdr".
#' @param feature.sig.level Numeric, significance level cutoff for highlighting features, default is 0.1.
#' @param feature.box.axis.transform A string indicating the transformation to apply to the y-axis of the feature's boxplot visualization before plotting. Options are:
#' - "identity": No transformation (default).
#' - "sqrt": Square root transformation.
#' - "log": Logarithmic transformation. Zeros are replaced with half of the minimum non-zero value for each taxon before log transformation.
#' @param base.size Base font size for the generated plots.
#' @param theme.choice
#' Plot theme choice. Specifies the visual style of the plot. Can be one of the following pre-defined themes:
#'   - "prism": Utilizes the ggprism::theme_prism() function from the ggprism package, offering a polished and visually appealing style.
#'   - "classic": Applies theme_classic() from ggplot2, providing a clean and traditional look with minimal styling.
#'   - "gray": Uses theme_gray() from ggplot2, which offers a simple and modern look with a light gray background.
#'   - "bw": Employs theme_bw() from ggplot2, creating a classic black and white plot, ideal for formal publications and situations where color is best minimized.
#'   - "light": Implements theme_light() from ggplot2, featuring a light theme with subtle grey lines and axes, suitable for a fresh, modern look.
#'   - "dark": Uses theme_dark() from ggplot2, offering a dark background, ideal for presentations or situations where a high-contrast theme is desired.
#'   - "minimal": Applies theme_minimal() from ggplot2, providing a minimalist theme with the least amount of background annotations and colors.
#'   - "void": Employs theme_void() from ggplot2, creating a blank canvas with no axes, gridlines, or background, ideal for custom, creative plots.
#' Each theme option adjusts various elements like background color, grid lines, and font styles to match the specified aesthetic.
#' Default is "bw", offering a universally compatible black and white theme suitable for a wide range of applications.
#' @param custom.theme
#' A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. Custom themes are useful for creating publication-ready figures with specific formatting requirements.
#'
#' To use a custom theme, create a theme object with ggplot2::theme(), including any desired customizations. Common customizations for publication-ready figures might include adjusting text size for readability, altering line sizes for clarity, and repositioning or formatting the legend. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=14, face="bold"),        # Bold axis titles with larger font
#'   axis.text = ggplot2::element_text(size=12),                      # Slightly larger axis text
#'   legend.position = "top",                                         # Move legend to the top
#'   legend.background = ggplot2::element_rect(fill="lightgray"),     # Light gray background for legend
#'   panel.background = ggplot2::element_rect(fill="white", colour="black"), # White panel background with black border
#'   panel.grid.major = ggplot2::element_line(colour = "grey90"),     # Lighter color for major grid lines
#'   panel.grid.minor = ggplot2::element_blank(),                     # Remove minor grid lines
#'   plot.title = ggplot2::element_text(size=16, hjust=0.5)           # Centered plot title with larger font
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. If `custom.theme` is NULL (the default), the theme is determined by `theme.choice`. This flexibility allows for both easy theme selection for general use and detailed customization for specific presentation or publication needs.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf Logical indicating whether to save plots as PDF files (default: TRUE).
#' @param file.ann Annotation text for the PDF file names.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param output.file A character string specifying the output file name for the report. The file extension
#' (.pdf or .html) will be automatically added based on the output.format if not already present.
#' This parameter also determines the title of the generated report. For example, if you set it to
#' "path_to_your_location/report", the title of the report will be "report".
#' @param output.format A character string specifying the desired output format of the report.
#' Must be either "pdf" or "html". Default is c("pdf", "html"), which will use the first value ("pdf")
#' if not explicitly specified. This parameter determines whether the report will be generated as a PDF
#' or HTML document.
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
