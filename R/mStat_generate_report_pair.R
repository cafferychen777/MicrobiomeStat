#' Generate a report for microbial ecology analysis of paired data
#'
#' This function generates a comprehensive report for microbial ecology analysis,
#' including changes in alpha diversity, beta diversity, and taxonomic features between paired data.
#' The function is specifically designed for analysis of paired data.
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
#' @param alpha.name Names of alpha diversity indices to include in the analysis.
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
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param pc.obj A list containing the results of dimension reduction/Principal Component Analysis.
#' This should be the output from functions like \code{\link[MicrobiomeStat]{mStat_calculate_PC}}, containing the PC coordinates and other metadata.
#' If NULL (default), dimension reduction will be automatically performed using metric multidimensional scaling (MDS) via \code{\link[MicrobiomeStat]{mStat_calculate_PC}}.
#' The pc.obj list structure should contain:
#' \itemize{
#'  \item{$points:}{A matrix with samples as rows and PCs as columns containing the coordinates.}
#'  \item{$eig:}{Eigenvalues for each PC dimension.}
#'  \item{$vectors:}{Loadings vectors for features onto each PC.}
#'  \item{Other metadata like $method, $dist.name, etc.}
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
#' @param theme.choice Plot theme choice. Can be one of:
#'   - "prism": ggprism::theme_prism()
#'   - "classic": theme_classic()
#'   - "gray": theme_gray()
#'   - "bw": theme_bw()
#' Default is "bw".
#' @param custom.theme A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. To use a custom theme, first create a theme object with ggplot2::theme(), then pass it to this argument. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=16, color="red"),
#'   legend.position = "none"
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. Default is NULL, which will use the default theme based on `theme.choice`.
#' @param palette Color palette used for the plots.
#' @param pdf Logical indicating whether to save plots as PDF files (default: TRUE).
#' @param file.ann Annotation text for the PDF file names.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param output.file Output file name for the report.
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
#'   test.adj.vars = c("sex"),
#'   vis.adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   dist.name = c("BC",'Jaccard'),
#'   change.base = "1",
#'   feature.change.func = "relative change",
#'   strata.var = "sex",
#'   vis.feature.level = c("Phylum","Family","Genus"),
#'   test.feature.level = c("Genus"),
#'   feature.dat.type = "count",
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
#'   theme.choice = "bw",
#'   base.size = 18,
#'   output.file = "path/to/report.pdf"
#' )
#'
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
#'   output.file = "/Users/apple/MicrobiomeStat/report.pdf"
#' )
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
                                       feature.change.func = "log fold change",
                                       feature.mt.method = c("fdr"),
                                       feature.sig.level = 0.1,
                                       feature.box.axis.transform = c("identity"),
                                       base.size = 16,
                                       theme.choice = "bw",
                                       custom.theme = NULL,
                                       palette = NULL,
                                       pdf = TRUE,
                                       file.ann = NULL,
                                       pdf.wid = 11,
                                       pdf.hei = 8.5,
                                       output.file) {
  template <- "
---
title: '`r sub(\".pdf$\", \"\", basename(output.file))`'
author: '[Powered by MicrobiomeStat (Ver 1.1.1)](http://www.microbiomestat.wiki)'
date: '`r Sys.Date()`'
output:
  pdf_document:
    toc: true
    toc_depth: 3
    latex_engine: lualatex
---

# 1. Data overview and summary statistics

## 1.1 Parameter setting

```{r input-parameters-summary, echo=FALSE, message=FALSE, results='asis'}

custom_depth_status <- ifelse(is.null(depth), 'NULL', toString(depth))

custom_theme_status <- ifelse(is.null(custom.theme), 'NULL', 'Not NULL')

custom_palette_status <- ifelse(is.null(palette), 'NULL', 'Not NULL')

custom_file.ann_status <- ifelse(is.null(file.ann), 'NULL', 'Not NULL')

custom_change.base_status <- ifelse(is.null(change.base), 'NULL', change.base)

custom_test.adj.vars_status <- ifelse(is.null(test.adj.vars), 'NULL', toString(test.adj.vars))

custom_vis.adj.vars_status <- ifelse(is.null(vis.adj.vars), 'NULL', toString(vis.adj.vars))

params_data <- data.frame(Parameter = c('data.obj',
                                        'feature.dat.type',
                                        'group.var',
                                        'test.adj.vars',
                                        'vis.adj.vars',
                                        'strata.var',
                                        'subject.var',
                                        'time.var',
                                        'change.base',
                                        'alpha.obj',
                                        'alpha.name',
                                        'alpha.change.func',
                                        'depth',
                                        'dist.obj',
                                        'dist.name',
                                        'pc.obj',
                                        'prev.filter',
                                        'abund.filter',
                                        'feature.analysis.rarafy',
                                        'feature.change.func',
                                        'bar.area.feature.no',
                                        'heatmap.feature.no',
                                        'dotplot.feature.no',
                                        'vis.feature.level',
                                        'test.feature.level',
                                        'feature.mt.method',
                                        'feature.sig.level',
                                        'feature.box.axis.transform',
                                        'base.size',
                                        'theme.choice',
                                        'custom.theme',
                                        'palette',
                                        'pdf',
                                        'file.ann',
                                        'pdf.wid',
                                        'pdf.hei'),
                          Value =       c(deparse(substitute(data.obj)),
                                        feature.dat.type,
                                        group.var,
                                        custom_test.adj.vars_status,
                                        custom_vis.adj.vars_status,
                                        strata.var,
                                        subject.var,
                                        time.var,
                                        custom_change.base_status,
                                        deparse(substitute(alpha.obj)),
                                        toString(alpha.name),
                                        alpha.change.func,
                                        custom_depth_status,
                                        deparse(substitute(dist.obj)),
                                        toString(dist.name),
                                        deparse(substitute(pc.obj)),
                                        prev.filter,
                                        abund.filter,
                                        feature.analysis.rarafy,
                                        feature.change.func,
                                        bar.area.feature.no,
                                        heatmap.feature.no,
                                        dotplot.feature.no,
                                        toString(vis.feature.level),
                                        toString(test.feature.level),
                                        feature.mt.method,
                                        feature.sig.level,
                                        feature.box.axis.transform,
                                        base.size,
                                        theme.choice,
                                        custom_theme_status,
                                        custom_palette_status,
                                        pdf,
                                        custom_file.ann_status,
                                        pdf.wid,
                                        pdf.hei))

pander::pander(params_data)
```

## 1.2 Summary statistics

```{r mStat-data-summary, message=FALSE}
mStat_results <- mStat_summarize_data_obj(data.obj = data.obj,
                                          time.var = time.var,
                                          group.var = group.var,
                                          palette = palette)
```

```{r mStat-data-summary-print, echo=FALSE, message=FALSE, results='asis'}
# Display the results
pander::pander(mStat_results)
```

```{r object-pre-calculation, echo=FALSE, message=FALSE, results='asis'}

original.data.obj <- data.obj

rarefy.data.obj <- mStat_normalize_data(data.obj = data.obj, method = 'Rarefy', depth = depth)$data.obj.norm

if (is.null(depth)){
  depth <- min(colSums(data.obj$feature.tab))
  cat(sprintf('No rarefaction depth is specified. The minimum depth, %d, is used as the rarefaction depth. ', depth))
}

unique_levels <- unique(c(vis.feature.level, test.feature.level))

if (!all(unique_levels %in% names(data.obj$feature.agg.list))) {
  original.data.obj <- mStat_aggregate_by_taxonomy(original.data.obj, feature.level = unique_levels)
  rarefy.data.obj <- mStat_aggregate_by_taxonomy(rarefy.data.obj, feature.level = unique_levels)
}

# Now, let's tell the user how many samples remain after rarefaction.
  remaining_samples <- ncol(rarefy.data.obj$feature.tab)
  cat(sprintf('After rarefaction, %d samples remain in the analysis. ', remaining_samples))

if (is.null(alpha.obj)){
  alpha.obj <- mStat_calculate_alpha_diversity(x = rarefy.data.obj$feature.tab, alpha.name = alpha.name)
  cat('alpha.obj is calculated based on the rarefied data.obj. ')
}

if (is.null(dist.obj)){
  dist.obj <- mStat_calculate_beta_diversity(data.obj = rarefy.data.obj, dist.name = dist.name)
  cat('dist.obj is calculated based on the rarefied data.obj.\n')
}

if (is.null(pc.obj)){
  pc.obj <- mStat_calculate_PC(dist.obj = dist.obj, dist.name = dist.name)
  cat('pc.obj is calculated based on the dist.obj using multi-dimensional scaling.\n')
}

```

## 1.3 Data visualization (overall)

```{r Check-and-Select-Rarefaction1, echo=FALSE, message=FALSE, results='asis'}

if (feature.analysis.rarafy) {
    cat('Rarefaction has been enabled for feature-level analysis and visualization.\n\n',
        'Reason: The observed abundance of rare/low-abundance features can be strongly influenced by the sequencing depth. ',
        'Rarefaction is an effective method to control the effect of sequencing depth variation. ',
        'By employing rarefaction, presence/absence status of the features are more comparable and we can potentially increase the power of detecting those rare/low-abundance features, ',
        'even though it introduces some variation due to under-sampling.\n',
        'In essence, this step improves comparability across samples across samples with varying sequencing depth.\n\n',
        'If you do not wish to perform rarefaction during feature-level analysis, please turn feature.analysis.rarafy to FALSE.\n')

    data.obj <- rarefy.data.obj
} else {
  data.obj <- original.data.obj
}

```

### 1.3.1 Feature heatmap

```{r taxa-heatmap-generation, message=FALSE, fig.align='center', fig.width = 16, fig.height = 12, results='asis'}
taxa_heatmap_pair_results <- generate_taxa_heatmap_pair(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = vis.feature.level,
  feature.dat.type = feature.dat.type,
  features.plot = NULL,
  top.k.plot = heatmap.feature.no,
  top.k.func = 'mean',
  prev.filter = prev.filter,
  abund.filter = abund.filter,
  base.size = base.size,
  palette = palette,
  cluster.rows = NULL,
  cluster.cols = NULL,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)

```

```{r taxa-heatmap-pair-avergae-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 12, warning = FALSE}
cat('The following plots display the average proportions for each time point, group, and stratum. \n\n')

indiv_list <- lapply(taxa_heatmap_pair_results, function(x) x$indiv)

average_list <- lapply(taxa_heatmap_pair_results, function(x) x$average)

average_list
```

```{r taxa-heatmap-pair-indiv-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 30, fig.height = 15, warning = FALSE}
cat('The following plots display the individual proportions for each sample. \n\n')
indiv_list
```

### 1.3.2 Feature barplot

```{r taxa-barplot-generation, message=FALSE, fig.align='center', fig.width = 30, fig.height = 15, warning = FALSE, results='asis'}
taxa_barplot_pair_results <- generate_taxa_barplot_pair(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = vis.feature.level,
  feature.dat.type = feature.dat.type,
  feature.number = bar.area.feature.no,
  base.size = base.size,
  theme.choice = theme.choice,
  custom.theme = custom.theme,
  palette = palette,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)
```

```{r taxa-barplot-pair-avergae-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 20, warning = FALSE}
cat('The following plots display the average proportions for each time point, group, and stratum. \n\n')

indiv_list <- lapply(taxa_barplot_pair_results, function(x) x$indiv)

average_list <- lapply(taxa_barplot_pair_results, function(x) x$average)

average_list
```

```{r taxa-barplot-pair-indiv-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 20, warning = FALSE}
cat('The following plots display the individual proportions for each sample. \n\n')
indiv_list
```

### 1.3.3 Feature dotplot

```{r taxa-dotplot-generation, message=FALSE, fig.align='center', fig.width = 25, fig.height = 12, results='asis'}
taxa_dotplot_results <- generate_taxa_dotplot_pair(
                                              data.obj = data.obj,
                                              subject.var = subject.var,
                                              time.var = time.var,
                                              group.var = group.var,
                                              strata.var = strata.var,
                                              feature.level = vis.feature.level,
                                              feature.dat.type = feature.dat.type,
                                              features.plot = NULL,
                                              top.k.plot = dotplot.feature.no,
                                              top.k.func = 'mean',
                                              prev.filter = prev.filter,
                                              abund.filter = abund.filter,
                                              base.size = base.size,
                                              theme.choice = theme.choice,
                                              custom.theme = custom.theme,
                                              palette = palette,
                                              pdf = pdf,
                                              file.ann = file.ann,
                                              pdf.wid = pdf.wid,
                                              pdf.hei = pdf.hei)
```

```{r taxa-dotplot-pair-avergae-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 12}
cat('The following plots display the average proportions for each time point, group, and stratum. \n\n')
taxa_dotplot_results
```

### 1.3.4 Feature change heatmap

```{r taxa-change-heatmap, message=FALSE, fig.align='center', fig.width = 25, fig.height = 20, results='asis'}
taxa_change_heatmap_results <- generate_taxa_change_heatmap_pair(
                                      data.obj = data.obj,
                                      subject.var = subject.var,
                                      time.var = time.var,
                                      change.base = change.base,
                                      group.var = group.var,
                                      strata.var = strata.var,
                                      feature.level = vis.feature.level,
                                      feature.change.func = feature.change.func,
                                      feature.dat.type = feature.dat.type,
                                      features.plot = NULL,
                                      top.k.plot = heatmap.feature.no,
                                      top.k.func = 'mean',
                                      prev.filter = prev.filter,
                                      abund.filter = abund.filter,
                                      base.size = base.size,
                                      palette = palette,
                                      cluster.cols = NULL,
                                      cluster.rows = NULL,
                                      pdf = pdf,
                                      file.ann = file.ann,
                                      pdf.wid = pdf.wid,
                                      pdf.hei = pdf.hei)
```

```{r taxa-change-heatmap-pair-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 20}
if (is.function(feature.change.func)) {
  cat('The changes were computed using a custom function provided by the user.')
} else if (feature.change.func == 'relative change') {
  cat('The changes were relative changes, which were computed as (after.abund - before.abund) / (after.abund + before.abund) so the values lie between [-1, 1].')
} else if (feature.change.func == 'absolute change') {
  cat('The changes were absolute changes, computed as the difference between after.abund and before.abund.')
} else if (feature.change.func == 'log fold change') {
  cat('The changes were log2 fold changes, computed as the logarithm of the ratio of after.abund to before.abund, with a small constant added to avoid taking the log of zero.')
}


cat(' The following plots display the change for each subject. \n\n')
indiv_list <- lapply(taxa_change_heatmap_results, function(x) x$indiv)

average_list <- lapply(taxa_change_heatmap_results, function(x) x$average)

indiv_list
```

```{r taxa-change-heatmap-pair-print-indiv, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 20}

cat(' The following plots display the average change for each time point, group, and stratum. \n\n')
average_list

```

### 1.3.5 Feature change dotplot

```{r taxa-change-dotplot, message=FALSE, fig.align='center', fig.width = 25, fig.height = 12, results='asis'}
taxa_change_dotplot_results <- generate_taxa_change_dotplot_pair(
                                             data.obj = data.obj,
                                             subject.var = subject.var,
                                             time.var = time.var,
                                             change.base = change.base,
                                             feature.change.func = feature.change.func,
                                             group.var = group.var,
                                             strata.var = strata.var,
                                             feature.level = vis.feature.level,
                                             feature.dat.type = feature.dat.type,
                                             features.plot = NULL,
                                             top.k.plot = dotplot.feature.no,
                                             top.k.func = 'mean',
                                             prev.filter = prev.filter,
                                             abund.filter = abund.filter,
                                             base.size = base.size,
                                             theme.choice = theme.choice,
                                             custom.theme = custom.theme,
                                             palette = palette,
                                             pdf = pdf,
                                             file.ann = file.ann,
                                             pdf.wid = pdf.wid,
                                             pdf.hei = pdf.hei)
```

```{r taxa-change-dotplot-pair-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 20}
if (is.function(feature.change.func)) {
  cat('The changes were computed using a custom function provided by the user.')
} else if (feature.change.func == 'relative change') {
  cat('The changes were relative changes, which were computed as (after.abund - before.abund) / (after.abund + before.abund) so the values lie between [-1, 1].')
} else if (feature.change.func == 'absolute change') {
  cat('The changes were absolute changes, computed as the difference between after.abund and before.abund.')
} else if (feature.change.func == 'log fold change') {
  cat('The changes were log2 fold changes, computed as the logarithm of the ratio of after.abund to before.abund, with a small constant added to avoid taking the log of zero.')
}

cat(' The following plots display the average change for each group, and stratum. \n\n')
taxa_change_dotplot_results
```

# 2. Alpha diversity analysis

## 2.1 Data visualization

### 2.1.1 Alpha diversity boxplot

```{r alpha-boxplot-long-generation, message=FALSE, fig.align='center', results='asis', fig.width=8, fig.height=3}
alpha_boxplot_results <- generate_alpha_boxplot_long(
                                                data.obj = data.obj,
                                                alpha.obj = alpha.obj,
                                                alpha.name = alpha.name,
                                                depth = depth,
                                                subject.var = subject.var,
                                                time.var = time.var,
                                                t0.level = change.base,
                                                ts.levels = NULL,
                                                group.var = group.var,
                                                strata.var = strata.var,
                                                base.size = base.size,
                                                theme.choice = theme.choice,
                                                custom.theme = custom.theme,
                                                palette = palette,
                                                pdf = pdf,
                                                file.ann = file.ann,
                                                pdf.wid = pdf.wid,
                                                pdf.hei = pdf.hei)
alpha_boxplot_results
```

### 2.1.2 Alpha diversity change boxplot

```{r alpha-diversity-change-boxplot, message=FALSE, fig.align='center', results='asis'}
alpha_change_boxplot_results <- generate_alpha_change_boxplot_pair(
                                               data.obj = data.obj,
                                               alpha.obj = alpha.obj,
                                               alpha.name = alpha.name,
                                               depth = depth,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               change.base = change.base,
                                               alpha.change.func = alpha.change.func,
                                               group.var = group.var,
                                               strata.var = strata.var,
                                               base.size = base.size,
                                               theme.choice = theme.choice,
                                               custom.theme = custom.theme,
                                               palette = palette,
                                               pdf = pdf,
                                               file.ann = file.ann,
                                               pdf.wid = pdf.wid,
                                               pdf.hei = pdf.hei)
```

```{r alpha-diversity-change-boxplot-print, message=FALSE, fig.align='center', results='asis', echo = FALSE, fig.width=8, fig.height=3}
if (is.function(alpha.change.func)) {
  cat('The changes from change.base were computed using a custom function provided by the user.\n\n')
} else if (alpha.change.func == 'log fold change') {
  cat('The changes from change.base were computed as the log2 fold change of alpha diversity at the current timepoint versus change.base.\n\n')
} else {
  cat('The changes from change.base were computed as the absolute difference in alpha diversity at the current timepoint versus change.base.\n\n')
}

alpha_change_boxplot_results
```

## 2.2 Alpha diversity association test based on LMM

```{r alpha-test-pair-generation, message=FALSE}
alpha_test_results <- generate_alpha_test_pair(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 depth = depth,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = test.adj.vars)
```

```{r alpha-test-results-analysis, echo=FALSE, message=FALSE, results='asis'}

group_levels <- data.obj$meta.dat %>% dplyr::select(!!sym(group.var)) %>% dplyr::pull() %>% as.factor() %>% levels

reference_level <- group_levels[1]

  # Extract p-value for a given term
  extract_p_val <- function(term) {
    return(data_frame[data_frame$Term == term,]$P.Value)
  }

if (!is.null(test.adj.vars)) {
    adj.vars_string <- paste(test.adj.vars, collapse = ', ')

    cat(sprintf('In this analysis, we utilized a linear mixed effects model with a random intercept and possibly a random slope for time to investigate a potential difference in alpha diversity across different levels of %s. The model includes an interaction term between %s and %s. Additionally, we included %s, %s and %s as covariates.', group.var, time.var, group.var, time.var, group.var, adj.vars_string))
} else {
    cat(sprintf('In this analysis, we utilized a linear mixed effects model with a random intercept and possibly a random slope for time to investigate a potential difference in alpha diversity across different levels of %s. The model includes an interaction term between %s and %s. Specifically, we included %s and %s as covariates.', group.var, time.var, group.var, time.var, group.var))
}

# Define a function to report the significance of interaction terms
report_significance <- function(data_frame, group.var) {

    if (length(group_levels) > 2) {
        p_val_group_var <- data_frame[data_frame$Term == paste0(group.var),]$P.Value
        cat(sprintf(
            '\n An ANOVA test of the null hypothesis of no group difference among the %s levels produces a p-value of %.3f.\n\n',
            length(group_levels),
            p_val_group_var
        ))

        p_val_group_var_time_var <- data_frame[data_frame$Term == paste0(group.var, ':', time.var),]$P.Value
        cat(sprintf(
            '\n An ANOVA test of the null hypothesis of no trend difference among the %s groups produces a p-value of %.3f.\n\n',
            length(group_levels),
            p_val_group_var
        ))

    } else {
        # Extracting interaction terms
    interaction_terms <- data_frame$Term[grepl(paste0('^', group.var, '[^:]*$'), data_frame$Term)]
        for (term in interaction_terms) {
            p_val <- data_frame[data_frame$Term == term,]$P.Value
            level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
            level <- gsub('_', '', level) # Remove any underscores if they exist

            # Describing interaction terms
            if (p_val < 0.05) {
                cat(sprintf(
                    '\n Based on the linear mixed effects model, a significant group difference was observed between %s and %s of the variable %s, with a p-value of %.3f.',
                    reference_level,
                    level,
                    group.var,
                    p_val
                ))
            } else {
                cat(sprintf(
                    '\n Based on the linear mixed effects model, no significant group difference was detected between %s and %s of the variable %s, with a p-value of %.3f.',
                    reference_level,
                    level,
                    group.var,
                    p_val
                ))
            }
        }
          # Extracting interaction terms
  interaction_terms <- grep(paste0(group.var, '.+:', time.var), data_frame$Term, value = TRUE)

    # Report significant trend for interaction term
  report_interaction_significance <- function(term, group.var) {
    p_val <- data_frame[data_frame$Term == term,]$P.Value
    level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
    level <- gsub('_', '', level) # Remove underscores

    if(p_val < 0.05) {
      cat(sprintf('\n Based on the linear mixed effects model, a significant trend difference was observed between %s and %s of the variable %s, with a p-value of %.3f.', reference_level, level, group.var, p_val))
    } else {
      cat(sprintf('\n Based on the linear mixed effects model, no significant trend difference was detected between %s and %s of the variable %s, with a p-value of %.3f.', reference_level, level, group.var, p_val))
    }
  }

    for(term in interaction_terms) {
    report_interaction_significance(term, group.var)
    }

    }
}


# Function to convert first letter to uppercase
firstToUpper <- function(s) {
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

# Initialize the sub-section counter
counter <- 1

# Report significance for each diversity index
for(index_name in names(alpha_test_results)) {
  # Print with updated counter and index name
  cat(sprintf('\n\n### 2.2.%d %s index \n\n', counter, firstToUpper(ifelse(index_name == 'observed_species', 'observed species', index_name))))
  cat('\n')

  # Report significance
  report_significance(data_frame = alpha_test_results[[index_name]], group.var = group.var)
  cat('\n\n')

  output <- pander::pander(alpha_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```

## 2.3 Alpha diversity association test based on changes

```{r alpha-diversity-change-test, message=FALSE}
alpha_change_test_results <- generate_alpha_change_test_pair(
                                                 data.obj = data.obj,
                                                 time.var = time.var,
                                                 change.base = change.base,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = test.adj.vars,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 depth = depth,
                                                 alpha.change.func = alpha.change.func)
```

```{r alpha-diversity-change-analysis, echo=FALSE, message=FALSE, results='asis'}

# Initial description for volatility
num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('\n In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on change.\n', group.var))
} else {
    cat(sprintf('\n In this analysis, we utilized a general linear model to examine the influence of the variable %s on change.\n', group.var))
}

if (is.function(alpha.change.func)) {
    cat('The alpha diversity change is calculated using a custom function supplied by the user to compute the rate of change between consecutive time points.\n\n')
} else {
    if (alpha.change.func == 'log fold change') {
        cat('The alpha diversity change is calculated by taking the logarithm of the fold change between consecutive time points.\n\n')
    } else {
        cat('The alpha diversity change is calculated by computing the direct difference in alpha diversity between consecutive time points.\n\n')
    }
}

# Define a function to report the significance of volatility based on group.var
report_change_significance <- function(data_frame, group.var) {

  # Extracting terms excluding ANOVA and intercept
  terms <- grep(group.var, data_frame$Term, value = TRUE)
  terms <- terms[!terms %in% c('(Intercept)', 'Residuals', group.var)]

  for(term in terms) {
      p_val <- data_frame[data_frame$Term == term,]$P.Value

      # Extract only the level part from the term by removing the group.var prefix and underscore
      level <- sub(group.var, '', term)

      # Describing significance based on lm model
      if(p_val < 0.05) {
          cat(sprintf('\n Based on the general linear model, the level %s of the variable %s significantly differs from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
      } else {
          cat(sprintf('\n Based on the general linear model, the level %s of the variable %s did not significantly differ from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
      }
  }

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(num_levels > 2) {
      cat(sprintf('An ANOVA test of the null hypothesis of no difference in change among %s groups produces a p-value of %.3f. ', num_levels, p_val_anova))
  }
}

# Initialize the sub-section counter
counter <- 1

# Report significance for each diversity index
for(index_name in names(alpha_change_test_results)) {
  # Print with updated counter and index name
  cat(sprintf('\n### 2.3.%d %s index \n\n', counter, firstToUpper(ifelse(index_name == 'observed_species', 'observed species', index_name))))
  cat('\n')

  # Report significance
  report_change_significance(data_frame = alpha_change_test_results[[index_name]], group.var = group.var)
  cat('\n')

  output <- pander::pander(alpha_change_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```

# 3. Beta diversity analysis

## 3.1 Data visualization

### 3.1.1 Beta diversity ordinationplot

```{r beta-ordination-pair-generation, message=FALSE, fig.align='center', warning = FALSE, fig.width = 18, fig.height = 8, results='asis'}
beta_ordination_results <- generate_beta_ordination_pair(
                                                    data.obj = data.obj,
                                                    dist.obj = dist.obj,
                                                    pc.obj = pc.obj,
                                                    subject.var = subject.var,
                                                    time.var = time.var,
                                                    group.var = group.var,
                                                    strata.var = NULL,
                                                    dist.name = dist.name,
                                                    base.size = base.size,
                                                    theme.choice = theme.choice,
                                                    custom.theme = custom.theme,
                                                    palette = palette,
                                                    pdf = pdf,
                                                    file.ann = file.ann,
                                                    pdf.wid = pdf.wid,
                                                    pdf.hei = pdf.hei)
beta_ordination_results

if (!is.null(strata.var)){
beta_ordination_stratified_results <- generate_beta_ordination_pair(
                                                    data.obj = data.obj,
                                                    dist.obj = dist.obj,
                                                    pc.obj = pc.obj,
                                                    subject.var = subject.var,
                                                    time.var = time.var,
                                                    group.var = group.var,
                                                    strata.var = strata.var,
                                                    dist.name = dist.name,
                                                    base.size = base.size,
                                                    theme.choice = theme.choice,
                                                    custom.theme = custom.theme,
                                                    palette = palette,
                                                    pdf = pdf,
                                                    file.ann = file.ann,
                                                    pdf.wid = pdf.wid,
                                                    pdf.hei = pdf.hei)
beta_ordination_stratified_results
}
```

### 3.1.2 Beta diversity change boxplot

```{r beta-diversity-change-boxplot, message=FALSE, fig.align='center', results='hide', results='asis', fig.width=8, fig.height=3}
beta_change_boxplot_results <- generate_beta_change_boxplot_pair(
                                                      data.obj = data.obj,
                                                      dist.obj = dist.obj,
                                                      subject.var = subject.var,
                                                      time.var = time.var,
                                                      group.var = group.var,
                                                      strata.var = strata.var,
                                                      change.base = change.base,
                                                      dist.name = dist.name,
                                                      base.size = base.size,
                                                      theme.choice = theme.choice,
                                                      custom.theme = custom.theme,
                                                      palette = palette,
                                                      pdf = pdf,
                                                      file.ann = file.ann,
                                                      pdf.wid = pdf.wid,
                                                      pdf.hei = pdf.hei)
```

```{r beta-diversity-change-boxplot-print, message=FALSE, fig.align='center', results='hide', results='asis', fig.width=8, fig.height=3, echo = FALSE}
cat(sprintf('\n Beta change represents the distance of each subject from their change.base.\n\n'))
beta_change_boxplot_results
```

### 3.1.3 Beta diversity PC boxplot

```{r beta-pc-boxplot-pair-generation, message=FALSE, fig.align='center', results='asis', fig.width=8, fig.height=3}
pc_boxplot_longitudinal_results <- generate_beta_pc_boxplot_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  pc.obj = pc.obj,
  pc.ind = c(1, 2),
  subject.var = subject.var,
  time.var = time.var,
  t0.level = change.base,
  ts.levels = NULL,
  group.var = group.var,
  strata.var = strata.var,
  dist.name = dist.name,
  base.size = base.size,
  theme.choice = theme.choice,
  custom.theme = custom.theme,
  palette = palette,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)

pc_boxplot_longitudinal_results
```

### 3.1.4 Beta diversity PC change boxplot

```{r pc-change-boxplot-pairs, message=FALSE, fig.align='center', results='hide', results='asis', fig.width=8, fig.height=3}
pc_change_boxplot_pairs <- generate_beta_pc_change_boxplot_pair(
  data.obj = data.obj,
  dist.obj = dist.obj,
  pc.obj = pc.obj,
  pc.ind = c(1, 2),
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  change.base = change.base,
  change.func = 'absolute change',
  dist.name = dist.name,
  base.size = base.size,
  theme.choice = theme.choice,
  custom.theme = custom.theme,
  palette = palette,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)

pc_change_boxplot_pairs
```

## 3.2 Beta diversity association test based on changes

```{r beta-diversity-change-test, message=FALSE, results='asis'}
beta_change_test_results <- generate_beta_change_test_pair(
                                               data.obj = data.obj,
                                               dist.obj = dist.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars,
                                               dist.name = dist.name,
                                               change.base = change.base)
```

```{r beta-diversity-change-analysis, echo=FALSE, message=FALSE, results='asis'}

# Initial description for volatility
num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on beta diversity change.', group.var))
} else {
    cat(sprintf('In this analysis, we utilized a general linear model to examine the influence of the variable %s on beta diversity change.', group.var))
}

cat(sprintf('\n Beta change represents the distance of each subject from their change.base.\n\n'))

# Define a function to report the significance of volatility based on group.var
report_beta_change_significance <- function(data_frame, group.var) {

  # Extracting terms excluding ANOVA and intercept
  terms <- grep(group.var, data_frame$Term, value = TRUE)
  terms <- terms[!terms %in% c('(Intercept)', 'Residuals', group.var)]

  for(term in terms) {
    p_val <- data_frame[data_frame$Term == term,]$P.Value

    # Extract only the level part from the term by removing the group.var prefix and underscore
    level <- sub(group.var, '', term)

    # Describing significance based on lm model
    if(p_val < 0.05) {
      cat(sprintf('\n Based on the general linear model, the level %s of the variable %s significantly differs from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
    } else {
      cat(sprintf('\n Based on the general linear model, the level %s of the variable %s did not significantly differ from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
    }
  }

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(num_levels > 2) {
    cat(sprintf('An ANOVA test of the null hypothesis of no difference in volatility among %s groups produces a p-value of %.3f. ', num_levels, p_val_anova))
  }

}

counter <- 1

for(index_name in names(beta_change_test_results)) {
  cat(sprintf('\n### 3.3.%d %s distance \n\n', counter, ifelse(index_name == 'BC', 'Bray-Curtis', index_name)))

  report_beta_change_significance(beta_change_test_results[[index_name]], group.var)

  cat('\n')
  output <- pander::pander(beta_change_test_results[[index_name]])
  cat(output)

  counter <- counter + 1
}

```

# 4. Feature-level Analysis

```{r Check-and-Select-Rarefaction2, echo=FALSE, message=FALSE, results='asis'}

if (feature.analysis.rarafy) {
    cat('Rarefaction has been enabled for feature-level analysis and visualization.\n\n',
        'Reason: The observed abundance of rare/low-abundance features can be strongly influenced by the sequencing depth. ',
        'Rarefaction is an effective method to control the effect of sequencing depth variation. ',
        'By employing rarefaction, presence/absence status of the features are more comparable and we can potentially increase the power of detecting those rare/low-abundance features, ',
        'even though it introduces some variation due to under-sampling.\n',
        'In essence, this step improves comparability across samples across samples with varying sequencing depth.\n\n',
        'If you do not wish to perform rarefaction during feature-level analysis, please turn feature.analysis.rarafy to FALSE.\n')

    data.obj <- rarefy.data.obj
} else {
  data.obj <- original.data.obj
}

```

## 4.1 Feature-level association test based on LinDA-LMM

```{r taxa-test-generation, message=FALSE, results='asis', fig.align='center'}
taxa_test_results <- generate_taxa_test_pair(data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = test.feature.level,
                                               feature.dat.type = feature.dat.type)
```

```{r taxa-test-results-display, echo=FALSE, message=FALSE, results='asis', warning = FALSE, fig.align='center', fig.width = 10, fig.height = 8}
volcano_plots <- generate_taxa_volcano_single(
                                  data.obj = data.obj,
                                  group.var = group.var,
                                  test.list = taxa_test_results,
                                  feature.sig.level = feature.sig.level,
                                  feature.mt.method = feature.mt.method
)
volcano_plots

cat(sprintf('In this analysis, we utilized the LinDA linear model to investigate potential differences in trend. Specifically, we tested the effect of the variable %s and the interaction between %s and %s for different taxa, while adjusting for other covariates.\n\n',
            group.var, group.var, time.var))

# Iterate over each taxonomic rank in taxa_test_results
for(taxon_rank in names(taxa_test_results)) {

  # Extract the specific taxonomic rank results
  taxon_results <- taxa_test_results[[taxon_rank]]

  # Iterate over each comparison in taxon_results
  for(comparison in names(taxon_results)) {

    # Extract specific comparison results data frame
    comparison_df <- taxon_results[[comparison]]

    # Filter interaction terms based on the selected multiple testing method
    if (feature.mt.method == 'fdr') {
        interaction_terms_results <- comparison_df %>%
            dplyr::filter(Adjusted.P.Value < feature.sig.level)
        p_value_str = 'adjusted p-value'
    } else if (feature.mt.method == 'none') {
        interaction_terms_results <- comparison_df %>%
            dplyr::filter(P.Value < feature.sig.level)
        p_value_str = 'p-value'
    }

    # Check if filtered results have rows
    if (nrow(interaction_terms_results) == 0) {
        cat(sprintf('For the taxa investigated under the %s category in comparison %s, no significant results were detected using the %s method for p-value adjustment, at a %s threshold of %s. ', taxon_rank, comparison, feature.mt.method, p_value_str, feature.sig.level))
    } else {
        cat(sprintf('For the taxon %s in comparison %s, significant results were identified using the %s method for p-value adjustment, based on a threshold of %s. ', taxon_rank, comparison, feature.mt.method, feature.sig.level))
        cat('\n')
        output <- pander::pander(interaction_terms_results)
        cat(output)
    }
  }
}

filename_prefix <- 'taxa_test_results_'
file_ext <- '.csv'

for(taxon_rank in names(taxa_test_results)) {

    comparisons <- names(taxa_test_results[[taxon_rank]])

    for(comparison in comparisons) {

        file_name <- paste0(filename_prefix, taxon_rank, '_', gsub(' ', '_', gsub('/', '_or_', comparison)), file_ext)


        write.csv(taxa_test_results[[taxon_rank]][[comparison]],
                  file = file_name,
                  row.names = FALSE)
    }
}

cat(sprintf('\n\nThe results for features have been saved in the current working directory.
Each taxa rank has its own file named in the format: %s followed by the taxon rank, the comparison, and the file extension %s.
Please refer to these files for more detailed results.', filename_prefix, file_ext))

```

## 4.2 Feature-level association test based on changes

```{r taxa-change-test-pair, message=FALSE, results='hide'}
taxa_change_test_results <- generate_taxa_change_test_pair(data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               change.base = change.base,
                                               feature.change.func = feature.change.func,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = test.feature.level,
                                               feature.dat.type = feature.dat.type)
```

```{r taxa-change-test-results-display, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 10, fig.height = 8}

change_volcano_plots <- generate_taxa_volcano_single(data.obj = data.obj,
                                                     group.var = group.var,
                                                     test.list = taxa_change_test_results,
                                                     feature.sig.level = feature.sig.level,
                                                     feature.mt.method = feature.mt.method)

change_volcano_plots

# Initial description for Feature Volatility
num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, a general linear model followed by ANOVA was employed to test the effect of %s on the change of various taxa abundances.', group.var))
} else {
    cat(sprintf('In this analysis, a general linear model was utilized to investigate the influence of the variable %s on the change of various taxa abundances.', group.var))
}

if (is.function(feature.change.func)) {
  cat('The changes were computed using a custom function provided by the user.')
} else if (feature.change.func == 'relative change') {
  cat('The changes were relative changes, which were computed as (after.abund - before.abund) / (after.abund + before.abund) so the values lie between [-1, 1].')
} else if (feature.change.func == 'absolute change') {
  cat('The changes were absolute changes, computed as the difference between after.abund and before.abund.')
} else if (feature.change.func == 'log fold change') {
  cat('The changes were log2 fold changes, computed as the logarithm of the ratio of after.abund to before.abund, with a small constant added to avoid taking the log of zero.')
}

# Iterate over each taxonomic rank in taxa_change_test_results
for(taxon_rank in names(taxa_change_test_results)) {

  # Extract the specific taxonomic rank results
  taxon_results <- taxa_change_test_results[[taxon_rank]]

  # Iterate over each comparison in taxon_results
  for(comparison in names(taxon_results)) {

    # Extract specific comparison results data frame
    comparison_df <- taxon_results[[comparison]]

    # Filter interaction terms based on the selected multiple testing method
    if (feature.mt.method == 'fdr') {
        interaction_terms_results <- comparison_df %>%
            dplyr::filter(Adjusted.P.Value < feature.sig.level)
        p_value_str = 'adjusted p-value'
    } else if (feature.mt.method == 'none') {
        interaction_terms_results <- comparison_df %>%
            dplyr::filter(P.Value < feature.sig.level)
        p_value_str = 'p-value'
    }

    # Check if filtered results have rows
    if (nrow(interaction_terms_results) == 0) {
        cat(sprintf('For the taxa investigated under the %s category in comparison %s, no significant results were detected using the %s method for p-value adjustment, at a %s threshold of %s. ', taxon_rank, comparison, feature.mt.method, p_value_str, feature.sig.level))
    } else {
        cat(sprintf('For the taxon %s in comparison %s, significant results were identified using the %s method for p-value adjustment, based on a threshold of %s. ', taxon_rank, comparison, feature.mt.method, feature.sig.level))
        cat('\n')
        output <- pander::pander(interaction_terms_results)
        cat(output)
    }
  }
}

filename_prefix <- 'taxa_change_test_results_'
file_ext <- '.csv'

for(taxon_rank in names(taxa_change_test_results)) {

    comparisons <- names(taxa_change_test_results[[taxon_rank]])

    for(comparison in comparisons) {

        file_name <- paste0(filename_prefix, taxon_rank, '_', gsub(' ', '_', gsub('/', '_or_', comparison)), file_ext)

        write.csv(taxa_change_test_results[[taxon_rank]][[comparison]],
                  file = file_name,
                  row.names = FALSE)
    }
}

cat(sprintf('\n\n The change test results for individual feature have been saved in the current working directory. Each taxa rank and its corresponding comparison have their own file named with the prefix: %s followed by the taxon rank, the comparison, and the file extension %s. Please refer to these files for more detailed data.', filename_prefix, file_ext))

```

```{r extract_significant_taxa, echo=FALSE, results='hide'}
# Extract Variables based on the provided condition
extract_significant_variables <- function(data_frame, p_value_column) {
    filtered_df <- data_frame %>%
        dplyr::filter(!!sym(p_value_column) < feature.sig.level)
    return(filtered_df$Variable)
}

# Decide which column to filter on based on feature.mt.method
if (feature.mt.method == 'fdr') {
    p_value_column <- 'Adjusted.P.Value'
} else if (feature.mt.method == 'none') {
    p_value_column <- 'P.Value'
} else {
    stop('Invalid feature.mt.method provided!')
}

# Function to process each list in the main list
process_list <- function(test_results) {
    significant_vars <- lapply(test_results, function(taxon_results) {
        lapply(taxon_results, function(comparison_df) {
            extract_significant_variables(comparison_df, p_value_column)
        })
    })
    return(unlist(significant_vars))
}

# Process each main list and extract significant variables
significant_vars <- process_list(taxa_test_results)
significant_vars_change <- process_list(taxa_change_test_results)

# Combine and de-duplicate the results
combined_significant_taxa <- unique(c(significant_vars, significant_vars_change))

```

## 4.3 Data visualization (significant features)

### 4.3.1 Significant features boxplot

```{r taxa-test-boxplot-pair-generation, message = FALSE, warning = FALSE, fig.width = 8, fig.height = 3, fig.align='center', results='asis'}

if (length(significant_vars) != 0){
taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_long(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = change.base,
                                   ts.levels = NULL,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = test.feature.level,
                                   features.plot = significant_vars,
                                   transform = feature.box.axis.transform,
                                   feature.dat.type = feature.dat.type,
                                   top.k.plot = NULL,
                                   top.k.func = NULL,
                                   prev.filter = prev.filter,
                                   abund.filter = abund.filter,
                                   base.size = base.size,
                                   theme.choice = theme.choice,
                                   custom.theme = custom.theme,
                                   palette = palette,
                                   pdf = TRUE,
                                   file.ann = file.ann,
                                   pdf.wid = pdf.wid,
                                   pdf.hei = pdf.hei)

taxa_indiv_boxplot_results
}

```

```{r boxplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

if (length(significant_vars) != 0){
  pdf_name <- paste0(
          'taxa_indiv_boxplot_long',
          '_',
          'subject_',
          subject.var,
          '_',
          'time_',
          time.var,
          '_',
          'feature_level_',
          test.feature.level,
          '_',
          'transform_',
          feature.box.axis.transform,
          '_',
          'prev_filter_',
          prev.filter,
          '_',
          'abund_filter_',
          abund.filter
        )
        if (!is.null(group.var)) {
          pdf_name <- paste0(pdf_name, '_', 'group_', group.var)
        }
        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, '_', 'strata_', strata.var)
        }
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, '_', file.ann)
        }

cat('\n')
cat(paste0('\n\n The boxplot results for individual features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
}
```

### 4.3.2 Significant features boxplot (change)

```{r taxa-change-boxplot-generation, message=FALSE, fig.align='center', fig.width = 8, fig.height = 3, results='asis'}

if (length(significant_vars_change) != 0){
taxa_indiv_change_boxplot_results <- generate_taxa_indiv_change_boxplot_pair(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   change.base = change.base,
                                   feature.change.func = feature.change.func,
                                   feature.level = test.feature.level,
                                   features.plot = significant_vars_change,
                                   feature.dat.type = feature.dat.type,
                                   top.k.plot = NULL,
                                   top.k.func = NULL,
                                   prev.filter = prev.filter,
                                   abund.filter = abund.filter,
                                   base.size = base.size,
                                   theme.choice = theme.choice,
                                   custom.theme = custom.theme,
                                   palette = palette,
                                   pdf = pdf,
                                   file.ann = file.ann,
                                   pdf.wid = pdf.wid,
                                   pdf.hei = pdf.hei)
taxa_indiv_change_boxplot_results
}

```

```{r change-boxplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

if (length(significant_vars_change) != 0){
  pdf_name <- paste0(
          'taxa_indiv_change_boxplot_pair',
          '_',
          'subject_',
          subject.var,
          '_',
          'time_',
          time.var,
          '_',
          'change_base_',
          change.base,
          '_',
          'feature_level_',
          test.feature.level,
          '_',
          'prev_filter_',
          prev.filter,
          '_',
          'abund_filter_',
          abund.filter
        )
        if (!is.null(group.var)) {
          pdf_name <- paste0(pdf_name, '_', 'group_', group.var)
        }
        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, '_', 'strata_', strata.var)
        }
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, '_', file.ann)
        }
        pdf_name <- paste0(pdf_name, '.pdf')

cat('\n')
cat(paste0('\n\n The change boxplot results for individual features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
}

```

"

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
