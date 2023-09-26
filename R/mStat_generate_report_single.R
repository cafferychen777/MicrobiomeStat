#' Generate a report for microbial ecology analysis of cross-sectional, longitudinal or paired data at a single time point
#'
#' This function generates a comprehensive report for microbial ecology analysis,
#' including alpha diversity, beta diversity, and taxonomic feature analyses. The function is designed
#' to perform analysis on cross-sectional data, a single time point from longitudinal or paired data.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var Variable name used for grouping samples.
#' @param test.adj.vars Character vector, names of columns in the metadata containing covariates to be adjusted for in statistical tests and models, such as linear mixed effects models for longitudinal data analysis. This allows the user to account for the effects of additional variables in assessing the effects of primary variables of interest such as time and groups. Default is NULL, which indicates no covariates are adjusted for in statistical testing.
#' @param vis.adj.vars Character vector, names of columns in the metadata containing covariates to visualize in plots, in addition to the primary variables of interest such as groups. For example, if sex is provided in vis.adj.vars, plots will display facets or colors for different sex groups. This allows visualization of effects across multiple covariates. Default is NULL, which indicates only the primary variables of interest will be visualized without additional covariates.
#' @param strata.var Variable to stratify the analysis by (optional).
#' @param subject.var Variable name used for subject identification.
#' @param time.var Variable name used for time points.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
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
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
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
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param feature.analysis.rarafy Logical, indicating whether to rarefy the data at the feature-level for analysis.
#' If TRUE, the feature data will be rarefied before analysis. Default is TRUE.
#' @param bar.area.feature.no A numeric value indicating the number of top abundant features to retain in both barplot and areaplot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 20.
#' @param heatmap.feature.no A numeric value indicating the number of top abundant features to retain in the heatmap. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 20.
#' @param dotplot.feature.no A numeric value indicating the number of top abundant features to retain in the dotplot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 40.
#' @param feature.mt.method Character, multiple testing method for features, "fdr" or "none", default is "fdr".
#' @param feature.sig.level Numeric, significance level cutoff for highlighting features, default is 0.1.
#' @param feature.box.axis.transform A string indicating the transformation to apply to the data before plotting. Options are:
#' - "identity": No transformation (default)
#' - "sqrt": Square root transformation
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
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A report file containing the microbial ecology analysis results.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   vis.adj.vars = c("sex"),
#'   test.adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "sex",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Microbiome/Longitudinal/MicrobiomeStat_Paper/Report/devtools_test_report_single2.pdf"
#' )
#' data(peerj32.obj)
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   vis.adj.vars = c("sex"),
#'   test.adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "sex",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Microbiome/Longitudinal/MicrobiomeStat_Paper/Report/devtools_test_report_single1.pdf"
#' )
#' }
#' data(peerj32.obj)
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   vis.adj.vars = c("sex"),
#'   test.adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "sex",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Microbiome/Longitudinal/MicrobiomeStat_Paper/Report/devtools_test_report_single2.pdf"
#' )
#' data(peerj32.obj)
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   vis.adj.vars = c("sex"),
#'   test.adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon", "observed_species"),
#'   depth = NULL,
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   feature.box.axis.transform = "sqrt",
#'   strata.var = "sex",
#'   vis.feature.level = c("Phylum", "Family", "Genus"),
#'   test.feature.level = "Family",
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 20,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.2,
#'   output.file = "/Users/apple/Microbiome/Longitudinal/MicrobiomeStat_Paper/Report/devtools_test_report_single1.pdf"
#' )
#' @export
mStat_generate_report_single <- function(data.obj,
                                         group.var,
                                         vis.adj.vars = NULL,
                                         test.adj.vars = NULL,
                                         strata.var = NULL,
                                         subject.var,
                                         time.var,
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
                                         feature.box.axis.transform = c("identity", "sqrt", "log"),
                                         base.size = 16,
                                         theme.choice = "bw",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         output.file,
                                         ...) {

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

custom_t.level_status <- ifelse(is.null(t.level), 'NULL', t.level)

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
                                        't.level',
                                        'alpha.obj',
                                        'alpha.name',
                                        'depth',
                                        'dist.obj',
                                        'dist.name',
                                        'pc.obj',
                                        'prev.filter',
                                        'abund.filter',
                                        'feature.analysis.rarafy',
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
                                        custom_t.level_status,
                                        deparse(substitute(alpha.obj)),
                                        toString(alpha.name),
                                        custom_depth_status,
                                        deparse(substitute(dist.obj)),
                                        toString(dist.name),
                                        deparse(substitute(pc.obj)),
                                        prev.filter,
                                        abund.filter,
                                        feature.analysis.rarafy,
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
pander::pander(mStat_results)
```

```{r object-pre-calculation, echo=FALSE, message=FALSE, results='asis'}

original.data.obj <- data.obj

rarefy.data.obj <- mStat_normalize_data(data.obj = data.obj, method = 'Rarefy', depth = depth)$data.obj.norm

if (is.null(depth)){
  depth <- min(colSums(data.obj$feature.tab))
  cat(sprintf('No rarefaction depth is specified. The minimum depth, %d, is used as the rarefaction depth. ', depth))
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
  cat('pc.obj is calculated based on the dist.obj.\n')
}

```

## 1.3 Data visualization(overall)

```{r Check-and-Select-Rarefaction1, echo=FALSE, message=FALSE, results='asis'}

if (feature.analysis.rarafy) {
    cat('Rarefaction has been enabled for feature-level analysis.\n\n',
        'Reason: The observed abundance of rare/low-abundance features can be strongly influenced by sequence depth. ',
        'Rarefaction is an effective method to control the effect of sequence depth variation. ',
        'By employing rarefaction, we can potentially increase the power of detecting those rare/low-abundance features, ',
        'even though it introduces some variation due to under-sampling.\n',
        'In essence, this step helps to ensure more accurate and consistent results across samples with varying sequence depths.\n\n',
        'If you do not wish to perform rarefaction during feature-level analysis, please turn feature.analysis.rarafy to FALSE.\n')

    data.obj <- rarefy.data.obj
} else {
  data.obj <- original.data.obj
}

```

### 1.3.1 Feature barplot

```{r taxa-barplot-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8}
taxa_barplot_results <- generate_taxa_barplot_single(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t.level = t.level,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = vis.feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     feature.number = bar.area.feature.no,
                                                     base.size = base.size,
                                                     theme.choice = theme.choice,
                                                     custom.theme = custom.theme,
                                                     palette = NULL,
                                                     pdf = pdf,
                                                     file.ann = file.ann,
                                                     pdf.wid = pdf.wid,
                                                     pdf.hei = pdf.hei)
```

```{r taxa-barplot-avergae-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 15}
cat('The following plots display the average proportions for each group, and stratum. \n\n')
indiv_list <- lapply(taxa_barplot_results, function(x) x$indiv)

average_list <- lapply(taxa_barplot_results, function(x) x$average)

average_list
```

```{r taxa-barplot-indiv-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 15}
cat('The following plots display the individual proportions for each group, and stratum. \n\n')
indiv_list
```

### 1.3.2 Feature dotplot

```{r taxa-dotplot-generation, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 8}
taxa_dotplot_results <- generate_taxa_dotplot_single(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t0.level = t.level,
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
taxa_dotplot_results
```

### 1.3.3 Feature heatmap

```{r taxa-heatmap-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_heatmap_results <- generate_taxa_heatmap_single(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t.level = t.level,
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
                                                     palette = NULL,
                                                     cluster.cols = NULL,
                                                     cluster.rows = NULL,
                                                     pdf = pdf,
                                                     file.ann = file.ann,
                                                     pdf.wid = pdf.wid,
                                                     pdf.hei = pdf.hei)
```

```{r taxa-heatmap-indiv-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 12}
cat('The following plots display the individual proportions for each sample. \n\n')
taxa_heatmap_results
```

# 2. Alpha diversity analysis

## 2.1 Data visualization

### 2.1.1 Alpha diversity boxplot

```{r alpha-boxplot-generation, message=FALSE, warning = FALSE, fig.align='center', fig.width = 16, fig.height = 8, results='asis'}
alpha_boxplot_results <- generate_alpha_boxplot_single(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       depth = depth,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t.level = t.level,
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

## 2.2 Anova test

```{r alpha-diversity-test-generation, message = FALSE, warning = FALSE}
alpha_test_results <- generate_alpha_test_single(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 depth = depth,
                                                 time.var = time.var,
                                                 t.level = t.level,
                                                 group.var = group.var,
                                                 adj.vars = test.adj.vars)
```

```{r alpha-index-analysis, echo=FALSE, message=FALSE, results='asis'}

# Function to convert first letter to uppercase
firstToUpper <- function(s) {
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

# Initial description for volatility
num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

group_levels <- data.obj$meta.dat %>% select(!!sym(group.var)) %>% pull() %>% as.factor() %>% levels

reference_level <- group_levels[1]

if (!is.null(test.adj.vars)) {
    adj.description <- sprintf(' while adjusting for covariates %s', paste(test.adj.vars, collapse=' and '))
} else {
    adj.description <- ''
}

if(num_levels > 2) {
    cat(sprintf('\n In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on alpha diversity%s.\n', group.var, adj.description))
} else {
    cat(sprintf('\n In this analysis, we utilized a general linear model to examine the influence of the variable %s on alpha diversity%s.\n', group.var, adj.description))
}

# Define a function to report the significance of volatility based on group.var
report_significance <- function(data_frame, group.var) {

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
      cat(sprintf('An ANOVA test of the null hypothesis of no difference in alpha diversity among %s groups produces a p-value of %.3f. ', num_levels, p_val_anova))
  }
}

# Initialize the sub-section counter
counter <- 1

# Report significance for each diversity index
for(index_name in names(alpha_test_results)) {
  # Print with updated counter and index name
  cat(sprintf('\n### 2.2.%d %s index \n\n', counter, firstToUpper(index_name)))
  cat('\n')

  # Report significance
  report_significance(data_frame = alpha_test_results[[index_name]], group.var = group.var)
  cat('\n')

  output <- pander::pander(alpha_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```

# 3. Beta diversity analysis

## 3.1 Data visualization

### 3.1.1 Beta diversity ordinationplot

```{r beta-ordination-generation, message=FALSE, fig.align='center', warning = FALSE, fig.width = 10, fig.height = 8, results='asis'}
beta_ordination_results <- generate_beta_ordination_single(data.obj = data.obj,
                                                           dist.obj = dist.obj,
                                                           pc.obj = pc.obj,
                                                           subject.var = subject.var,
                                                           time.var = time.var,
                                                           t.level = t.level,
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
  beta_ordination_stratified_results <- generate_beta_ordination_single(data.obj = data.obj,
                                                           dist.obj = dist.obj,
                                                           pc.obj = pc.obj,
                                                           subject.var = subject.var,
                                                           time.var = time.var,
                                                           t.level = t.level,
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

## 3.2 Permanova test

```{r beta-diversity-analysis, message=FALSE, results='asis'}
beta_test_results <- generate_beta_test_single(data.obj = data.obj,
                                               dist.obj = dist.obj,
                                               dist.name = dist.name,
                                               time.var = time.var,
                                               t.level = t.level,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars)
```



```{r p-tab-results-and-permanova-analysis, echo=FALSE, message=FALSE, results='asis'}

if (!is.null(test.adj.vars)) {
    adj.description <- sprintf(' while adjusting for covariates %s', paste(test.adj.vars, collapse=' and '))
} else {
    adj.description <- ''
}

if (length(dist.name) == 1) {
    cat(sprintf('\n In this analysis, we employed the PermanovaG2 function from the GUniFrac package to assess the impact of %s on beta diversity using the %s distance metric%s.\n', group.var, dist.name[1], adj.description))
} else {
    cat(sprintf('\n In this analysis, we utilized the PermanovaG2 function from the GUniFrac package to evaluate the effect of %s on beta diversity leveraging multiple distance matrices, specifically %s %s. Additionally, an omnibus test was conducted to combine the power from these matrices.\n', group.var, paste(dist.name, collapse=' and '), adj.description))
}

# Define a function to report the significance of beta diversity based on group.var
report_beta_significance <- function(data_frame, distance) {

  # Extracting rows related to group.var
  group_data <- data_frame[data_frame$Distance == distance & data_frame$Variable == group.var,]

  for(row in 1:nrow(group_data)) {

    # Fetch the p-value
    p_val <- as.numeric(group_data[row, 'P.Value'])

    # Describing significance
    if(p_val < 0.05) {
      cat(sprintf('\n Based on the %s distance metric, the variable %s has a statistically significant effect on beta diversity, with a p-value of %.3f. ', distance, group.var, p_val))
    } else {
      cat(sprintf('\n Based on the %s distance metric, the variable %s did not have a statistically significant effect on beta diversity, with a p-value of %.3f. ', distance, group.var, p_val))
    }
  }

}

# Initialize the sub-section counter for beta diversity
counter <- 1

# Report significance for each distance metric
distance_metrics <- unique(beta_test_results$aov.tab$Distance)
for(distance in distance_metrics) {

  # Skip distance named 'Total'
  if(distance != 'Total') {
    # Print with updated counter and distance name
    cat(sprintf('\n### 3.2.%d %s distance \n\n', counter, ifelse(distance == 'BC', 'Bray-Curtis', distance)))

    # Report significance
    report_beta_significance(data_frame = beta_test_results$aov.tab, distance = distance)
    cat('\n')

    output <- pander::pander(beta_test_results$aov.tab %>% filter(Distance == distance) %>% select(-all_of(c('Distance'))))
    cat(output)


    # Increment the counter
    counter <- counter + 1
  }
}

# Reporting omnibus test significance ONLY if length of dist.name is more than 1
if(length(dist.name) > 1) {
  omni_p_val <- beta_test_results$p.tab[beta_test_results$p.tab$Term == group.var, 'omni.p.value']
  if(omni_p_val < 0.05) {
    cat(sprintf('\n### 3.2.%d Omnibus distance \n\n', counter))
    cat(sprintf('\n The omnibus test indicates that the variable %s has a statistically significant effect on beta diversity across the combined distance matrices, with a p-value of %.3f.\n\n', group.var, omni_p_val))
    output <- pander::pander(beta_test_results$p.tab)
    cat(output)
  } else {
    cat(sprintf('\n### 3.2.%d Omnibus distance \n\n', counter))
    cat(sprintf('\n The omnibus test indicates that the variable %s did not have a statistically significant effect on beta diversity across the combined distance matrices, with a p-value of %.3f.\n\n', group.var, omni_p_val))
    output <- pander::pander(beta_test_results$p.tab)
    cat(output)
  }
}

```

# 4. Feature-level Analysis

```{r Check-and-Select-Rarefaction2, echo=FALSE, message=FALSE, results='asis'}

if (feature.analysis.rarafy) {
    cat('Rarefaction has been enabled for feature-level analysis.\n\n',
        'Reason: The observed abundance of rare/low-abundance features can be strongly influenced by sequence depth. ',
        'Rarefaction is an effective method to control the effect of sequence depth variation. ',
        'By employing rarefaction, we can potentially increase the power of detecting those rare/low-abundance features, ',
        'even though it introduces some variation due to under-sampling.\n',
        'In essence, this step helps to ensure more accurate and consistent results across samples with varying sequence depths.\n\n',
        'If you do not wish to perform rarefaction during feature-level analysis, please turn feature.analysis.rarafy to FALSE.\n')

    data.obj <- rarefy.data.obj
} else {
  data.obj <- original.data.obj
}

```

## 4.1 Differential abundance test

```{r taxa-test-execution, message=FALSE, results='asis', fig.align='center'}
taxa_test_results <- generate_taxa_test_single(data.obj = data.obj,
                                               time.var = time.var,
                                               t.level = t.level,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = test.feature.level,
                                               feature.dat.type = feature.dat.type,
                                               feature.sig.level = feature.sig.level,
                                               feature.mt.method = feature.mt.method,
                                               ...)
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

cat(sprintf('In this analysis, we utilized the LinDA linear model to investigate potential differences in abundance. Specifically, we tested the effect of variables %s for different taxa, while adjusting for other covariates.\n\n', group.var))

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

cat(sprintf('\n\nThe differential abundance test results for features have been saved in the current working directory. Each taxa rank and its corresponding comparison have their own file named with the prefix: %s followed by the taxon rank, the comparison, and the file extension %s. Please refer to these files for more detailed data.', filename_prefix, file_ext))

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
significant_vars <- as.vector(significant_vars)

```

## 4.2 Data visualization(significant features)

### 4.2.1 Significant features boxplot

```{r taxa-significant-boxplot-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 15, results='asis'}

taxa_boxplot_results <- generate_taxa_boxplot_single(
                                              data.obj = data.obj,
                                              subject.var = subject.var,
                                              time.var = time.var,
                                              t.level = t.level,
                                              group.var = group.var,
                                              strata.var = strata.var,
                                              feature.level = test.feature.level,
                                              feature.dat.type = feature.dat.type,
                                              features.plot = significant_vars,
                                              top.k.plot = NULL,
                                              top.k.func = NULL,
                                              transform = feature.box.axis.transform,
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
taxa_boxplot_results

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_single(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t.level = t.level,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = test.feature.level,
                                   features.plot = significant_vars,
                                   feature.dat.type = feature.dat.type,
                                   top.k.plot = NULL,
                                   top.k.func = NULL,
                                   transform = feature.box.axis.transform,
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

```{r taxa-indiv-boxplot-pdf-name, echo=FALSE, message=FALSE, results='asis'}
pdf_name <- paste0(
  'taxa_indiv_boxplot_single',
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

pdf_name <- paste0(pdf_name, '_', test.feature.level, '.pdf')

cat(paste0('The boxplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
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
