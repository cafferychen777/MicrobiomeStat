#' Generate Comprehensive Longitudinal Microbiome Analysis Report for MicrobiomeStat
#'
#' This function performs extensive and in-depth longitudinal analysis of microbiome data using MicrobiomeStat package and generates a comprehensive PDF report.
#' The report contains thorough statistical analysis along with data visualizations, statistical graphics, and result tables for microbial alpha diversity, beta diversity, taxonomic composition, and their temporal dynamics with MicrobiomeStat.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list).
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
#' @param group.var Character, column name in metadata containing grouping variable, e.g. "treatment". Required if present in metadata.
#' @param strata.var Character, column name in metadata containing stratification variable, e.g "sex". Optional.
#' @param test.adj.vars Character vector, names of columns in the metadata containing covariates to be adjusted for in statistical tests and models, such as linear mixed effects models for longitudinal data analysis. This allows the user to account for the effects of additional variables in assessing the effects of primary variables of interest such as time and groups. Default is NULL, which indicates no covariates are adjusted for in statistical testing.
#' @param vis.adj.vars Character vector, names of columns in the metadata containing covariates to visualize in plots, in addition to the primary variables of interest such as groups. For example, if sex is provided in vis.adj.vars, plots will display facets or colors for different sex groups. This allows visualization of effects across multiple covariates. Default is NULL, which indicates only the primary variables of interest will be visualized without additional covariates.
#' @param subject.var Character, column name in metadata containing subject/sample IDs, e.g. "subject_id". Required.
#' @param time.var Character, column name in metadata containing time variable, e.g. "week". Required.
#' @param t0.level Character or numeric, baseline time point for longitudinal analysis, e.g. "week_0" or 0. Required.
#' @param ts.levels Character vector, names of follow-up time points, e.g. c("week_4", "week_8"). Required.
#' @param feature.change.func A function or character string specifying how to calculate
#' the change from baseline value. This allows flexible options:
#' - If a function is provided, it will be applied to each row to calculate change.
#'   The function should take 2 arguments: value at timepoint t and value at baseline t0.
#' - If a character string is provided, following options are supported:
#'   - 'relative change': (value_t - value_t0) / (value_t + value_t0)
#'   - 'difference': value_t - value_t0
#'   - 'lfc': log2(value_t + 1e-5) - log2(value_t0 + 1e-5)
#' - Default is 'relative change'.
#'
#' If none of the above options are matched, an error will be thrown indicating
#' the acceptable options or prompting the user to provide a custom function.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Feature with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Feature with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param bar.area.feature.no A numeric value indicating the number of top abundant features to retain in both barplot and areaplot. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 20.
#' @param heatmap.feature.no A numeric value indicating the number of top abundant features to retain in the heatmap. Features with average relative abundance ranked below this number will be grouped into 'Other'. Default 20.
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
#' @param feature.mt.method Character, multiple testing method for features, "fdr" or "none", default is "fdr".
#' @param feature.sig.level Numeric, significance level cutoff for highlighting features, default is 0.1.
#' @param feature.box.axis.transform A string indicating the transformation to apply to the y-axis of the feature's boxplot visualization before plotting. Options are:
#' - "identity": No transformation (default).
#' - "sqrt": Square root transformation.
#' - "log": Logarithmic transformation. Zeros are replaced with half of the minimum non-zero value for each taxon before log transformation.
#' @param base.size Numeric, base font size for all plots, default is 16.
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
#' @param palette Character vector or function defining color palette for plots, default is NULL.
#' @param pdf Logical, if TRUE save plots as PDF files, default is TRUE.
#' @param file.ann Character, annotation text to add to PDF plot filenames, default is NULL.
#' @param pdf.wid Numeric, width of PDF plots in inches, default is 11.
#' @param pdf.hei Numeric, height of PDF plots in inches, default is 8.5.
#' @param output.file Character, output PDF report filename (required).
#' @param ... Additional arguments passed to generate_taxa_trend_test_long().
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
#' # Assuming peerj32.obj, dist.obj, alpha.obj are pre-defined objects
#' data(peerj32.obj)
#' mStat_generate_report_long(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   test.adj.vars = c("sex"),
#'   vis.test.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon","observed_species"),
#'   dist.name = c("BC",'Jaccard'),
#'   t0.level = "1",
#'   ts.levels = "2",
#'   strata.var = "sex",
#'   vis.feature.level = c("Phylum"),
#'   test.feature.level = c("Phylum"),
#'   feature.dat.type = "count",
#'   feature.box.axis.transform = "log",
#'   theme.choice = "bw",
#'   base.size = 12,
#'   output.file = "path/report.pdf"
#' )
#' data(subset_T2D.obj)
#' mStat_generate_report_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   pc.obj = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   test.adj.vars = NULL,
#'   vis.adj.vars = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number_num",
#'   alpha.name = c("shannon","observed_species"),
#'   dist.name = c("BC",'Jaccard'),
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   feature.mt.method = "none",
#'   feature.sig.level = 0.1,
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
#'   output.file = "/Users/apple/Microbiome/Longitudinal/MicrobiomeStat_Paper/报告/Omics Analysis Report.pdf"
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
                                       t0.level,
                                       ts.levels,
                                       alpha.obj = NULL,
                                       alpha.name = c("shannon", "observed_species"),
                                       depth = NULL,
                                       dist.obj = NULL,
                                       dist.name = c('BC', 'Jaccard'),
                                       pc.obj = NULL,
                                       prev.filter = 0,
                                       abund.filter = 0,
                                       bar.area.feature.no = 15,
                                       heatmap.feature.no = 20,
                                       vis.feature.level = NULL,
                                       test.feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       feature.analysis.rarafy = TRUE,
                                       feature.change.func = "relative change",
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
author: 'Powered by MicrobiomeStat (Ver 1.1.1)'
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

# 判断 custom.theme 是否为 NULL
custom_theme_status <- ifelse(is.null(custom.theme), 'NULL', 'Not NULL')

custom_palette_status <- ifelse(is.null(palette), 'NULL', 'Not NULL')

custom_file.ann_status <- ifelse(is.null(file.ann), 'NULL', 'Not NULL')

custom_t0.level_status <- ifelse(is.null(t0.level), 'NULL', t0.level)

custom_ts.levels_status <- ifelse(is.null(ts.levels), 'NULL', toString(ts.levels))

custom_test.adj.vars_status <- ifelse(is.null(test.adj.vars), 'NULL', toString(test.adj.vars))

custom_vis.adj.vars_status <- ifelse(is.null(vis.adj.vars), 'NULL', toString(vis.adj.vars))

# 创建一个数据框，其中包含参数的名称和对应的值
params_data <- data.frame(Parameter = c('data.obj',
                                        'feature.dat.type',
                                        'group.var',
                                        'test.adj.vars',
                                        'vis.adj.vars',
                                        'strata.var',
                                        'subject.var',
                                        'time.var',
                                        't0.level',
                                        'ts.levels',
                                        'alpha.obj',
                                        'alpha.name',
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
                                        custom_t0.level_status,
                                        custom_ts.levels_status,
                                        deparse(substitute(alpha.obj)),
                                        toString(alpha.name),
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

# 使用pander来渲染数据框
pander::pander(params_data)
cat('## 1.2 Summary statistics \n')
```

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

# 2. Alpha diversity analysis

## 2.1 Data visualization

### 2.1.1 Alpha diversity boxplot

```{r alpha-boxplot-generation, message=FALSE, warning = FALSE, fig.align='center', fig.width = 20, fig.height = 8, results='asis'}
alpha_boxplot_results <- generate_alpha_boxplot_long(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       depth = depth,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t0.level = t0.level,
                                                       ts.levels = ts.levels,
                                                       group.var = group.var,
                                                       strata.var = strata.var,
                                                       adj.vars = vis.adj.vars,
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

### 2.1.2 Alpha diversity spaghettiplot

```{r alpha-spaghettiplot-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8, results='asis'}

alpha_spaghettiplot_results <- generate_alpha_spaghettiplot_long(
                                                       data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       depth = depth,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t0.level = t0.level,
                                                       ts.levels = ts.levels,
                                                       group.var = group.var,
                                                       strata.var = strata.var,
                                                       adj.vars = vis.adj.vars,
                                                       base.size = base.size,
                                                       theme.choice = theme.choice,
                                                       custom.theme = custom.theme,
                                                       palette = palette,
                                                       pdf = pdf,
                                                       file.ann = file.ann,
                                                       pdf.wid = pdf.wid,
                                                       pdf.hei = pdf.hei)
alpha_spaghettiplot_results
```

## 2.2 Trend test

```{r alpha-trend-test-generation, message=FALSE}
alpha_trend_test_results <- generate_alpha_trend_test_long(
                                                 data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 depth = depth,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = test.adj.vars)
```

```{r alpha-trend-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# 如果adj.vars是一个向量，将其转换为逗号分隔的字符串，否则设置为空字符串
if (!is.null(adj.vars)) {
    adj_vars_string <- paste0(', including ', paste(adj.vars, collapse = ', '), ' as covariates')
} else {
    adj_vars_string <- ''
}

if (!is.null(group.var)) {
    cat(sprintf('In this analysis, we utilized a linear mixed effects model with both a random intercept and a random slope%s. Specifically, we tested the interaction between the variables %s and %s.\n\n', adj_vars_string, group.var, time.var))
} else {
    cat(sprintf('In this analysis, we utilized a linear mixed effects model with both a random intercept and a random slope%s. Since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of %s.\n\n', adj_vars_string, time.var))
}

# Define a function to report the significance of interaction terms
report_significance <- function(data_frame, group.var, time.var) {
  if (!is.null(group.var)) {
    # Extracting interaction terms
    interaction_terms <- grep(paste0(group.var, '.+:', time.var), data_frame$Term, value = TRUE)

    if (length(interaction_terms) > 1){
      p_val <- data_frame[data_frame$Term == paste0(group.var, ':', time.var),]$P.Value

      if(p_val < 0.05) {
        cat(sprintf('\n Based on the linear mixed effects model, a significant interaction was observed between %s and %s, with a p-value of %.3f.\n\n', time.var, group.var, p_val))
      } else {
        cat(sprintf('\n Based on the linear mixed effects model, no significant interaction was detected between %s and %s, with a p-value of %.3f.\n\n', time.var, group.var, p_val))
      }

    } else {
      for(term in interaction_terms) {
      p_val <- data_frame[data_frame$Term == term,]$P.Value

      level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
      level <- gsub('_', '', level) # Remove any underscores if they exist

      # Describing interaction terms
      if(p_val < 0.05) {
        cat(sprintf('\n Based on the linear mixed effects model, a significant interaction was observed between %s and the level %s of the variable %s, with a p-value of %.3f.\n\n', time.var, level, group.var, p_val))
      } else {
        cat(sprintf('\n Based on the linear mixed effects model, no significant interaction was detected between %s and the level %s of the variable %s, with a p-value of %.3f.\n\n', time.var, level, group.var, p_val))
      }
    }
    }

  } else {
    p_val <- data_frame[data_frame$Term == time.var,]$P.Value

    # Describing the linear trend
    if(p_val < 0.05) {
      cat(sprintf('\n Based on the linear mixed effects model, a significant linear trend with respect to %s was identified, with a p-value of %.3f.\n\n', time.var, p_val))
    } else {
      cat(sprintf('\n Based on the linear mixed effects model, no significant linear trend was observed with respect to %s, with a p-value of %.3f.\n\n', time.var, p_val))
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
for(index_name in names(alpha_trend_test_results)) {
  # Print with updated counter and index name
  cat(sprintf('\n### 2.2.%d %s index \n\n', counter, firstToUpper(index_name)))
  cat('\n')

  # Report significance
  report_significance(alpha_trend_test_results[[index_name]], group.var, time.var)

  output <- pander::pander(alpha_trend_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```

## 2.3 Volatility test

```{r alpha-volatility-test-generation, message=FALSE, results='asis'}
alpha_volatility_test_results <- generate_alpha_volatility_test_long(
                                                 data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 depth = depth,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = test.adj.vars)
```

```{r alpha-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description for volatility
    num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

    if(num_levels > 2) {
        cat(sprintf('\n In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on volatility.\n', group.var))
    } else {
        cat(sprintf('\n In this analysis, we utilized a general linear model to examine the influence of the variable %s on volatility.\n', group.var))
    }

cat('The alpha diversity volatility is calculated by averaging the rate of change in alpha diversity across consecutive time points. Specifically, for each pair of adjacent time points, we compute the difference in alpha diversity, normalize it by the time difference, and then take the average over all such pairs.\n\n')

# Define a function to report the significance of volatility based on group.var
report_volatility_significance <- function(data_frame, group.var) {

  # Extracting terms excluding ANOVA and intercept
terms <- grep(group.var, data_frame$Term, value = TRUE)
terms <- terms[!terms %in% c('(Intercept)', 'Residuals', group.var)]

for(term in terms) {
    p_val <- data_frame[data_frame$Term == term,]$P.Value

    # Extract only the level part from the term by removing the group.var prefix and underscore
    level <- sub(group.var, '', term)

    # Describing significance based on lm model
    if(p_val < 0.05) {
      cat(sprintf('\n Based on the general linear model, the level %s of the variable %s significantly affected the alpha diversity volatility, with a p-value of %.3f.', level, group.var, p_val))
    } else {
      cat(sprintf('\n Based on the general linear model, the level %s of the variable %s did not significantly influence the alpha diversity volatility, with a p-value of %.3f.', level, group.var, p_val))
    }
}

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(num_levels > 2) {
    if(p_val_anova < 0.05) {
      cat(sprintf('The ANOVA test indicated a significant effect of the variable %s on alpha diversity volatility, with a p-value of %.3f.', group.var, p_val_anova))
    } else {
      cat(sprintf('The ANOVA test showed no significant effect of the variable %s on alpha diversity volatility, with a p-value of %.3f.', group.var, p_val_anova))
    }
  }

}

# Initialize the sub-section counter
counter <- 1

# Report significance for each diversity index
for(index_name in names(alpha_volatility_test_results)) {
  # Print with updated counter and index name
  cat(sprintf('\n### 2.3.%d %s index \n\n', counter, firstToUpper(index_name)))
  cat('\n')

  # Report significance
  report_volatility_significance(data_frame = alpha_volatility_test_results[[index_name]], group.var = group.var)
  cat('\n')

  output <- pander::pander(alpha_volatility_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```

# 3. Beta diversity analysis

## 3.1 Data visualization

### 3.1.1 Beta diversity ordinationplot

```{r beta-ordination-generation, message=FALSE, fig.align='center', warning = FALSE, fig.width = 18, fig.height = 8, results='asis'}
beta_ordination_results <- generate_beta_ordination_long(data.obj = data.obj,
                                                           dist.obj = dist.obj,
                                                           pc.obj = pc.obj,
                                                           subject.var = subject.var,
                                                           time.var = time.var,
                                                           t0.level = t0.level,
                                                           ts.levels = ts.levels,
                                                           group.var = group.var,
                                                           strata.var = strata.var,
                                                           adj.vars = vis.adj.vars,
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
```

### 3.1.2 Beta diversity principal coordinate spaghettiplot

```{r pc-boxplot-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8, results='asis'}
pc_boxplot_longitudinal_results <- generate_beta_pc_spaghettiplot_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  pc.obj = pc.obj,
  pc.ind = c(1, 2),
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  group.var = group.var,
  strata.var = strata.var,
  adj.vars = vis.adj.vars,
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

### 3.1.3 Beta diversity change spaghettiplot

```{r spaghettiplot-longitudinal-generation, message=FALSE, fig.align='center', results='asis', fig.width = 20, fig.height = 8}
spaghettiplot_longitudinal_results <- generate_beta_change_spaghettiplot_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  group.var = group.var,
  strata.var = strata.var,
  adj.vars = vis.adj.vars,
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
```

```{r spaghettiplot-longitudinal-print, echo = FALSE, message=FALSE, fig.align='center', results='asis', fig.width = 20, fig.height = 8}
cat(sprintf('\n In this visualization, the beta change represents the distance of each subject from their first/reference time point.\n\n'))

spaghettiplot_longitudinal_results
```

## 3.2 Distance-based trend test

```{r beta-trend-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_trend_test_longitudinal_results <- generate_beta_trend_test_long(
                                                  data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = test.adj.vars,
                                                  dist.name = dist.name)
```

```{r beta-trend-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description
if (!is.null(group.var)) {
    cat(sprintf('\n In this analysis, we utilized a linear mixed effects model to investigate potential interactions. Specifically, we tested the interaction between the variables %s and %s, while considering the distances to the first/reference time point.\n\n', group.var, time.var))
} else {
    cat('\n In this analysis, we utilized a linear mixed effects model. Since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of', time.var, 'only, with respect to the distances to the first/reference time point.\n\n')
}

# Define a function to report the significance of interaction terms for Beta diversity
report_beta_significance <- function(data_frame, group.var, time.var) {
  if (!is.null(group.var)) {
    # Extracting interaction terms
    interaction_terms <- grep(paste0(group.var, '.+:', time.var), data_frame$Term, value = TRUE)

    if (length(interaction_terms) > 1) {
      p_val <- data_frame[data_frame$Term == paste0(group.var, ':', time.var),]$P.Value

      # Describing interaction terms
      if (p_val < 0.05) {
        cat(sprintf(
          '\n Based on the linear mixed effects model, a significant interaction was observed between %s and %s, with regards to the distances to the first/reference time point, with a p-value of %.3f.\n\n',
          time.var, group.var, p_val
        ))
      } else {
        cat(sprintf(
          '\n Based on the linear mixed effects model, no significant interaction was detected between %s and %s, in terms of the distances to the first/reference time point, with a p-value of %.3f.\n\n',
          time.var, group.var, p_val
        ))
      }
    } else {
      for (term in interaction_terms) {
        p_val <- data_frame[data_frame$Term == term,]$P.Value

        level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
        level <- gsub('_', '', level)  # Remove any underscores if they exist

        # Describing interaction terms
        if (p_val < 0.05) {
          cat(sprintf(
            '\n Based on the linear mixed effects model, a significant interaction was observed between %s and the level %s of the variable %s, with regards to the distances to the first/reference time point, with a p-value of %.3f.\n\n',
            time.var, level, group.var, p_val
          ))
        } else {
          cat(sprintf(
            '\n Based on the linear mixed effects model, no significant interaction was detected between %s and the level %s of the variable %s, in terms of the distances to the first/reference time point, with a p-value of %.3f.\n\n',
            time.var, level, group.var, p_val
          ))
        }
      }
    }
  } else {
    p_val <- data_frame[data_frame$Term == time.var,]$P.Value

    # Describing the linear trend
    if (p_val < 0.05) {
      cat(sprintf(
        '\n Based on the linear mixed effects model, a significant linear trend with respect to %s was identified, concerning the distances to the first/reference time point, with a p-value of %.3f.\n\n',
        time.var, p_val
      ))
    } else {
      cat(sprintf(
        '\n Based on the linear mixed effects model, no significant linear trend was observed with respect to %s, when considering the distances to the first/reference time point, with a p-value of %.3f.\n\n',
        time.var, p_val
      ))
    }
  }
}

counter <- 1

# Report significance for each Beta diversity index in beta_trend_test_longitudinal_results
for(index_name in names(beta_trend_test_longitudinal_results)) {
  cat(sprintf('\n### 3.2.%d %s distance \n\n', counter, ifelse(index_name == 'BC', 'Bray–Curtis', index_name)))
  cat('\n')

  report_beta_significance(beta_trend_test_longitudinal_results[[index_name]], group.var, time.var)
  output <- pander::pander(beta_trend_test_longitudinal_results[[index_name]])
  cat(output)
  counter <- counter + 1
}

```

## 3.3 Distance-based volatility test

```{r beta-volatility-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_volatility_test_longitudinal_results <- generate_beta_volatility_test_long(
                                                  data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = test.adj.vars,
                                                  dist.name = dist.name)
```

```{r beta-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description for volatility
num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on beta diversity volatility.\n', group.var))
} else {
    cat(sprintf('In this analysis, we utilized a general linear model to examine the influence of the variable %s on beta diversity volatility.\n', group.var))
}

cat('The beta diversity volatility is calculated by averaging the rate of change in beta diversity across consecutive time points. Specifically, for each pair of adjacent time points, we compute the difference in beta diversity, normalize it by the time difference, and then take the average over all such pairs.\n\n')

# Define a function to report the significance of volatility based on group.var
report_beta_volatility_significance <- function(data_frame, group.var) {

  # Extracting terms excluding ANOVA and intercept
  terms <- grep(group.var, data_frame$Term, value = TRUE)
  terms <- terms[!terms %in% c('(Intercept)', 'Residuals', group.var)]

  for(term in terms) {
    p_val <- data_frame[data_frame$Term == term,]$P.Value

    # Extract only the level part from the term by removing the group.var prefix and underscore
    level <- sub(group.var, '', term)

    # Describing significance based on lm model
    if(p_val < 0.05) {
      cat(sprintf('\n Based on the general linear model, the level %s of the variable %s significantly affected the beta diversity volatility, with a p-value of %.3f.', level, group.var, p_val))
    } else {
      cat(sprintf('\n Based on the general linear model, the level %s of the variable %s did not significantly influence the beta diversity volatility, with a p-value of %.3f.', level, group.var, p_val))
    }
  }

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(num_levels > 2) {
    if(p_val_anova < 0.05) {
      cat(sprintf('The ANOVA test indicated a significant effect of the variable %s on beta diversity volatility, with a p-value of %.3f.', group.var, p_val_anova))
    } else {
      cat(sprintf('The ANOVA test showed no significant effect of the variable %s on beta diversity volatility, with a p-value of %.3f.', group.var, p_val_anova))
    }
  }

}

counter <- 1

for(index_name in names(beta_volatility_test_longitudinal_results)) {
  cat(sprintf('\n### 3.3.%d %s distance \n\n', counter, ifelse(index_name == 'BC', 'Bray–Curtis', index_name)))

  report_beta_volatility_significance(beta_volatility_test_longitudinal_results[[index_name]], group.var)

  cat('\n')
  output <- pander::pander(beta_volatility_test_longitudinal_results[[index_name]])
  cat(output)

  counter <- counter + 1
}

```

# 4. Feature-level Analysis

```{r 4, echo=FALSE, message=FALSE, results='asis'}

if (feature.analysis.rarafy) {
    cat('Rarefaction has been enabled for feature-level analysis.\n\n',
        'Reason: The observed abundance of rare/low-abundance features can be strongly influenced by sequence depth. ',
        'Rarefaction is an effective method to control the effect of sequence depth variation. ',
        'By employing rarefaction, we can potentially increase the power of detecting those rare/low-abundance features, ',
        'even though it introduces some variation due to under-sampling.\n',
        'In essence, this step helps to ensure more accurate and consistent results across samples with varying sequence depths.\n\n',
        'If you do not wish to perform rarefaction during feature-level analysis, please turn feature.analysis.rarafy to FALSE.\n')

    data.obj <- mStat_rarefy_data(data.obj, depth = depth)
}

```

## 4.1 Data visualization(overall)

### 4.1.1 Feature areaplot

```{r taxa-areaplot-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8}
taxa_areaplot_long_results <- generate_taxa_areaplot_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
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

```{r taxa-areaplot-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 25, fig.height = 15}
cat('The following plots display the average proportions for each time point, group, and stratum. \n\n')
taxa_areaplot_long_results
```

### 4.1.2 Feature heatmap

```{r taxa-heatmap-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8, results='hide', warning = FALSE}
taxa_heatmap_long_results <- generate_taxa_heatmap_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
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
  cluster.cols = NULL,
  cluster.rows = NULL,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)
```

```{r taxa-heatmap-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 12}
cat('The following plots display the average proportions for each time point, group, and stratum. \n\n')
taxa_heatmap_long_results
```

### 4.1.3 Feature change heatmap

```{r taxa-change-heatmap-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 25, fig.height = 12, results='hide', warning = FALSE}
taxa_change_heatmap_long_results <- generate_taxa_change_heatmap_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = vis.feature.level,
  feature.dat.type = feature.dat.type,
  features.plot = NULL,
  top.k.plot = heatmap.feature.no,
  top.k.func = 'mean',
  feature.change.func = feature.change.func,
  prev.filter = prev.filter,
  abund.filter = abund.filter,
  base.size = base.size,
  palette = palette,
  cluster.cols = NULL,
  cluster.rows = NULL,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)
```

```{r taxa-change-heatmap-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 12}
if (is.function(feature.change.func)) {
  cat('The changes from t0.level were computed using a custom function provided by the user.')
} else if (feature.change.func == 'relative change') {
  cat('The changes from t0.level were computed as the difference between the current value and t0.level divided by the sum of the two.')
} else if (feature.change.func == 'difference') {
  cat('The changes from t0.level were computed as the difference between the current value and t0.level.')
} else if (feature.change.func == 'lfc') {
  cat('The changes from t0.level were computed as the log2 difference between the current value and t0.level, with a small constant added to avoid taking log of zero.')
}

cat('The following plots display the average changes for each time point, group, and stratum. \n\n')
taxa_change_heatmap_long_results
```

### 4.1.4 Feature barplot

```{r taxa-barplot-longitudinal-generation, message=FALSE, warning = FALSE}
taxa_barplot_long_results <- generate_taxa_barplot_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = vis.feature.level,
  feature.dat.type = feature.dat.type,
  feature.number = bar.area.feature.no,
  t0.level = t0.level,
  ts.levels = ts.levels,
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

```{r taxa-barplot-longitudinal-print, echo=FALSE, message=FALSE, warning = FALSE, results='asis', fig.width = 25, fig.height = 15, fig.align='center'}
cat('The following plots display the average proportions for each time point, group, and stratum. \n\n')
taxa_barplot_long_results
```

## 4.2 Trend test

```{r taxa-trend-test-longitudinal-generation, message=FALSE, results='hide', warning = FALSE}
taxa_trend_test_results <- generate_taxa_trend_test_long(
                                               data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = test.feature.level,
                                               feature.dat.type = feature.dat.type,
                                               ...)
```

```{r taxa-trend-test-results-print, echo=FALSE, message=FALSE, results='asis', warning = FALSE, fig.align='center', fig.width = 5, fig.height = 5}

trend_volcano_plots <- generate_taxa_trend_volcano_long(
                                               data.obj = data.obj,
                                               group.var = group.var,
                                               time.var = time.var,
                                               test.list = taxa_trend_test_results,
                                               feature.sig.level = feature.sig.level,
                                               feature.mt.method = feature.mt.method
                                                  )

trend_volcano_plots

# Initial description
if (!is.null(group.var)) {
    cat(sprintf('In this analysis, we utilized the LinDA linear mixed effects model to investigate potential interactions in the context of trend test. Specifically, we tested the interaction between the variables %s and %s, for different taxa, while adjusting for other covariates.\n\n', group.var, time.var))
} else {
    cat(sprintf('In this analysis, we utilized the LinDA linear mixed effects model for the trend test. For different taxa, since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of %s only, while adjusting for other covariates.\n\n', time.var))
}

# Iterate over each taxonomic rank in taxa_trend_test_results
for(taxon_rank in names(taxa_trend_test_results)) {

    # Filter interaction terms based on the selected multiple testing method
    if (feature.mt.method == 'fdr') {
        interaction_terms_results <- taxa_trend_test_results[[taxon_rank]] %>%
            dplyr::filter(grepl(paste0('^', group.var, '.*', time.var, '$'), Output.Element)) %>%
            filter(Adjusted.P.Value < feature.sig.level)
        p_value_str = 'adjusted p-value'
    } else if (feature.mt.method == 'none') {
        interaction_terms_results <- taxa_trend_test_results[[taxon_rank]] %>%
            dplyr::filter(grepl(paste0('^', group.var, '.*', time.var, '$'), Output.Element)) %>%
            filter(P.Value < feature.sig.level)
        p_value_str = 'p-value'
    }

    # Check if filtered results have rows
    if (nrow(interaction_terms_results) == 0) {
        cat(sprintf('For the taxa investigated under the %s category, no significant interactions were detected between %s and %s using the %s method for p-value adjustment, at a %s threshold of %s.\n\n', taxon_rank, group.var, time.var, feature.mt.method, p_value_str, feature.sig.level))
    } else {
        cat('- Significant features in trend test results \n\n')
        cat(sprintf('For the taxon %s, significant interactions were identified in the trend test Results using the %s method for p-value adjustment, based on a threshold of %s:\n\n', taxon_rank, feature.mt.method, feature.sig.level))
    }
}

pander::pander(interaction_terms_results)

# 指定文件名前缀和后缀
filename_prefix <- 'taxa_trend_test_results_'
file_ext <- '.csv'

# 遍历列表并保存每个data.frame
for(taxon_rank in names(taxa_trend_test_results)) {
    write.csv(taxa_trend_test_results[[taxon_rank]],
              file = paste0(filename_prefix, taxon_rank, file_ext),
              row.names = FALSE)
}

# 通知用户
cat(sprintf('The trend test results for individual taxa or features have been saved in the current working directory. Each taxa rank has its own file named with the prefix: %s followed by the taxon rank and the file extension %s. Please refer to these files for more detailed data.', filename_prefix, file_ext))

```

## 4.3 Volatility test

```{r taxa-volatility-test-longitudinal-generation, message=FALSE, results='hide', warning = FALSE}
taxa_volatility_test_results <- generate_taxa_volatility_test_long(
                                               data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = test.adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = test.feature.level,
                                               feature.dat.type = feature.dat.type
                                               )
```

```{r taxa-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 7, fig.height = 5, warning = FALSE}

volatility_volcano_plots <- generate_taxa_volatility_volcano_long(data.obj = data.obj,
                                                                  group.var = group.var,
                                                                  test.list = taxa_volatility_test_results,
                                                                  feature.sig.level = feature.sig.level,
                                                                  feature.mt.method = feature.mt.method)

volatility_volcano_plots

# Initial description for Feature Volatility
num_levels <- length(unique(data.obj$meta.dat[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, a general linear model followed by ANOVA was employed to test the effect of %s on the volatility of various taxa abundances.\n\n', group.var))
} else {
    cat(sprintf('In this analysis, a general linear model was utilized to investigate the influence of the variable %s on the volatility of various taxa abundances.\n\n', group.var))
}

cat('Feature abundances were transformed using the centered log-ratio (CLR) transformation. For count data, 0.5 was added to all counts before performing the CLR. For proportion data, zeros were replaced by half the minimum non-zero proportion for each taxon.\n\n')

# Function to check and report significance for taxa
report_taxa_volatility_significance <- function(data_frame, group_var) {
    significant_group_var <- data_frame %>%
    filter(
            ifelse(
            num_levels > 2,
            Term == group_var,
            grepl(paste0('^', group_var), Term)
        )
    ) %>%
    filter(P.Value < feature.sig.level)

    if(nrow(significant_group_var) > 0) {
        return(data_frame)
    } else {
        return(NULL)
    }
}

# Iterate over each taxonomic rank in taxa_volatility_test_results
for(taxon_rank in names(taxa_volatility_test_results)) {
    # Obtain significant results for this rank
    significant_results_list <- lapply(taxa_volatility_test_results[[taxon_rank]], function(df) {
        report_taxa_volatility_significance(df, group.var)
    })

    # Filter out NULL entries
    significant_results_list <- Filter(Negate(is.null), significant_results_list)

    # Report significant results for each taxon under the current rank
    if(length(significant_results_list) > 0) {
      cat('- Significant features in volatility test results \n\n')
      output <- pander::pander(significant_results_list)
      cat(output)
    } else {
        cat(sprintf('No significant results were detected for the volatility test at a p-value threshold of %s for %s.\n\n', feature.sig.level, taxon_rank))
    }
}

# 指定文件名前缀和后缀
filename_prefix <- 'feature_volatility_test_results_'
file_ext <- '.csv'

# 遍历主列表的每个元素
for(main_taxon in names(taxa_volatility_test_results)) {
    # 获取次级列表
    sub_list <- taxa_volatility_test_results[[main_taxon]]

    # 遍历次级列表的每个元素
    for(sub_taxon in names(sub_list)) {
        # 获取完整的文件名
        full_filename <- paste0(filename_prefix, main_taxon, '_', sub_taxon, file_ext)

        # 保存data.frame
        write.csv(sub_list[[sub_taxon]], file = full_filename, row.names = FALSE)
    }
}

# 通知用户
cat(sprintf('The volatility test results for individual feature or features have been saved in the current working directory. Each feature rank and sub-rank combination has its own file named with the prefix: %s followed by the main taxon, sub-taxon, and the file extension %s. Please refer to these files for more detailed data.', filename_prefix, file_ext))

```

```{r extract_significant_taxa, echo=FALSE, results='hide'}
# 从taxa_trend_test_results提取具有统计学意义的taxon
significant_taxa_from_trend <- interaction_terms_results$Variable

# 从taxa_volatility_test_results提取具有统计学意义的taxon
significant_taxa_from_volatility <- names(significant_results_list)[sapply(significant_results_list, nrow) > 0]

# 结合并去重
combined_significant_taxa <- unique(c(significant_taxa_from_trend, significant_taxa_from_volatility))
```

## 4.4 Data visualization(significant features)

### 4.4.1 Significant features boxplot

```{r taxa-test-boxplot-longitudinal-generation, message=FALSE, fig.height=10, fig.width=10, fig.align='center', results='asis'}

if (!is.null(combined_significant_taxa)){
  taxa_boxplot_results <- generate_taxa_boxplot_long(
                                              data.obj = data.obj,
                                              subject.var = subject.var,
                                              time.var = time.var,
                                              t0.level = t0.level,
                                              ts.levels = ts.levels,
                                              group.var = group.var,
                                              strata.var = strata.var,
                                              feature.level = test.feature.level,
                                              feature.dat.type = feature.dat.type,
                                              features.plot = combined_significant_taxa,
                                              transform = feature.box.axis.transform,
                                              top.k.plot = top.k.plot,
                                              top.k.func = top.k.func,
                                              prev.filter = prev.filter,
                                              abund.filter = abund.filter,
                                              base.size = 10,
                                              theme.choice = theme.choice,
                                              custom.theme = custom.theme,
                                              palette = palette,
                                              pdf = pdf,
                                              file.ann = file.ann,
                                              pdf.wid = pdf.wid,
                                              pdf.hei = pdf.hei)

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_long(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = t0.level,
                                   ts.levels = ts.levels,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = test.feature.level,
                                   features.plot = combined_significant_taxa,
                                   transform = feature.box.axis.transform,
                                   feature.dat.type = feature.dat.type,
                                   top.k.plot = top.k.plot,
                                   top.k.func = top.k.func,
                                   prev.filter = prev.filter,
                                   abund.filter = abund.filter,
                                   base.size = 10,
                                   theme.choice = theme.choice,
                                   custom.theme = custom.theme,
                                   palette = palette,
                                   pdf = TRUE,
                                   file.ann = file.ann,
                                   pdf.wid = pdf.wid,
                                   pdf.hei = pdf.hei)
}

```

```{r taxa-test-boxplot-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 12, fig.height = 12}
#taxa_boxplot_results
```

```{r boxplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

if (!is.null(combined_significant_taxa)){
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
cat(paste0('\n\n The boxplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
}
```

### 4.4.2 Significant features spaghettiplot

```{r taxa-spaghettiplot-longitudinal-generation, message=FALSE, fig.height=15, fig.width=10, fig.align='center', results='asis'}
if (!is.null(combined_significant_taxa)){
  taxa_spaghettiplot_results <- generate_taxa_spaghettiplot_long(
                                          data.obj = data.obj,
                                          subject.var = subject.var,
                                          time.var = time.var,
                                          group.var = group.var,
                                          strata.var = strata.var,
                                          t0.level = t0.level,
                                          ts.levels = ts.levels,
                                          feature.level = test.feature.level,
                                          feature.dat.type = feature.dat.type,
                                          features.plot = combined_significant_taxa,
                                          top.k.plot = top.k.plot,
                                          top.k.func = top.k.func,
                                          prev.filter = prev.filter,
                                          abund.filter = abund.filter,
                                          base.size = 10,
                                          theme.choice = theme.choice,
                                          custom.theme = custom.theme,
                                          palette = palette,
                                          pdf = pdf,
                                          file.ann = file.ann,
                                          pdf.wid = pdf.wid,
                                          pdf.hei = pdf.hei)

taxa_indiv_spaghettiplot_results <- generate_taxa_indiv_spaghettiplot_long(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = t0.level,
                                   ts.levels = ts.levels,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   change.base = change.base,
                                   feature.change.func = feature.change.func,
                                   feature.level = test.feature.level,
                                   features.plot = combined_significant_taxa,
                                   feature.dat.type = feature.dat.type,
                                   top.k.plot = top.k.plot,
                                   top.k.func = top.k.func,
                                   prev.filter = prev.filter,
                                   abund.filter = abund.filter,
                                   base.size = 10,
                                   theme.choice = theme.choice,
                                   custom.theme = custom.theme,
                                   palette = palette,
                                   pdf = TRUE,
                                   file.ann = file.ann,
                                   pdf.wid = pdf.wid,
                                   pdf.hei = pdf.hei)
}

```

```{r taxa-spaghettiplot-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 10, fig.height = 8}
taxa_spaghettiplot_results
```

```{r spaghettiplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

if (!is.null(combined_significant_taxa)){
  pdf_name <- paste0(
          'taxa_indiv_spaghettiplot_long',
          '_',
          'subject_',
          subject.var,
          '_',
          'time_',
          time.var,
          '_',
          'group_',
          group.var,
          '_',
          'strata_',
          strata.var,
          '_',
          'feature_level_',
          test.feature.level,
          '_',
          'prev_filter_',
          prev.filter,
          '_',
          'abund_filter_',
          abund.filter,
          '_',
          'base_size_',
          base.size,
          '_',
          'theme_choice_',
          theme.choice,
          '_',
          'pdf_wid_',
          pdf.wid,
          '_',
          'pdf_hei_',
          pdf.hei
        )

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, '_', file.ann)
        }

        pdf_name <- paste0(pdf_name, '.pdf')

cat(paste0('The spaghettiplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations. \n\n'))
}

```

"

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
                    quiet = FALSE)

return(report_file)
}
