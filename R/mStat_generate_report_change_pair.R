#' Generate a report for microbial ecology analysis of paired data
#'
#' This function generates a comprehensive report for microbial ecology analysis,
#' including changes in alpha diversity, beta diversity, and taxonomic features between paired data.
#' The function is specifically designed for analysis of paired data.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param group.var Variable name used for grouping samples.
#' @param adj.vars Variables to adjust for in the analysis.
#' @param subject.var Variable name used for subject identification.
#' @param time.var Variable name used for time points in paired data.
#' @param alpha.name Names of alpha diversity indices to include in the analysis.
#' @param dist.name Names of beta diversity distance metrics to include in the analysis.
#' @param change.base The base level for calculating changes in paired data.
#' @param change.func The function for calculating changes in paired data.
#' @param strata.var Variable to stratify the analysis by (optional).
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
#' @param prev.filter Prevalence filter for feature analysis.
#' @param abund.filter Abundance filter for feature analysis.
#' @param feature.level Taxonomic level for feature analysis.
#' @param feature.dat.type Data type for feature analysis (count, proportion, or other).
#' @param output.file Output file name for the report.
#' @param top.k.plot Number of top taxa to plot for taxa change analysis.
#' @param top.k.func Function for selecting top taxa (default is based on adjusted p-value).
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A report file containing the microbial ecology analysis results for paired data.
#' @examples
#' \dontrun{
#' library(GUniFrac)
#' library(pheatmap)
#' library(vegan)
#' library(ggh4x)
#'
#'
#' data(peerj32.obj)
#'
#' # Generate a report for microbial ecology analysis
#' mStat_generate_report_change_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon","simpson"),
#'   dist.name = c("BC",'Jaccard'),
#'   change.base = "1",
#'   change.func = "relative difference",
#'   strata.var = "sex",
#'   features.plot = NULL,
#'   feature.level = c("Phylum","Family"),
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 12,
#'   palette = NULL,
#'   output.file = "mStat_generate_report_change_pair_example.pdf"
#' )
#' }
#' @export
mStat_generate_report_change_pair <- function(data.obj,
                                         dist.obj = NULL,
                                         alpha.obj = NULL,
                                         group.var,
                                         adj.vars,
                                         subject.var,
                                         time.var,
                                         alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
                                         dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'),
                                         change.base,
                                         change.func,
                                         top.k.plot = NULL,
                                         top.k.func = NULL,
                                         strata.var = NULL,
                                         base.size = 16,
                                         theme.choice = "prism",
                                         custom.theme = NULL,
                                         palette = NULL,
                                         pdf = TRUE,
                                         file.ann = NULL,
                                         pdf.wid = 11,
                                         pdf.hei = 8.5,
                                         prev.filter = 0,
                                         abund.filter = 0,
                                         features.plot = NULL,
                                         feature.level = NULL,
                                         feature.dat.type = c("count", "proportion", "other"),
                                         output.file,
                                         ...) {

  template <- "
---
title: 'Microbial Ecology Analysis Report'
output:
  pdf_document:
    toc: true
    latex_engine: lualatex
---

## 1. Data Summary and Preparation

```{r mStat-data-summary, message=FALSE}
mStat_results <- mStat_summarize_data_obj(data.obj = data.obj,
                                          time.var = time.var,
                                          group.var = group.var,
                                          palette = palette)

# Display the results
cat('## mStat Results \n')
pander::pander(mStat_results)
```

## 2. Alpha Diversity Change Analysis

### 2.1 Alpha Diversity Change Boxplots

```{r alpha-diversity-change-boxplot, message=FALSE, fig.align='center'}
alpha_change_boxplot_results <- generate_alpha_change_boxplot_pair(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       change.base = change.base,
                                                       change.func = change.func,
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
alpha_change_boxplot_results
```

### 2.2 Alpha Diversity Change Test Results

```{r alpha-diversity-change-test, message=FALSE}
alpha_change_test_results <- generate_alpha_change_test_pair(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 time.var = time.var,
                                                 alpha.name = alpha.name,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars,
                                                 change.base = change.base,
                                                 change.func = change.func)
```

```{r alpha-diversity-change-analysis, echo=FALSE, message=FALSE, results='asis'}

indices <- names(alpha_change_test_results)

for (index in indices) {

  cat(paste0('\n## ', index, ' Index Change\n\n'))

  cat(as.character(pander::pander(alpha_change_test_results[[index]])), '\n')

  results <- alpha_change_test_results[[index]]

  cat(paste0('\n### ', index, ' Index Change Analysis\n\n'))

  for (i in 1:nrow(results)) {

    term <- results$Term[i]
    estimate <- results$Estimate[i]
    p_value <- results$P.Value[i]

    if (p_value < 0.05) {
      message <- paste0('The variable ', term, ' has a statistically significant impact on the ',
                        index, ' diversity index change with an estimate of ', round(estimate, 2), '.\n')
    } else {
      message <- paste0('The variable ', term, ' does not appear to have a statistically significant effect on the ',
                        index, ' diversity index change. The estimate of its effect is ', round(estimate, 2), '.\n')
    }

    message_lines <- strwrap(message, width = 100)

    cat(paste(message_lines, collapse = '\n'), '\n\n')
  }
}

```

## 3. Beta Diversity Change Analysis

### 3.1 Beta Diversity PC Change Boxplot Pairs

```{r pc-change-boxplot-pairs, message=FALSE, fig.align='center', results='hide'}
pc_change_boxplot_pairs <- generate_beta_pc_change_boxplot_pair(
  data.obj = data.obj,
  dist.obj = dist.obj,
  pc.obj = NULL,
  pc.ind = c(1, 2),
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  change.base = change.base,
  change.func = change.func,
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

### 3.2 Beta Diversity Change Boxplot

```{r beta-diversity-change-boxplot, message=FALSE, fig.align='center', results='hide'}
beta_change_boxplot_results <- generate_beta_change_boxplot_pair(data.obj = data.obj,
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
beta_change_boxplot_results
```

### 3.3 Beta Diversity Change Test Results

```{r beta-diversity-change-test, message=FALSE, results='asis'}
beta_change_test_results <- generate_beta_change_test_pair(data.obj = data.obj,
                                               dist.obj = dist.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               dist.name = dist.name,
                                               change.base = change.base)
```

```{r beta-diversity-change-analysis, echo=FALSE, message=FALSE, results='asis'}
distances <- names(beta_change_test_results)

for (distance in distances) {

  cat(paste0('\n## ', distance, ' Results\n\n'))

  cat(as.character(pander::pander(beta_change_test_results[[distance]])), '\n')

  results <- beta_change_test_results[[distance]]

  cat(paste0('\n### ', distance, ' Distance Analysis\n\n'))

  for (i in 1:nrow(results)) {

    term <- results$Term[i]
    estimate <- results$Estimate[i]
    p_value <- results$P.Value[i]

    if (p_value < 0.05) {
      message <- paste0('The variable ', term, ' has a statistically significant impact on the ',
                        distance, ' beta diversity change with an estimate of ', round(estimate, 2), '.\n')
    } else {
      message <- paste0('The variable ', term, ' does not appear to have a statistically significant effect on the ',
                        distance, ' beta diversity change. The estimate of its effect is ', round(estimate, 2), '.\n')
    }

    message_lines <- strwrap(message, width = 100)

    cat(paste(message_lines, collapse = '\n'), '\n\n')
  }
}
```

## 4. Taxonomic Feature Analysis

### 4.1 Taxa Change Dotplot

```{r taxa-change-dotplot, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_change_dotplot_results <- generate_taxa_change_dotplot_pair(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     change.base = change.base,
                                                     change.func = change.func,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     features.plot = features.plot,
                                                     top.k.plot = top.k.plot,
                                                     top.k.func = top.k.func,
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
taxa_change_dotplot_results
```

### 4.2 Taxa Change Heatmap

```{r taxa-change-heatmap, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_heatmap_results <- generate_taxa_change_heatmap_pair(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     change.base = change.base,
                                                     change.func = change.func,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     features.plot = features.plot,
                                                     top.k.plot = top.k.plot,
                                                     top.k.func = top.k.func,
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

### 4.3 Taxa Change Scatterplot

```{r taxa-change-scatterplot, message=FALSE, fig.align='center'}
if (is_continuous_numeric(data.obj[[group.var]])) {
  taxa_change_scatterplot_results <- generate_taxa_indiv_change_scatterplot_pair(data.obj = data.obj,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       change.base = change.base,
                                                       change.func = change.func,
                                                       group.var = group.var,
                                                       strata.var = strata.var,
                                                       feature.level = feature.level,
                                                       features.plot = features.plot,
                                                       feature.dat.type = feature.dat.type,
                                                       top.k.plot = top.k.plot,
                                                       top.k.func = top.k.func,
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
  taxa_change_scatterplot_results
} else {
  message(paste0('The variable ', group.var, ' is not a continuous numeric. The function generate_taxa_indiv_change_scatterplot_pair was not executed.'))
}
```

### 4.4 Taxa Change Test

```{r taxa-change-test-pair, message=FALSE, results='hide'}
taxa_change_test_results <- generate_taxa_change_test_pair(data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               change.base = change.base,
                                               change.func = change.func,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = feature.level,
                                               feature.dat.type = feature.dat.type)
```

```{r taxa-change-test-results-display, echo=FALSE, message=FALSE, results='asis'}
cat('## Taxa Change Test Results \n')
pander::pander(taxa_change_test_results)
```

### 4.5 Taxa Boxplot for Significant Taxa

```{r taxa-change-boxplot, message=FALSE, fig.align='center', fig.width = 8, fig.height = 16}
combined_df <- do.call('rbind', taxa_change_test_results)

significant_taxa <- combined_df$Variable[combined_df$Adjusted.P.Value < 1]

if (!is.null(significant_taxa)){
taxa_change_boxplot_results <- generate_taxa_change_boxplot_pair(data.obj = data.obj,
                                                               subject.var = subject.var,
                                                               time.var = time.var,
                                                               group.var = group.var,
                                                               strata.var = strata.var,
                                                               change.base = change.base,
                                                               change.func = change.func,
                                                               feature.level = feature.level,
                                                               feature.dat.type = feature.dat.type,
                                                               features.plot = significant_taxa,
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
taxa_change_boxplot_results

taxa_indiv_change_boxplot_results <- generate_taxa_indiv_change_boxplot_pair(data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   change.base = change.base,
                                   change.func = change.func,
                                   feature.level = feature.level,
                                   features.plot = significant_taxa,
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

```{r taxa-change-boxplot-pdf, echo=FALSE, message=FALSE, results='asis'}
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
  feature.level,
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

if (!is.null(significant_taxa)){
  cat(paste0('The boxplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
}
```

"

rmd_code <- knitr::knit_expand(text = template, data.obj = data.obj,
                        dist.obj = dist.obj, alpha.obj = alpha.obj,
                        group.var = group.var,
                        adj.vars = adj.vars, subject.var = subject.var,
                        time.var = time.var, alpha.name = alpha.name,
                        dist.name = dist.name, top.k.plot = top.k.plot,
                        top.k.func = top.k.func, change.base = change.base, change.func = change.func,
                        strata.var = strata.var, base.size = base.size,
                        theme.choice = theme.choice, custom.theme = custom.theme,
                        palette = palette, pdf = pdf, file.ann = file.ann,
                        pdf.wid = pdf.wid, pdf.hei = pdf.hei,
                        prev.filter = prev.filter, abund.filter = abund.filter,
                        features.plot = features.plot,
                        feature.level = feature.level,
                        feature.dat.type = feature.dat.type)

rmd_file <- tempfile(fileext = ".Rmd")
writeLines(rmd_code, con = rmd_file)

report_file <- rmarkdown::render(input = rmd_file, output_file = output.file, quiet = FALSE)

return(report_file)
}
