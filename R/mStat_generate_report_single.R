#' Generate a report for microbial ecology analysis of cross-sectional, longitudinal or paired data at a single time point
#'
#' This function generates a comprehensive report for microbial ecology analysis,
#' including alpha diversity, beta diversity, and taxonomic feature analyses. The function is designed
#' to perform analysis on cross-sectional data, a single time point from longitudinal or paired data.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param group.var Variable name used for grouping samples.
#' @param adj.vars Variables to adjust for in the longitudinal analysis.
#' @param subject.var Variable name used for subject identification.
#' @param time.var Variable name used for time points.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou".
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
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
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param transform A string indicating the transformation to apply to the data before plotting. Options are:
#' - "identity": No transformation (default)
#' - "sqrt": Square root transformation
#' - "log": Logarithmic transformation. Zeros are replaced with half of the minimum non-zero value for each taxon before log transformation.
#' @param output.file Output file name for the report.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A report file containing the microbial ecology analysis results.
#'
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
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#'
#' # Generate a report for microbial ecology analysis
#' mStat_generate_report_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon","simpson"),
#'   dist.name = c("BC",'Jaccard'),
#'   t.level = "1",
#'   transform = "identity",
#'   strata.var = "sex",
#'   feature.level = c("Phylum","Family"),
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 12,
#'   output.file = "mStat_generate_report_single_example.pdf"
#' )
#' }
#' @export
mStat_generate_report_single <- function(data.obj,
                                         dist.obj = NULL,
                                         alpha.obj = NULL,
                                         group.var,
                                         adj.vars,
                                         subject.var,
                                         time.var,
                                         transform,
                                         alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
                                         dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'),
                                         t.level = NULL,
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
author: 'Powered by MicrobiomeStat (Ver 1.1.1)'
date: '`r Sys.Date()`'
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

## 2. Alpha Diversity Analysis

### 2.1 Alpha Diversity Boxplots

```{r alpha-diversity-boxplot-generation, message=FALSE}
alpha_boxplot_results <- generate_alpha_boxplot_single(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
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

### 2.2 Alpha Diversity Test Results

```{r alpha-diversity-test-generation, message=FALSE}
alpha_test_results <- generate_alpha_test_single(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 time.var = time.var,
                                                 t.level = t.level,
                                                 alpha.name = alpha.name,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

```{r alpha-index-analysis, echo=FALSE, message=FALSE, results='asis'}

# 遍历alpha指数并生成解析报告
for (index in names(alpha_test_results)) {

  # 打印指数名称的标题
  cat(paste0('## ', index, ' Index \n'))

  # 提取当前指数的分析结果并输出表格
  cat(as.character(pander::pander(alpha_test_results[[index]])))

  cat(paste0('\n### ', index, ' Index Analysis\n'))

  # 提取并解释每个变量的结果
  results <- alpha_test_results[[index]]
  for (i in 1:nrow(results)) {

    # 提取变量名称和相关统计结果
    term <- results$term[i]
    estimate <- results$Estimate[i]
    p_value <- results$P.Value[i]

    if (is.na(p_value)) next

    # 根据P值决定输出的消息
    if (p_value < 0.05) {
      message <- paste0('The variable ', term, ' has a statistically significant impact on the ',
                        index, ' diversity index with an estimate of ', round(estimate, 2), '.\n')
    } else {
      message <- paste0('The variable ', term, ' does not appear to have a statistically significant effect on the ',
                        index, ' diversity index. The estimate of its effect is ', round(estimate, 2), '.\n')
    }

    # 使用strwrap函数将消息断行
    message_lines <- strwrap(message, width = 100)

    # 打印分析结果
    cat(paste0(message_lines, collapse = '\n\n'), '\n\n')

  }
}

```

## 3. Beta Diversity Analysis

### 3.1 Beta Diversity Ordination

```{r beta-ordination-generation, message=FALSE}
beta_ordination_results <- generate_beta_ordination_single(data.obj = data.obj,
                                                           dist.obj = dist.obj,
                                                           pc.obj = NULL,
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
beta_ordination_results
```

### 3.2 Beta Diversity Test Results

```{r beta-diversity-analysis, message=FALSE, results='asis'}
beta_test_results <- generate_beta_test_single(data.obj = data.obj,
                                               dist.obj = dist.obj,
                                               time.var = time.var,
                                               t.level = t.level,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               dist.name = dist.name)
```

```{r p-tab-results-and-permanova-analysis, echo=FALSE, message=FALSE, results='asis'}
cat('## P-Tab Results \n')
pander::pander(beta_test_results$p.tab)

for (i in 1:nrow(beta_test_results$p.tab)) {

  # 提取变量名称
  term <- beta_test_results$p.tab$Term[i]

  cat(paste0('\n### Beta Diversity PERMANOVA Analysis for Variable: ', term, '\n'))

  # 提取并解释每个距离矩阵的结果
  for (j in 1:length(dist.name)) {

    # 提取相关统计结果
    p_value <- beta_test_results$p.tab[paste0('D', j, '.p.value')][i,]

    # 根据P值决定输出的消息
    if (p_value < 0.05) {
      message <- paste0('The variable ', term, ' has a statistically significant impact on the ',
                        'beta diversity according to the PERMANOVA test with the ', dist.name[j], ' distance matrix.')
    } else {
      message <- paste0('The variable ', term, ' does not appear to have a statistically significant effect on the ',
                        'beta diversity according to the PERMANOVA test with the ', dist.name[j], ' distance matrix.')
    }

    # 使用strwrap函数将消息断行，并在每行末尾添加两个空格
    message_lines <- paste0(strwrap(message, width = 100), '  ')

    # 打印分析结果
    cat(paste(message_lines, collapse = '\n'), '\n')
  }

  # 提取omnibus检验的结果
  p_value <- beta_test_results$p.tab['omni.p.value'][i,]

  # 根据P值决定输出的消息
  if (p_value < 0.05) {
    message <- paste0('The variable ', term, ' has a statistically significant impact on the ',
                      'beta diversity according to the omnibus PERMANOVA test.')
  } else {
    message <- paste0('The variable ', term, ' does not appear to have a statistically significant effect on the ',
                      'beta diversity according to the omnibus PERMANOVA test.')
  }

  # 使用strwrap函数将消息断行，并在每行末尾添加两个空格
  message_lines <- paste0(strwrap(message, width = 100), '  ')

  # 打印分析结果
  cat(paste(message_lines, collapse = '\n'), '\n')
}

cat('\n## AOV-Tab Results \n\n')
pander::pander(beta_test_results$aov.tab)

# 遍历aov.tab并生成解析报告
for (variable in unique(beta_test_results$aov.tab$Variable)) {

  # 提取当前变量的分析结果
  aov_results <- subset(beta_test_results$aov.tab, Variable == variable)

  # 打印变量名称的标题
  cat(paste0('\n### ', variable, ' Variable Analysis\n'))

  # 提取并解释每个距离矩阵的结果
  for (i in 1:nrow(aov_results)) {

    # 提取距离矩阵名称和相关统计结果
    distance <- aov_results$Distance[i]
    p_value <- as.numeric(aov_results$P_Value[i])

    if (is.na(p_value)){
    next
    }

    # 根据P值决定输出的消息
    if (p_value < 0.05) {
      message <- paste0('The variable ', variable, ' has a statistically significant impact on the beta diversity ',
                        'according to the PERMANOVA test with the ', distance, ' distance matrix.\n')
    } else {
      message <- paste0('The variable ', variable, ' does not appear to have a statistically significant effect on the beta diversity ',
                        'according to the PERMANOVA test with the ', distance, ' distance matrix.\n')
    }

    # 使用strwrap函数将消息断行，并在每行末尾添加两个空格
    message_lines <- paste0(strwrap(message, width = 100), '  ')

    # 打印分析结果
    cat(paste(message_lines, collapse = '\n'), '\n')
  }
}

```

## 4. Taxonomic Feature Analysis

### 4.1 Taxa Barplot

```{r taxa-barplot-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_barplot_results <- generate_taxa_barplot_single(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t.level = t.level,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     feature.number = 20,
                                                     base.size = base.size,
                                                     theme.choice = theme.choice,
                                                     custom.theme = custom.theme,
                                                     palette = NULL,
                                                     pdf = pdf,
                                                     file.ann = file.ann,
                                                     pdf.wid = pdf.wid,
                                                     pdf.hei = pdf.hei)
taxa_barplot_results
```

### 4.2 Taxa Dotplot

```{r taxa-dotplot-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_dotplot_results <- generate_taxa_dotplot_single(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t0.level = t.level,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     features.plot = features.plot,
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
taxa_dotplot_results
```

### 4.3 Taxa Heatmap

```{r taxa-heatmap-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_heatmap_results <- generate_taxa_heatmap_single(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t.level = t.level,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     features.plot = features.plot,
                                                     top.k.plot = NULL,
                                                     top.k.func = NULL,
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

### 4.4 Taxa Test

```{r taxa-test-execution, message=FALSE}
taxa_test_results <- generate_taxa_test_single(data.obj = data.obj,
                                               time.var = time.var,
                                               t.level = t.level,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = feature.level,
                                               feature.dat.type = feature.dat.type,
                                               ...)
```

```{r taxa-test-results-display, echo=FALSE, message=FALSE}
cat('## Taxa Test Results \n')
pander::pander(taxa_test_results)
```

### 4.5 Taxa Boxplot for Significant Taxa

```{r taxa-significant-boxplot-generation, message=FALSE, fig.align='center', fig.width = 8, fig.height = 16}
significant_taxa <- taxa_test_results$Taxa[taxa_test_results$Adjusted.P.Value < 0.05]

taxa_boxplot_results <- generate_taxa_boxplot_single(data.obj = data.obj,
                                                               subject.var = subject.var,
                                                               time.var = time.var,
                                                               t.level = t.level,
                                                               group.var = group.var,
                                                               strata.var = strata.var,
                                                               feature.level = feature.level,
                                                               feature.dat.type = feature.dat.type,
                                                               features.plot = significant_taxa,
                                                               top.k.plot = NULL,
                                                               top.k.func = NULL,
                                                               transform = transform,
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
names(taxa_boxplot_results) <- feature.level
taxa_boxplot_results

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_single(data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t.level = t.level,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = feature.level,
                                   features.plot = significant_taxa,
                                   feature.dat.type = feature.dat.type,
                                   top.k.plot = NULL,
                                   top.k.func = NULL,
                                   transform = transform,
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
  feature.level,
  '_',
  'transform_',
  transform,
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

pdf_name <- paste0(pdf_name, '_', feature.level, '.pdf')

cat(paste0('The boxplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
```

"

  rmd_code <- knitr::knit_expand(text = template, data.obj = data.obj,
                          dist.obj = dist.obj, alpha.obj = alpha.obj,
                          group.var = group.var,
                          adj.vars = adj.vars, subject.var = subject.var,
                          time.var = time.var, alpha.name = alpha.name,
                          dist.name = dist.name, t.level = t.level,
                          transform = transform,
                          strata.var = strata.var, base.size = base.size,
                          theme.choice = theme.choice, custom.theme = custom.theme,
                          palette = palette, pdf = pdf, file.ann = file.ann,
                          pdf.wid = pdf.wid, pdf.hei = pdf.hei,
                          prev.filter = prev.filter, abund.filter = abund.filter,
                          feature.level = feature.level, features.plot = features.plot,
                          feature.dat.type = feature.dat.type)

  # 将生成的R Markdown代码写入.Rmd文件
  rmd_file <- tempfile(fileext = ".Rmd")
  writeLines(rmd_code, con = rmd_file)

  # 使用rmarkdown::render()函数将.Rmd文件转化为报告
  report_file <- rmarkdown::render(input = rmd_file, output_file = output.file, quiet = FALSE)

  return(report_file)
}
