#' @title Generate a Comprehensive Report for Microbial Ecology Analysis of Longitudinal Data
#'
#' @description This function generates a comprehensive report for microbial ecology analysis,
#'              including changes in alpha diversity, beta diversity, and taxonomic features in longitudinal data.
#'
#' @param data.obj A data object created by mStat_convert_phyloseq_to_data_obj.
#' @param dist.obj A distance object created by mStat_calculate_beta_diversity.
#' @param alpha.obj An alpha diversity object (optional).
#' @param group.var Variable name used for grouping samples.
#' @param change.func A function to calculate change between time points, e.g. log2fold change. Default is NULL.
#' @param adj.vars Variables to adjust for in the analysis.
#' @param subject.var Variable name used for subject identification.
#' @param time.var Variable name used for time points in longitudinal data.
#' @param alpha.name Names of alpha diversity indices to include in the analysis.
#' @param dist.name Names of beta diversity distance metrics to include in the analysis.
#' @param t0.level The base level for time points in longitudinal data.
#' @param ts.levels The levels for time points in longitudinal data.
#' @param strata.var Variable to stratify the analysis by (optional).
#' @param base.size Base font size for the generated plots.
#' @param theme.choice Plot theme choice (default: "prism").
#' @param custom.theme Custom ggplot2 theme (optional).
#' @param palette Color palette used for the plots.
#' @param pdf Logical indicating whether to save plots as PDF files (default: TRUE).
#' @param file.ann Annotation text for the PDF file names.
#' @param pdf.wid Width of the PDF plots.
#' @param pdf.hei Height of the PDF plots.
#' @param prev.filter Prevalence filter for feature analysis.
#' @param abund.filter Abundance filter for feature analysis.
#' @param feature.level Taxonomic level for feature analysis.
#' @param feature.dat.type Data type for feature analysis (count, proportion, or other).
#' @param Transform The transformation function applied to the data (default: "log").
#' @param output.file Output file name for the report.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A PDF report file containing the microbial ecology analysis results for longitudinal data.
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
#'   adj.vars = c("sex"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   alpha.name = c("shannon","simpson"),
#'   dist.name = c("BC",'Jaccard'),
#'   t0.level = "1",
#'   ts.levels = "2",
#'   strata.var = "sex",
#'   feature.level = c("Phylum"),
#'   feature.dat.type = "count",
#'   Transform = "log",
#'   theme.choice = "bw",
#'   base.size = 12,
#'   output.file = "path/report.pdf"
#' )
#'
#' mStat_generate_report_long(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = NULL,
#'   alpha.obj = NULL,
#'   pc.obj = NULL,
#'   group.var = "subject_gender",
#'   adj.vars = c("subject_race"),
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   alpha.name = c("shannon","simpson"),
#'   dist.name = c("BC",'Jaccard'),
#'   t0.level = unique(sort(subset_T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = unique(sort(subset_T2D.obj$meta.dat$visit_number))[-1],
#'   strata.var = "subject_race",
#'   feature.level = c("Phylum"),
#'   change.func = "relative difference",
#'   feature.dat.type = "count",
#'   prev.filter = 0.0001,
#'   abund.filter = 0.0001,
#'   Transform = "log",
#'   theme.choice = "bw",
#'   base.size = 12,
#'   output.file = "/Users/apple/Microbiome/Longitudinal/MicrobiomeStat_Paper/报告/mStat_generate_report_long_example.pdf"
#' )
#' }
#'
#' @export
mStat_generate_report_long <- function(data.obj,
                                       alpha.obj = NULL,
                                       dist.obj = NULL,
                                       pc.obj = NULL,
                                       group.var,
                                       strata.var = NULL,
                                       adj.vars = NULL,
                                       subject.var,
                                       time.var,
                                       alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
                                       dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'),
                                       t0.level,
                                       ts.levels,
                                       change.func = "relative difference",
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
                                       feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       Transform = c("identity", "sqrt", "log"),
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

```{r object-pre-calculation, echo=FALSE, message=FALSE, results='asis'}

if (is.null(alpha.obj)){
  alpha.obj <- mStat_calculate_alpha_diversity(x = data.obj$feature.tab, alpha.name = alpha.name)
}

if (is.null(dist.obj)){
  dist.obj <- mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
}

if (is.null(pc.obj)){
  pc.obj <- mStat_calculate_PC(dist.obj = dist.obj, dist.name = dist.name)
}

```

## 2. Alpha Diversity Analysis

### 2.1 Alpha Diversity Boxplots

```{r alpha-boxplot-generation, message=FALSE, fig.align='center'}
alpha_boxplot_results <- generate_alpha_boxplot_long(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t0.level = t0.level,
                                                       ts.levels = ts.levels,
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

### 2.2 Alpha Diversity Spaghettiplots

```{r alpha-spaghettiplot-generation, message=FALSE, fig.align='center'}
alpha_spaghettiplot_results <- generate_alpha_spaghettiplot_long(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t0.level = t0.level,
                                                       ts.levels = ts.levels,
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
alpha_spaghettiplot_results
```

### 2.3 Alpha Diversity Test Results

```{r alpha-test-generation, message=FALSE}
alpha_test_results <- generate_alpha_test_long(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 time.var = time.var,
                                                 t0.level = t0.level,
                                                 ts.levels = ts.levels,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

```{r alpha-diversity-index-analysis, echo=FALSE, message=FALSE, results='asis'}

indices <- names(alpha_test_results)

for (index in indices) {

  # 打印指数名称的标题
  cat(paste0('\n## ', index, ' Index \n\n'))

  # 打印pander函数的结果
  cat(as.character(pander::pander(alpha_test_results[[index]])), '\n')

  # 提取当前指数的分析结果
  results <- alpha_test_results[[index]]

  # 打印指数分析的标题
  cat(paste0('\n### ', index, ' Index Analysis\n\n'))

  # 提取并解释每个变量的结果
  for (i in 1:nrow(results)) {

    # 提取变量名称和相关统计结果
    term <- results$Term[i]
    estimate <- results$Estimate[i]
    p_value <- results$P.Value[i]

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
    cat(paste(message_lines, collapse = '\n'), '\n\n')
  }
}

```

### 2.4 Alpha Diversity Trend Test Results

```{r alpha-trend-test-generation, message=FALSE}
alpha_trend_test_results <- generate_alpha_trend_test_long(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

### 2.5 Alpha Diversity Volatility Test Results

```{r alpha-volatility-test-generation, message=FALSE}
alpha_volatility_test_results <- generate_alpha_volatility_test_long(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

## 3. Beta Diversity Analysis

### 3.1 Beta Diversity Ordination

```{r beta-ordination-generation, message=FALSE, fig.align='center', warning = FALSE}
beta_ordination_results <- generate_beta_ordination_long(data.obj = data.obj,
                                                           dist.obj = dist.obj,
                                                           pc.obj = pc.obj,
                                                           subject.var = subject.var,
                                                           time.var = time.var,
                                                           t0.level = t0.level,
                                                           ts.levels = ts.levels,
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

### 3.2 Beta Diversity PC Boxplot

```{r pc-boxplot-longitudinal-generation, message=FALSE, fig.align='center'}
pc_boxplot_longitudinal_results <- generate_beta_pc_boxplot_long(
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

### 3.3 Beta Diversity Change Spaghetti Plot

```{r spaghettiplot-longitudinal-generation, message=FALSE, fig.align='center'}
spaghettiplot_longitudinal_results <- generate_beta_change_spaghettiplot_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
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

spaghettiplot_longitudinal_results
```



### 3.5 Beta Diversity Trend Test Results

```{r beta-trend-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_trend_test_longitudinal_results <- generate_beta_trend_test_long(data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

### 3.6 Beta Diversity PC Trend Test Results

```{r beta-pc-trend-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_pc_trend_test_longitudinal_results <- generate_beta_pc_trend_test_long(data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  pc.obj = pc.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

### 3.7 Beta Diversity Volatility Test Results

```{r beta-volatility-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_volatility_test_longitudinal_results <- generate_beta_volatility_test_long(data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

### 3.8 Beta Diversity PC Volatility Test Results

```{r beta-pc-volatility-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_pc_volatility_test_longitudinal_results <- generate_beta_pc_volatility_test_long(data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  pc.obj = pc.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

## 4. Taxonomic Feature Analysis

### 4.1 Taxa Areaplot Longitudinal

```{r taxa-areaplot-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_areaplot_long_results <- generate_taxa_areaplot_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = feature.level,
  feature.dat.type = feature.dat.type,
  feature.number = 20,
  base.size = base.size,
  theme.choice = theme.choice,
  custom.theme = custom.theme,
  palette = palette,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)

taxa_areaplot_long_results
```

### 4.2 Taxa Heatmap Longitudinal

```{r taxa-heatmap-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_heatmap_long_results <- generate_taxa_heatmap_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = feature.level,
  feature.dat.type = feature.dat.type,
  features.plot = NULL,
  top.k.plot = NULL,
  top.k.func = NULL,
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

### 4.3 Taxa Change Heatmap Longitudinal

```{r taxa-change-heatmap-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_change_heatmap_long_results <- generate_taxa_change_heatmap_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  t0.level = t0.level,
  ts.levels = ts.levels,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = feature.level,
  feature.dat.type = feature.dat.type,
  features.plot = NULL,
  top.k.plot = NULL,
  top.k.func = NULL,
  change.func = change.func,
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

### 4.4 Taxa Barplot Longitudinal

```{r taxa-barplot-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8, warning = FALSE}
taxa_barplot_long_results <- generate_taxa_barplot_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = feature.level,
  feature.dat.type = feature.dat.type,
  feature.number = 20,
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

taxa_barplot_long_results
```

### 4.5 Taxa Test

```{r taxa-test-longitudinal-generation, message=FALSE, results='asis', warning = FALSE}
taxa_test_results <- generate_taxa_test_long(data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               t0.level = t0.level,
                                               ts.levels = ts.levels,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = feature.level,
                                               feature.dat.type = feature.dat.type,
                                               ...)
```

```{r taxa-test-results-print, echo=FALSE, message=FALSE}
cat('## Taxa Test Results \n')
pander::pander(taxa_test_results)
```

### 4.6 Taxa Trend Test

```{r taxa-trend-test-longitudinal-generation, message=FALSE, results='asis', warning = FALSE}
taxa_trend_test_results <- generate_taxa_trend_test_long(data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = feature.level,
                                               feature.dat.type = feature.dat.type,
                                               ...)
```

```{r taxa-trend-test-results-print, echo=FALSE, message=FALSE}
cat('## Taxa Trend Test Results \n')
pander::pander(taxa_trend_test_results)
```

### 4.7 Taxa Volatility Test

```{r taxa-volatility-test-longitudinal-generation, message=FALSE, results='asis', warning = FALSE}
taxa_volatility_test_results <- generate_taxa_volatility_test_long(data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = feature.level,
                                               feature.dat.type = feature.dat.type,
                                               ...)
```

```{r taxa-volatility-test-results-print, echo=FALSE, message=FALSE}
cat('## Taxa Volatility Test Results \n')
pander::pander(taxa_volatility_test_results)
```

### 4.8 Taxa Boxplot for Significant Taxa

```{r taxa-test-boxplot-longitudinal-generation, message=FALSE, fig.height=20, fig.width=15, fig.align='center'}
taxa_test_results <- do.call('rbind', taxa_test_results)
significant_taxa <- taxa_test_results$Variable[taxa_test_results$Adjusted.P.Value < 0.05]

if (!is.null(significant_taxa)){
  taxa_boxplot_results <- generate_taxa_boxplot_long(data.obj = data.obj,
                                                               subject.var = subject.var,
                                                               time.var = time.var,
                                                               t0.level = t0.level,
                                                               ts.levels = ts.levels,
                                                               group.var = group.var,
                                                               strata.var = strata.var,
                                                               feature.level = feature.level,
                                                               feature.dat.type = feature.dat.type,
                                                               features.plot = significant_taxa,
                                                               Transform = Transform,
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
taxa_boxplot_results

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_long(data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = t0.level,
                                   ts.levels = ts.levels,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = feature.level,
                                   features.plot = significant_taxa,
                                   Transform = Transform,
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

```{r boxplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

if (!is.null(significant_taxa)){
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
          feature.level,
          '_',
          'transform_',
          Transform,
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
        pdf_name <- paste0(pdf_name,'_', feature.level, '.pdf')

cat(paste0('The boxplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
}
```

### 4.9 Taxa Spaghettiplot for Significant Taxa

```{r taxa-spaghettiplot-longitudinal-generation, message=FALSE, fig.height=20, fig.width=15, fig.align='center'}

if (!is.null(significant_taxa)){
  taxa_spaghettiplot_results <- generate_taxa_spaghettiplot_long(data.obj = data.obj,
                                                               subject.var = subject.var,
                                                               time.var = time.var,
                                                               group.var = group.var,
                                                               strata.var = strata.var,
                                                               t0.level = t0.level,
                                                               ts.levels = ts.levels,
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
taxa_spaghettiplot_results

taxa_indiv_spaghettiplot_results <- generate_taxa_indiv_spaghettiplot_long(data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = t0.level,
                                   ts.levels = ts.levels,
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

```{r spaghettiplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

if (!is.null(significant_taxa)){
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
          feature.level,
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

cat(paste0('The spaghettiplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
}

```

"

rmd_code <- knitr::knit_expand(text = template, data.obj = data.obj,
                        dist.obj = dist.obj, alpha.obj = alpha.obj,
                        pc.obj = pc.obj, group.var = group.var,
                        adj.vars = adj.vars, subject.var = subject.var,
                        time.var = time.var, alpha.name = alpha.name,
                        dist.name = dist.name, t0.level = t0.level, ts.levels = ts.levels,
                        strata.var = strata.var, base.size = base.size,
                        theme.choice = theme.choice, custom.theme = custom.theme,
                        palette = palette, pdf = pdf, file.ann = file.ann,
                        pdf.wid = pdf.wid, pdf.hei = pdf.hei, change.func = change.func,
                        prev.filter = prev.filter, abund.filter = abund.filter,
                        feature.level = feature.level,
                        feature.dat.type = feature.dat.type)

rmd_file <- tempfile(fileext = ".Rmd")
writeLines(rmd_code, con = rmd_file)

report_file <- rmarkdown::render(input = rmd_file, output_file = output.file, quiet = FALSE)

return(report_file)
}
