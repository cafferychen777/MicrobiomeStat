#' Generate a report for microbial ecology analysis of paired data
#'
#' This function generates a comprehensive report for microbial ecology analysis,
#' including changes in alpha diversity, beta diversity, and taxonomic features between paired data.
#' The function is specifically designed for analysis of paired data.
#'
#' @param data.obj A MicrobiomeStat Data Object.
#' @param dist.obj A distance object created by mStat_calculate_beta_diversity.
#' @param alpha.obj An alpha diversity object (optional).
#' @param depth Sampling depth (optional).
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
#' @param output.file Output file name for the report.
#' @param features.plot Custom list of taxa names to plot for taxa change analysis.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A report file containing the microbial ecology analysis results for paired data.
#'
#' @export
#'
#' @examples
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
#' mStat_generate_report_pair(
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
#'   feature.level = c("Phylum","Family"),
#'   feature.dat.type = "count",
#'   theme.choice = "bw",
#'   base.size = 12,
#'   output.file = "mStat_generate_report_pair_example.pdf"
#' )
#' @export
mStat_generate_report_pair <- function(data.obj,
                                              dist.obj = NULL,
                                              alpha.obj = NULL,
                                              depth = NULL,
                                              group.var,
                                              adj.vars,
                                              subject.var,
                                              time.var,
                                              alpha.name = c("shannon", "simpson", "observed_species", "chao1", "ace", "pielou"),
                                              dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS'),
                                              change.base,
                                              change.func,
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
                                              feature.level = NULL,
                                              features.plot = NULL,
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

## 1. Alpha Diversity Analysis

### 1.1 Alpha Diversity Boxplots

```{r, message=FALSE, fig.align='center'}
alpha_boxplot_results <- generate_alpha_boxplot_long(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
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

### 1.2 Alpha Diversity Test Results

```{r, message=FALSE}
alpha_test_results <- generate_alpha_test_pair(data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 time.var = time.var,
                                                 alpha.name = alpha.name,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

```{r echo=FALSE, message=FALSE, results='asis'}

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

## 2. Beta Diversity Change Analysis

### 2.1 Beta Diversity Ordination

```{r, message=FALSE, fig.align='center', warning = FALSE}
beta_ordination_results <- generate_beta_ordination_pair(data.obj = data.obj,
                                                           dist.obj = dist.obj,
                                                           pc.obj = NULL,
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
beta_ordination_results
```

### 2.2 Beta Diversity Test Pair

```{r, message=FALSE, fig.align='center'}
beta_test_pair_results <- generate_beta_test_pair(data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  time.var = time.var,
                                                  subject.var = subject.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

```{r echo=FALSE, message=FALSE, results='asis'}
cat('## P-Tab Results \n')
pander::pander(beta_test_pair_results$p.tab)

for (i in 1:nrow(beta_test_pair_results$p.tab)) {

  # 提取变量名称
  term <- beta_test_pair_results$p.tab$Term[i]

  cat(paste0('\n### Beta Diversity PERMANOVA Analysis for Variable: ', term, '\n'))

  # 提取并解释每个距离矩阵的结果
  for (j in 1:length(dist.name)) {

    # 提取相关统计结果
    p_value <- beta_test_pair_results$p.tab[paste0('D', j, '.p.value')][i,]

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
  p_value <- beta_test_pair_results$p.tab['omni.p.value'][i,]

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
pander::pander(beta_test_pair_results$aov.tab)

# 遍历aov.tab并生成解析报告
for (variable in unique(beta_test_pair_results$aov.tab$Variable)) {

  # 提取当前变量的分析结果
  aov_results <- subset(beta_test_pair_results$aov.tab, Variable == variable)

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

### 2.3 Beta Diversity PC Boxplot

```{r, message=FALSE, fig.align='center'}
pc_boxplot_longitudinal_results <- generate_beta_pc_boxplot_long(
  data.obj = data.obj,
  dist.obj = dist.obj,
  pc.obj = NULL,
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

## 3. Taxonomic Feature Analysis

### 3.1 Taxa Change Dotplot

```{r, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_dotplot_results <- generate_taxa_dotplot_pair(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
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

### 3.2 Taxa Heatmap Pair

```{r, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8}
taxa_heatmap_pair_results <- generate_taxa_heatmap_pair(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
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
  palette = palette,
  cluster.rows = NULL,
  cluster.cols = NULL,
  pdf = pdf,
  file.ann = file.ann,
  pdf.wid = pdf.wid,
  pdf.hei = pdf.hei
)
```

### 3.3 Taxa Barplot Pair

```{r, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8, warning = FALSE}
taxa_barplot_pair_results <- generate_taxa_barplot_pair(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
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

taxa_barplot_pair_results
```

### 3.4 Taxa Change Test

```{r, message=FALSE, results='asis'}
taxa_test_results <- generate_taxa_test_pair(data.obj = data.obj,
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

```{r echo=FALSE, message=FALSE}
cat('## Taxa Test Results \n')
pander::pander(taxa_test_results)
```

### 3.5 Taxa Boxplot for Significant Taxa

```{r, message=FALSE, fig.align='center', fig.width = 8, fig.height = 16}
taxa_test_results <- do.call('rbind', taxa_test_results)
significant_taxa <- taxa_test_results$Variable[taxa_test_results$Adjusted.P.Value < 1]

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

```

```{r echo=FALSE, message=FALSE, results='asis'}
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

cat(paste0('The boxplot results for individual taxa or features can be found in the current working directory. The relevant file is named: ', pdf_name, '. Please refer to this file for more detailed visualizations.'))
```

"

rmd_code <- knitr::knit_expand(text = template, data.obj = data.obj,
                        dist.obj = dist.obj, alpha.obj = alpha.obj,
                        depth = depth, group.var = group.var,
                        adj.vars = adj.vars, subject.var = subject.var,
                        time.var = time.var, alpha.name = alpha.name,
                        dist.name = dist.name, change.base = change.base, change.func = change.func,
                        strata.var = strata.var, base.size = base.size,
                        theme.choice = theme.choice, custom.theme = custom.theme,
                        palette = palette, pdf = pdf, file.ann = file.ann,
                        pdf.wid = pdf.wid, pdf.hei = pdf.hei,
                        prev.filter = prev.filter, abund.filter = abund.filter,
                        feature.level = feature.level, features.plot = features.plot,
                        feature.dat.type = feature.dat.type)

rmd_file <- tempfile(fileext = ".Rmd")
writeLines(rmd_code, con = rmd_file)

report_file <- rmarkdown::render(input = rmd_file, output_file = output.file, quiet = TRUE)

return(report_file)
}
