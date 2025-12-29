#' Report Section Generators for Pair Analysis
#'
#' Internal functions that generate R Markdown template sections
#' for the pair (two time points) report.
#'
#' @name pair_report_sections
#' @keywords internal
NULL

#' Generate YAML Header for Pair Report
#'
#' @param yaml_output The YAML output format string
#' @return YAML header template string
#' @noRd
generate_pair_report_yaml_header <- function(yaml_output) {
  paste0("
---
title: '`r sub(\".pdf$|.html$\", \"\", basename(output.file))`'
author: '[Powered by MicrobiomeStat (Ver 1.2.1)](http://www.microbiomestat.wiki)'
date: '`r Sys.Date()`'
", yaml_output, "
---
")
}

#' Generate Section 1: Data Overview for Pair Report
#'
#' @return Template string for data overview section
#' @noRd
generate_pair_report_section_overview <- function() {
  "
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

custom_strata_status <- ifelse(is.null(strata.var), 'NULL', toString(strata.var))

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
                                        custom_strata_status,
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

```{r mStat-data-summary, message=FALSE, fig.width = 3, fig.height = 4}
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
  depth <- round(min(colSums(data.obj$feature.tab)))
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
  # Extract tree if faith_pd is requested
  tree <- NULL
  if ('faith_pd' %in% alpha.name) {
    tree <- rarefy.data.obj$tree
  }

  alpha.obj <- mStat_calculate_alpha_diversity(x = rarefy.data.obj$feature.tab, alpha.name = alpha.name, tree = tree)
  cat('alpha.obj is calculated based on the rarefied data.obj. ')
}

if (is.null(dist.obj)){
  dist.obj <- mStat_calculate_beta_diversity(data.obj = rarefy.data.obj, dist.name = dist.name)
  cat('dist.obj is calculated based on the rarefied data.obj.\\n')
}

if (is.null(pc.obj)){
  pc.obj <- mStat_calculate_PC(dist.obj = dist.obj, dist.name = dist.name)
  cat('pc.obj is calculated based on the dist.obj using multi-dimensional scaling.\\n')
}

```

## 1.3 Data visualization (overall)

```{r Check-and-Select-Rarefaction1, echo=FALSE, message=FALSE, results='asis'}

if (feature.analysis.rarafy) {
    cat('Rarefaction has been enabled for feature-level analysis and visualization.\\n\\n',
        'Reason: The observed abundance of rare/low-abundance features can be strongly influenced by the sequencing depth. ',
        'Rarefaction is an effective method to control the effect of sequencing depth variation. ',
        'By employing rarefaction, presence/absence status of the features are more comparable and we can potentially increase the power of detecting those rare/low-abundance features, ',
        'even though it introduces some variation due to under-sampling.\\n',
        'In essence, this step improves comparability across samples across samples with varying sequencing depth.\\n\\n',
        'If you do not wish to perform rarefaction during feature-level analysis, please turn feature.analysis.rarafy to FALSE.\\n')

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

```{r taxa-heatmap-pair-avergae-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 20, fig.height = 12, warning = FALSE}
cat('The following plots display the average proportions for each time point, group, and stratum. \\n\\n')

indiv_list <- lapply(taxa_heatmap_pair_results, function(x) x$indiv)

average_list <- lapply(taxa_heatmap_pair_results, function(x) x$average)

average_list
```

```{r taxa-heatmap-pair-indiv-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 30, fig.height = 15, warning = FALSE}
cat('The following plots display the individual proportions for each sample. \\n\\n')
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

```{r taxa-barplot-pair-avergae-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 25, fig.height = 20, warning = FALSE}
cat('The following plots display the average proportions for each time point, group, and stratum. \\n\\n')

indiv_list <- lapply(taxa_barplot_pair_results, function(x) x$indiv)

average_list <- lapply(taxa_barplot_pair_results, function(x) x$average)

average_list
```

```{r taxa-barplot-pair-indiv-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 25, fig.height = 20, warning = FALSE}
cat('The following plots display the individual proportions for each sample. \\n\\n')
indiv_list
```

### 1.3.3 Feature dotplot

```{r taxa-dotplot-generation, message=FALSE, results='asis'}
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

```{r taxa-dotplot-pair-avergae-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 25, fig.height = 12}
cat('The following plots display the average proportions for each time point, group, and stratum. \\n\\n')
taxa_dotplot_results
```

### 1.3.4 Feature change heatmap

```{r taxa-change-heatmap, message=FALSE, fig.align='center', fig.width = 25, fig.height = 20, results=result.output}
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


cat(' The following plots display the change for each subject. \\n\\n')
indiv_list <- lapply(taxa_change_heatmap_results, function(x) x$indiv)

average_list <- lapply(taxa_change_heatmap_results, function(x) x$average)

indiv_list
```

```{r taxa-change-heatmap-pair-print-indiv, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 25, fig.height = 20}

cat(' The following plots display the average change for each time point, group, and stratum. \\n\\n')
average_list

```

### 1.3.5 Feature change dotplot

```{r taxa-change-dotplot, message=FALSE, results='asis'}
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

```{r taxa-change-dotplot-pair-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 25, fig.height = 20}
if (is.function(feature.change.func)) {
  cat('The changes were computed using a custom function provided by the user.')
} else if (feature.change.func == 'relative change') {
  cat('The changes were relative changes, which were computed as (after.abund - before.abund) / (after.abund + before.abund) so the values lie between [-1, 1].')
} else if (feature.change.func == 'absolute change') {
  cat('The changes were absolute changes, computed as the difference between after.abund and before.abund.')
} else if (feature.change.func == 'log fold change') {
  cat('The changes were log2 fold changes, computed as the logarithm of the ratio of after.abund to before.abund, with a small constant added to avoid taking the log of zero.')
}

cat(' The following plots display the average change for each group, and stratum. \\n\\n')
taxa_change_dotplot_results
```
"
}
