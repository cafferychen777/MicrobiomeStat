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

#' Generate Section 2: Alpha Diversity Analysis for Pair Report
#'
#' @return Template string for alpha diversity section
#' @noRd
generate_pair_report_section_alpha <- function() {
  "
# 2. Alpha diversity analysis

## 2.1 Data visualization

### 2.1.1 Alpha diversity boxplot

```{r alpha-boxplot-long-generation, message=FALSE, fig.align='center', results=result.output, fig.width=8, fig.height=3}
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

```{r alpha-diversity-change-boxplot, message=FALSE, fig.align='center', results=result.output}
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

```{r alpha-diversity-change-boxplot-print, message=FALSE, fig.align='center', results=result.output, echo = FALSE, fig.width=8, fig.height=3}
if (is.function(alpha.change.func)) {
  cat('The changes from change.base were computed using a custom function provided by the user.\\n\\n')
} else if (alpha.change.func == 'log fold change') {
  cat('The changes from change.base were computed as the log2 fold change of alpha diversity at the current timepoint versus change.base.\\n\\n')
} else {
  cat('The changes from change.base were computed as the absolute difference in alpha diversity at the current timepoint versus change.base.\\n\\n')
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
            '\\n An ANOVA test of the null hypothesis of no group difference among the %s levels produces a p-value of %.3f.\\n\\n',
            length(group_levels),
            p_val_group_var
        ))

        p_val_group_var_time_var <- data_frame[data_frame$Term == paste0(group.var, ':', time.var),]$P.Value
        cat(sprintf(
            '\\n An ANOVA test of the null hypothesis of no trend difference among the %s groups produces a p-value of %.3f.\\n\\n',
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
                    '\\n Based on the linear mixed effects model, a significant group difference was observed between %s and %s of the variable %s, with a p-value of %.3f.',
                    reference_level,
                    level,
                    group.var,
                    p_val
                ))
            } else {
                cat(sprintf(
                    '\\n Based on the linear mixed effects model, no significant group difference was detected between %s and %s of the variable %s, with a p-value of %.3f.',
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
      cat(sprintf('\\n Based on the linear mixed effects model, a significant trend difference was observed between %s and %s of the variable %s, with a p-value of %.3f.', reference_level, level, group.var, p_val))
    } else {
      cat(sprintf('\\n Based on the linear mixed effects model, no significant trend difference was detected between %s and %s of the variable %s, with a p-value of %.3f.', reference_level, level, group.var, p_val))
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
  cat(sprintf('\\n\\n### 2.2.%d %s index \\n\\n', counter, firstToUpper(ifelse(index_name == 'observed_species', 'observed species', index_name))))
  cat('\\n')

  # Report significance
  report_significance(data_frame = alpha_test_results[[index_name]], group.var = group.var)
  cat('\\n\\n')

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
    cat(sprintf('\\n In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on change.\\n', group.var))
} else {
    cat(sprintf('\\n In this analysis, we utilized a general linear model to examine the influence of the variable %s on change.\\n', group.var))
}

if (is.function(alpha.change.func)) {
    cat('The alpha diversity change is calculated using a custom function supplied by the user to compute the rate of change between consecutive time points.\\n\\n')
} else {
    if (alpha.change.func == 'log fold change') {
        cat('The alpha diversity change is calculated by taking the logarithm of the fold change between consecutive time points.\\n\\n')
    } else {
        cat('The alpha diversity change is calculated by computing the direct difference in alpha diversity between consecutive time points.\\n\\n')
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
          cat(sprintf('\\n Based on the general linear model, the level %s of the variable %s significantly differs from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
      } else {
          cat(sprintf('\\n Based on the general linear model, the level %s of the variable %s did not significantly differ from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
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
  cat(sprintf('\\n### 2.3.%d %s index \\n\\n', counter, firstToUpper(ifelse(index_name == 'observed_species', 'observed species', index_name))))
  cat('\\n')

  # Report significance
  report_change_significance(data_frame = alpha_change_test_results[[index_name]], group.var = group.var)
  cat('\\n')

  output <- pander::pander(alpha_change_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```
"
}

#' Generate Section 3: Beta Diversity Analysis for Pair Report
#'
#' @return Template string for beta diversity section
#' @noRd
generate_pair_report_section_beta <- function() {
  "
# 3. Beta diversity analysis

## 3.1 Data visualization

### 3.1.1 Beta diversity ordinationplot

```{r beta-ordination-pair-generation, message=FALSE, fig.align='center', warning = FALSE, fig.width = 14, fig.height = 12, results=result.output}
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
```

```{r beta-ordination-pair-generation-2, message=FALSE, fig.align='center', warning = FALSE, fig.width = 16, fig.height = 8, results=result.output}
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

```{r beta-diversity-change-boxplot-print, message=FALSE, fig.align='center', results='hide', results=result.output, fig.width=8, fig.height=3, echo = FALSE}
cat(sprintf('\\n Beta change represents the distance of each subject from their change.base.\\n\\n'))
beta_change_boxplot_results
```

### 3.1.3 Beta diversity PC boxplot

```{r beta-pc-boxplot-pair-generation, message=FALSE, fig.align='center', results=result.output, fig.width=8, fig.height=3}
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

```{r pc-change-boxplot-pairs, message=FALSE, fig.align='center', results='hide', results=result.output, fig.width=8, fig.height=3}
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

cat(sprintf('\\n Beta change represents the distance of each subject from their change.base.\\n\\n'))

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
      cat(sprintf('\\n Based on the general linear model, the level %s of the variable %s significantly differs from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
    } else {
      cat(sprintf('\\n Based on the general linear model, the level %s of the variable %s did not significantly differ from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
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
  cat(sprintf('\\n### 3.3.%d %s distance \\n\\n', counter, ifelse(index_name == 'BC', 'Bray-Curtis', index_name)))

  report_beta_change_significance(beta_change_test_results[[index_name]], group.var)

  cat('\\n')
  output <- pander::pander(beta_change_test_results[[index_name]])
  cat(output)

  counter <- counter + 1
}

```
"
}

#' Generate Section 4: Feature-level Analysis for Pair Report
#'
#' @return Template string for feature-level analysis section
#' @noRd
generate_pair_report_section_taxa <- function() {
  "
# 4. Feature-level Analysis

```{r Check-and-Select-Rarefaction2, echo=FALSE, message=FALSE, results='asis'}

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

```{r taxa-test-results-display, echo=TRUE, message=FALSE, results=result.output, warning = FALSE, fig.align='center', fig.width = 6.5, fig.height = 6.5}
volcano_plots <- generate_taxa_volcano_single(
                                  data.obj = data.obj,
                                  group.var = group.var,
                                  test.list = taxa_test_results,
                                  feature.sig.level = feature.sig.level,
                                  feature.mt.method = feature.mt.method
)
volcano_plots
```

```{r taxa-test-results-display2, echo=FALSE, message=FALSE, results='asis', warning = FALSE, fig.align='center', fig.width = 6.5, fig.height = 6.5}
cat(sprintf('In this analysis, we utilized the LinDA linear model to investigate potential differences in trend. Specifically, we tested the effect of the variable %s and the interaction between %s and %s for different taxa, while adjusting for other covariates.\\n\\n',
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
        cat('\\n')
        output <- pander::pander(interaction_terms_results)
        cat(output)
    }
  }
}

filename_prefix <- 'taxa_test_results_'
file_ext <- '.csv'

# Extract the directory path from output.file
output_dir <- dirname(output.file)
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

for(taxon_rank in names(taxa_test_results)) {

    comparisons <- names(taxa_test_results[[taxon_rank]])

    for(comparison in comparisons) {

        file_name <- paste0(filename_prefix, taxon_rank, '_', gsub(' ', '_', gsub('/', '_or_', comparison)), file_ext)

        # Include the output directory in the file path
        file_path <- file.path(output_dir, file_name)

        write.csv(taxa_test_results[[taxon_rank]][[comparison]],
                  file = file_path,
                  row.names = FALSE)
    }
}

cat(sprintf('\\n\\nThe differential abundance test results for features have been saved in the directory: %s. Each taxa rank and its corresponding comparison have their own file named with the prefix: %s followed by the taxon rank, the comparison, and the file extension %s. Please refer to these files for more detailed data.', output_dir, filename_prefix, file_ext))

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

```{r taxa-change-test-results-display, echo=TRUE, message=FALSE, results=result.output, fig.align='center', fig.width = 6.5, fig.height = 6.5}

change_volcano_plots <- generate_taxa_volcano_single(data.obj = data.obj,
                                                     group.var = group.var,
                                                     test.list = taxa_change_test_results,
                                                     feature.sig.level = feature.sig.level,
                                                     feature.mt.method = feature.mt.method)

change_volcano_plots
```


```{r taxa-change-test-results-display2, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 6.5, fig.height = 6.5}
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
        cat('\\n')
        output <- pander::pander(interaction_terms_results)
        cat(output)
    }
  }
}

filename_prefix <- 'taxa_change_test_results_'
file_ext <- '.csv'

# Extract the directory path from output.file
output_dir <- dirname(output.file)
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

for(taxon_rank in names(taxa_change_test_results)) {

    comparisons <- names(taxa_change_test_results[[taxon_rank]])

    for(comparison in comparisons) {

        file_name <- paste0(filename_prefix, taxon_rank, '_', gsub(' ', '_', gsub('/', '_or_', comparison)), file_ext)

        # Include the output directory in the file path
        file_path <- file.path(output_dir, file_name)

        write.csv(taxa_change_test_results[[taxon_rank]][[comparison]],
                  file = file_path,
                  row.names = FALSE)
    }
}

cat(sprintf('\\n\\nThe change test results for individual feature have been saved in the directory: %s. Each taxa rank and its corresponding comparison have their own file named with the prefix: %s followed by the taxon rank, the comparison, and the file extension %s. Please refer to these files for more detailed data.', output_dir, filename_prefix, file_ext))

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
taxa_indiv_boxplot_results_sig_features <- generate_taxa_indiv_boxplot_long(
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
}

```

```{r taxa-boxplot-pair-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 8, fig.height = 4}
if (length(significant_vars) != 0){
taxa_indiv_boxplot_results_sig_features
}
```

```{r boxplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}
output_dir <- dirname(output.file) # Extract the directory path from output.file

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_long(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = change.base,
                                   ts.levels = NULL,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = test.feature.level,
                                   transform = feature.box.axis.transform,
                                   feature.dat.type = feature.dat.type,
                                   features.plot = NULL,
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

# Loop over each feature level in taxa_indiv_boxplot_results
for (feature_level in names(taxa_indiv_boxplot_results)) {
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
    feature_level,
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

  pdf_name <- paste0(pdf_name, '.pdf')

  full_file_path <- file.path(output_dir, pdf_name)

  pdf(full_file_path, width = pdf.wid, height = pdf.hei)
  lapply(taxa_indiv_boxplot_results[[feature_level]], print)
  dev.off()

  cat(paste0('The boxplot results for individual features at the ', feature_level, ' level can be found at: ', full_file_path, '. Please refer to this file for more detailed visualizations.\\n'))
}
```

### 4.3.2 Significant features boxplot (change)

```{r taxa-change-boxplot-generation, message=FALSE, fig.align='center', fig.width = 8, fig.height = 3, results='asis'}

if (length(significant_vars_change) != 0){
taxa_indiv_change_boxplot_results_sig_features <- generate_taxa_indiv_change_boxplot_pair(
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

}

```

```{r taxa-change-boxplot-pair-print, echo=FALSE, message=FALSE, results=result.output, fig.align='center', fig.width = 8, fig.height = 4}
if (length(significant_vars_change) != 0){
taxa_indiv_change_boxplot_results_sig_features
}
```

```{r change-boxplot-pdf-name-creation, echo=FALSE, message=FALSE, results='asis'}

output_dir <- dirname(output.file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

taxa_indiv_change_boxplot_results <- generate_taxa_indiv_change_boxplot_pair(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   change.base = change.base,
                                   feature.change.func = feature.change.func,
                                   feature.level = test.feature.level,
                                   feature.dat.type = feature.dat.type,
                                   features.plot = NULL,
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

for (feature_level in names(taxa_indiv_change_boxplot_results)) {
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
    feature_level,
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

  full_file_path <- file.path(output_dir, pdf_name)

  pdf(full_file_path, width = pdf.wid, height = pdf.hei)
  lapply(taxa_indiv_change_boxplot_results[[feature_level]], print)
  dev.off()

  cat(paste0('The change boxplot results for individual features at the ', feature_level, ' level can be found at: ', full_file_path, '. Please refer to this file for more detailed visualizations.\\n'))
}

```
"
}
