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
#' @param feature.number The number of features to plot. Default is 15.
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
#'   adj.vars = NULL,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   alpha.name = c("shannon","simpson"),
#'   dist.name = c("BC",'Jaccard'),
#'   t0.level = unique(sort(subset_T2D.obj$meta.dat$visit_number))[1],
#'   ts.levels = unique(sort(subset_T2D.obj$meta.dat$visit_number))[-1],
#'   strata.var = "subject_race",
#'   feature.level = c("Class"),
#'   change.func = "relative difference",
#'   feature.dat.type = "count",
#'   prev.filter = 1e-17,
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
                                       feature.number = 15,
                                       feature.level = NULL,
                                       feature.dat.type = c("count", "proportion", "other"),
                                       Transform = c("identity", "sqrt", "log"),
                                       output.file,
                                       ...) {

  template <- "
---
title: 'Microbial Ecology Analysis Report'
author: 'Powered by MicrobiomeStat'
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
```

```{r mStat-data-summary-print, echo=FALSE, message=FALSE, results='asis'}
# Display the results
cat('## Summary Statistics and Data Overview \n')
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

```{r alpha-boxplot-generation, message=FALSE, warning = FALSE, fig.align='center', fig.width = 20, fig.height = 8, results='asis'}
alpha_boxplot_results <- generate_alpha_boxplot_long(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t0.level = t0.level,
                                                       ts.levels = ts.levels,
                                                       group.var = group.var,
                                                       strata.var = strata.var,
                                                       adj.vars = adj.vars,
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

```{r alpha-spaghettiplot-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8, results='asis'}
alpha_spaghettiplot_results <- generate_alpha_spaghettiplot_long(
                                                       data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       subject.var = subject.var,
                                                       time.var = time.var,
                                                       t0.level = t0.level,
                                                       ts.levels = ts.levels,
                                                       group.var = group.var,
                                                       strata.var = strata.var,
                                                       adj.vars = adj.vars,
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

### 2.3 Alpha Diversity Trend Test Results

```{r alpha-trend-test-generation, message=FALSE}
alpha_trend_test_results <- generate_alpha_trend_test_long(
                                                 data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

```{r alpha-trend-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description
if (!is.null(group.var)) {
    cat(sprintf('In this analysis, we utilized a linear mixed effects model to investigate potential interactions. Specifically, we tested the interaction between the variables %s and %s.\n\n', group.var, time.var))
} else {
    cat('In this analysis, we utilized a linear mixed effects model. Since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of', time.var, 'only.\n\n')
}

# Define a function to report the significance of interaction terms
report_significance <- function(data_frame, group.var, time.var) {
  if (!is.null(group.var)) {
    # Extracting interaction terms
    interaction_terms <- grep(paste0(group.var, '.+:', time.var), data_frame$Term, value = TRUE)

    for(term in interaction_terms) {
      p_val <- data_frame[data_frame$Term == term,]$P.Value

      level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
      level <- gsub('_', '', level) # Remove any underscores if they exist

      # Describing interaction terms
      if(p_val < 0.05) {
        cat(sprintf('Based on the linear mixed effects model, a significant interaction was observed between %s and the level %s of the variable %s, with a p-value of %.3f.\n\n', time.var, level, group.var, p_val))
      } else {
        cat(sprintf('Based on the linear mixed effects model, no significant interaction was detected between %s and the level %s of the variable %s, with a p-value of %.3f.\n\n', time.var, level, group.var, p_val))
      }
    }
  } else {
    p_val <- data_frame[data_frame$Term == time.var,]$P.Value

    # Describing the linear trend
    if(p_val < 0.05) {
      cat(sprintf('Based on the linear mixed effects model, a significant linear trend with respect to %s was identified, with a p-value of %.3f.\n\n', time.var, p_val))
    } else {
      cat(sprintf('Based on the linear mixed effects model, no significant linear trend was observed with respect to %s, with a p-value of %.3f.\n\n', time.var, p_val))
    }
  }
}

# Report significance for each diversity index
for(index_name in names(alpha_trend_test_results)) {
  cat(sprintf('\n## Results for %s diversity index: \n', index_name))
  report_significance(alpha_trend_test_results[[index_name]], group.var, time.var)
}

# Display detailed results
cat('\n## Detailed Results for Alpha Diversity Trend Test: \n')
pander::pander(alpha_trend_test_results)

```

### 2.4 Alpha Diversity Volatility Test Results

```{r alpha-volatility-test-generation, message=FALSE, results='asis'}
alpha_volatility_test_results <- generate_alpha_volatility_test_long(
                                                 data.obj = data.obj,
                                                 alpha.obj = alpha.obj,
                                                 alpha.name = alpha.name,
                                                 time.var = time.var,
                                                 subject.var = subject.var,
                                                 group.var = group.var,
                                                 adj.vars = adj.vars)
```

```{r alpha-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description for volatility
    num_levels <- length(unique(data.obj[[group.var]]))

    if(num_levels > 2) {
        cat(sprintf('In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on volatility.\n\n', group.var))
    } else {
        cat(sprintf('In this analysis, we utilized a general linear model to examine the influence of the variable %s on volatility.\n\n', group.var))
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
      cat(sprintf('Based on the general linear model, the level %s of the variable %s significantly affected the alpha diversity volatility, with a p-value of %.3f.\n\n', level, group.var, p_val))
    } else {
      cat(sprintf('Based on the general linear model, the level %s of the variable %s did not significantly influence the alpha diversity volatility, with a p-value of %.3f.\n\n', level, group.var, p_val))
    }
}

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(!is.na(p_val_anova)) {
    if(p_val_anova < 0.05) {
      cat(sprintf('The ANOVA test indicated a significant effect of the variable %s on alpha diversity volatility, with a p-value of %.3f.\n\n', group.var, p_val_anova))
    } else {
      cat(sprintf('The ANOVA test showed no significant effect of the variable %s on alpha diversity volatility, with a p-value of %.3f.\n\n', group.var, p_val_anova))
    }
  }

}

# Report significance for each diversity index
for(index_name in names(alpha_volatility_test_results)) {
  cat(sprintf('\n## Results for %s diversity index: \n', index_name))
  report_volatility_significance(alpha_volatility_test_results[[index_name]], group.var)
}

cat('## Alpha Diversity Volatility Test Results \n')
pander::pander(alpha_volatility_test_results)
```

## 3. Beta Diversity Analysis

### 3.1 Beta Diversity Ordination

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

```{r pc-boxplot-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8, results='asis'}
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
beta_trend_test_longitudinal_results <- generate_beta_trend_test_long(
                                                  data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

```{r beta-trend-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description
if (!is.null(group.var)) {
    cat(sprintf('In this analysis, we utilized a linear mixed effects model to investigate potential interactions. Specifically, we tested the interaction between the variables %s and %s, while considering the distances to the first/reference time point.\n\n', group.var, time.var))
} else {
    cat('In this analysis, we utilized a linear mixed effects model. Since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of', time.var, 'only, with respect to the distances to the first/reference time point.\n\n')
}

# Define a function to report the significance of interaction terms for Beta diversity
report_beta_significance <- function(data_frame, group.var, time.var) {
  if (!is.null(group.var)) {
    # Extracting interaction terms
    interaction_terms <- grep(paste0(group.var, '.+:', time.var), data_frame$Term, value = TRUE)

    for(term in interaction_terms) {
      p_val <- data_frame[data_frame$Term == term,]$P.Value

      level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
      level <- gsub('_', '', level) # Remove any underscores if they exist

      # Describing interaction terms
      if(p_val < 0.05) {
        cat(sprintf('Based on the linear mixed effects model, a significant interaction was observed between %s and the level %s of the variable %s, with regards to the distances to the first/reference time point, with a p-value of %.3f.\n\n', time.var, level, group.var, p_val))
      } else {
        cat(sprintf('Based on the linear mixed effects model, no significant interaction was detected between %s and the level %s of the variable %s, in terms of the distances to the first/reference time point, with a p-value of %.3f.\n\n', time.var, level, group.var, p_val))
      }
    }
  } else {
    p_val <- data_frame[data_frame$Term == time.var,]$P.Value

    # Describing the linear trend
    if(p_val < 0.05) {
      cat(sprintf('Based on the linear mixed effects model, a significant linear trend with respect to %s was identified, concerning the distances to the first/reference time point, with a p-value of %.3f.\n\n', time.var, p_val))
    } else {
      cat(sprintf('Based on the linear mixed effects model, no significant linear trend was observed with respect to %s, when considering the distances to the first/reference time point, with a p-value of %.3f.\n\n', time.var, p_val))
    }
  }
}

# Report significance for each Beta diversity index in beta_trend_test_longitudinal_results
for(index_name in names(beta_trend_test_longitudinal_results)) {
  cat(sprintf('\n## Results for %s Beta Diversity Index: \n', index_name))
  report_beta_significance(beta_trend_test_longitudinal_results[[index_name]], group.var, time.var)
}

# Display detailed results
cat('\n## Detailed Results for Beta Diversity Trend Test: \n')
  pander::pander(beta_trend_test_longitudinal_results)

```

### 3.6 Beta Diversity PC Trend Test Results

```{r beta-pc-trend-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_pc_trend_test_longitudinal_results <- generate_beta_pc_trend_test_long(
                                                  data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  pc.obj = pc.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

```{r beta-pc-trend-test-results-print, echo=FALSE, message=FALSE, results='asis'}
# Initial description
cat('## Beta Diversity PC Trend Test Results \n')
if (!is.null(group.var)) {
    cat(sprintf('In this analysis, we utilized a linear mixed effects model to investigate potential interactions. Specifically, we tested the interaction between the variables %s and %s, while considering individual principal components.\n\n', group.var, time.var))
} else {
    cat('In this analysis, we utilized a linear mixed effects model. Since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of', time.var, 'only, for each individual principal component.\n\n')
}

report_beta_pc_significance <- function(data_frame, group.var, time.var) {
  if (!is.null(group.var)) {
    # Extracting interaction terms
    interaction_terms <- grep(paste0(group.var, '.+:', time.var), data_frame$Term, value = TRUE)

    for(term in interaction_terms) {
      p_val <- data_frame[data_frame$Term == term,]$P.Value

      level <- gsub(group.var, '', strsplit(term, ':')[[1]][1])
      level <- gsub('_', '', level) # Remove any underscores if they exist

      # Describing interaction terms
      if(p_val < 0.05) {
        cat(sprintf('Significant interaction observed between %s and the level %s of the variable %s, p-value = %.3f.\n\n', time.var, level, group.var, p_val))
      } else {
        cat(sprintf('No significant interaction detected between %s and the level %s of the variable %s, p-value = %.3f.\n\n', time.var, level, group.var, p_val))
      }
    }
  } else {
    p_val <- data_frame[data_frame$Term == time.var,]$P.Value

    # Describing the linear trend for this PC
    if(p_val < 0.05) {
      cat(sprintf('Significant linear trend observed with respect to %s, p-value = %.3f.\n\n', time.var, p_val))
    } else {
      cat(sprintf('No significant linear trend observed with respect to %s, p-value = %.3f.\n\n', time.var, p_val))
    }
  }
}

for(index_name in names(beta_pc_trend_test_longitudinal_results)) {
  cat(sprintf('\n## Results for %s diversity index: \n', index_name))

  # For each principal component in the index
  for(pc_name in names(beta_pc_trend_test_longitudinal_results[[index_name]])) {
    cat(sprintf('\n### Results for principal component: %s\n', pc_name))

    df <- beta_pc_trend_test_longitudinal_results[[index_name]][[pc_name]]
    report_beta_pc_significance(df, group.var, time.var)
  }
}

# Display detailed results
cat('\n## Detailed Results for Beta Diversity Principal Coordinate Trend Test: \n')
pander::pander(beta_pc_trend_test_longitudinal_results)
```

### 3.7 Beta Diversity Volatility Test Results

```{r beta-volatility-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_volatility_test_longitudinal_results <- generate_beta_volatility_test_long(
                                                  data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

```{r beta-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description for volatility
num_levels <- length(unique(data.obj[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on beta diversity volatility.\n\n', group.var))
} else {
    cat(sprintf('In this analysis, we utilized a general linear model to examine the influence of the variable %s on beta diversity volatility.\n\n', group.var))
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
      cat(sprintf('Based on the general linear model, the level %s of the variable %s significantly affected the beta diversity volatility, with a p-value of %.3f.\n\n', level, group.var, p_val))
    } else {
      cat(sprintf('Based on the general linear model, the level %s of the variable %s did not significantly influence the beta diversity volatility, with a p-value of %.3f.\n\n', level, group.var, p_val))
    }
  }

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(!is.na(p_val_anova)) {
    if(p_val_anova < 0.05) {
      cat(sprintf('The ANOVA test indicated a significant effect of the variable %s on beta diversity volatility, with a p-value of %.3f.\n\n', group.var, p_val_anova))
    } else {
      cat(sprintf('The ANOVA test showed no significant effect of the variable %s on beta diversity volatility, with a p-value of %.3f.\n\n', group.var, p_val_anova))
    }
  }

}

# Report significance for each diversity index
for(index_name in names(beta_volatility_test_longitudinal_results)) {
  cat(sprintf('\n## Results for %s diversity index: \n', index_name))
  report_beta_volatility_significance(beta_volatility_test_longitudinal_results[[index_name]], group.var)
}

cat('## Beta Diversity Volatility Test Results \n')
pander::pander(beta_volatility_test_longitudinal_results)

```

### 3.8 Beta Diversity PC Volatility Test Results

```{r beta-pc-volatility-test-longitudinal-generation, message=FALSE, fig.align='center'}
beta_pc_volatility_test_longitudinal_results <- generate_beta_pc_volatility_test_long(
                                                  data.obj = data.obj,
                                                  dist.obj = dist.obj,
                                                  pc.obj = pc.obj,
                                                  subject.var = subject.var,
                                                  time.var = time.var,
                                                  group.var = group.var,
                                                  adj.vars = adj.vars,
                                                  dist.name = dist.name)
```

```{r beta-pc-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description for PC volatility
num_levels <- length(unique(data.obj[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on beta diversity PC volatility.\n\n', group.var))
} else {
    cat(sprintf('In this analysis, we utilized a general linear model to examine the influence of the variable %s on beta diversity PC volatility.\n\n', group.var))
}

cat('The beta diversity PC volatility is calculated by averaging the rate of change in principal components of beta diversity across consecutive time points. Specifically, for each pair of adjacent time points, we compute the difference in principal components, normalize it by the time difference, and then take the average over all such pairs.\n\n')

# Define a function to report the significance of PC volatility based on group.var
report_pc_volatility_significance <- function(data_frame, group.var) {

  # Extracting terms excluding ANOVA and intercept
  terms <- grep(group.var, data_frame$Term, value = TRUE)
  terms <- terms[!terms %in% c('(Intercept)', 'Residuals', group.var)]

  for(term in terms) {
    p_val <- data_frame[data_frame$Term == term,]$P.Value

    # Extract only the level part from the term by removing the group.var prefix and underscore
    level <- sub(group.var, '', term)

    # Describing significance based on lm model
    if(p_val < 0.05) {
      cat(sprintf('Based on the general linear model, the level %s of the variable %s significantly affected the beta diversity PC volatility, with a p-value of %.3f.\n\n', level, group.var, p_val))
    } else {
      cat(sprintf('Based on the general linear model, the level %s of the variable %s did not significantly influence the beta diversity PC volatility, with a p-value of %.3f.\n\n', level, group.var, p_val))
    }
  }

  # Reporting significance for ANOVA
  p_val_anova <- data_frame[data_frame$Term == group.var,]$P.Value
  if(!is.na(p_val_anova)) {
    if(p_val_anova < 0.05) {
      cat(sprintf('The ANOVA test indicated a significant effect of the variable %s on beta diversity PC volatility, with a p-value of %.3f.\n\n', group.var, p_val_anova))
    } else {
      cat(sprintf('The ANOVA test showed no significant effect of the variable %s on beta diversity PC volatility, with a p-value of %.3f.\n\n', group.var, p_val_anova))
    }
  }

}

# Report significance for each distance and its PCs
for(dist_name in names(beta_pc_volatility_test_longitudinal_results)) {
  cat(sprintf('\n## Results for %s distance: \n', dist_name))
  for(pc_name in names(beta_pc_volatility_test_longitudinal_results[[dist_name]])) {
    cat(sprintf('\n### Results for %s PC: \n', pc_name))
    report_pc_volatility_significance(beta_pc_volatility_test_longitudinal_results[[dist_name]][[pc_name]], group.var)
  }
}

cat('## Beta Diversity PC Volatility Test Results \n')
pander::pander(beta_pc_volatility_test_longitudinal_results)

```

## 4. Taxonomic Feature Analysis

### 4.1 Taxa Areaplot Longitudinal

```{r taxa-areaplot-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8}
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
  feature.number = feature.number,
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

```{r taxa-areaplot-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 8}
cat('### Average Version: This plot displays the average proportions for each time point, group, and strata. \n')
taxa_areaplot_long_results
```

### 4.2 Taxa Heatmap Longitudinal

```{r taxa-heatmap-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8, results='hide', warning = FALSE}
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

```{r taxa-heatmap-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 12}
cat('### Average Version: This plot displays the average proportions for each time point, group, and strata. \n')
taxa_heatmap_long_results
```

### 4.3 Taxa Change Heatmap Longitudinal

```{r taxa-change-heatmap-longitudinal-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 12, results='hide', warning = FALSE}
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

```{r taxa-change-heatmap-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 20, fig.height = 12}
if (is.function(change.func)) {
  cat('### Change Calculation: Custom Function\n')
  cat('The changes from t0.level were computed using a custom function provided by the user.\n\n')
} else if (change.func == 'relative difference') {
  cat('### Change Calculation: Relative Difference\n')
  cat('The changes from t0.level were computed as the difference between the current value and t0.level divided by the sum of the two.\n\n')
} else if (change.func == 'difference') {
  cat('### Change Calculation: Difference\n')
  cat('The changes from t0.level were computed as the difference between the current value and t0.level.\n\n')
} else if (change.func == 'lfc') {
  cat('### Change Calculation: Log2 Fold Change (lfc)\n')
  cat('The changes from t0.level were computed as the log2 difference between the current value and t0.level, with a small constant added to avoid taking log of zero.\n\n')
}

cat('### Average Version: This plot displays the average proportions for each time point, group, and strata. \n')
taxa_change_heatmap_long_results
```

### 4.4 Taxa Barplot Longitudinal

```{r taxa-barplot-longitudinal-generation, message=FALSE, warning = FALSE}
taxa_barplot_long_results <- generate_taxa_barplot_long(
  data.obj = data.obj,
  subject.var = subject.var,
  time.var = time.var,
  group.var = group.var,
  strata.var = strata.var,
  feature.level = feature.level,
  feature.dat.type = feature.dat.type,
  feature.number = feature.number,
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

```{r taxa-barplot-longitudinal-print, echo=FALSE, message=FALSE, warning = FALSE, results='asis', fig.width = 15, fig.height = 8, fig.align='center'}
cat('### Average Version: This plot displays the average proportions for each time point, group, and strata. \n')
taxa_barplot_long_results
```

### 4.5 Taxa Trend Test

```{r taxa-trend-test-longitudinal-generation, message=FALSE, results='asis', warning = FALSE}
taxa_trend_test_results <- generate_taxa_trend_test_long(
                                               data.obj = data.obj,
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

```{r taxa-trend-test-results-print, echo=FALSE, message=FALSE, results='asis'}
# Initial description
if (!is.null(group.var)) {
    cat(sprintf('In this analysis, we utilized the LinDA linear mixed effects model to investigate potential interactions in the context of Taxa Trend Test. Specifically, we tested the interaction between the variables %s and %s, for different taxa, while adjusting for other covariates.\n\n', group.var, time.var))
} else {
    cat(sprintf('In this analysis, we utilized the LinDA linear mixed effects model for the Taxa Trend Test. For different taxa, since no group variable (group.var) was provided, we tested the slope, i.e., the linear trend, of %s only, while adjusting for other covariates.\n\n', time.var))
}

# Iterate over each taxonomic rank in taxa_trend_test_results
for(taxon_rank in names(taxa_trend_test_results)) {
    # Filter interaction terms
    interaction_terms_results <- taxa_trend_test_results[[taxon_rank]] %>%
        filter(grepl(paste0(group.var, ':', time.var), Output.Element)) %>%
        filter(Adjusted.P.Value < 0.05)

    # Check if filtered results have rows
    if (nrow(interaction_terms_results) == 0) {
        cat(sprintf('For the investigated taxa under %s, no significant interactions between %s and %s were detected at an adjusted p-value threshold of 0.05.\n\n', taxon_rank, group.var, time.var))
    } else {
        cat(sprintf('## Significant Interactions for %s in Taxa Trend Test Results \n', taxon_rank))
        pander::pander(interaction_terms_results)
    }
}

```

### 4.6 Taxa Volatility Test

```{r taxa-volatility-test-longitudinal-generation, message=FALSE, results='asis', warning = FALSE}
taxa_volatility_test_results <- generate_taxa_volatility_test_long(
                                               data.obj = data.obj,
                                               subject.var = subject.var,
                                               time.var = time.var,
                                               group.var = group.var,
                                               adj.vars = adj.vars,
                                               prev.filter = prev.filter,
                                               abund.filter = abund.filter,
                                               feature.level = feature.level,
                                               feature.dat.type = feature.dat.type
                                               )
```

```{r taxa-volatility-test-results-print, echo=FALSE, message=FALSE, results='asis'}

# Initial description for Taxa Volatility
num_levels <- length(unique(data.obj[[group.var]]))

if(num_levels > 2) {
    cat(sprintf('In this analysis, a general linear model followed by ANOVA was employed to test the effect of %s on the volatility of various taxa abundances.\n\n', group.var))
} else {
    cat(sprintf('In this analysis, a general linear model was utilized to investigate the influence of the variable %s on the volatility of various taxa abundances.\n\n', group.var))
}

cat('Taxa abundances were transformed using the centered log-ratio (CLR) transformation. For count data, 0.5 was added to all counts before performing the CLR. For proportion data, zeros were replaced by half the minimum non-zero proportion for each taxon.\n\n')

# Function to check and report significance for taxa
report_taxa_volatility_significance <- function(data_frame, group_var) {
    significant_group_var <- data_frame %>%
                             filter(Term == group_var) %>%
                             filter(P.Value < 0.05)

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
      cat('## Significant Interactions for Different Taxa in Taxa Volatility Test Results \n')
        for(taxon_name in names(significant_results_list)) {
            cat(sprintf('\n### Results for the %s %s: \n', taxon_rank, taxon_name))
            df_to_display <- as.data.frame(significant_results_list[[taxon_name]])
            pander::pander(df_to_display)
        }
    } else {
        cat(sprintf('No significant results were detected for the taxa volatility at a p-value threshold of 0.05 for %s.\n\n', taxon_rank))
    }
}

```

```{r extract_significant_taxa, echo=FALSE, results='hide'}
# 从taxa_trend_test_results提取具有统计学意义的taxon
significant_taxa_from_trend <- rownames(interaction_terms_results)

# 从taxa_volatility_test_results提取具有统计学意义的taxon
significant_taxa_from_volatility <- names(significant_results_list)[sapply(significant_results_list, nrow) > 0]

# 结合并去重
combined_significant_taxa <- unique(c(significant_taxa_from_trend, significant_taxa_from_volatility))
combined_significant_taxa
```

### 4.7 Taxa Boxplot for Significant Taxa

```{r taxa-test-boxplot-longitudinal-generation, message=FALSE, fig.height=15, fig.width=10, fig.align='center', results='asis'}

if (!is.null(combined_significant_taxa)){
  taxa_boxplot_results <- generate_taxa_boxplot_long(data.obj = data.obj,
                                                     subject.var = subject.var,
                                                     time.var = time.var,
                                                     t0.level = t0.level,
                                                     ts.levels = ts.levels,
                                                     group.var = group.var,
                                                     strata.var = strata.var,
                                                     feature.level = feature.level,
                                                     feature.dat.type = feature.dat.type,
                                                     features.plot = combined_significant_taxa,
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

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_long(
                                   data.obj = data.obj,
                                   subject.var = subject.var,
                                   time.var = time.var,
                                   t0.level = t0.level,
                                   ts.levels = ts.levels,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = feature.level,
                                   features.plot = combined_significant_taxa,
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

```{r taxa-test-boxplot-longitudinal-print, echo=FALSE, message=FALSE, results='asis', fig.align='center', fig.width = 10, fig.height = 8}
taxa_boxplot_results
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

### 4.8 Taxa Spaghettiplot for Significant Taxa

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
                                          feature.level = feature.level,
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
                                   change.func = change.func,
                                   feature.level = feature.level,
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
                        feature.number = feature.number,
                        feature.level = feature.level,
                        feature.dat.type = feature.dat.type)

rmd_file <- tempfile(fileext = ".Rmd")
writeLines(rmd_code, con = rmd_file)

report_file <- rmarkdown::render(input = rmd_file, output_file = output.file, quiet = FALSE)

return(report_file)
}
