#' Report Section Generators for Single Time Point Analysis
#'
#' Internal functions that generate R Markdown template sections
#' for the single time point (cross-sectional) report.
#'
#' @name single_report_sections
#' @keywords internal
NULL

#' Generate YAML Header for Single Report
#'
#' @param yaml_output The YAML output format string
#' @return YAML header template string
#' @noRd
generate_single_report_yaml_header <- function(yaml_output) {
  paste0("
---
title: '`r sub(\".pdf$|.html$\", \"\", basename(output.file))`'
author: '[Powered by MicrobiomeStat (Ver 1.2.1)](http://www.microbiomestat.wiki)'
date: '`r Sys.Date()`'
", yaml_output, "
---
")
}

#' Generate Section 1: Data Overview for Single Report
#'
#' @return Template string for data overview section
#' @noRd
generate_single_report_section_overview <- function() {
  "
# 1. Data overview and summary statistics

## 1.1 Parameter setting

```{r input-parameters-summary, echo=FALSE, message=FALSE, results='asis'}

custom_depth_status <- ifelse(is.null(depth), 'NULL', toString(depth))

custom_theme_status <- ifelse(is.null(custom.theme), 'NULL', 'Not NULL')

custom_palette_status <- ifelse(is.null(palette), 'NULL', 'Not NULL')

custom_file.ann_status <- ifelse(is.null(file.ann), 'NULL', 'Not NULL')

custom_time.var_status <- ifelse(is.null(time.var), 'NULL', time.var)

custom_t.level_status <- ifelse(is.null(t.level), 'NULL', t.level)

custom_test.adj.vars_status <- ifelse(is.null(test.adj.vars), 'NULL', toString(test.adj.vars))

custom_vis.adj.vars_status <- ifelse(is.null(vis.adj.vars), 'NULL', toString(vis.adj.vars))

custom_strata_status <- ifelse(is.null(strata.var), 'NULL', toString(strata.var))

params_data <- data.frame(Parameter = c('data.obj',
                                        'feature.dat.type',
                                        'group.var',
                                        'test.adj.vars',
                                        'vis.adj.vars',
                                        'strata.var',
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
                                        custom_strata_status,
                                        custom_time.var_status,
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

```{r mStat-data-summary, message=FALSE, fig.width = width, fig.height = 4}
mStat_results <- mStat_summarize_data_obj(data.obj = data.obj,
                                          time.var = time.var,
                                          group.var = group.var,
                                          palette = palette)
```

```{r mStat-data-summary-print, echo=FALSE, message=FALSE, results='asis'}
pander::pander(mStat_results)
```

```{r object-pre-calculation, echo=FALSE, message=FALSE, results='asis'}

if (!is.null(time.var) & !is.null(t.level)){
  condition <- paste(time.var, '==', t.level, sep = ' ')
  data.obj <- mStat_subset_data(data.obj, condition = condition)
}

original.data.obj <- data.obj

rarefy.data.obj <- mStat_normalize_data(data.obj = data.obj, method = 'Rarefy', depth = depth)$data.obj.norm

unique_levels <- unique(c(vis.feature.level, test.feature.level))

if (!all(unique_levels %in% names(data.obj$feature.agg.list))) {
  original.data.obj <- mStat_aggregate_by_taxonomy(original.data.obj, feature.level = unique_levels)
  rarefy.data.obj <- mStat_aggregate_by_taxonomy(rarefy.data.obj, feature.level = unique_levels)
}

if (is.null(depth)){
  depth <- round(min(colSums(data.obj$feature.tab)))
  cat(sprintf('No rarefaction depth is specified. The minimum depth, %.0f, is used as the rarefaction depth. ', depth))
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

### 1.3.1 Feature barplot

```{r taxa-barplot-generation, message=FALSE, fig.align='center', fig.width = 20, fig.height = 8, warning = FALSE}
taxa_barplot_results <- generate_taxa_barplot_single(data.obj = data.obj,
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

```{r taxa-barplot-average-print, echo=FALSE, message=FALSE, result = result.output, fig.align='center', fig.width = 25, fig.height = 15, warning = FALSE}
cat('The following plots display the average proportions for each group, and stratum. \\n\\n')
indiv_list <- lapply(taxa_barplot_results, function(x) x$indiv)

average_list <- lapply(taxa_barplot_results, function(x) x$average)

average_list
```

```{r taxa-barplot-indiv-print, echo=FALSE, message=FALSE, result = result.output, fig.align='center', fig.width = 25, fig.height = 15, warning = FALSE}
cat('The following plots display the individual proportions for each group, and stratum. \\n\\n')
indiv_list
```

### 1.3.2 Feature dotplot

```{r taxa-dotplot-generation, message=FALSE, result = result.output, fig.align='center', fig.width = 20, fig.height = 8, warning = FALSE}
taxa_dotplot_results <- generate_taxa_dotplot_single(data.obj = data.obj,
                                                     time.var = time.var,
                                                     t.level = t.level,
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

```{r taxa-heatmap-generation, message=FALSE, fig.align='center', fig.width = 15, fig.height = 8, warning = FALSE}
taxa_heatmap_results <- generate_taxa_heatmap_single(data.obj = data.obj,
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

```{r taxa-heatmap-indiv-print, echo=FALSE, message=FALSE, result = result.output, fig.align='center', fig.width = 20, fig.height = 12, warning = FALSE}
cat('The following plots display the individual proportions for each sample. \\n\\n')
taxa_heatmap_results
```
"
}

#' Generate Section 2: Alpha Diversity Analysis for Single Report
#'
#' @return Template string for alpha diversity section
#' @noRd
generate_single_report_section_alpha <- function() {
  "
# 2. Alpha diversity analysis

## 2.1 Data visualization

### 2.1.1 Alpha diversity boxplot

```{r alpha-boxplot-generation, message=FALSE, warning = FALSE, fig.align='center', fig.width = 16, fig.height = 8, result = result.output}
alpha_boxplot_results <- generate_alpha_boxplot_single(data.obj = data.obj,
                                                       alpha.obj = alpha.obj,
                                                       alpha.name = alpha.name,
                                                       depth = depth,
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

group_levels <- data.obj$meta.dat %>% dplyr::select(!!sym(group.var)) %>% dplyr::pull() %>% as.factor() %>% levels

reference_level <- group_levels[1]

if (!is.null(test.adj.vars)) {
    adj.description <- sprintf(' while adjusting for covariates %s', paste(test.adj.vars, collapse=' and '))
} else {
    adj.description <- ''
}

if(num_levels > 2) {
    cat(sprintf('\\n In this analysis, we employed a general linear model followed by ANOVA to test the effect of %s on alpha diversity%s.\\n', group.var, adj.description))
} else {
    cat(sprintf('\\n In this analysis, we utilized a general linear model to examine the influence of the variable %s on alpha diversity%s.\\n', group.var, adj.description))
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
          cat(sprintf('\\n Based on the general linear model, the level %s of the variable %s significantly differs from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
      } else {
          cat(sprintf('\\n Based on the general linear model, the level %s of the variable %s did not significantly differ from level %s, with a p-value of %.3f. ', level, group.var, reference_level, p_val))
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
  cat(sprintf('\\n### 2.2.%d %s index \\n\\n', counter, firstToUpper(index_name)))
  cat('\\n')

  # Report significance
  report_significance(data_frame = alpha_test_results[[index_name]], group.var = group.var)
  cat('\\n')

  output <- pander::pander(alpha_test_results[[index_name]])
  cat(output)

  # Increment the counter
  counter <- counter + 1
}

```
"
}

#' Generate Section 3: Beta Diversity Analysis for Single Report
#'
#' @return Template string for beta diversity section
#' @noRd
generate_single_report_section_beta <- function() {
  "
# 3. Beta diversity analysis

## 3.1 Data visualization

### 3.1.1 Beta diversity ordinationplot

```{r beta-ordination-generation, message=FALSE, fig.align='center', warning = FALSE, fig.width = 10, fig.height = 8, result = result.output}
beta_ordination_results <- generate_beta_ordination_single(data.obj = data.obj,
                                                           time.var = time.var,
                                                           t.level = t.level,
                                                           group.var = group.var,
                                                           adj.vars = vis.adj.vars,
                                                           strata.var = NULL,
                                                           dist.obj = dist.obj,
                                                           dist.name = dist.name,
                                                           pc.obj = pc.obj,
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
  beta_ordination_stratified_results <- generate_beta_ordination_single(
                                                           data.obj = data.obj,
                                                           time.var = time.var,
                                                           t.level = t.level,
                                                           group.var = group.var,
                                                           adj.vars = vis.adj.vars,
                                                           strata.var = strata.var,
                                                           dist.obj = dist.obj,
                                                           dist.name = dist.name,
                                                           pc.obj = pc.obj,
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
    cat(sprintf('\\n In this analysis, we employed the PermanovaG2 function from the GUniFrac package to assess the impact of %s on beta diversity using the %s distance metric%s.\\n', group.var, dist.name[1], adj.description))
} else {
    cat(sprintf('\\n In this analysis, we utilized the PermanovaG2 function from the GUniFrac package to evaluate the effect of %s on beta diversity leveraging multiple distance matrices, specifically %s %s. Additionally, an omnibus test was conducted to combine the power from these matrices.\\n', group.var, paste(dist.name, collapse=' and '), adj.description))
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
      cat(sprintf('\\n Based on the %s distance metric, the variable %s has a statistically significant effect on beta diversity, with a p-value of %.3f. ', distance, group.var, p_val))
    } else {
      cat(sprintf('\\n Based on the %s distance metric, the variable %s did not have a statistically significant effect on beta diversity, with a p-value of %.3f. ', distance, group.var, p_val))
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
    cat(sprintf('\\n### 3.2.%d %s distance \\n\\n', counter, ifelse(distance == 'BC', 'Bray-Curtis', distance)))

    # Report significance
    report_beta_significance(data_frame = beta_test_results$aov.tab, distance = distance)
    cat('\\n')

    output <- pander::pander(beta_test_results$aov.tab %>% filter(Distance == distance) %>% dplyr::select(-all_of(c('Distance'))))
    cat(output)


    # Increment the counter
    counter <- counter + 1
  }
}

# Reporting omnibus test significance ONLY if length of dist.name is more than 1
if(length(dist.name) > 1) {
  omni_p_val <- beta_test_results$p.tab[beta_test_results$p.tab$Term == group.var, 'omni.p.value']
  if(omni_p_val < 0.05) {
    cat(sprintf('\\n### 3.2.%d Omnibus distance \\n\\n', counter))
    cat(sprintf('\\n The omnibus test indicates that the variable %s has a statistically significant effect on beta diversity across the combined distance matrices, with a p-value of %.3f.\\n\\n', group.var, omni_p_val))
    output <- pander::pander(beta_test_results$p.tab)
    cat(output)
  } else {
    cat(sprintf('\\n### 3.2.%d Omnibus distance \\n\\n', counter))
    cat(sprintf('\\n The omnibus test indicates that the variable %s did not have a statistically significant effect on beta diversity across the combined distance matrices, with a p-value of %.3f.\\n\\n', group.var, omni_p_val))
    output <- pander::pander(beta_test_results$p.tab)
    cat(output)
  }
}

```
"
}

#' Generate Section 4: Feature-level Analysis for Single Report
#'
#' @return Template string for feature-level analysis section
#' @noRd
generate_single_report_section_taxa <- function() {
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
                                               feature.dat.type = feature.dat.type)
```

```{r taxa-test-description, echo=FALSE, results='asis'}
cat(sprintf('In this analysis, we utilized the LinDA linear model to investigate potential differences in abundance. Specifically, we tested the effect of variables %s for different taxa, while adjusting for other covariates.\\n\\n', group.var))
```

```{r taxa-cladogram, message=FALSE, warning=FALSE, fig.align='center', fig.width=12, fig.height=12, result = result.output}
cladogram_plots <- generate_taxa_cladogram_single(
  data.obj = data.obj,
  test.list = taxa_test_results,
  group.var = group.var,
  feature.level = test.feature.level,
  feature.mt.method = feature.mt.method,
  cutoff = 1,
  color.group.level = test.feature.level[length(test.feature.level)-1]
)

cladogram_plots

```

```{r taxa-volcano , message = FALSE, warning = FALSE, fig.align = 'center', fig.width = 6.5, fig.height = 6.5, result = result.output}
volcano_plots <- generate_taxa_volcano_single(
                                  data.obj = data.obj,
                                  group.var = group.var,
                                  test.list = taxa_test_results,
                                  feature.sig.level = feature.sig.level,
                                  feature.mt.method = feature.mt.method
)
volcano_plots
```

```{r taxa-test-results-output, echo=FALSE, results='asis'}
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

# extract the directory from the full file path
output_dir <- dirname(output.file)
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

for(taxon_rank in names(taxa_test_results)) {

    comparisons <- names(taxa_test_results[[taxon_rank]])

    for(comparison in comparisons) {

        file_name <- paste0(filename_prefix, taxon_rank, '_', gsub(' ', '_', gsub('/', '_or_', comparison)), file_ext)

        # include the output directory in the file path
        file_path <- file.path(output_dir, file_name)

        write.csv(taxa_test_results[[taxon_rank]][[comparison]],
                  file = file_path,
                  row.names = FALSE)
    }
}

cat(sprintf('\\n\\nThe differential abundance test results for features have been saved in the directory: %s. Each taxa rank and its corresponding comparison have their own file named with the prefix: %s followed by the taxon rank, the comparison, and the file extension %s. Please refer to these files for more detailed data.', output_dir, filename_prefix, file_ext))
```

```{r extract_significant_taxa, echo=FALSE, results='hide', warning = FALSE}
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

## 4.2 Data visualization (significant features)

### 4.2.1 Significant features boxplot

```{r taxa-significant-boxplot-generation, message=FALSE, fig.align='center', fig.width = 8, fig.height = 4, results='asis', warning = FALSE}
if (length(significant_vars) != 0) {
  taxa_indiv_boxplot_results_sig_features <- generate_taxa_indiv_boxplot_single(
    data.obj = data.obj,
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
    pdf.hei = pdf.hei
  )
}
```

```{r taxa-boxplot-single-print, echo=FALSE, message=FALSE, result = result.output, fig.align='center', fig.width = 8, fig.height = 4, warning = FALSE}
if (length(significant_vars) != 0){
taxa_indiv_boxplot_results_sig_features
}
```

```{r taxa-indiv-boxplot-pdf-name, echo=FALSE, message=FALSE, results='asis', warning = FALSE}
output_dir <- dirname(output.file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

taxa_indiv_boxplot_results <- generate_taxa_indiv_boxplot_single(
                                   data.obj = data.obj,
                                   time.var = time.var,
                                   t.level = t.level,
                                   group.var = group.var,
                                   strata.var = strata.var,
                                   feature.level = test.feature.level,
                                   features.plot = NULL,
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

for (feature_level in names(taxa_indiv_boxplot_results)) {
  pdf_name <- paste0(
    'taxa_indiv_boxplot_single',
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

  cat(paste0('The boxplot results for all individual taxa or features at the ', feature_level, ' level can be found at: ', full_file_path, '. Please refer to this file for more detailed visualizations.\\n\\n'))
}

```
"
}
