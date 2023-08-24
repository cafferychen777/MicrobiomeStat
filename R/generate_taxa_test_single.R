#' Conduct Differential Abundance Testing Using ZicoSeq Method in MicrobiomeStat Package
#'
#' This function applies a differential abundance analysis using ZicoSeq on a data set. The function filters taxa based on prevalence and abundance, then it aggregates and applies the ZicoSeq method. Finally, it creates a report of significant taxa with relevant statistics.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param time.var Character string specifying the column name in metadata containing time variable.
#'                Used to subset data to a single timepoint if provided. Default NULL does not subset.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var Character string specifying the column name in metadata containing grouping
#'                 categories. This will be used as the predictor in differential abundance testing.
#' @param adj.vars Character vector specifying column names in metadata containing covariates.
#'                These will be used for adjustment in differential abundance testing.
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
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param ... Additional arguments to be passed to the ZicoSeq function.
#'
#' @return A list of tibble(s) containing information about significant taxa, including R.Squared, F.Statistic, Estimate, P.Value, Adjusted.P.Value, Mean.Proportion, Mean.Prevalence, SD.Abundance and SD.Prevalence.
#'
#' @details
#' This function performs differential abundance analysis using ZicoSeq method:
#'
#' 1. Subset data to specific time point if time variable and level provided.
#'
#' 2. Extract OTU table, taxonomy table, and sample metadata.
#'
#' 3. Filter OTUs by prevalence and abundance thresholds if counting data.
#'
#' 4. Aggregate OTU table to taxonomic levels specified by \code{feature.level}.
#'
#' 5. Run ZicoSeq on aggregated table, with grouping and adjustment variables.
#'
#' 6. Extract and compile statistics for significant taxa into results tables,
#' including R-squared, F statistic, coefficients, p-values, mean proportions, etc.
#'
#' 7. Return list of tables, with each element corresponding to one taxonomic level.
#'
#' In summary, it applies preprocessing steps, runs ZicoSeq differential abundance testing,
#' and compiles informative results tables for significant taxa. Adjustment for covariates is
#' supported. Customizable taxonomic aggregation allows testing at different levels.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' generate_taxa_test_single(
#'     data.obj = peerj32.obj,
#'     time.var = "time",
#'     t.level = "1",
#'     group.var = "group",
#'     adj.vars = "sex",
#'     feature.dat.type = "count",
#'     feature.level = "Genus",
#'     prev.filter = 0,
#'     abund.filter = 0,
#'     is.winsor = TRUE,
#'     outlier.pct = 0.001,
#'     winsor.end = 'top',
#'     is.post.sample = TRUE,
#'     post.sample.no = 25,
#'     list(function (x) x^0.5, function (x) x^0.25),
#'     stats.combine.func = max,
#'     perm.no = 99,
#'     strata = NULL,
#'     ref.pct = 0.5,
#'     stage.no = 6,
#'     excl.pct = 0.2,
#'     is.fwer = TRUE,
#'     verbose = TRUE
#' )
#' }
#' @export
generate_taxa_test_single <- function(data.obj,
                                      time.var = NULL,
                                      t.level = NULL,
                                      group.var,
                                      adj.vars,
                                      prev.filter = 0,
                                      abund.filter = 0,
                                      feature.level,
                                      feature.dat.type = c("count", "proportion", "other"),
                                      ...) {
  # Extract data
  mStat_validate_data(data.obj)

  if (!is.null(time.var)) {
    if (!is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <- mStat_subset_data(data.obj, condition = condition)
    }
  }

  if (feature.dat.type == "count") {
    message(
      "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
    )
    otu_tab <-
      load_data_obj_count(mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm)
  } else{
    otu_tab <- load_data_obj_count(data.obj)
  }

  tax_tab <- load_data_obj_taxonomy(data.obj) %>%
    as.data.frame() %>%
    {
      if ("original" %in% feature.level)
        dplyr::mutate(., original = rownames(.))
      else
        .
    } %>%
    select(all_of(feature.level))

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% select(all_of(c(
      time.var, group.var, adj.vars
    )))

  # 将 OTU 表与分类表合并
  otu_tax <-
    cbind(otu_tab, tax_tab)

  if (feature.dat.type == "other") {
    prev.filter <- 0
    abund.filter <- 0
  }

  test.list <- lapply(feature.level, function(feature.level) {
    # Filter taxa based on prevalence and abundance
    otu_tax_filtered <- otu_tax %>%
      tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(!!sym(feature.level))) %>%
      dplyr::summarise(total_count = mean(count),
                       prevalence = sum(count > 0) / dplyr::n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count, -prevalence) %>%
      dplyr::left_join(otu_tax, by = feature.level)

    # 聚合 OTU 表
    otu_tax_agg <- otu_tax_filtered %>%
      tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      dplyr::group_by_at(vars(sample, !!sym(feature.level))) %>%
      dplyr::summarise(count = sum(count)) %>%
      tidyr::spread(key = "sample", value = "count")

    # 转换计数为数值类型
    otu_tax_agg_numeric <-
      dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    otu_tax_agg_numeric <- otu_tax_agg_numeric %>%
      dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(feature.level) %>%
      as.matrix()

    # Remove rows that are all zeros
    otu_tax_agg_numeric <-
      otu_tax_agg_numeric[rowSums(otu_tax_agg_numeric != 0) > 0,]

    # Run ZicoSeq
    zico.obj <- GUniFrac::ZicoSeq(
      meta.dat = meta_tab,
      feature.dat = otu_tax_agg_numeric,
      grp.name = group.var,
      adj.name = adj.vars,
      feature.dat.type = "other",
      prev.filter = prev.filter,
      max.abund.filter = abund.filter,
      ...
    )

    # Extract relevant information
    significant_taxa <- names(which(zico.obj$p.adj.fdr <= 1))

    # If no significant taxa were found, stop the function
    if (length(significant_taxa) == 0) {
      stop("No significant taxa were found. Stopping the function.")
    }

    # Initialize results table
    results <- data.frame()

    prop_prev_data <- tidyr::gather(
      otu_tax_agg_numeric %>%
        as.data.frame() %>% rownames_to_column(feature.level),
      key = "sample",
      value = "count",-feature.level
    )  %>%
      dplyr::left_join(data.obj$meta.dat %>% select(all_of(c(
        group.var, adj.vars
      ))) %>%
        rownames_to_column("sample"),
      by = "sample") %>%
      dplyr::group_by_at(vars(!!sym(feature.level),!!sym(group.var))) %>%
      dplyr::summarise(
        mean_proportion = mean(count),
        sdev_count = sd(count),
        prevalence = sum(count > 0) / dplyr::n(),
        sdev_prevalence = sd(ifelse(count > 0, 1, 0))
      )

    is_categorical <- function(x) {
      if (is.factor(x) || is.character(x) || is.logical(x)) {
        return(TRUE)
      } else if (is.numeric(x) &&
                 length(unique(x)) < 10) {
        # 这里的10你可以根据需要调整
        return(TRUE)
      } else {
        return(FALSE)
      }
    }

    for (taxa in significant_taxa) {
      R.Squared <- zico.obj$R2[taxa, 1]
      F.Statistic <- zico.obj$F0[taxa, 1]

      Estimate <-
        zico.obj$coef.list[[1]][startsWith(rownames(zico.obj$coef.list[[1]]), group.var), taxa]  # 选择 group.var 的估计值

      P.Value <- zico.obj$p.raw[taxa]
      Adjusted.P.Value <- zico.obj$p.adj.fdr[taxa]

      # 检查group.var是不是因子
      if (is_categorical(data.obj$meta.dat[[group.var]])) {
        for (group in unique(data.obj$meta.dat[[group.var]])) {
          group_data <-
            prop_prev_data[Matrix::which(prop_prev_data[[feature.level]] == taxa &
                                   prop_prev_data[[group.var]] == group),]
          mean_prop <- group_data$mean_proportion
          mean_prev <- group_data$prevalence
          sd_abundance <- group_data$sdev_count
          sd_prevalence <- group_data$sdev_prevalence

          results <- rbind(
            results,
            data.frame(
              Variable = taxa,
              Group = group,
              R.Squared = R.Squared,
              F.Statistic = F.Statistic,
              Estimate = toString(Estimate),
              P.Value = P.Value,
              Adjusted.P.Value = Adjusted.P.Value,
              Mean.Proportion = mean_prop,
              Mean.Prevalence = mean_prev,
              SD.Abundance = sd_abundance,
              SD.Prevalence = sd_prevalence
            )
          )
        }
      } else {
        total_data <-
          prop_prev_data[which(prop_prev_data[[feature.level]] == taxa),]
        mean_prop <- total_data$mean_proportion
        mean_prev <- total_data$prevalence
        sd_abundance <- total_data$sdev_count
        sd_prevalence <- total_data$sdev_prevalence

        results <- rbind(
          results,
          data.frame(
            Variable = taxa,
            R.Squared = R.Squared,
            F.Statistic = F.Statistic,
            Estimate = toString(Estimate),
            P.Value = P.Value,
            Adjusted.P.Value = Adjusted.P.Value,
            Mean.Proportion = mean_prop,
            Mean.Prevalence = mean_prev,
            SD.Abundance = sd_abundance,
            SD.Prevalence = sd_prevalence
          )
        )
      }
    }
    return(as_tibble(results))
  })

  # Assign names to the elements of test.list
  names(test.list) <- feature.level

  # Return the results table
  return(test.list)
}
