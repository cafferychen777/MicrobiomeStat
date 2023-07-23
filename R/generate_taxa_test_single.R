#' Conduct Differential Abundance Testing Using ZicoSeq Method in MicrobiomeStat Package
#'
#' This function applies a differential abundance analysis using ZicoSeq on a data set. The function filters taxa based on prevalence and abundance, then it aggregates and applies the ZicoSeq method. Finally, it creates a report of significant taxa with relevant statistics.
#'
#' @param data.obj A list containing the metadata and feature data.
#' @param time.var A string representing the time variable. Default is NULL.
#' @param t.level A string representing the time level. Default is NULL.
#' @param group.var A string indicating the variable used for grouping.
#' @param adj.vars A string indicating the adjustment variables.
#' @param prev.filter A numeric value indicating the prevalence filter threshold. Default is 0.
#' @param abund.filter A numeric value indicating the abundance filter threshold. Default is 0.
#' @param feature.level A character vector indicating the feature level(s).
#' @param feature.dat.type A character string representing the type of feature data. Choices are "count", "proportion", "other". Default is "count".
#' @param ... Additional arguments to be passed to the ZicoSeq function.
#'
#' @return A list of tibble(s) containing information about significant taxa, including R.Squared, F.Statistic, Estimate, P.Value, Adjusted.P.Value, Mean.Proportion, Mean.Prevalence, SD.Abundance and SD.Prevalence.
#'
#' @export
#'
#' @examples
#'
#' data(peerj32.obj)
#' da_report <- generate_taxa_test_single(
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
#'     verbose = TRUE,
#'     return.feature.dat = T
#' )
#'
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
      load_data_obj_count(mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm)
  } else{
    otu_tab <- load_data_obj_count(data.obj)
  }

  tax_tab <- load_data_obj_taxonomy(data.obj) %>%
    as.data.frame() %>%
    {
      if ("original" %in% feature.level)
        mutate(., original = rownames(.))
      else
        .
    } %>%
    select(all_of(feature.level))

  meta_tab <-
    load_data_obj_metadata(data.obj) %>% select(all_of(c(
      time.var, group.var, adj.vars, subject.var
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
      group_by_at(vars(!!sym(feature.level))) %>%
      dplyr::summarise(total_count = mean(count),
                       prevalence = sum(count > 0) / dplyr::n()) %>%
      filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
      select(-total_count, -prevalence) %>%
      dplyr::left_join(otu_tax, by = feature.level)

    # 聚合 OTU 表
    otu_tax_agg <- otu_tax_filtered %>%
      tidyr::gather(key = "sample", value = "count", -one_of(colnames(tax_tab))) %>%
      group_by_at(vars(sample, !!sym(feature.level))) %>%
      dplyr::summarise(count = sum(count)) %>%
      spread(key = "sample", value = "count")

    # 转换计数为数值类型
    otu_tax_agg_numeric <-
      mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

    otu_tax_agg_numeric <- otu_tax_agg_numeric %>%
      mutate(!!sym(feature.level) := replace_na(!!sym(feature.level), "Unclassified")) %>%
      column_to_rownames(feature.level) %>%
      as.matrix()

    # Remove rows that are all zeros
    otu_tax_agg_numeric <-
      otu_tax_agg_numeric[rowSums(otu_tax_agg_numeric != 0) > 0,]

    # Run ZicoSeq
    zico.obj <- ZicoSeq(
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
      group_by_at(vars(!!sym(feature.level),!!sym(group.var))) %>%
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
            prop_prev_data[which(prop_prev_data[[feature.level]] == taxa &
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
