#' @title Longitudinal Taxa Test Generation
#' @description Generate longitudinal taxa test statistics based on the provided MicrobiomeStat data object and user-specified variables.
#'
#' @param data.obj A MicrobiomeStat formatted list object.
#' @param subject.var A string specifying the subject variable column in the metadata.
#' @param time.var A string specifying the time variable column in the metadata. Default is NULL.
#' @param t0.level A numeric specifying the baseline time level.
#' @param ts.levels A vector of numerics specifying the time points to be considered in the analysis.
#' @param group.var A string specifying the grouping variable column in the metadata.
#' @param adj.vars A string vector specifying the additional covariates in the metadata.
#' @param feature.level A string specifying the taxonomic level for analysis.
#' @param prev.filter A numeric specifying the minimum prevalence filter for taxa. Default is 0.
#' @param abund.filter A numeric specifying the minimum abundance filter for taxa. Default is 0.
#' @param feature.dat.type A string specifying the feature data type. Either "count" or "proportion". Default is "count".
#' @param ... Additional arguments to pass to the function.
#' @return A list of dataframes, each corresponding to a feature level. Each dataframe contains calculated taxa changes, p-values and other statistics from the linear model.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' generate_taxa_trend_test_long(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = "sex",
#'   feature.level = "Phylum",
#'   prev.filter = 0.001,
#'   abund.filter = 0.001,
#'   feature.dat.type = "count"
#' )
#'
#' data(ecam.obj)
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   feature.level = c("Phylum","Class"),
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   feature.dat.type = "proportion"
#' )
#' }
#' @export
generate_taxa_trend_test_long <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion"),
           ...) {
    # Extract data
    mStat_validate_data(data.obj)

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
          dplyr::mutate(., original = rownames(.))
        else
          .
      } %>%
      select(all_of(feature.level))

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    if (is.null(group.var)) {
      fixed_effects <- paste(paste(adj.vars, collapse = " + "), time.var, sep = " + ")
    } else {
      fixed_effects <- paste(paste(adj.vars, collapse = " + "), group.var, time.var, sep = " + ")
    }

    random_effects <- paste("(1|", subject.var, ")", sep = "")
    formula <- paste(fixed_effects, random_effects, sep = " + ")

    test.list <- lapply(feature.level, function(feature.level) {

      otu_tax <-
        cbind(otu_tab,
              tax_tab %>% select(all_of(feature.level)))

      otu_tax_filtered <- otu_tax %>%
        tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        dplyr::group_by_at(vars(!!sym(feature.level))) %>%
        dplyr::summarise(total_count = mean(value),
                         prevalence = sum(value > 0) / dplyr::n(),
                         , na.rm = TRUE) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count,-prevalence) %>%
        dplyr::left_join(otu_tax, by = feature.level)

      otu_tax_agg <- otu_tax_filtered %>%
        tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        dplyr::group_by_at(vars(sample,!!sym(feature.level))) %>%
        dplyr::summarise(value = sum(value)) %>%
        tidyr::spread(key = "sample", value = "value")

      otu_tax_agg_numeric <-
        otu_tax_agg %>%
        dplyr::mutate_at(vars(-!!sym(feature.level)), as.numeric) %>%
        dplyr::mutate(!!sym(feature.level) := tidyr::replace_na(!!sym(feature.level), "unclassified")) %>%
        column_to_rownames(feature.level)

      linda.obj <- linda(feature.dat = otu_tax_agg_numeric,
                         meta.dat = meta_tab,
                         formula = paste("~", formula),
                         feature.dat.type = "proportion",
                         prev.filter = prev.filter,
                         mean.abund.filter = abund.filter,
                         ...)

      # Extract relevant information
      significant_taxa <- rownames(linda.obj$feature.dat.use)

      if (is.null(group.var)){
        linda.obj$meta.dat.use$group <- "No Group"
        group.var <- "group"
      }

      # 计算每个分组的平均丰度
      prop_prev_data <-
        linda.obj$feature.dat.use %>% as.data.frame() %>% rownames_to_column(feature.level) %>%
        tidyr::gather(-!!sym(feature.level),
                      key = "sample",
                      value = "count") %>%
        dplyr::inner_join(linda.obj$meta.dat.use %>% rownames_to_column("sample"), by = "sample", relationship = "many-to-many") %>%
        dplyr::group_by(!!sym(group.var), !!sym(feature.level)) %>%
        dplyr::summarise(mean_abundance = mean(count),
                         prevalence = sum(count > 0) / dplyr::n())

      # Initialize results table
      results <- data.frame()

      # Iterate over each element in linda.obj$output
      for (output_element in names(linda.obj$output)) {
        current_output <- linda.obj$output[[output_element]]

        for (taxa in significant_taxa) {
          baseMean <- current_output[taxa, "baseMean"]
          log2FoldChange <- current_output[taxa, "log2FoldChange"]
          lfcSE <- current_output[taxa, "lfcSE"]
          stat <- current_output[taxa, "stat"]
          pvalue <- current_output[taxa, "pvalue"]
          padj <- current_output[taxa, "padj"]

          # 检查group.var是不是因子
          if (is_categorical(linda.obj$meta.dat.use[[group.var]])) {
            for (group in unique(linda.obj$meta.dat.use[[group.var]])) {
              group_data <-
                prop_prev_data[which(prop_prev_data[[feature.level]] == taxa &
                                       prop_prev_data[[group.var]] == group), ]
              mean_prop <- group_data$mean_abundance
              mean_prev <- group_data$prevalence

              results <- rbind(
                results,
                data.frame(
                  Variable = taxa,
                  Group = group,
                  Base.Mean = baseMean,
                  Log2.Fold.Change = log2FoldChange,
                  LFC.SE = lfcSE,
                  Stat = stat,
                  P.Value = pvalue,
                  Adjusted.P.Value = padj,
                  Mean.Abundance = mean_prop,
                  Mean.Prevalence = mean_prev,
                  Output.Element = output_element
                )
              )
            }
          } else {
            total_data <-
              prop_prev_data[which(prop_prev_data[[feature.level]] == taxa), ]
            mean_prop <- total_data$mean_abundance
            mean_prev <- total_data$prevalence

            results <- rbind(
              results,
              data.frame(
                Variable = taxa,
                Base.Mean = baseMean,
                Log2.Fold.Change = log2FoldChange,
                LFC.SE = lfcSE,
                Stat = stat,
                P.Value = pvalue,
                Adjusted.P.Value = padj,
                Mean.Abundance = mean_prop,
                Mean.Prevalence = mean_prev,
                Output.Element = output_element
              )
            )
          }
        }
      }

      return(results)

    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)

  }
