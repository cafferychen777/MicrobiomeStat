#' Generate Taxa Test Pair
#'
#' This function takes as input a MicrobiomeStat data object and various specifications regarding the variables,
#' and returns an analysis of changes in taxa over time dplyr::across different groups. It also performs linear modelling on these changes,
#' utilizing specified group variables and covariates.
#'
#' @param data.obj A list object in MicrobiomeStat format.
#' @param subject.var A string that specifies the name of the subject variable column in the metadata.
#' @param time.var A string that specifies the name of the time variable column in the metadata. If not provided, it's NULL by default.
#' @param group.var A string that specifies the name of the grouping variable column in the metadata for linear modelling.
#' @param adj.vars A vector of strings that specify the names of additional variables to be used as covariates in the analysis.
#' @param feature.level A string specifying the taxonomic level at which to perform the analysis.
#' @param prev.filter A numeric value specifying the minimum prevalence filter for the taxa. It's set to 0 by default.
#' @param abund.filter A numeric value specifying the minimum abundance filter for the taxa. It's set to 0 by default.
#' @param feature.dat.type A string that specifies the type of the feature data. Options include: "count", "proportion". "count" is the default.
#' @param ... Additional parameters to be passed to the function.
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(GUniFrac)
#' library(ape)
#' library(philentropy)
#' library(MicrobiomeStat)
#'
#' data(peerj32.obj)
#' results <- generate_taxa_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   feature.level = "Family",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   feature.dat.type = "count"
#' )
#' }
#'
#' @return A named list where each element corresponds to a feature level. Each list element contains a dataframe with the calculated taxa changes, their corresponding p-values,
#' and other statistics obtained from the linear model.
#' @export
#' @name generate_taxa_test_pair
generate_taxa_test_pair <-
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

    # 将 OTU 表与分类表合并
    otu_tax <-
      cbind(otu_tab, tax_tab)

    # Create the formula
    fixed_effects <- paste(paste(adj.vars, collapse = " + "), group.var, time.var, sep = " + ")
    random_effects <- paste("(1|", subject.var, ")", sep = "")
    formula <- paste(fixed_effects, random_effects, sep = " + ")

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
        dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric) %>% column_to_rownames(feature.level)

      linda.obj <- linda(feature.dat = otu_tax_agg_numeric, meta.dat = meta_tab,
                         formula = paste("~",formula), feature.dat.type = "proportion")

      # Extract relevant information
      significant_taxa <- rownames(linda.obj$feature.dat.use)

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
