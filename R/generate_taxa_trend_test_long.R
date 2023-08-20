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
#'
#' data("ecam.obj")
#' generate_taxa_trend_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month_num",
#'   group.var = "diet",
#'   adj.vars = NULL,
#'   feature.level = c("Phylum","Class"),
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   feature.dat.type = "proportion"
#' )
#'
#' data("subset_T2D.obj")
#' a <- generate_taxa_trend_test_long(
#'   data.obj = subset_T2D.obj,
#'   subject.var = "subject_id",
#'   time.var = "visit_number",
#'   group.var = "subject_gender",
#'   adj.vars = NULL,
#'   feature.level = "Phylum",
#'   prev.filter = 1e-4,
#'   abund.filter = 1e-5,
#'   feature.dat.type = "count",
#'   feature.sig.level = 0.1,
#'   feature.mt.method = "fdr"
#' )
#'
#' }
#' @export
generate_taxa_trend_test_long <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var,
           adj.vars = NULL,
           feature.level,
           prev.filter = 0,
           abund.filter = 0,
           feature.dat.type = c("count", "proportion"),
           feature.sig.level = 0.1,
           feature.mt.method = "fdr",
           ...) {
    # Extract data
    mStat_validate_data(data.obj)

    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      otu_tab <-
        load_data_obj_count(mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm)
    } else{
      otu_tab <- load_data_obj_count(data.obj)
    }

    message(
      "The volatility calculation in generate_taxa_volatility_test_long relies on a numeric time variable.\n",
      "Please check that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the results of the volatility test.\n",
      "You can ensure the time variable is numeric by mutating it in the metadata."
    )

    data.obj$meta.dat <-
      data.obj$meta.dat %>% dplyr::mutate(!!sym(time.var) := as.numeric(!!sym(time.var)))

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

    # Adjusting the fixed effects for interaction and checking adj.vars
    if (is.null(group.var)) {
      if (!is.null(adj.vars)) {
        fixed_effects <- paste(adj.vars, "+", time.var)
      } else {
        fixed_effects <- time.var
      }
    } else {
      if (!is.null(adj.vars)) {
        fixed_effects <- paste(adj.vars, "+", group.var, "*", time.var)
      } else {
        fixed_effects <- paste(group.var, "*", time.var)
      }
    }

    random_effects <- paste("(1 +", time.var, "|", subject.var, ")")
    formula <- paste(fixed_effects, random_effects, sep = " + ")

    test.list <- lapply(feature.level, function(feature.level) {

      otu_tax <-
        cbind(otu_tab,
              tax_tab %>% select(all_of(feature.level)))

      otu_tax_filtered <- otu_tax %>%
        tidyr::gather(key = "sample", value = "value",-one_of(feature.level)) %>%
        dplyr::group_by_at(vars(!!sym(feature.level))) %>%
        dplyr::summarise(total_count = mean(value),
                         prevalence = sum(value > 0) / dplyr::n()) %>%
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
                         mean.abund.filter = abund.filter)

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

    # Filtering and generating volcano plots for the selected test results
    plots.list <- lapply(test.list, function(test.result) {
      filtered_result <- dplyr::filter(test.result, grepl(paste0("^", group.var, ".*", time.var, "$"), Output.Element))

      levels <- unique(filtered_result$Output.Element)

      sub_plots.list <- lapply(levels, function(level){
        sub_filtered_result <- filtered_result %>% filter(Output.Element == level)
        if (nrow(filtered_result) > 0) { # Check if there are any rows after filtering
          return(generate_taxa_trend_volcano_long(sub_filtered_result))
        } else {
          return(NULL) # If no rows after filtering, return NULL
        }
      })

    })

    print(plots.list)

    return(test.list)

  }


generate_taxa_trend_volcano_long <- function(data.obj, group.var, time.var, test.list, feature.sig.level = 0.1, feature.mt.method = c("fdr","none")) {

  meta_tab <- load_data_obj_metadata(data.obj) %>%
    dplyr::select(all_of(c(
      group.var
    ))) %>% rownames_to_column("sample")

  feature.level <- names(test.list)

  group_level <- meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% unique()

  reference_level <- group_level[1]

  # Find max absolute log2FoldChange for symmetric x-axis
  max_abs_log2FC <- max(abs(test.result$Log2.Fold.Change), na.rm = TRUE)

  # Define the custom color palette
  color_palette <- c("#2A9D8F", "#F9F871", "#F4A261", "#FF6347")

  p <- ggplot(test.result, aes(x = Log2.Fold.Change, y = -log10(P.Value), color = -log10(P.Value))) +
    geom_point(aes(shape = Adjusted.P.Value < 0.05), size = 7) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.5, color = "grey") +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed", size = 1.5, color = "grey") +
    geom_text(aes(label = ifelse(Adjusted.P.Value < feature.sig.level, as.character(Variable), '')), vjust = -0.5, hjust = 0.5, size = 3.5) +
    scale_shape_manual(values = c(16, 17)) +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(p-value)", shape = "Significant", color = "-log10(p-value)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_line(color = "grey", linetype = "dashed"),
      panel.grid.minor = element_line(color = "grey", linetype = "dotted"),
      legend.position = "bottom"
    ) +
    scale_color_gradientn(colors = color_palette) +
    coord_cartesian(xlim = c(-max_abs_log2FC, max_abs_log2FC))



  return(p)
}



