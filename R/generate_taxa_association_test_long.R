#' Longitudinal Taxa Association Test Generation
#'
#' This function performs association testing between taxa abundances and a grouping variable in longitudinal microbiome data. It discerns how taxa abundances differ between experimental or observational groups over time.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string that indicates the column name in the metadata which uniquely identifies each subject or sample.
#' @param group.var A character string specifying the grouping variable column in the metadata. This variable differentiates between different experimental or observational groups.
#' @param adj.vars A vector of character strings. Each string should denote a column name in the metadata that will serve as a covariate in the analysis. These variables might account for potential confounding influences. Default is NULL.
#' @param feature.level A character string indicating the taxonomic resolution for analysis (e.g., "Phylum", "Class"). This choice will determine the granularity of the analysis.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param feature.dat.type A character string, either "count" or "proportion", indicating the nature of the data in the `data.obj`. This helps the function to determine if normalization is required. Default is "count".
#' @param feature.sig.level A numeric threshold, usually between 0 and 1, for assessing the significance of individual taxa. Default is 0.1.
#' @param feature.mt.method A character string specifying the method employed for multiple testing correction (e.g., "fdr" for False Discovery Rate). Default is "fdr".
#' @param ... Additional arguments to cater to any specialized requirements. For now, these are placeholder and not used.
#' @details
#' Based on whether group.var and adj.vars are NULL, the formula tests:
#'
#' - When group.var is NULL and adj.vars is NOT NULL:
#'   - Tests adj.vars main effects only.
#'   - Adjusted for adj.vars but not group.var.
#'
#' - When group.var is NOT NULL and adj.vars is NOT NULL:
#'   - Tests adj.vars and group.var main effects.
#'   - Adjusted for adj.vars.
#'
#' - When group.var is NOT NULL and adj.vars is NULL:
#'   - Tests group.var main effects only.
#'   - Unadjusted analysis.
#'
#' - When both group.var and adj.vars are NULL:
#'   - Tests the intercept only.
#'   - Unadjusted analysis.
#'
#' The formula combines the appropriate terms based on which variables are NULL.
#' Subject variability is accounted for through random effects.
#'
#' When group.var and adj.vars are NULL, the intercept is tested without adjusting for any covariates.
#'
#' @return
#' A list of dataframes, with each dataframe representing a specific taxonomic level (as specified in `feature.level`). These dataframes contain essential statistics, including taxa changes, p-values, and other metrics derived from the linear model.
#'
#' @examples
#' \dontrun{
#' # Example 1
#' data("ecam.obj")
#' generate_taxa_association_test_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   group.var = "delivery",
#'   feature.level = c("Phylum","Class"),
#'   feature.dat.type = c("count")
#' )
#'
#' # Example 2
#' data("subset_T2D.obj")
#' generate_taxa_association_test_long(
#'   data.obj = subset_T2D.obj2,
#'   subject.var = "subject_id",
#'   feature.level = "Genus",
#'   group.var = "subject_race",
#'   feature.dat.type = c("count"),
#'   prev.filter = 0.1,
#'   abund.filter = 0.001
#' )
#' }
#' @export
generate_taxa_association_test_long <-
  function(data.obj,
           subject.var,
           group.var = NULL,
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

    message(
      "The trend test calculation relies on a numeric time variable.\n",
      "Please check that your time variable is coded as numeric.\n",
      "If the time variable is not numeric, it may cause issues in computing the test results.\n",
      "You can ensure the time variable is numeric by mutating it in the metadata."
    )

    meta_tab <-
      load_data_obj_metadata(data.obj) %>% select(all_of(c(
        group.var, adj.vars, subject.var
      )))

    # Function to generate the formula
    generate_formula <- function(group.var=NULL, adj.vars=NULL, time.var=NULL, subject.var="Subject") {

      # Initialize the fixed_effects and random_effects variables
      fixed_effects <- NULL
      random_effects <- NULL

      # Combine multiple adj.vars into a single string if they exist
      if (!is.null(adj.vars)) {
        adj.vars_str <- paste(adj.vars, collapse = " + ")
      } else {
        adj.vars_str <- NULL
      }

      # Case where time.var is NULL
      if (is.null(time.var)) {
        if (is.null(group.var)) {
          fixed_effects <- adj.vars_str
          if (is.null(fixed_effects)) {
            fixed_effects <- "1"  # Intercept-only model
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var)
          } else {
            fixed_effects <- group.var
          }
        }
        random_effects <- paste("(1 |", subject.var, ")")

        # Case where time.var is NOT NULL
      } else {
        if (is.null(group.var)) {
          fixed_effects <- paste(adj.vars_str, "+", time.var)
          if (is.null(adj.vars_str)) {
            fixed_effects <- time.var
          }
        } else {
          if (!is.null(adj.vars_str)) {
            fixed_effects <- paste(adj.vars_str, "+", group.var, "*", time.var)
          } else {
            fixed_effects <- paste(group.var, "*", time.var)
          }
        }
        random_effects <- paste("(1 +", time.var, "|", subject.var, ")")
      }

      # Generate the full formula
      formula <- paste(fixed_effects, random_effects, sep = " + ")
      return(formula)
    }

    formula <- generate_formula(group.var = group.var,
                                adj.vars = adj.vars,
                                time.var = NULL,
                                subject.var = subject.var)

    test.list <- lapply(feature.level, function(feature.level) {

      if (feature.dat.type == "count"){
        message(
          "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
        )
        data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
      }

      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- load_data_obj_count(data.obj)
      }

      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter)

      linda.obj <- linda(
        feature.dat = otu_tax_agg,
        meta.dat = meta_tab,
        formula = paste("~", formula),
        feature.dat.type = "proportion",
        prev.filter = prev.filter,
        mean.abund.filter = abund.filter,
        ...
      )

      # Extract relevant information
      significant_taxa <- rownames(linda.obj$feature.dat.use)

      if (is.null(group.var)) {
        linda.obj$meta.dat.use$group <- "No Group"
        group.var <- "group"
      }

      # 计算每个分组的平均丰度
      prop_prev_data <-
        linda.obj$feature.dat.use %>% as.data.frame() %>% rownames_to_column(feature.level) %>%
        tidyr::gather(-!!sym(feature.level),
                      key = "sample",
                      value = "count") %>%
        dplyr::inner_join(
          linda.obj$meta.dat.use %>% rownames_to_column("sample"),
          by = "sample",
          relationship = "many-to-many"
        ) %>%
        dplyr::group_by(!!sym(group.var),!!sym(feature.level)) %>%
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
                                       prop_prev_data[[group.var]] == group),]
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
              prop_prev_data[which(prop_prev_data[[feature.level]] == taxa),]
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

    plot.list <-
      generate_taxa_association_volcano_long(
        data.obj = data.obj,
        group.var = group.var,
        test.list = test.list,
        feature.sig.level = feature.sig.level,
        feature.mt.method = feature.mt.method
      )
    print(plot.list)

    return(test.list)
  }

#' Generate Volcano Plots for Taxa Association Test
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var The grouping variable tested, found in metadata
#' @param test.list The list of test results returned by generate_taxa_trend_test_long
#' @param feature.sig.level The significance level cutoff for highlighting taxa
#' @param feature.mt.method Multiple testing correction method, "fdr" or "none"
#'
#' @return A list of ggplot objects of volcano plots for each taxonomic level
#'
#' @examples
#' # Generate test results
#' # test_list <- generate_taxa_association_test_long(...)

#' # Create volcano plots
#' # volcano_plots <- generate_taxa_association_volcano_long(
#' #   data.obj, group.var, test_list,
#' #   feature.sig.level = 0.05, feature.mt.method = "fdr"
#' # )
#'
#' @importFrom dplyr distinct pull
#' @export
generate_taxa_association_volcano_long <-
  function(data.obj,
           group.var = NULL,
           test.list,
           feature.sig.level = 0.1,
           feature.mt.method = "fdr") {
    meta_tab <- load_data_obj_metadata(data.obj) %>%
      dplyr::select(all_of(c(group.var))) %>% rownames_to_column("sample")

    # Define the custom color palette
    color_palette <- c("#2A9D8F", "#F9F871", "#F4A261", "#FF6347")

    feature.level <- names(test.list)

    # 使用条件表达式设置要使用的p值变量
    p_val_var <-
      ifelse(feature.mt.method == "fdr",
             "Adjusted.P.Value",
             "P.Value")

    plot.list <- lapply(feature.level, function(feature.level) {
      sub_test.list <- test.list[[feature.level]]

        group_level <-
          meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

        reference_level <- group_level[1]

        test.result <- sub_test.list %>%
          dplyr::filter(grepl(
            paste0('^', group.var, '.*', NULL, '$'),
            Output.Element
          ))

        sub_plot.list <-
          lapply(group_level[-1], function(group_level) {

              sub_test.result <- test.result %>%
                dplyr::filter(Output.Element == paste0(group.var, group_level))

            # Find max absolute log2FoldChange for symmetric x-axis
            max_abs_log2FC <-
              max(abs(sub_test.result$Log2.Fold.Change), na.rm = TRUE)

            p <-
              ggplot(sub_test.result,
                     aes(
                       x = Log2.Fold.Change,
                       y = -log10(get(p_val_var)),
                       color = -log10(get(p_val_var))
                     )) +
              geom_point(aes(shape = get(p_val_var) < feature.sig.level), size = 7) +
              geom_vline(
                aes(xintercept = 0),
                linetype = "dashed",
                linewidth = 1.5,
                color = "grey"
              ) +
              geom_hline(
                aes(yintercept = -log10(feature.sig.level)),
                linetype = "dashed",
                linewidth = 1.5,
                color = "grey"
              ) +
              geom_text(
                aes(label = ifelse(
                  get(p_val_var) < feature.sig.level,
                  as.character(Variable),
                  ''
                )),
                vjust = -0.5,
                hjust = 0.5,
                size = 3.5
              ) +
              scale_shape_manual(values = c(16, 17)) +
              labs(
                title = paste(group_level, "vs", reference_level, "(Reference)"),
                x = "Log2 Fold Change",
                y = "-log10(p-value)",
                shape = "Significant",
                color = "-log10(p-value)"
              ) +
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
          })

        names(sub_plot.list) <-
          paste(group_level[-1], "vs", reference_level, "(Reference)")

      return(sub_plot.list)
    })


    names(plot.list) <- feature.level
    return(plot.list)
  }
