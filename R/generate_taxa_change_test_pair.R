#' Compute taxa changes and analyze differential abundance
#'
#' This function calculates taxa abundance changes between two time points and performs differential abundance analysis between groups using linear models or ANOVA.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var The name of the subject variable column in the metadata.
#' @param time.var The name of the time variable column in the metadata (optional).
#' @param group.var The name of the grouping variable column for linear modeling in the metadata.
#' @param adj.vars Names of additional variables to be used as covariates in the analysis.
#' @param change.base The baseline time point for detecting changes in taxa. If NULL, the first unique value from the time.var column will be used (optional).
#' @param feature.change.func Specifies the method or function used to compute the change between two time points. Options include:
#'
#' - "absolute change" (default): Computes the absolute difference between the values at the two time points (`value_time_2` and `value_time_1`).
#'
#' - "log fold change": Computes the log2 fold change between the two time points. For zero values, imputation is performed using half of the minimum nonzero value for each feature level at the respective time point before taking the logarithm.
#'
#' - "relative change": Computes the relative change as `(value_time_2 - value_time_1) / (value_time_2 + value_time_1)`. If both time points have a value of 0, the change is defined as 0.
#'
#' - A custom function: If a user-defined function is provided, it should take two numeric vectors as input corresponding to the values at the two time points (`value_time_1` and `value_time_2`) and return a numeric vector of the computed change. This custom function will be applied directly to calculate the difference.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#'
#' @examples
#' \dontrun{
#' data(peerj32.obj)
#' generate_taxa_change_test_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   feature.dat.type = "count"
#' )
#'
#' data(subset_pairs.obj)
#' generate_taxa_change_test_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = c("Genus"),
#'   prev.filter = 0.1,
#'   abund.filter = 1e-4,
#'   feature.dat.type = "count"
#' )
#' }
#'
#' @return A named list where each element corresponds to a feature level and contains a dataframe with the calculated taxa changes, their corresponding p-values, and other statistics from the linear model.
#' @export
#' @name generate_taxa_change_test_pair
generate_taxa_change_test_pair <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           group.var = NULL,
           adj.vars = NULL,
           change.base,
           feature.change.func = "relative change",
           feature.level,
           prev.filter = 0.1,
           abund.filter = 1e-4,
           feature.dat.type = c("count", "proportion", "other")) {
    # Extract data
    mStat_validate_data(data.obj)

    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    if (!is.null(adj.vars)){
      # Use the modified mStat_identify_time_varying_vars function
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # Check if there are any time-varying variables
      if (length(time_varying_info$time_varying_vars) > 0) {
        stop("Feature-level analysis does not yet support adjustment for time-varying variables. Found time-varying variables: ",
             paste(time_varying_info$time_varying_vars, collapse = ", "),
             ". Future versions will support this feature.")
      }
    }

    group_level <-
      meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

    reference_level <- group_level[1]

    # Create a formula including the group variable and adjustment variables (if any)
    formula_str <- paste("value ~", group.var)

    if (!is.null(adj.vars)) {
      formula_str <-
        paste(formula_str, "+", paste(adj.vars, collapse = " + "))
    }

    formula <- as.formula(formula_str)

    if (is.null(change.base)) {
      change.base <- unique(meta_tab %>% select(all_of(c(time.var))))[1,]
      message(
        "The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'meta.dat' data frame: ",
        change.base
      )
    }

    if (feature.dat.type == "other") {
      prev.filter <- 0
      abund.filter <- 0
    }

    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      data.obj <-
        mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    test.list <- lapply(feature.level, function(feature.level) {
      if (is.null(data.obj$feature.agg.list[[feature.level]]) &
          feature.level != "original") {
        data.obj <-
          mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      if (feature.level != "original") {
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      otu_tax_agg_filter <- otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        rownames_to_column(feature.level)

      # Convert the otu_tax_agg_filter from wide format to long format.
      otu_tax_long <- otu_tax_agg_filter %>%
        tidyr::gather(key = "sample", value = "value", -feature.level)

      # Connect otu_tax_long and meta_tab by the sample column.
      merged_data <- otu_tax_long %>%
        dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      # Based on the time column
      grouped_data <- merged_data %>%
        dplyr::group_by(!!sym(time.var))

      change.after <-
        unique(grouped_data %>% select(all_of(c(time.var))))[unique(grouped_data %>% select(all_of(c(time.var)))) != change.base]

      # Split into a list, with each time value having an independent tibble.
      split_data <-
        split(merged_data, f = grouped_data %>% select(all_of(c(time.var))))

      # Extract the first and second tables from split_data.
      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # Join these two tables together to calculate the difference
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      if (is.function(feature.change.func)) {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = feature.change.func(value_time_2, value_time_1))
      } else if (feature.change.func == "log fold change") {
        half_nonzero_min_time_2 <- combined_data %>%
          filter(value_time_2 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(value_time_2) / 2,
                           .groups = "drop")
        half_nonzero_min_time_1 <- combined_data %>%
          filter(value_time_1 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(value_time_1) / 2,
                           .groups = "drop")

        combined_data <-
          dplyr::left_join(
            combined_data,
            half_nonzero_min_time_2,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data <-
          dplyr::left_join(
            combined_data,
            half_nonzero_min_time_1,
            by = feature.level,
            suffix = c("_time_1", "_time_2")
          )
        combined_data$value_time_2[combined_data$value_time_2 == 0] <-
          combined_data$half_nonzero_min_time_2[combined_data$value_time_2 == 0]
        combined_data$value_time_1[combined_data$value_time_1 == 0] <-
          combined_data$half_nonzero_min_time_1[combined_data$value_time_1 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message(
          "Imputation was performed using half the minimum nonzero proportion for each taxon at different time points."
        )

        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (feature.change.func == "relative change") {
        combined_data <- combined_data %>%
          dplyr::mutate(value_diff = dplyr::case_when(
            value_time_2 == 0 & value_time_1 == 0 ~ 0,
            TRUE ~ (value_time_2 - value_time_1) / (value_time_2 + value_time_1)
          ))
      } else if (feature.change.func == "absolute change"){
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      } else {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      }

      message("Note: For repeated measurements of the same subject at the same time point, the average will be taken.")

      value_diff_matrix <- combined_data %>%
        select(feature.level, !!sym(subject.var), value_diff) %>%
        dplyr::group_by(!!sym(feature.level), !!sym(subject.var)) %>%
        dplyr::summarise(value_diff = mean(value_diff, na.rm = TRUE)) %>%
        tidyr::spread(key = !!sym(subject.var), value = value_diff) %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      unique_meta_tab <- meta_tab %>%
        filter(!!sym(subject.var) %in% colnames(value_diff_matrix)) %>%
        select(all_of(c(subject.var, group.var, adj.vars))) %>%
        dplyr::distinct(!!sym(subject.var), .keep_all = TRUE) %>% as_tibble()

      cols_order <- colnames(na.omit(value_diff_matrix))

      unique_meta_tab <-
        unique_meta_tab %>% column_to_rownames(subject.var)

      sorted_unique_meta_tab <- unique_meta_tab %>%
        dplyr::slice(match(cols_order, rownames(unique_meta_tab)))

      prop_prev_data <-
        otu_tax_agg %>%
        as.matrix() %>%
        as.table() %>%
        as.data.frame() %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(avg_abundance = mean(Freq),
                         prevalence = sum(Freq > 0) / dplyr::n()) %>% column_to_rownames("Var1") %>%
        rownames_to_column(feature.level)

      value_diff_long <- value_diff_matrix %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::gather(key = !!sym(subject.var), value = "value", -feature.level)

      sub_test.list <-
        lapply(value_diff_long %>% select(all_of(feature.level)) %>% pull() %>% unique(), function(taxon) {
          test_df <- value_diff_long %>%
            dplyr::filter(!!sym(feature.level) == taxon) %>%
            dplyr::left_join(sorted_unique_meta_tab %>%
                               as.data.frame() %>%
                               rownames_to_column(subject.var),
                             by = subject.var)

          # Run the linear model
          test_result <- lm(formula, data = test_df)

          coef.tab <- extract_coef(test_result)

          # Run ANOVA on the model if group.var is multi-categorical
          if (length(unique(test_df[[group.var]])) > 2) {
            anova <- anova(test_result)
            anova.tab <- anova %>%
              as.data.frame() %>%
              rownames_to_column("Term") %>%
              select(
                Term,
                Statistic = `F value`,
                P.Value = `Pr(>F)`
              ) %>%
              as_tibble() %>%
              dplyr::mutate(Estimate = NA, Std.Error = NA)

            # Reorder the columns to match coef.tab
            anova.tab <- anova.tab %>%
              select(
                Term,
                Estimate,
                Std.Error,
                Statistic,
                P.Value
              )

            coef.tab <-
              rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
          }
          return(as_tibble(coef.tab))
        })

      # Assign names to the elements of test.list
      names(sub_test.list) <- value_diff_long %>%
        select(all_of(feature.level)) %>%
        pull() %>%
        unique()

      unique_terms <-
        grep(paste0("^", group.var, "$|^", group.var, ".*"),
             unique(unlist(
               lapply(sub_test.list, function(df)
                 unique(df$Term))
             )),
             value = TRUE)

      result_list <- lapply(unique_terms, function(term) {
        do.call(rbind, lapply(sub_test.list, function(df) {
          df %>% dplyr::filter(Term == term)
        })) %>%
          dplyr::mutate(!!sym(feature.level) := names(sub_test.list)) %>%
          dplyr::left_join(prop_prev_data, by = feature.level) %>%
          dplyr::mutate(Adjusted.P.Value = p.adjust(P.Value, method = "fdr")) %>%
          dplyr::select(all_of(
            c(
              feature.level,
              "Estimate",
              "Std.Error",
              "P.Value",
              "Adjusted.P.Value",
              "avg_abundance",
              "prevalence"
            )
          )) %>%
          dplyr::rename(
            Coefficient = Estimate,
            SE = Std.Error,
            Variable = feature.level,
            Mean.Abundance = avg_abundance,
            Prevalence = prevalence
          )
      })

      names(result_list) <- unique_terms

      new_names <- sapply(names(result_list), function(name) {

        if (grepl(paste0("^", group.var), name) &&
            !grepl(paste0("^", group.var, "$"), name)) {
          sub_name <- sub(paste0(group.var), "", name)
          return(paste(sub_name, "vs", reference_level, "(Reference)"))
        }
        return(name)
      })

      names(result_list) <- new_names

      return(result_list)

    })

    # Assign names to the elements of test.list
    names(test.list) <- feature.level

    return(test.list)

  }
