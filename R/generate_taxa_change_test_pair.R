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
#' @param winsor.qt A numeric value between 0 and 1, specifying the quantile for winsorization (default: 0.97).
#'   Winsorization is a data preprocessing method used to limit extreme values or outliers in the data.
#'   The `winsor.qt` parameter determines the upper and lower quantiles for winsorization.
#'   For example, if `winsor.qt` is set to 0.97, the lower quantile will be (1 - 0.97) / 2 = 0.015,
#'   and the upper quantile will be 1 - (1 - 0.97) / 2 = 0.985.
#'   Values below the lower quantile will be replaced with the lower quantile,
#'   and values above the upper quantile will be replaced with the upper quantile.
#'   This helps to reduce the impact of extreme values or outliers on subsequent analyses.
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
           feature.dat.type = c("count", "proportion", "other"),
           winsor.qt = 0.97) {
    # Validate the input data object
    mStat_validate_data(data.obj)

    # Match the feature data type argument
    feature.dat.type <- match.arg(feature.dat.type)

    # Extract relevant columns from the metadata
    meta_tab <-
      data.obj$meta.dat %>% select(all_of(c(
        time.var, group.var, adj.vars, subject.var
      )))

    # Check for time-varying covariates
    if (!is.null(adj.vars)){
      # Identify time-varying variables in the metadata
      time_varying_info <- mStat_identify_time_varying_vars(meta.dat = meta_tab, adj.vars = adj.vars, subject.var = subject.var)

      # If time-varying variables are found, stop the analysis
      if (length(time_varying_info$time_varying_vars) > 0) {
        stop("Feature-level analysis does not yet support adjustment for time-varying variables. Found time-varying variables: ",
             paste(time_varying_info$time_varying_vars, collapse = ", "),
             ". Future versions will support this feature.")
      }
    }

    # Extract levels of the grouping variable
    group_level <-
      meta_tab %>% select(all_of(c(group.var))) %>% pull() %>% as.factor() %>% levels

    # Set the reference level for the grouping variable
    reference_level <- group_level[1]

    # Construct the formula for the linear model
    formula_str <- paste("value ~", group.var)

    # Add adjustment variables to the formula if provided
    if (!is.null(adj.vars)) {
      formula_str <-
        paste(formula_str, "+", paste(adj.vars, collapse = " + "))
    }

    # Convert the formula string to a formula object
    formula <- as.formula(formula_str)

    # Set the baseline time point if not provided
    if (is.null(change.base)) {
      change.base <- unique(meta_tab %>% select(all_of(c(time.var))))[1,]
      message(
        "The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'meta.dat' data frame: ",
        change.base
      )
    }

    # Adjust filters for 'other' data type
    if (feature.dat.type == "other") {
      prev.filter <- 0
      abund.filter <- 0
    }

    # Normalize count data if necessary
    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      data.obj <-
        mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
    }

    # Perform analysis for each feature level
    test.list <- lapply(feature.level, function(feature.level) {
      # Aggregate data by taxonomy if necessary
      if (is.null(data.obj$feature.agg.list[[feature.level]]) &
          feature.level != "original") {
        data.obj <-
          mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      # Select the appropriate feature table
      if (feature.level != "original") {
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      # Filter the feature table based on prevalence and abundance
      otu_tax_agg_filter <- otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        rownames_to_column(feature.level)

      # Perform data imputation and winsorization for count or proportion data
      if (feature.dat.type %in% c("count", "proportion")) {
        # Calculate half of the minimum non-zero value for each taxon across all samples and time points
        # This ensures the same pseudocount is used for each taxon at both time points,
        # eliminating time-point specific bias in log fold change calculations
        half_nonzero_min <- apply(otu_tax_agg_filter[, -1], 1, function(x) {
          nonzero_values <- x[x > 0]
          if (length(nonzero_values) > 0) {
            min(nonzero_values) / 2
          } else {
            1e-10  # Fallback for taxa with all zeros
          }
        })

        # Create a logical matrix identifying zero values
        zero_matrix <- otu_tax_agg_filter[, -1] == 0

        # Create a matrix of imputation values (per-taxon pseudocount applied across all samples)
        half_nonzero_min_matrix <- matrix(half_nonzero_min, nrow = nrow(zero_matrix),
                                          ncol = ncol(zero_matrix), byrow = FALSE)

        # Impute zero values with half of the minimum non-zero value
        otu_tax_agg_filter[, -1][zero_matrix] <- half_nonzero_min_matrix[zero_matrix]

        message("Zero-handling: Per-taxon half-minimum pseudocount calculated across ALL samples and time points combined. ",
                "This ensures unbiased change calculations by using the same pseudocount at both time points.")

        # Apply winsorization to limit extreme values
        otu_tax_agg_filter[, -1] <- apply(otu_tax_agg_filter[, -1], 2, function(x) {
          qt <- quantile(x, probs = c((1 - winsor.qt) / 2, 1 - (1 - winsor.qt) / 2))
          x[x < qt[1]] <- qt[1]
          x[x > qt[2]] <- qt[2]
          return(x)
        })

      } else if (feature.dat.type == "other") {
        # Apply winsorization for 'other' data type
        otu_tax_agg_filter[, -1] <- apply(otu_tax_agg_filter[, -1], 2, function(x) {
          qt <- quantile(x, probs = c((1 - winsor.qt) / 2, 1 - (1 - winsor.qt) / 2))
          x[x < qt[1]] <- qt[1]
          x[x > qt[2]] <- qt[2]
          return(x)
        })
      }

      # Convert the filtered data from wide to long format
      otu_tax_long <- otu_tax_agg_filter %>%
        tidyr::gather(key = "sample", value = "value", -feature.level)

      # Merge the long-format data with metadata
      merged_data <- otu_tax_long %>%
        dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      # Group the data by time
      grouped_data <- merged_data %>%
        dplyr::group_by(!!sym(time.var))

      # Identify the time point after the baseline
      change.after <-
        unique(grouped_data %>% select(all_of(c(time.var))))[unique(grouped_data %>% select(all_of(c(time.var)))) != change.base]

      # Split the data into separate time points
      split_data <-
        split(merged_data, f = grouped_data %>% select(all_of(c(time.var))))

      # Extract data for the baseline and follow-up time points
      data_time_1 <- split_data[[change.base]]
      data_time_2 <- split_data[[change.after]]

      # Join the baseline and follow-up data
      combined_data <- data_time_1 %>%
        dplyr::inner_join(
          data_time_2,
          by = c(feature.level, subject.var),
          suffix = c("_time_1", "_time_2")
        )

      # Calculate the change in feature values based on the specified method
      if (is.function(feature.change.func)) {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = feature.change.func(value_time_2, value_time_1))
      } else if (feature.change.func == "log fold change") {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = log2(value_time_2) - log2(value_time_1))
      } else if (feature.change.func == "relative change") {
        combined_data <- combined_data %>%
          dplyr::mutate(value_diff = (value_time_2 - value_time_1) / (value_time_2 + value_time_1))
      } else if (feature.change.func == "absolute change"){
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      } else {
        combined_data <-
          combined_data %>% dplyr::mutate(value_diff = value_time_2 - value_time_1)
      }

      message("Note: For repeated measurements of the same subject at the same time point, the average will be taken.")

      # Create a matrix of value differences
      value_diff_matrix <- combined_data %>%
        select(feature.level, !!sym(subject.var), value_diff) %>%
        dplyr::group_by(!!sym(feature.level), !!sym(subject.var)) %>%
        dplyr::summarise(value_diff = mean(value_diff, na.rm = TRUE)) %>%
        tidyr::spread(key = !!sym(subject.var), value = value_diff) %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      # Prepare metadata for analysis
      unique_meta_tab <- meta_tab %>%
        filter(!!sym(subject.var) %in% colnames(value_diff_matrix)) %>%
        select(all_of(c(subject.var, group.var, adj.vars))) %>%
        dplyr::distinct(!!sym(subject.var), .keep_all = TRUE) %>% as_tibble()

      cols_order <- colnames(na.omit(value_diff_matrix))

      unique_meta_tab <-
        unique_meta_tab %>% column_to_rownames(subject.var)

      sorted_unique_meta_tab <- unique_meta_tab %>%
        dplyr::slice(match(cols_order, rownames(unique_meta_tab)))

      # Calculate average abundance and prevalence for each feature
      prop_prev_data <-
        otu_tax_agg %>%
        as.matrix() %>%
        as.table() %>%
        as.data.frame() %>%
        dplyr::group_by(Var1) %>%
        dplyr::summarise(avg_abundance = mean(Freq),
                         prevalence = sum(Freq > 0) / dplyr::n()) %>% column_to_rownames("Var1") %>%
        rownames_to_column(feature.level)

      # Convert the value difference matrix to long format
      value_diff_long <- value_diff_matrix %>%
        as.data.frame() %>%
        rownames_to_column(feature.level) %>%
        tidyr::gather(key = !!sym(subject.var), value = "value", -feature.level)

      # Perform statistical tests for each feature
      sub_test.list <-
        lapply(value_diff_long %>% select(all_of(feature.level)) %>% pull() %>% unique(), function(taxon) {
          # Prepare data for the current feature
          test_df <- value_diff_long %>%
            dplyr::filter(!!sym(feature.level) == taxon) %>%
            dplyr::left_join(sorted_unique_meta_tab %>%
                               as.data.frame() %>%
                               rownames_to_column(subject.var),
                             by = subject.var)

          # Fit the linear model
          test_result <- lm(formula, data = test_df)

          # Extract coefficients from the linear model
          coef.tab <- extract_coef(test_result)

          # Perform ANOVA if the grouping variable has more than two levels
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

            # Combine the coefficient table and ANOVA table
            coef.tab <-
              rbind(coef.tab, anova.tab)
          }
          return(as_tibble(coef.tab))
        })

      # Assign names to the elements of sub_test.list
      names(sub_test.list) <- value_diff_long %>%
        select(all_of(feature.level)) %>%
        pull() %>%
        unique()

      # Extract unique terms related to the grouping variable
      unique_terms <-
        grep(paste0("^", group.var, "$|^", group.var, ".*"),
             unique(unlist(
               lapply(sub_test.list, function(df)
                 unique(df$Term))
             )),
             value = TRUE)

      # Compile results for each term
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

      # Assign names to the result list
      names(result_list) <- unique_terms

      # Modify result names to include reference level information
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