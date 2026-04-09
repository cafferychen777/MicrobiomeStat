#' @title Taxa Data Utilities
#' @description Helper functions for retrieving, aggregating, and filtering
#'   feature (OTU/ASV/taxa) tables. Centralizes the aggregation + extraction +
#'   filtering pipeline that was previously duplicated across 30+ files.
#' @name taxa_utils
#' @keywords internal
NULL


#' Retrieve a filtered feature table for a given taxonomic level
#'
#' Performs the three-step pipeline that appears in virtually every
#' \code{generate_taxa_*} function:
#' \enumerate{
#'   \item Aggregate features to \code{feature.level} if not already cached.
#'   \item Extract the aggregated table (or the raw table for \code{"original"}).
#'   \item Optionally filter by prevalence / abundance, then convert rownames
#'         to a column named \code{feature.level}.
#' }
#'
#' @param data.obj A MicrobiomeStat data object containing at least
#'   \code{feature.tab} and optionally \code{feature.agg.list}.
#' @param feature.level Character string specifying the taxonomic level to
#'   extract. Use \code{"original"} to retrieve \code{data.obj$feature.tab}
#'   without aggregation.
#' @param prev.filter Numeric prevalence threshold passed to
#'   \code{\link{mStat_filter}}. Set to 0 to skip filtering. Default 0.
#' @param abund.filter Numeric abundance threshold passed to
#'   \code{\link{mStat_filter}}. Set to 0 to skip filtering. Default 0.
#' @param feature.col Logical; if \code{TRUE} (default), rownames are converted
#'   to a leading column named \code{feature.level} via
#'   \code{tibble::rownames_to_column}. Set to \code{FALSE} when the downstream
#'   consumer expects a matrix-like object with rownames (e.g. \code{linda()}).
#'
#' @return When \code{feature.col = TRUE}, a data.frame with a leading column
#'   named \code{feature.level}. When \code{FALSE}, a data.frame with rownames.
#'
#' @keywords internal
get_taxa_data <- function(data.obj,
                          feature.level,
                          prev.filter = 0,
                          abund.filter = 0,
                          feature.col = TRUE) {

  # Step 1: Aggregate if needed (idempotent, cached in data.obj$feature.agg.list)
  if (is.null(data.obj$feature.agg.list[[feature.level]]) &
      feature.level != "original") {
    data.obj <-
      mStat_aggregate_by_taxonomy(data.obj = data.obj,
                                  feature.level = feature.level)
  }

  # Step 2: Extract the appropriate table
  if (feature.level != "original") {
    otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
  } else {
    otu_tax_agg <- data.obj$feature.tab
  }

  # Step 3: Filter
  result <- otu_tax_agg %>%
    as.data.frame() %>%
    mStat_filter(prev.filter = prev.filter,
                 abund.filter = abund.filter)

  if (feature.col) {
    result <- tibble::rownames_to_column(result, feature.level)
  }

  result
}


#' @keywords internal
mStat_as_taxa_feature_matrix <- function(feature.dat,
                                         feature.level,
                                         feature_in_column = FALSE) {
  feature.df <- as.data.frame(
    feature.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (feature_in_column) {
    if (!feature.level %in% colnames(feature.df)) {
      stop(
        "The feature column `", feature.level, "` is missing from the taxa data.",
        call. = FALSE
      )
    }
    feature.df <- tibble::column_to_rownames(feature.df, var = feature.level)
  }

  feature_ids <- rownames(feature.df)
  if (is.null(feature_ids) || anyNA(feature_ids) || anyDuplicated(feature_ids)) {
    stop(
      "Taxa data must have unique, non-missing feature row names.",
      call. = FALSE
    )
  }

  data.matrix(feature.df)
}


#' @keywords internal
mStat_as_taxa_composition_matrix <- function(feature.dat,
                                             feature.level,
                                             feature_in_column = TRUE) {
  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.dat,
    feature.level = feature.level,
    feature_in_column = feature_in_column
  )

  sample_totals <- colSums(feature.mat, na.rm = TRUE)
  normalized <- matrix(
    0,
    nrow = nrow(feature.mat),
    ncol = ncol(feature.mat),
    dimnames = dimnames(feature.mat)
  )

  nonzero_samples <- is.finite(sample_totals) & sample_totals != 0
  if (any(nonzero_samples)) {
    normalized[, nonzero_samples] <- sweep(
      feature.mat[, nonzero_samples, drop = FALSE],
      2,
      sample_totals[nonzero_samples],
      "/"
    )
  }

  normalized
}


#' @keywords internal
mStat_meta_with_sample <- function(meta.dat,
                                   sample_col = "sample") {
  meta.df <- as.data.frame(
    meta.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (sample_col %in% colnames(meta.df)) {
    return(meta.df)
  }

  tibble::rownames_to_column(meta.df, var = sample_col)
}


#' @keywords internal
mStat_resolve_optional_flag <- function(flag,
                                        default,
                                        name) {
  if (is.null(flag)) {
    return(default)
  }

  if (!is.logical(flag) || length(flag) != 1 || is.na(flag)) {
    stop("`", name, "` should be a single TRUE/FALSE value or NULL.", call. = FALSE)
  }

  flag
}


#' @keywords internal
mStat_prepare_taxa_long_data <- function(feature.dat,
                                         feature.level,
                                         value_col = "value",
                                         meta.dat = NULL,
                                         sample_col = "sample",
                                         feature_in_column = TRUE,
                                         join = c("left", "inner")) {
  join <- match.arg(join)

  feature.df <- as.data.frame(
    feature.dat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!feature_in_column) {
    feature.df <- tibble::rownames_to_column(feature.df, var = feature.level)
  }

  long.df <- feature.df %>%
    tidyr::pivot_longer(
      cols = -all_of(feature.level),
      names_to = sample_col,
      values_to = value_col
    ) %>%
    dplyr::mutate(!!value_col := as.numeric(.data[[value_col]]))

  if (is.null(meta.dat)) {
    return(long.df)
  }

  meta.df <- mStat_meta_with_sample(meta.dat = meta.dat, sample_col = sample_col)

  if (join == "inner") {
    dplyr::inner_join(long.df, meta.df, by = sample_col)
  } else {
    dplyr::left_join(long.df, meta.df, by = sample_col)
  }
}


#' @keywords internal
mStat_summarize_grouped_taxa_long <- function(long.df,
                                              feature.level,
                                              group_vars = NULL,
                                              value_col = "value",
                                              mean_col = "mean_abundance",
                                              prevalence_col = "prevalence",
                                              mean_transform = NULL) {
  summary.df <- long.df %>%
    dplyr::group_by(dplyr::across(all_of(c(group_vars, feature.level)))) %>%
    dplyr::summarise(
      .mean_value = mean(.data[[value_col]], na.rm = TRUE),
      .prevalence = mean(.data[[value_col]] != 0, na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(mean_transform)) {
    summary.df$.mean_value <- mean_transform(summary.df$.mean_value)
  }

  summary.df %>%
    dplyr::rename(
      !!mean_col := .mean_value,
      !!prevalence_col := .prevalence
    )
}


#' @keywords internal
mStat_transform_taxa_long_values <- function(long.df,
                                             feature.level,
                                             value_col = "value",
                                             feature.dat.type = c("count", "proportion", "other"),
                                             transform = c("identity", "sqrt", "log")) {
  feature.dat.type <- match.arg(feature.dat.type)
  transform <- match.arg(transform)

  if (!feature.dat.type %in% c("count", "proportion") || transform == "identity") {
    return(long.df)
  }

  if (transform == "sqrt") {
    long.df[[value_col]] <- sqrt(long.df[[value_col]])
    return(long.df)
  }

  min_half_nonzero <- long.df %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::filter(sum(.data[[value_col]], na.rm = TRUE) != 0) %>%
    dplyr::summarise(
      min_half_value = min(.data[[value_col]][.data[[value_col]] > 0], na.rm = TRUE) / 2,
      .groups = "drop"
    )

  long.df %>%
    dplyr::group_by(!!rlang::sym(feature.level)) %>%
    dplyr::filter(sum(.data[[value_col]], na.rm = TRUE) != 0) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(min_half_nonzero, by = feature.level) %>%
    dplyr::mutate(
      !!value_col := ifelse(
        .data[[value_col]] == 0,
        log10(min_half_value),
        log10(.data[[value_col]])
      )
    ) %>%
    dplyr::select(-min_half_value)
}


#' @keywords internal
mStat_prepare_taxa_clr_long_data <- function(feature.dat,
                                             feature.level,
                                             prev.filter = 0,
                                             abund.filter = -Inf,
                                             sample_col = "sample",
                                             value_col = "value",
                                             feature_in_column = FALSE) {
  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.dat,
    feature.level = feature.level,
    feature_in_column = feature_in_column
  )

  sample_totals <- colSums(feature.mat, na.rm = TRUE)
  valid_samples <- is.finite(sample_totals) & sample_totals > 0

  if (!any(valid_samples)) {
    return(tibble::tibble(
      !!feature.level := character(),
      !!sample_col := character(),
      !!value_col := numeric()
    ))
  }

  if (any(!valid_samples)) {
    warning(
      "Found ", sum(!valid_samples),
      " samples with zero total abundance during CLR preparation; these samples were excluded.",
      call. = FALSE
    )
  }

  sample_ids <- colnames(feature.mat)[valid_samples]
  transformed_samples <- lapply(sample_ids, function(sample_id) {
    sample_values <- feature.mat[, sample_id]
    positive_values <- sample_values[sample_values > 0 & is.finite(sample_values)]

    pseudocount <- min(positive_values, na.rm = TRUE) / 2
    imputed_values <- ifelse(sample_values == 0, pseudocount, sample_values)
    geometric_mean <- exp(mean(log(imputed_values), na.rm = TRUE))

    log(imputed_values / geometric_mean)
  })

  clr.mat <- do.call(cbind, transformed_samples)
  rownames(clr.mat) <- rownames(feature.mat)
  colnames(clr.mat) <- sample_ids

  clr.filtered <- mStat_filter(
    clr.mat,
    prev.filter = prev.filter,
    abund.filter = abund.filter
  )

  mStat_prepare_taxa_long_data(
    feature.dat = clr.filtered,
    feature.level = feature.level,
    value_col = value_col,
    sample_col = sample_col,
    feature_in_column = FALSE
  )
}


#' @keywords internal
mStat_summarize_taxa_features <- function(feature.dat,
                                          feature.level,
                                          feature_in_column = FALSE) {
  feature.mat <- mStat_as_taxa_feature_matrix(
    feature.dat = feature.dat,
    feature.level = feature.level,
    feature_in_column = feature_in_column
  )
  feature_ids <- rownames(feature.mat)

  tibble::tibble(
    !!feature.level := feature_ids,
    avg_abundance = unname(rowMeans(feature.mat, na.rm = TRUE)),
    prevalence = unname(rowMeans(feature.mat != 0, na.rm = TRUE))
  )
}


#' @keywords internal
mStat_format_linda_feature_results <- function(result.df,
                                               feature.level,
                                               feature.stats,
                                               include_significant = FALSE) {
  required_cols <- c("log2FoldChange", "lfcSE", "pvalue", "padj")
  missing_cols <- setdiff(required_cols, colnames(result.df))
  if (length(missing_cols) != 0) {
    stop(
      "LinDA-style results are missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (include_significant && !"reject" %in% colnames(result.df)) {
    stop(
      "LinDA-style results must include `reject` when `include_significant = TRUE`.",
      call. = FALSE
    )
  }

  formatted <- tibble::rownames_to_column(
    as.data.frame(result.df, stringsAsFactors = FALSE),
    var = feature.level
  ) %>%
    dplyr::left_join(feature.stats, by = feature.level)

  select_cols <- c(
    feature.level,
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj",
    if (include_significant) "reject",
    "avg_abundance",
    "prevalence"
  )

  formatted %>%
    dplyr::select(all_of(select_cols)) %>%
    dplyr::rename(
      Variable = all_of(feature.level),
      Coefficient = log2FoldChange,
      SE = lfcSE,
      P.Value = pvalue,
      Adjusted.P.Value = padj,
      Mean.Abundance = avg_abundance,
      Prevalence = prevalence
    ) %>%
    {
      if (include_significant) {
        dplyr::rename(., Significant = reject)
      } else {
        .
      }
    }
}
