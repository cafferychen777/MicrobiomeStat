#' mStat_combine_data
#'
#' This function is designed specifically for combining two MicrobiomeStat data objects,
#' assuming that they are both in 'Raw' format. Each MicrobiomeStat data object should be a list that contains
#' 'feature.tab' (a matrix of features),
#' 'meta.dat' (a data frame of metadata), and 'feature.ann' (a matrix of feature annotations).
#'
#' @param data.obj1 The first MicrobiomeStat data object to be combined.
#' @param data.obj2 The second MicrobiomeStat data object to be combined.
#'
#' @return A list with the combined 'feature.tab', 'meta.dat' and 'feature.ann'.
#'
#' @examples
#' # combined_data <- mStat_combine_data(data_obj1, data_obj2)
#'
#' @export
mStat_combine_data <- function(data.obj1, data.obj2) {

  message("Combining two MicrobiomeStat data objects in Raw format...")

  if (is.null(data.obj1$feature.tab) || is.null(data.obj2$feature.tab)) {
    stop("Both data objects must contain 'feature.tab'.")
  }
  if (is.null(data.obj1$meta.dat) || is.null(data.obj2$meta.dat)) {
    stop("Both data objects must contain 'meta.dat'.")
  }

  combine_feature_table <- function(tab1, tab2) {
    tab1 <- as.matrix(tab1)
    tab2 <- as.matrix(tab2)

    if (is.null(rownames(tab1)) || is.null(colnames(tab1)) ||
        is.null(rownames(tab2)) || is.null(colnames(tab2))) {
      stop("feature.tab in both data objects must have row and column names.")
    }
    if (anyDuplicated(rownames(tab1)) != 0 || anyDuplicated(colnames(tab1)) != 0 ||
        anyDuplicated(rownames(tab2)) != 0 || anyDuplicated(colnames(tab2)) != 0) {
      stop("feature.tab row and column names must be unique in both data objects.")
    }
    if (!is.numeric(tab1) || !is.numeric(tab2)) {
      stop("feature.tab in both data objects must be numeric.")
    }

    common_features <- intersect(rownames(tab1), rownames(tab2))
    common_samples <- intersect(colnames(tab1), colnames(tab2))

    if (length(common_features) > 0 && length(common_samples) > 0) {
      tab1_common <- tab1[common_features, common_samples, drop = FALSE]
      tab2_common <- tab2[common_features, common_samples, drop = FALSE]
      conflict <- !is.na(tab1_common) & !is.na(tab2_common) & (tab1_common != tab2_common)
      if (any(conflict)) {
        stop("Inconsistent values found in overlapping feature.tab entries.")
      }
    }

    all_features <- union(rownames(tab1), rownames(tab2))
    all_samples <- union(colnames(tab1), colnames(tab2))

    combined <- matrix(
      NA_real_,
      nrow = length(all_features),
      ncol = length(all_samples),
      dimnames = list(all_features, all_samples)
    )

    combined[rownames(tab1), colnames(tab1)] <- tab1

    block <- combined[rownames(tab2), colnames(tab2), drop = FALSE]
    block_conflict <- !is.na(block) & !is.na(tab2) & (block != tab2)
    if (any(block_conflict)) {
      stop("Inconsistent values found while merging feature.tab.")
    }
    block[is.na(block)] <- tab2[is.na(block)]
    combined[rownames(tab2), colnames(tab2)] <- block

    combined[is.na(combined)] <- 0
    combined
  }

  combine_row_annotated <- function(df1, df2, object_name, as_matrix = FALSE) {
    df1 <- as.data.frame(df1, stringsAsFactors = FALSE)
    df2 <- as.data.frame(df2, stringsAsFactors = FALSE)

    if (is.null(rownames(df1)) || is.null(rownames(df2))) {
      stop(object_name, " in both data objects must have row names.")
    }
    if (anyDuplicated(rownames(df1)) != 0 || anyDuplicated(rownames(df2)) != 0) {
      stop(object_name, " row names must be unique in both data objects.")
    }

    id_col <- ".row_id"
    left <- tibble::rownames_to_column(df1, id_col)
    right <- tibble::rownames_to_column(df2, id_col)

    merged <- dplyr::full_join(left, right, by = id_col, suffix = c(".x", ".y"))
    shared_cols <- intersect(colnames(df1), colnames(df2))

    for (col in shared_cols) {
      x_col <- paste0(col, ".x")
      y_col <- paste0(col, ".y")

      x <- merged[[x_col]]
      y <- merged[[y_col]]

      conflict <- !is.na(x) & !is.na(y) & (as.character(x) != as.character(y))
      if (any(conflict)) {
        stop(
          "Inconsistent values found in ", object_name,
          " for overlapping row IDs in column '", col, "'."
        )
      }

      merged[[col]] <- dplyr::coalesce(x, y)
      merged[[x_col]] <- NULL
      merged[[y_col]] <- NULL
    }

    ordered_cols <- unique(c(id_col, colnames(df1), colnames(df2)))
    merged <- merged[, ordered_cols, drop = FALSE]

    rownames(merged) <- merged[[id_col]]
    merged[[id_col]] <- NULL

    if (as_matrix) {
      return(as.matrix(merged))
    }

    as.data.frame(merged, stringsAsFactors = FALSE)
  }

  combined_feature_tab <- combine_feature_table(data.obj1$feature.tab, data.obj2$feature.tab)

  combined_feature_ann <- NULL
  if (!is.null(data.obj1$feature.ann) || !is.null(data.obj2$feature.ann)) {
    if (is.null(data.obj1$feature.ann) || is.null(data.obj2$feature.ann)) {
      stop("Both data objects must contain 'feature.ann' when either one provides it.")
    }
    combined_feature_ann <- combine_row_annotated(
      data.obj1$feature.ann,
      data.obj2$feature.ann,
      object_name = "feature.ann",
      as_matrix = TRUE
    )
  }

  combined_meta_dat <- combine_row_annotated(
    data.obj1$meta.dat,
    data.obj2$meta.dat,
    object_name = "meta.dat",
    as_matrix = FALSE
  )

  combined <- list(
    feature.tab = combined_feature_tab,
    meta.dat = combined_meta_dat,
    feature.ann = combined_feature_ann
  )

  message("Returning combined data object.")
  mStat_validate_data(combined)
}
