#' Generate Single Time Point Taxa Boxplots
#'
#' This function generates boxplots to visualize the taxonomic composition of samples for a single time point in a longitudinal study.
#' It provides options for grouping and stratifying data, and selecting the top k features based on a user-defined function.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A string specifying the variable for subjects.
#' @param time.var A string specifying the variable for time. If NULL, the function assumes that data for a single time point is provided.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var Optional string specifying the variable for groups.
#' @param strata.var Optional string specifying the variable for strata.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param top.k.plot An integer specifying the top k features to plot based on the function specified in `top.k.func`.
#' @param top.k.func A function to determine the top k features to plot.
#' @param transform A string indicating the transformation to apply to the data before plotting. Options are:
#' - "identity": No transformation (default)
#' - "sqrt": Square root transformation
#' - "log": Logarithmic transformation. Zeros are replaced with half of the minimum non-zero value for each taxon before log transformation.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param base.size A numeric value specifying the base font size for the plot.
#' @param theme.choice Plot theme choice. Can be one of:
#'   - "prism": ggprism::theme_prism()
#'   - "classic": theme_classic()
#'   - "gray": theme_gray()
#'   - "bw": theme_bw()
#' Default is "bw".
#' @param custom.theme A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. To use a custom theme, first create a theme object with ggplot2::theme(), then pass it to this argument. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=16, color="red"),
#'   legend.position = "none"
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. Default is NULL, which will use the default theme based on `theme.choice`.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf A logical value indicating whether to save the plot as a PDF. Default is TRUE.
#' @param file.ann A string for additional annotation to the file name. Default is NULL.
#' @param pdf.wid A numeric value specifying the width of the PDF file. Default is 11.
#' @param pdf.hei A numeric value specifying the height of the PDF file. Default is 8.5.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of ggplot objects, one for each taxonomic level.
#' @details
#' This function generates a boxplot of taxa abundances for a single time point in a longitudinal study.
#' The boxplot can be stratified by a group variable and/or other variables.
#' It also allows for different taxonomic levels to be used and a specific number of features to be included in the plot.
#' The function also has options to customize the size, theme, and color palette of the plot, and to save the plot as a PDF.
#'
#' @examples
#' \dontrun{
#' # Generate the boxplot pair
#' data(ecam.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t.level = "1",
#'   group.var = "diet",
#'   strata.var = NULL,
#'   feature.level = c("Phylum"),
#'   features.plot = sample(unique(ecam.obj$feature.ann[,"Phylum"]),3),
#'   feature.dat.type = "proportion",
#'   transform = "log",
#'   prev.filter = 0,
#'   abund.filter = 0,
#'   base.size = 12,
#'   theme.choice = "classic",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test",
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data(peerj32.obj)
#' generate_taxa_boxplot_single(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   t.level = "1",
#'   group.var = "group",
#'   strata.var = "sex",
#'   feature.level = c("Family"),
#'   feature.dat.type = "count",
#'   features.plot = NULL,
#'   top.k.plot = NULL,
#'   top.k.func = NULL,
#'   transform = "log",
#'   prev.filter = 0.1,
#'   abund.filter = 0.0001,
#'   base.size = 12,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_boxplot_single <-
  function(data.obj,
           subject.var,
           time.var = NULL,
           t.level = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           transform = c("identity", "sqrt", "log"),
           prev.filter = 0.05,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "prism",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {
    feature.dat.type <- match.arg(feature.dat.type)

    mStat_validate_data(data.obj)

    if (!is.character(subject.var))
      stop("`subject.var` should be a character string.")
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    if (!is.null(time.var) & !is.null(t.level)) {
      condition <- paste(time.var, "== '", t.level, "'", sep = "")
      data.obj <-
        mStat_subset_data(data.obj, condition = condition)
    }

    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        subject.var, group.var, time.var, strata.var
      )))

    aes_function <- aes(
      x = !!sym(group.var),
      y = value,
      fill = !!sym(group.var)
    )

    col <- mStat_get_palette(palette)

    theme_function <- switch(
      theme.choice,
      prism = ggprism::theme_prism(),
      classic = theme_classic(),
      gray = theme_gray(),
      bw = theme_bw(),
      ggprism::theme_prism()
    ) # 根据用户选择设置主题

    # 使用用户自定义主题（如果提供），否则使用默认主题
    theme_to_use <-
      if (!is.null(custom.theme))
        custom.theme
    else
      theme_function

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    if (feature.dat.type == "count") {
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation."
      )
      data.obj <-
        mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm
    }

    plot_list <- lapply(feature.level, function(feature.level) {
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

      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        tibble::rownames_to_column(feature.level)

      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      otu_tax_agg_numeric <- otu_tax_agg %>%
        tidyr::gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        dplyr::mutate(value = as.numeric(value))

      otu_tax_agg_merged <-
        dplyr::left_join(otu_tax_agg_numeric,
                         meta_tab %>% rownames_to_column("sample"),
                         by = "sample") %>%
        select(one_of(
          c(
            "sample",
            feature.level,
            subject.var,
            time.var,
            group.var,
            strata.var,
            "value"
          )
        ))

      if (feature.dat.type %in% c("count", "proportion")) {
        # Apply transformation
        if (transform %in% c("identity", "sqrt", "log")) {
          if (transform == "identity") {
            # No transformation needed
          } else if (transform == "sqrt") {
            otu_tax_agg_merged$value <- sqrt(otu_tax_agg_merged$value)
          } else if (transform == "log") {
            # Find the half of the minimum non-zero proportion for each taxon
            min_half_nonzero <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::summarise(min_half_value = min(value[value > 0]) / 2) %>%
              dplyr::ungroup()
            # Replace zeros with the log of the half minimum non-zero proportion
            otu_tax_agg_merged <- otu_tax_agg_merged %>%
              dplyr::group_by(!!sym(feature.level)) %>%
              filter(sum(value) != 0) %>%
              dplyr::ungroup() %>%
              dplyr::left_join(min_half_nonzero, by = feature.level) %>%
              dplyr::mutate(value = ifelse(value == 0, log10(min_half_value), log10(value))) %>%
              select(-min_half_value)
          }
        }
      }

      taxa.levels <-
        otu_tax_agg_merged %>% select(feature.level) %>% dplyr::distinct() %>% dplyr::pull()

      n_subjects <-
        length(unique(otu_tax_agg_merged[[subject.var]]))

      sub_otu_tax_agg_merged <- otu_tax_agg_merged

      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 6) {
          features.plot <- taxa.levels[1:6]
        } else {
          features.plot <- taxa.levels
        }
      }

      if (!is.null(group.var)){
        group.levels <- sub_otu_tax_agg_merged %>% select(!!sym(group.var)) %>% pull() %>% unique() %>% length()
      }

      boxplot <-
        ggplot(sub_otu_tax_agg_merged %>% filter(!!sym(feature.level) %in% features.plot),
               aes_function) +
        geom_violin(trim = FALSE, alpha = 0.8) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
          fill = "white"
        ) +
        geom_jitter(width = 0.1,
                    alpha = 0.1,
                    size = 1) +
        scale_fill_manual(values = col) +
        {
          if (feature.dat.type == "other") {
            labs(y = "Abundance")
          } else {
            labs(y = paste("Relative Abundance(", transform, ")"))
          }
        } +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          strip.text.x = element_text(size = base.size, color = "black"),
          strip.text.y = element_text(size = base.size, color = "black"),
          axis.text = element_text(color = "black"),
          axis.text.x =  element_blank(),
          axis.text.y = element_text(color = "black", size = (base.size -
                                                                2)),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      if (!is.null(group.var)) {
        if (is.null(strata.var)) {
          boxplot <-
            boxplot + ggh4x::facet_nested_wrap(as.formula(paste(
              "~", feature.level, "+", group.var
            )), scales = "free_x", ncol = group.levels*4)
        } else {
          boxplot <- boxplot + ggh4x::facet_nested_wrap(
            as.formula(paste('~', feature.level, '+', strata.var, '+', group.var)),
            scales = "free_x", ncol = group.levels*4
          )
        }
      }

      if (feature.dat.type != "other") {
        # 添加对Y轴刻度的修改
        if (transform == "sqrt") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  a ^ b, list(a = i, b = 2)
                )))
          )
        } else if (transform == "log") {
          boxplot <- boxplot + scale_y_continuous(
            labels = function(x)
              sapply(x, function(i)
                as.expression(substitute(
                  10 ^ a, list(a = i)
                )))
          )
        }
      }

      return(boxplot)
    })

    # Save the plots as a PDF file
    if (pdf) {
      pdf_name <- paste0(
        "taxa_boxplot_single",
        "_",
        "subject_",
        subject.var,
        "_",
        "feature_level_",
        feature.level,
        "_",
        "transform_",
        transform,
        "_",
        "prev_filter_",
        prev.filter,
        "_",
        "abund_filter_",
        abund.filter
      )
      if (!is.null(group.var)) {
        pdf_name <- paste0(pdf_name, "_", "group_", group.var)
      }
      if (!is.null(strata.var)) {
        pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
      }
      if (!is.null(file.ann)) {
        pdf_name <- paste0(pdf_name, "_", file.ann)
      }
      pdf_name <- paste0(pdf_name, "_", feature.level, ".pdf")
      # Create a multi-page PDF file
      pdf(pdf_name, width = pdf.wid, height = pdf.hei)
      # Use lapply to print each ggplot object in the list to a new PDF page
      lapply(plot_list, print)
      # Close the PDF device
      dev.off()
    }

    names(plot_list) <- feature.level

    return(plot_list)
  }
