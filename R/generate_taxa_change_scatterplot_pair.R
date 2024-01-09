#' Generate taxa-level difference scatterplot pairs
#'
#' This function generates scatterplot pairs showing the difference in abundance of taxa

#' at a specified taxonomic level between two time points, with grouping and stratification.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var Character string specifying the subject variable in metadata.
#' @param time.var Character string specifying the time variable in metadata.
#' @param group.var Character string specifying the grouping variable in metadata. Default NULL.
#' @param strata.var Character string specifying the stratification variable in metadata. Default NULL.
#' @param change.base Character string specifying the baseline time point. This should match one of the time points present in the metadata for the 'time.var' variable. The change will be calculated by comparing the other time points to this baseline time point.
#' @param feature.change.func Specifies the method or function used to compute the change between two time points.
#' The following options are supported:
#'
#' - "absolute change": Computes the difference between the counts at the two time points (`count_ts` and `count_t0`). This is the default behavior if no value or an unrecognized value is provided.
#'
#' - "log fold change": Computes the log2 fold change between the two time points. For zero counts, imputation is performed using half of the minimum nonzero count for each taxa level at the respective time point before taking the logarithm.
#'
#' - "relative change": Computes the relative change as `(count_ts - count_t0) / (count_ts + count_t0)`. If both time points have a count of 0, the change is defined as 0.
#'
#' - A custom function: If a user-defined function is provided, it should take two numeric vectors as input corresponding to the counts at the two time points (`count_t0` and `count_ts`) and return a numeric vector of the computed change. This custom function will be applied directly to calculate the new counts.
#' @param feature.level The column name in the feature annotation matrix (feature.ann) of data.obj
#' to use for summarization and plotting. This can be the taxonomic level like "Phylum", or any other
#' annotation columns like "Genus" or "OTU_ID". Should be a character vector specifying one or more
#' column names in feature.ann. Multiple columns can be provided, and data will be plotted separately
#' for each column. Default is NULL, which defaults to all columns in feature.ann if `features.plot`
#' is also NULL.
#' @param feature.dat.type The type of the feature data, which determines how the data is handled in downstream analyses.
#' Should be one of:
#' - "count": Raw count data, will be normalized by the function.
#' - "proportion": Data that has already been normalized to proportions/percentages.
#' - "other": Custom abundance data that has unknown scaling. No normalization applied.
#' The choice affects preprocessing steps as well as plot axis labels.
#' Default is "count", which assumes raw OTU table input.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering
#' taxa before analysis. Taxa with prevalence below this value will be removed.
#' Prevalence is calculated as the proportion of samples where the taxon is present.
#' Default 0 removes no taxa by prevalence filtering.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering
#' taxa before analysis. Taxa with mean abundance below this value will be removed.
#' Abundance refers to counts or proportions depending on \code{feature.dat.type}.
#' Default 0 removes no taxa by abundance filtering.
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL, in which case features will be selected based on `top.k.plot` and `top.k.func`.
#' @param top.k.plot Integer specifying number of top k features to plot, when `features.plot` is NULL.
#' Default is NULL, in which case all features passing filters will be plotted.
#' @param top.k.func Function to use for selecting top k features, when `features.plot` is NULL.
#' Options include inbuilt functions like "mean", "sd", or a custom function. Default is NULL, in which
#' case features will be selected by abundance.
#' @param base.size Numeric specifying base font size for plot text elements. Default 16.
#' @param theme.choice
#' Plot theme choice. Specifies the visual style of the plot. Can be one of the following pre-defined themes:
#'   - "prism": Utilizes the ggprism::theme_prism() function from the ggprism package, offering a polished and visually appealing style.
#'   - "classic": Applies theme_classic() from ggplot2, providing a clean and traditional look with minimal styling.
#'   - "gray": Uses theme_gray() from ggplot2, which offers a simple and modern look with a light gray background.
#'   - "bw": Employs theme_bw() from ggplot2, creating a classic black and white plot, ideal for formal publications and situations where color is best minimized.
#'   - "light": Implements theme_light() from ggplot2, featuring a light theme with subtle grey lines and axes, suitable for a fresh, modern look.
#'   - "dark": Uses theme_dark() from ggplot2, offering a dark background, ideal for presentations or situations where a high-contrast theme is desired.
#'   - "minimal": Applies theme_minimal() from ggplot2, providing a minimalist theme with the least amount of background annotations and colors.
#'   - "void": Employs theme_void() from ggplot2, creating a blank canvas with no axes, gridlines, or background, ideal for custom, creative plots.
#' Each theme option adjusts various elements like background color, grid lines, and font styles to match the specified aesthetic.
#' Default is "bw", offering a universally compatible black and white theme suitable for a wide range of applications.
#' @param custom.theme
#' A custom ggplot theme provided as a ggplot2 theme object. This allows users to override the default theme and provide their own theme for plotting. Custom themes are useful for creating publication-ready figures with specific formatting requirements.
#'
#' To use a custom theme, create a theme object with ggplot2::theme(), including any desired customizations. Common customizations for publication-ready figures might include adjusting text size for readability, altering line sizes for clarity, and repositioning or formatting the legend. For example:
#'
#' ```r
#' my_theme <- ggplot2::theme(
#'   axis.title = ggplot2::element_text(size=14, face="bold"),        # Bold axis titles with larger font
#'   axis.text = ggplot2::element_text(size=12),                      # Slightly larger axis text
#'   legend.position = "top",                                         # Move legend to the top
#'   legend.background = ggplot2::element_rect(fill="lightgray"),     # Light gray background for legend
#'   panel.background = ggplot2::element_rect(fill="white", colour="black"), # White panel background with black border
#'   panel.grid.major = ggplot2::element_line(colour = "grey90"),     # Lighter color for major grid lines
#'   panel.grid.minor = ggplot2::element_blank(),                     # Remove minor grid lines
#'   plot.title = ggplot2::element_text(size=16, hjust=0.5)           # Centered plot title with larger font
#' )
#' ```
#'
#' Then pass `my_theme` to `custom.theme`. If `custom.theme` is NULL (the default), the theme is determined by `theme.choice`. This flexibility allows for both easy theme selection for general use and detailed customization for specific presentation or publication needs.
#' @param palette An optional parameter specifying the color palette to be used for the plot.
#'                It can be either a character string specifying the name of a predefined
#'                palette or a vector of color codes in a format accepted by ggplot2
#'                (e.g., hexadecimal color codes). Available predefined palettes include
#'                'npg', 'aaas', 'nejm', 'lancet', 'jama', 'jco', and 'ucscgb', inspired
#'                by various scientific publications and the `ggsci` package. If `palette`
#'                is not provided or an unrecognized palette name is given, a default color
#'                palette will be used. Ensure the number of colors in the palette is at
#'                least as large as the number of groups being plotted.
#' @param pdf Logical, if TRUE, save plot(s) as PDF file(s). Default TRUE.
#' @param file.ann Character string specifying file annotation to add to PDF file name. Default NULL.
#' @param pdf.wid,pdf.hei Numeric specifying PDF width and height in inches. Default 11 x 8.5.
#' @param ... Additional arguments passed to ggplot2::ggsave().
#'
#' @return A ggplot object containing the scatterplot pair(s).
#'
#' @examples
#' \dontrun{
#' # Load required libraries and data
#' library(vegan)
#' data(peerj32.obj)
#' peerj32.obj$meta.dat <- peerj32.obj$meta.dat %>%
#' dplyr::select(all_of("subject")) %>% dplyr::distinct() %>%
#' dplyr::mutate(cons = runif(dplyr::n(),0,5)) %>%
#' dplyr::left_join(peerj32.obj$meta.dat %>% rownames_to_column("sample"),by = "subject") %>%
#' tibble::column_to_rownames("sample")
#' # Generate the boxplot pair
#' generate_taxa_change_scatterplot_pair(
#'   data.obj = peerj32.obj,
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "cons",
#'   strata.var = "sex",
#'   change.base = "1",
#'   feature.change.func = "log fold change",
#'   feature.level = "Genus",
#'   feature.dat.type = "other",
#'   features.plot = NULL,
#'   top.k.plot = 8,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_pairs.obj")
#' subset_pairs.obj$meta.dat <- subset_pairs.obj$meta.dat %>%
#' dplyr::select(all_of("MouseID")) %>% dplyr::distinct() %>%
#' dplyr::mutate(cons = runif(dplyr::n(),0,5)) %>%
#' dplyr::left_join(subset_pairs.obj$meta.dat %>% rownames_to_column("sample"),by = "MouseID") %>%
#' tibble::column_to_rownames("sample")
#' # Generate the boxplot pair
#' generate_taxa_change_scatterplot_pair(
#'   data.obj = subset_pairs.obj,
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "cons",
#'   strata.var = NULL,
#'   change.base = "Baseline",
#'   feature.change.func = "log fold change",
#'   feature.level = "Genus",
#'   feature.dat.type = "other",
#'   features.plot = NULL,
#'   top.k.plot = 8,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_taxa_change_scatterplot_pair <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           change.base = NULL,
           feature.change.func = "absolute change",
           feature.level = NULL,
           feature.dat.type = c("count", "proportion", "other"),
           features.plot = NULL,
           top.k.plot = NULL,
           top.k.func = NULL,
           prev.filter = 0.01,
           abund.filter = 0.01,
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    feature.dat.type <- match.arg(feature.dat.type)

    mStat_validate_data(data.obj)

    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% select(all_of(c(
        time.var, group.var, strata.var, subject.var
      )))

    if (is.null(group.var)) {
      group.var = "ALL"
      meta_tab$ALL <- ""
    }

    if (is.null(strata.var)) {
      strata.var = "ALL2"
      meta_tab$ALL2 <- ""
    }

    # Define the colors
    colors <- mStat_get_palette(palette)

    # Assuming mStat_get_theme function is already defined
    # Replace the existing theme selection code with this:
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    aes_function <- if (!is.null(strata.var)){
      aes(shape = !!sym(strata.var), color = !!sym(strata.var))
    } else {
      aes(color = !!sym(time.var))
    }

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    if (feature.dat.type == "count"){
      message(
        "Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'TSS' transformation."
      )
      data.obj <- mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm
    }

    plot_list <- lapply(feature.level, function(feature.level) {

      if (is.null(data.obj$feature.agg.list[[feature.level]]) & feature.level != "original"){
        data.obj <- mStat_aggregate_by_taxonomy(data.obj = data.obj, feature.level = feature.level)
      }

      if (feature.level != "original"){
        otu_tax_agg <- data.obj$feature.agg.list[[feature.level]]
      } else {
        otu_tax_agg <- data.obj$feature.tab
      }

      otu_tax_agg <-  otu_tax_agg %>%
        as.data.frame() %>%
        mStat_filter(prev.filter = prev.filter,
                     abund.filter = abund.filter) %>%
        rownames_to_column(feature.level)

      if (is.null(features.plot) && !is.null(top.k.plot) && !is.null(top.k.func)) {
      computed_values <- compute_function(top.k.func, otu_tax_agg, feature.level)
      features.plot <- names(sort(computed_values, decreasing = TRUE)[1:top.k.plot])
      }

      # 转换计数为数值类型
      otu_tax_agg_numeric <-
        dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      otu_tab_norm <- otu_tax_agg_numeric %>%
        column_to_rownames(var = feature.level) %>%
        as.matrix()

      otu_tab_norm_agg <- otu_tax_agg_numeric %>%
        tidyr::gather(-!!sym(feature.level), key = "sample", value = "count") %>%
        dplyr::inner_join(meta_tab %>% rownames_to_column("sample"), by = "sample")

      taxa.levels <-
        otu_tab_norm_agg %>% select(all_of(feature.level)) %>% dplyr::distinct() %>% dplyr::pull()

      # 首先，把数据分为两个子集，一个为change.base，一个为change.after
      df_t0 <- otu_tab_norm_agg %>% filter(!!sym(time.var) == change.base)
      df_ts <- otu_tab_norm_agg %>% filter(!!sym(time.var) != change.base)

      # 然后，使用dplyr::inner_join合并这两个子集，基于Phylum、subject和sex
      df <- dplyr::inner_join(df_ts, df_t0, by = c(feature.level, subject.var), suffix = c("_ts", "_t0"), relationship = "many-to-many")

      # 最后，计算新的count值
      if (is.function(feature.change.func)) {
        df <- df %>% dplyr::mutate(new_count = feature.change.func(count_ts, count_t0))
      } else if (feature.change.func == "log fold change") {
        half_nonzero_min_time_2 <- df %>%
          filter(count_ts > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(count_ts) / 2,
                    .groups = "drop")
        half_nonzero_min_time_1 <- df %>%
          filter(count_t0 > 0) %>%
          dplyr::group_by(!!sym(feature.level)) %>%
          dplyr::summarize(half_nonzero_min = min(count_t0) / 2,
                    .groups = "drop")

        df <- dplyr::left_join(df, half_nonzero_min_time_2, by = feature.level, suffix = c("_t0", "_ts"))
        df <- dplyr::left_join(df, half_nonzero_min_time_1, by = feature.level, suffix = c("_t0", "_ts"))
        df$count_ts[df$count_ts == 0] <- df$half_nonzero_min_ts[df$count_ts == 0]
        df$count_t0[df$count_t0 == 0] <- df$half_nonzero_min_t0[df$count_t0 == 0]

        # Add a message to inform users that an imputation operation was performed.
        message("Imputation was performed using half the minimum nonzero count for each taxa at different time points.")

        df <- df %>% dplyr::mutate(new_count = log2(count_ts) - log2(count_t0))
      } else if (feature.change.func == "relative change"){
        df <- df %>%
          dplyr::mutate(new_count = dplyr::case_when(
            count_ts == 0 & count_t0 == 0 ~ 0,
            TRUE ~ (count_ts - count_t0) / (count_ts + count_t0)
          ))
      } else if (feature.change.func == "absolute change"){
        df <- df %>% dplyr::mutate(new_count = count_ts - count_t0)
      } else {
        df <- df %>% dplyr::mutate(new_count = count_ts - count_t0)
      }

      df <- df %>% dplyr::left_join(meta_tab %>% select(-all_of(time.var)) %>% dplyr::distinct(), by = c(subject.var))

      df <- df %>% setNames(ifelse(names(.) == paste0(time.var,"_ts"), time.var, names(.)))

      # 提前判断feature.change.func的类型，如果是自定义函数则给出特定的标签，否则保持原样
      ylab_label <- if (feature.dat.type != "other") {
        if (is.function(feature.change.func)) {
          paste0("Change in Relative Abundance", " (custom function)")
        } else {
          paste0("Change in Relative Abundance", " (", feature.change.func, ")")
        }
      } else {
        if (is.function(feature.change.func)) {
          paste0("Change in Abundance", " (custom function)")
        } else {
          paste0("Change in Abundance", " (", feature.change.func, ")")
        }
      }

      if (!is.null(features.plot)) {

      } else {
        if (length(taxa.levels) >= 5) {
          features.plot <- taxa.levels[1:4]
        } else {
          features.plot <- taxa.levels
        }
      }

        scatterplot <-
          ggplot(df %>% filter(!!sym(feature.level) %in% features.plot),
                 aes(
                   x = !!sym(group.var),
                   y = new_count,
                   fill = !!sym(strata.var),
                   color = !!sym(strata.var),  # Add this line to tidyr::separate color by strata.var
                   group = !!sym(strata.var)   # Add this line to tidyr::separate smoothing line by strata.var
                 )) +
          geom_smooth(se = TRUE, method = 'lm') +
          geom_point(aes_function,data = df %>% filter(!!sym(feature.level) %in% features.plot),
                     size = 4) +
          scale_shape_manual(values = c(21, 22, 24, 25)) +
          scale_fill_manual(values = colors) +
          ylab(ylab_label) +
          #scale_linetype_manual(values = c("solid", "dashed")) + # 设置曲线类型
          scale_color_manual(values = colors, guide = guide_legend(override.aes = list(size = 0))) +
          theme_to_use +
          theme(
            axis.text.x = element_text(
              vjust = 0.5,
              hjust = 1,
              size = base.size
            ),
            axis.text.y = element_text(color = "black", size = base.size),
            axis.title.y = element_text(size = base.size),
            axis.title.x = element_text(size = base.size),
            legend.position = "right",
            legend.direction = "vertical",
            legend.box = "vertical",
            plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
            axis.text = ggplot2::element_text(color = "black",
                                              size = base.size),
            legend.text = ggplot2::element_text(size = 16),
            legend.title = ggplot2::element_text(size = 16),
            plot.title = element_text(hjust = 0.5, size = 20)
          ) +
          ggh4x::facet_nested_wrap(as.formula(paste(".~",feature.level)), scales = "fixed")


        if (group.var == "ALL"){
          scatterplot <- scatterplot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
        }

        if (strata.var == "ALL2") {
          scatterplot <- scatterplot + theme(
            legend.title = element_blank()
          )
        }

        # Save the stacked dotplot as a PDF file
        if (pdf) {
          pdf_name <- paste0(
            "taxa_change_scatterplot",
            "_",
            "subject_",
            subject.var,
            "_",
            "time_",
            time.var,
            "_",
            "feature_level_",
            feature.level,
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
          pdf_name <- paste0(pdf_name, ".pdf")
          pdf(pdf_name, width = pdf.wid, height = pdf.hei)
          # Use lapply to print each ggplot object in the list to a new PDF page
          print(scatterplot)
          # Close the PDF device
          dev.off()
        }

      return(scatterplot)
    })

    names(plot_list) <- feature.level
    return(plot_list)
  }
