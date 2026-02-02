#' @title Generate Taxa Change Scatterplot for Paired Data
#'
#' @description Creates scatterplots showing taxa abundance changes between paired time points
#' plotted against a continuous grouping variable. Supports stratification for comparative visualization.
#'
#' @inheritParams mStat_data_obj_doc
#' @inheritParams mStat_plot_params_doc
#' @param change.base Character or numeric specifying the baseline time point.
#' @param feature.change.func Method for calculating change: "relative change",
#'   "log fold change", "absolute change", or a custom function.
#' @param features.plot Character vector of specific feature IDs to plot.
#' @param top.k.plot Integer specifying number of top features to plot.
#' @param top.k.func Function for selecting top features (e.g., "mean", "sd").
#'
#' @return A list of ggplot scatterplot objects.
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
           feature.change.func = "relative change",
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
      data.obj <- mStat_normalize_data(data.obj, method = "TSS")$data.obj.norm
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

      # Convert counts to numeric type
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

      # First, divide the data into two subsets, one for change.base and one for change.after
      df_t0 <- otu_tab_norm_agg %>% filter(!!sym(time.var) == change.base)
      df_ts <- otu_tab_norm_agg %>% filter(!!sym(time.var) != change.base)

      # Then, use dplyr::inner_join to merge these two subsets based on Phylum, subject and sex
      df <- dplyr::inner_join(df_ts, df_t0, by = c(feature.level, subject.var), suffix = c("_ts", "_t0"), relationship = "many-to-many")

      # Calculate the change in abundance based on the specified function
      df <- df %>%
        dplyr::mutate(new_count = compute_taxa_change(
          value_after  = count_ts,
          value_before = count_t0,
          method       = feature.change.func,
          feature_id   = .data[[feature.level]]
        ))

      df <- df %>% dplyr::left_join(meta_tab %>% select(-all_of(time.var)) %>% dplyr::distinct(), by = c(subject.var))

      df <- df %>% setNames(ifelse(names(.) == paste0(time.var,"_ts"), time.var, names(.)))

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
          #scale_linetype_manual(values = c("solid", "dashed")) + # Set the line type
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
