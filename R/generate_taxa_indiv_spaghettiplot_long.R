#' Generate longitudinal line plots of taxonomic composition
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param subject.var A character string defining subject variable in meta_tab
#' @param time.var A character string defining time variable in meta_tab
#' @param group.var A character string defining group variable in meta_tab used for sorting and facetting
#' @param strata.var (Optional) A character string defining strata variable in meta_tab used for sorting and facetting
#' @param taxa.level A character string defining the taxonomic level to analyze ('Phylum', 'Family', or 'Genus')
#' @param prev.filter A numeric value defining the prevalence threshold to filter taxa, between 0 and 1
#' @param abund.filter A numeric value defining the abundance threshold to filter taxa
#' @param pdf A logical value. If TRUE (default), saves the dotplot as a PDF file. If FALSE, the dotplot will be displayed interactively without creating a PDF
#' @param file.ann (Optional) A character string specifying a file annotation to include in the generated PDF file's name
#' @param ... Additional parameters to be passed
#' @return If the `pdf` parameter is set to TRUE, the function will save a PDF file and return the final ggplot object. If `pdf` is set to FALSE, the function will return the final ggplot object without creating a PDF file.
#' @examples
#'
#' plot_list_all <- generate_taxa_indiv_spaghettiplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   t0.level = NULL,
#'   ts.levels = NULL,
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = "Phylum",
#'   feature.dat.type = "other",
#'   features.plot = NULL,
#'   top.k.plot = 5,
#'   top.k.func = "mean",
#'   prev.filter = 0.01,
#'   abund.filter = 0.01,
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = c(ggsci::pal_npg()(9),ggsci::pal_jama()(7),ggsci::pal_lancet()(9)),
#'   pdf = TRUE,
#'   file.ann = "test"
#' )
#'
#' @export
generate_taxa_indiv_spaghettiplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           t0.level = NULL,
           ts.levels = NULL,
           group.var = NULL,
           strata.var = NULL,
           feature.level = NULL,
           features.plot = NULL,
           feature.dat.type = c("count", "proportion", "other"),
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
    # Data validation
    mStat_validate_data(data.obj)

    feature.dat.type <- match.arg(feature.dat.type)

    if (!is.character(subject.var))
      stop("`subject.var` should be a character string.")
    if (!is.character(time.var))
      stop("`time.var` should be a character string.")
    if (!is.null(group.var) &&
        !is.character(group.var))
      stop("`group.var` should be a character string or NULL.")
    if (!is.null(strata.var) &&
        !is.character(strata.var))
      stop("`strata.var` should be a character string or NULL.")

    # 提取数据
    data.obj <- mStat_process_time_variable(data.obj, time.var, t0.level, ts.levels)

    meta_tab <- load_data_obj_metadata(data.obj) %>% as.data.frame() %>% select(all_of(c(subject.var,group.var,time.var,strata.var)))

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
          mutate(., original = rownames(.))
        else
          .
      } %>%
      select(all_of(feature.level))

    if (feature.dat.type == "other" || !is.null(features.plot) ||
        (!is.null(top.k.func) && !is.null(top.k.plot))) {
      prev.filter <- 0
      abund.filter <- 0
    }

    plot_list_all <- lapply(feature.level, function(feature.level) {
      # 将 OTU 表与分类表合并
      otu_tax <-
        cbind(otu_tab, tax_tab %>% select(all_of(feature.level)))

      # Filter taxa based on prevalence and abundance
      otu_tax_filtered <- otu_tax %>%
        gather(key = "sample", value = "count",-one_of(feature.level)) %>%
        group_by_at(vars(!!sym(feature.level))) %>%
        summarise(total_count = mean(count),
                  prevalence = sum(count > 0) / n()) %>%
        filter(prevalence >= prev.filter, total_count >= abund.filter) %>%
        select(-total_count,-prevalence) %>%
        left_join(otu_tax, by = feature.level)

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax_filtered %>%
        gather(key = "sample",  value = "count",-one_of(feature.level)) %>%
        group_by_at(vars(sample,!!sym(feature.level))) %>%
        summarise(count = sum(count)) %>%
        spread(key = "sample", value = "count")

      compute_function <- function(top.k.func) {
        if (is.function(top.k.func)) {
          results <-
            top.k.func(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
        } else {
          switch(top.k.func,
                 "mean" = {
                   results <-
                     rowMeans(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                              na.rm = TRUE)
                 },
                 "sd" = {
                   results <-
                     rowSds(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix(),
                            na.rm = TRUE)
                   names(results) <- rownames(otu_tax_agg %>% column_to_rownames(feature.level) %>% as.matrix())
                 },
                 stop("Invalid function specified"))
        }

        return(results)
      }

      if (is.null(features.plot) &&
          !is.null(top.k.plot) && !is.null(top.k.func)) {
        features.plot <- names(sort(compute_function(top.k.func), decreasing = TRUE)[1:top.k.plot])
      }

      # 转换计数为数值类型
      otu_tax_agg_numeric <-
        mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      df <- otu_tax_agg_numeric %>%
        gather(key = "sample", value = "count",-one_of(feature.level)) %>%
        left_join(meta_tab %>% rownames_to_column(var = "sample"), by = "sample")

      if (is.null(group.var)) {
        df <- df %>% mutate("ALL" = "ALL")
        group.var = "ALL"
      }

      if (!is.null(strata.var)) {
        mean_df <-
          df %>% group_by(!!sym(feature.level),!!sym(time.var),!!sym(group.var),!!sym(strata.var)) %>%
          summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          left_join(df,
                    mean_df,
                    by = c(feature.level, time.var, group.var, strata.var))
      } else {
        mean_df <-
          df %>% group_by(!!sym(feature.level),
                          !!sym(time.var),
                          !!sym(group.var)) %>%
          summarize(mean_count = mean(count), na.rm = TRUE)
        df <-
          left_join(df, mean_df, by = c(feature.level, time.var, group.var))
      }

      if (is.null(palette)) {
        col <-
          c(
            "#E31A1C",
            "#1F78B4",
            "#FB9A99",
            "#33A02C",
            "#FDBF6F",
            "#B2DF8A",
            "#A6CEE3",
            "#BA7A70",
            "#9D4E3F",
            "#829BAB"
          )
      } else{
        col <- palette
      }

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
          custom.theme else
        theme_function

      # Calculate new sizes based on base.size
      title.size = base.size * 1.25
      axis.title.size = base.size * 0.75
      axis.text.size = base.size * 0.5
      legend.title.size = base.size * 1
      legend.text.size = base.size * 0.75

      if (is.null(features.plot)){
        taxa.levels <- df %>% select(feature.level) %>% distinct() %>% pull()
      } else {
        taxa.levels <- df %>% filter(!!sym(feature.level) %in% features.plot) %>% select(feature.level) %>% distinct() %>% pull()
      }

      plot_list <- lapply(taxa.levels, function(tax) {
        sub_df <- df %>% filter(!!sym(feature.level) == tax)
        lineplot <- ggplot() +
          geom_line(
            data = sub_df,
            aes_string(
              x = time.var,
              y = "count",
              group = subject.var,
              color = group.var
            ),
            alpha = 0.5
          ) +
          geom_line(
            data = sub_df,
            aes_string(
              x = time.var,
              y = "mean_count",
              group = group.var,
              color = group.var
            ),
            size = 2
          ) +
          scale_color_manual(values = col) +
          {
            if (feature.dat.type != "other") {
              labs(
                x = time.var,
                y = "Relative Abundance",
                color = group.var,
                title = tax
              )
            } else {
              labs(
                x = time.var,
                y = "Abundance",
                color = group.var,
                title = tax
              )
            }
          } +
          {
            if (!is.null(strata.var)) {
              facet_wrap(as.formula(paste('~', strata.var)))  # Use facet_wrap with strata.var as the faceting variable
            }
          } +
          theme_to_use +
          theme(
            plot.title = element_text(
              size = title.size,
              face = "bold",
              hjust = 0.5
            ),
            axis.title.x = element_text(size = axis.title.size, face = "bold"),
            axis.title.y = element_text(size = axis.title.size, face = "bold"),
            axis.text.x = element_text(size = axis.text.size),
            axis.text.y = element_text(size = axis.text.size),
            legend.title = element_text(size = legend.title.size, face = "bold"),
            legend.text = element_text(size = legend.text.size)
          )

        if (group.var == "ALL") {
          lineplot <- lineplot + theme(legend.position = "none")
        }
        return(lineplot)
      })

      # Save the plots as a PDF file
      if (pdf) {
        pdf_name <- paste0(
          "taxa_indiv_spaghettiplot_long",
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "group_",
          group.var,
          "_",
          "strata_",
          strata.var,
          "_",
          "feature_level_",
          feature.level,
          "_",
          "prev_filter_",
          prev.filter,
          "_",
          "abund_filter_",
          abund.filter,
          "_",
          "base_size_",
          base.size,
          "_",
          "theme_choice_",
          theme.choice,
          "_",
          "pdf_wid_",
          pdf.wid,
          "_",
          "pdf_hei_",
          pdf.hei
        )

        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }

        pdf_name <- paste0(pdf_name, ".pdf")

        # Create a multi-page PDF file
        pdf(pdf_name, width = pdf.wid, height = pdf.hei)
        # Use lapply to print each ggplot object in the list to a new PDF page
        lapply(plot_list, print)
        # Close the PDF device
        dev.off()
      }
      return(plot_list)
    })

    return(plot_list_all)
  }
