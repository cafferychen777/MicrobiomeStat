#' Generate taxa area plots over time
#'
#' This function generates taxa area plots for a given data object. The plots will show the relative abundance of
#' different taxa over time. Raw count data will be automatically normalized using rarefaction and total sum scaling (TSS).
#' The function also supports the generation of plots for grouped data and stratified data.
#'
#' @param data.obj A data object containing the taxa count data.
#' @param subject.var Name of the subject variable.
#' @param time.var Name of the time variable.
#' @param group.var Optional, name of the group variable. Default is NULL.
#' @param strata.var Optional, name of the stratification variable. Default is NULL.
#' @param feature.level The taxonomic level to plot. Default is "original".
#' @param feature.dat.type The type of features data; can be "count", "proportion", or "other". Default is "count".
#' @param feature.number The number of features to plot. Default is 20.
#' @param t0.level The initial time level. Default is NULL.
#' @param ts.levels The time series levels. Default is NULL.
#' @param base.size The base size for the ggplot2 theme. Default is 10.
#' @param theme.choice The theme choice for the ggplot2 theme. Default is "bw".
#' @param custom.theme Optional, a custom ggplot2 theme. Default is NULL.
#' @param palette Optional, a palette to use for the plot. Default is NULL.
#' @param pdf Logical indicating if the plot should be saved as a PDF. Default is TRUE.
#' @param file.ann Optional, a file annotation. Default is NULL.
#' @param pdf.wid Width of the output PDF. Default is 11.
#' @param pdf.hei Height of the output PDF. Default is 8.5.
#' @param ... Additional arguments to pass to the function.
#'
#' @return A list of ggplot objects, each representing a taxa area plot for the specified feature level.
#'
#' @examples
#' library(ggh4x)
#' library(vegan)
#' data(ecam.obj)
#' plot_list_all <- generate_taxa_areaplot_long(
#'   data.obj = ecam.obj,
#'   subject.var = "studyid",
#'   time.var = "month",
#'   group.var = "diet",
#'   strata.var = "antiexposedall",
#'   feature.level = "Family",
#'   feature.dat.type = "proportion",
#'   feature.number = 8,
#'   t0.level = unique(ecam.obj$meta.dat$month)[1],
#'   ts.levels = unique(ecam.obj$meta.dat$month)[-1],
#'   base.size = 10,
#'   theme.choice = "bw",
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = "test"
#' )
#' @import rlang
#' @import tibble
#' @export
generate_taxa_areaplot_long <-
  function(data.obj,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           feature.level = "original",
           feature.dat.type = c("count", "proportion", "other"),
           feature.number = 20,
           t0.level = NULL,
           ts.levels = NULL,
           base.size = 10,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    feature.dat.type <- match.arg(feature.dat.type)
    # Data validation
    mStat_validate_data(data.obj)

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

    if (feature.dat.type == "count"){
      message("Your data is in raw format ('Raw'). Normalization is crucial for further analyses. Now, 'mStat_normalize_data' function is automatically applying 'Rarefy-TSS' transformation.")
      otu_tab <- load_data_obj_count(mStat_normalize_data(data.obj, method = "Rarefy-TSS")$data.obj.norm)
    } else if (feature.dat.type == "other"){
      stop("The 'other' type is suitable for situations where the user has analyzed the data using a method not provided in 'mStat_normalize_data' method, and the 'barplot' is only applicable to raw data that has not undergone any processing or proportion data that adds up to 1. If you believe your data falls into these two categories, please modify 'feature.dat.type'.")
    } else if (feature.dat.type == "proportion"){
      otu_tab <- load_data_obj_count(data.obj)
    }

    tax_tab <- load_data_obj_taxonomy(data.obj) %>%
      as.data.frame() %>%
      {if("original" %in% feature.level) dplyr::mutate(., original = rownames(.)) else .} %>%
      select(all_of(feature.level))

    theme_function <- switch(theme.choice,
                             prism = ggprism::theme_prism(),
                             classic = theme_classic(),
                             gray = theme_gray(),
                             bw = theme_bw(),
                             ggprism::theme_prism()) # 根据用户选择设置主题

    # 使用用户自定义主题（如果提供），否则使用默认主题
    theme_to_use <- if (!is.null(custom.theme)) custom.theme else theme_function

    if (is.null(palette)){
      pal <- rep(c("#E41A1C","#1E90FF","#FF8C00","#4DAF4A","#984EA3","#40E0D0","#FFC0CB",
                   "#00BFFF","#FFDEAD","#90EE90","#EE82EE","#00FFFF","#F0A3FF", "#0075DC",
                   "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00",
                   "#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010",
                   "#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",
                   "#FF6A6A", "#FF8247", "#FFE7BA", "#87CEFA", "#B0E0E6", "#48D1CC", "#5F9EA0", "#66CDAA", "#458B00",
                   "#BCEE68", "#FFF68F", "#EEEE00", "#FFFFE0", "#8B8682", "#FFB6C1", "#9370DB", "#FFDAB9", "#FA8072",
                   "#90EE90", "#00FFFF", "#00BFFF", "#FFA07A", "#F08080", "#FFD700", "#ADFF2F", "#9ACD32", "#9400D3",
                   "#7B68EE", "#BA55D3", "#FFC0CB", "#FAEBD7", "#F0E68C", "#FFFACD", "#D2B48C", "#C0C0C0", "#696969",
                   "#CD5C5C", "#F08080", "#FA8072", "#E9967A", "#FFA07A", "#FF7F50", "#FF6347", "#FF4500", "#FF8C00",
                   "#FFA500", "#FFDAB9", "#FFE4B5", "#FFE4C4", "#FFDEAD", "#EEE8AA", "#F0E68C", "#BDB76B", "#FFD700",
                   "#DAA520", "#808000", "#7CFC00", "#00FF00", "#32CD32", "#00FA9A", "#90EE90", "#98FB98", "#8FBC8F",
                   "#3CB371", "#2E8B57", "#228B22", "#008000", "#006400", "#66CDAA", "#00FFFF", "#40E0D0", "#48D1CC",
                   "#AFEEEE", "#7FFFD4", "#B0C4DE", "#ADD8E6", "#87CEEB", "#87CEFA", "#6495ED", "#00BFFF", "#1E90FF",
                   "#4169E1", "#0000FF", "#00008B", "#000080", "#191970", "#6A5ACD", "#483D8B", "#9400D3", "#8A2BE2",
                   "#4B0082", "#FF00FF", "#FF69B4", "#FF1493", "#C71585", "#DB7093", "#FFC0CB", "#FFB6C1", "#FF69B4",
                   "#FF5F5F", "#DC143C", "#C0C0C0", "#A9A9A9", "#808080", "#696969", "#000000", "#FF1493", "#FF69B4",
                   "#FFB6C1", "#FFC0CB", "#FFDAB9", "#F4A460", "#FFA07A", "#FF7F50", "#FF6347", "#FF4500", "#FF8C00",
                   "#FFA500", "#FFFF00", "#9ACD32", "#32CD32", "#00FF00", "#7FFF00", "#7CFC00", "#00FA9A", "#90EE90",
                   "#98FB98", "#8FBC8F", "#3CB371", "#2E8B57", "#228B22", "#008000", "#006400"
      ),5)
    } else{

      if (feature.number > length(palette)) {
        stop("The number of unique features exceeds the length of the provided palette. Please provide a larger palette.")
      }

      pal = palette
    }

    plot_list_all <- lapply(feature.level,function(feature.level){
      # 将 OTU 表与分类表合并
      otu_tax <- cbind(otu_tab, tax_tab %>% select(all_of(feature.level)))

      # 聚合 OTU 表
      otu_tax_agg <- otu_tax %>%
        tidyr::gather(key = "sample", value = "value", -one_of(feature.level)) %>%
        dplyr::group_by_at(vars(sample, !!sym(feature.level))) %>%
        dplyr::summarise(value = sum(value)) %>%
        tidyr::spread(key = "sample", value = "value")

      # 转换计数为数值类型
      otu_tax_agg_numeric <- dplyr::mutate_at(otu_tax_agg, vars(-!!sym(feature.level)), as.numeric)

      # 标准化数据(Have been dropped out)
      otu_tab_norm <- apply(t(otu_tax_agg_numeric %>% select(-feature.level)), 1, function(x) x)

      rownames(otu_tab_norm) <- as.matrix(otu_tax_agg_numeric[, feature.level])

      meta_tab_sorted <- meta_tab[colnames(otu_tab_norm), ]

      # 计算每个taxon的平均相对丰度
      avg_abund <- rowMeans(otu_tab_norm)

      # 将相对丰度低于阈值的taxon替换为"Other"
      otu_tab_other <- otu_tab_norm %>%
        as.data.frame() %>%
        rownames_to_column(feature.level)

      # 将feature.number之后以下的相对丰度定为阈值
      other.abund.cutoff <- sort(avg_abund, decreasing=TRUE)[feature.number]

      if (!is.na(other.abund.cutoff)){
        otu_tab_other[, feature.level][avg_abund < other.abund.cutoff] <- "Other"
      }

      # 转换数据框为长格式
      otu_tab_long <- otu_tab_other %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::summarize_all(sum) %>%
        tidyr::gather(key = "sample", value = "value", -feature.level)

      # 将 otu_tab_long 和 meta_tab_sorted 合并
      merged_long_df <- otu_tab_long %>%
        dplyr::inner_join(meta_tab_sorted  %>% rownames_to_column("sample"), by = "sample")

      sorted_merged_long_df <- merged_long_df %>%
        dplyr::arrange(!!sym(subject.var), !!sym(time.var))

      last_sample_ids <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(subject.var)) %>%
        dplyr::summarize(last_sample_id = dplyr::last(sample))

      sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(feature.level) := as.factor(!!sym(feature.level)))

      original_levels <- levels(sorted_merged_long_df[[feature.level]])

      if (!is.na(other.abund.cutoff)){
        new_levels <- c("Other", setdiff(original_levels, "Other"))
      } else {
        new_levels <- original_levels
      }

      sorted_merged_long_df <- sorted_merged_long_df %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = original_levels)) %>%
        dplyr::mutate(!!sym(feature.level) := forcats::fct_relevel(!!sym(feature.level), new_levels))

      df <- sorted_merged_long_df %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = original_levels)) %>%
        dplyr::mutate(!!sym(feature.level) := forcats::fct_relevel(!!sym(feature.level), new_levels)) %>%
        dplyr::arrange(match(!!sym(feature.level), new_levels)) %>%
        dplyr::mutate(cumulative_value = (1-cumsum(value))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::mutate(next_cumulative_value = dplyr::if_else(sample %in% last_sample_ids$last_sample_id, NA_real_, dplyr::lead(cumulative_value))) %>%
        dplyr::ungroup()

      color_pal <- setNames(pal, as.matrix(unique(df %>% select(!!sym(feature.level)))))

      bar_width <- 0.6
      bar_spacing <- bar_width / 2

      # 以下为average barplot的绘制
      last_time_ids <- dplyr::last(levels(meta_tab[,time.var]))

      if (!is.null(strata.var)){
        if (!is.null(group.var)){
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := interaction(!!sym(group.var),!!sym(strata.var)))
        } else {
          group.var = ""
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := "")
        }
      } else {
        if (!is.null(group.var)){
        } else {
          group.var = ""
          sorted_merged_long_df <- sorted_merged_long_df %>% dplyr::mutate(!!sym(group.var) := "")
        }
      }

      df_average <- sorted_merged_long_df %>%
        dplyr::group_by(!!sym(feature.level),!!sym(group.var),!!sym(time.var)) %>%
        dplyr::summarise(mean_value  = mean(value)) %>%
        dplyr::mutate(!!sym(feature.level) := factor(!!sym(feature.level), levels = original_levels)) %>%
        dplyr::mutate(!!sym(feature.level) := forcats::fct_relevel(!!sym(feature.level), new_levels)) %>%
        dplyr::arrange(match(!!sym(feature.level), new_levels),!!sym(group.var),!!sym(time.var)) %>%
        dplyr::group_by(!!sym(group.var),!!sym(time.var)) %>%
        dplyr::mutate(cumulative_mean_value = (1-cumsum(mean_value))) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(!!sym(feature.level)) %>%
        dplyr::mutate(next_cumulative_mean_value = dplyr::if_else(!!sym(time.var) %in% last_time_ids, NA_real_, dplyr::lead(cumulative_mean_value))) %>%
        dplyr::ungroup()

      if (group.var == ""){
        df_average <- df_average %>% dplyr::mutate(!!sym(group.var) := "")
        df_average <- df_average %>%
          dplyr::mutate(joint_factor = interaction(!!sym(time.var), ""))
      } else {
        df_average <- df_average %>%
          dplyr::mutate(joint_factor = interaction(!!sym(time.var), !!sym(group.var)))
      }

      df_average <- df_average %>%
        dplyr::mutate(x_offset = ifelse(cumulative_mean_value == 0, (bar_width + bar_spacing) / 2, -(bar_width + bar_spacing) / 2))

      df_average$joint_factor <- droplevels(df_average$joint_factor)

      df_average$joint_factor_numeric <- match(df_average$joint_factor, levels(df_average$joint_factor))

      labels <- sub("\\..*", "", levels(df_average$joint_factor))

      if(!is.null(strata.var)){
        df_average <- df_average %>%
          tidyr::separate(!!sym(group.var), into = c(group.var, strata.var), sep = "\\.")
      }

      stack_areaplot_average  <- # Main plot code
        df_average %>%
        ggplot(aes(x = joint_factor_numeric, y = mean_value, fill = !!sym(feature.level))) +
        geom_area(stat = "identity", position = "fill") +
        {
          if (all(round(apply(otu_tab, 2, sum),2) == 1)){
            scale_y_continuous(expand = c(0, 0), labels = scales::percent)
          } else {
            scale_y_continuous(expand = c(0, 0))
          }
        } +
        scale_x_continuous(expand = c(0.01, 0.01), breaks = unique(df_average$joint_factor_numeric), labels = labels) +
        {
          if (!is.null(group.var)){
            if (group.var == ""){
            } else {
              if (!is.null(strata.var)){
                ggh4x::facet_nested(as.formula(paste(". ~", group.var, "+", strata.var)), drop = T, scale = "free", space = "free")
              } else {
                ggh4x::facet_nested(as.formula(paste(". ~", group.var)), drop = T, scale = "free", space = "free")
              }
            }
          }
        } +
        labs(fill = feature.level, y = "", x = "") +
        scale_fill_manual(values = color_pal) +
        scale_color_manual(values = color_pal) +
        theme_to_use +
        theme(strip.background = element_rect(fill="white",color="black"),
              panel.spacing = unit(0,"lines"),
              strip.text.x = element_text(size= base.size,color="black"),
              axis.text.y=element_text(size= base.size,color="black"),
              axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, size = base.size),
              axis.title.y = element_text(size= base.size,color="black"),
              legend.key=element_blank(),
              legend.title = element_text(size= base.size+2), # Add this line to change the legend title size
              legend.text = element_text(color="black",size= base.size),
              legend.spacing.x=unit(0.1,'cm'),
              legend.spacing.y=unit(0.1,'cm'),
              legend.key.width=unit(0.4,'cm'),
              legend.key.height=unit(0.4,'cm'),
              legend.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())


      # Save the stacked barplots as a PDF file
      if (pdf) {
        pdf_name <- paste0("taxa_areaplot_pair",
                           "_",
                           "subject_", subject.var,
                           "_",
                           "time_", time.var,
                           "_",
                           "feature_level_", feature.level,
                           "_",
                           "feature_number_", feature.number)
        if (!is.null(group.var)) {
          pdf_name <- paste0(pdf_name, "_", "group_", group.var)
        }
        if (!is.null(strata.var)) {
          pdf_name <- paste0(pdf_name, "_", "strata_", strata.var)
        }
        if (!is.null(file.ann)) {
          pdf_name <- paste0(pdf_name, "_", file.ann)
        }
        pdf_name <- paste0(pdf_name,"_avergae", ".pdf")
        ggsave(filename = pdf_name, plot = stack_areaplot_average, width = pdf.wid, height = pdf.hei)
      }

      # 返回堆叠条形图以进行显示
      return(stack_areaplot_average)
    })
    return(plot_list_all)
  }
