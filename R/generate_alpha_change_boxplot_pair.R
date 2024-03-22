#' Generate boxplot comparing change in specified alpha diversity index
#'
#' This function generates a boxplot comparing the change in a specified alpha diversity index between two time points. The change can be calculated as the log or the absolute value. Several optional arguments are available for customizing the output, such as strata or group variables.
#' @name generate_alpha_change_boxplot_pair
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param alpha.obj An optional list containing pre-calculated alpha diversity indices. If NULL (default), alpha diversity indices will be calculated using mStat_calculate_alpha_diversity function from MicrobiomeStat package.
#' @param alpha.name The alpha diversity index to be plotted. Supported indices include "shannon", "simpson", "observed_species", "chao1", "ace", and "pielou". Previously named as `alpha.index`.
#' @param depth An integer. The sequencing depth to be used for the "Rarefy" and "Rarefy-TSS" methods. If NULL, the smallest total count dplyr::across samples is used as the rarefaction depth.
#' @param time.var The variable in the metadata table that represents the time.
#' @param group.var (Optional) The variable in the metadata table that represents the grouping factor.
#' @param strata.var (Optional) The variable in the metadata table that represents the stratification factor.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param change.base The base time for calculating the change in alpha diversity.
#' @param alpha.change.func Function or method for calculating change in alpha diversity
#'   between two timepoints. This allows flexible options to quantify change:
#'
#'   - If a function is provided: The function will be applied to compare alpha diversity
#'     at timepoint t vs baseline t0. The function should take two arguments
#'     representing the alpha diversity values at t and t0. For instance, a custom function to
#'     calculate the percentage change might look like:
#'     \preformatted{
#'       percentage_change <- function(t, t0) {
#'         return ((t - t0) / t0) * 100
#'       }
#'     }
#'     You can then pass this function as the value for `alpha.change.func`.
#'
#'   - If a string is provided, the following options are supported:
#'     - 'log fold change': Calculates the log2 fold change of alpha diversity at t compared to t0.
#'     - 'absolute change': Calculates the absolute difference in alpha diversity at t compared to t0.
#'     - Any other value: A warning will be given that the provided method is not recognized,
#'       and the default method ('absolute change') will be used.
#'
#'   - Default behavior (if no recognized string or function is provided) is to compute the absolute difference between t and t0.
#' @details This parameter allows flexible quantification of how alpha diversity
#'   changes from baseline. Log-ratio is commonly used to compare relative
#'   difference. Absolute difference indicates the magnitude of change.
#'   Custom functions can also be supplied to calculate change as needed.
#' @param base.size The base font size for the plot.
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
#' @param pdf.wid The width of the output PDF file. Default is 11.
#' @param pdf.hei The height of the output PDF file. Default is 8.5.
#' @param subject.var The variable in the metadata table that represents the subject.
#' @param pdf (Optional) A boolean indicating whether to save the output as a PDF file.
#' @param file.ann (Optional) A string for annotating the output file name.
#' @param ... (Optional) Additional arguments to pass to the plotting function.
#'
#' @return A boxplot displaying the change in the specified alpha diversity index between two time points, stratified by the specified grouping and/or strata variables (if provided). The boxplot will be saved as a PDF if `pdf` is set to `TRUE`.
#' @details This function extracts metadata, calculates alpha diversity for specified
#'   indices, compares values between two timepoints, applies log2 fold change by default,
#'   visualizes the change using boxplots and optional faceting and saving as PDF.
#' @examples
#' \dontrun{
#' library(vegan)
#' data(peerj32.obj)
#' generate_alpha_change_boxplot_pair(
#'   data.obj = peerj32.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("simpson"),
#'   subject.var = "subject",
#'   time.var = "time",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   change.base = "1",
#'   alpha.change.func = "absolute change",
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' data("subset_pairs.obj")
#' generate_alpha_change_boxplot_pair(
#'   data.obj = subset_pairs.obj,
#'   alpha.obj = NULL,
#'   alpha.name = c("simpson"),
#'   subject.var = "MouseID",
#'   time.var = "Antibiotic",
#'   group.var = "Sex",
#'   strata.var = NULL,
#'   adj.vars = NULL,
#'   change.base = "Baseline",
#'   alpha.change.func = "log fold change",
#'   base.size = 16,
#'   theme.choice = "bw",
#'   palette = "lancet",
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#' }
#' @export
generate_alpha_change_boxplot_pair <-
  function(data.obj,
           alpha.obj = NULL,
           alpha.name = c("shannon", "observed_species"),
           depth = NULL,
           subject.var,
           time.var,
           group.var = NULL,
           strata.var = NULL,
           adj.vars = NULL,
           change.base = NULL,
           alpha.change.func = c("log fold change"),
           base.size = 16,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    if (is.null(alpha.name)){
      return()
    }

    if (is.null(alpha.obj)) {
      if (!is_rarefied(data.obj)) {
        message(
          "Diversity analysis needs rarefaction! Call 'mStat_rarefy_data' to rarefy the data!"
        )
        data.obj <- mStat_rarefy_data(data.obj, depth = depth)
      }
      otu_tab <- data.obj$feature.tab
      alpha.obj <-
        mStat_calculate_alpha_diversity(x = otu_tab, alpha.name = alpha.name)
    }

    meta_tab <-
      data.obj$meta.dat %>% as.data.frame() %>% dplyr::select(all_of(c(
        subject.var, group.var, time.var, strata.var, adj.vars
      )))

    # Convert the alpha.obj list to a data frame
    alpha_df <-
      dplyr::bind_cols(alpha.obj) %>% rownames_to_column("sample") %>%
      dplyr::inner_join(meta_tab %>% rownames_to_column("sample"),
                        by = c("sample"))

    if (is.null(change.base)) {
      change.base <-
        unique(alpha_df %>% dplyr::select(all_of(c(time.var))))[1, ]
      message(
        "The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'alpha_df' data frame: ",
        change.base
      )
    }

    change.after <-
      unique(alpha_df %>% dplyr::select(all_of(c(time.var))))[unique(alpha_df %>% dplyr::select(all_of(c(time.var)))) != change.base]

    alpha_grouped <- alpha_df %>% dplyr::group_by(!!sym(time.var))
    alpha_split <- split(alpha_df, f = alpha_grouped[[time.var]])

    alpha_time_1 <- alpha_split[[change.base]]
    alpha_time_2 <- alpha_split[[change.after]]

    combined_alpha <- alpha_time_1 %>%
      dplyr::inner_join(
        alpha_time_2,
        by = c(subject.var, group.var),
        suffix = c("_time_1", "_time_2")
      )

    diff_columns <- lapply(alpha.name, function(index) {
      diff_col_name <- paste0(index, "_diff")

      if (is.function(alpha.change.func)) {
        combined_alpha <- combined_alpha %>%
          dplyr::mutate(!!diff_col_name := alpha.change.func(!!sym(paste0(
            index, "_time_2"
          )),!!sym(paste0(
            index, "_time_1"
          )))) %>%
          dplyr::select(all_of(diff_col_name))
      } else {
        if (alpha.change.func == "log fold change") {
          combined_alpha <- combined_alpha %>%
            dplyr::mutate(!!sym(diff_col_name) := log2(!!sym(paste0(
              index, "_time_2"
            )) / !!sym(paste0(
              index, "_time_1"
            )))) %>%
            dplyr::select(all_of(c(diff_col_name)))
        } else
          if (alpha.change.func == "absolute change") {
            combined_alpha <- combined_alpha %>%
              dplyr::mutate(!!diff_col_name := !!sym(paste0(index, "_time_2")) -!!sym(paste0(index, "_time_1"))) %>%
              dplyr::select(all_of(diff_col_name))
          } else {
            message(paste("No valid alpha.change.func provided for", index, ". Defaulting to 'absolute change'."))
            combined_alpha <- combined_alpha %>%
              dplyr::mutate(!!diff_col_name := !!sym(paste0(index, "_time_2")) -!!sym(paste0(index, "_time_1"))) %>%
              dplyr::select(all_of(diff_col_name))
          }
      }
    })

    combined_alpha <- dplyr::bind_cols(combined_alpha, diff_columns)

    col <- mStat_get_palette(palette)

    facet_formula <-
      if (!is.null(strata.var)) {
        paste(". ~", strata.var)
      } else {
        ". ~ 1"
      }

    if (is.null(group.var)) {
      combined_alpha$group <- "All"
    } else{
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_df %>% dplyr::select(all_of(c(
          subject.var, group.var
        )))
        ,
        by = c(subject.var, group.var)) %>% dplyr::rename(group = group.var)
    }

    if (!is.null(strata.var)) {
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_time_1 %>% dplyr::select(all_of(c(
          subject.var, strata.var
        )))
        , by = c(subject.var))
    }

    if (!is.null(adj.vars) &&
        (is.null(strata.var) || strata.var != adj.vars)) {
      combined_alpha <-
        combined_alpha %>% dplyr::left_join(alpha_time_1 %>% dplyr::select(all_of(c(
          subject.var, adj.vars
        )))
        , by = c(subject.var))
    }

    # Assuming mStat_get_theme function is already defined
    # Replace the existing theme selection code with this:
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    ylab_label <- if (is.function(alpha.change.func)) {
      base_label <-
        paste0("Change from ", change.base, " (custom function)")
    } else {
      base_label <-
        paste0("Change from ", change.base, " (", alpha.change.func, ")")
    }

    if (!is.null(adj.vars)) {
      covariates <- paste(adj.vars, collapse = ", ")
      ylab_label <-
        paste0(base_label, " (adjusted by: ", covariates, ")")
    } else {
      ylab_label <- base_label
    }

    plot_list <- lapply(alpha.name, function(index) {
      if (!is.null(adj.vars)) {
        # 对非数值型协变量进行因子转换
        data_subset <- combined_alpha %>%
          dplyr::select(all_of(adj.vars)) %>%
          dplyr::mutate(dplyr::across(where(~ is.character(.) & !is.factor(.)), factor))

        # 创建模型矩阵，并为非数值型协变量设定对比度
        M <- model.matrix(
          ~ 0 + .,
          data = data_subset,
          contrasts.arg = lapply(data_subset, stats::contrasts, contrasts = FALSE)
        )

        # Center the covariates (no scaling)
        M_centered <- scale(M, scale = FALSE)

        # Fit the regression model
        fit <-
          lm(combined_alpha[[paste0(index, "_diff")]] ~ M_centered)

        # 计算调整后的alpha多样性值
        adjusted_value <- fit$coefficients[1] + residuals(fit)

        # 在combined_alpha中更新alpha多样性值
        combined_alpha[[paste0(index, "_diff")]] <- adjusted_value

        # 显示消息，表示已经为特定的协变量调整了alpha多样性
        message(
          "Alpha diversity Change has been adjusted for the following covariates: ",
          paste(adj.vars, collapse = ", "),
          "."
        )
      }

      plot <-
        ggplot(combined_alpha, aes(
          x = group,
          y = !!sym(paste0(index, "_diff")),
          fill = group
        )) +
        geom_violin(trim = F, alpha = 0.8) +
        stat_boxplot(geom = "errorbar",
                     position = position_dodge(width = 0.2),
                     width = 0.1) +
        geom_boxplot(
          position = position_dodge(width = 0.8),
          width = 0.1,
          fill = "white"
        ) +
        geom_jitter(width = 0.1,
                    alpha = 0.5,
                    size = 1.5) +
        scale_fill_manual(values = col) +
        ylab(ylab_label) +
        theme_to_use +
        theme(
          panel.spacing.x = unit(0, "cm"),
          panel.spacing.y = unit(0.1, "cm"),
          strip.text.x = element_text(size = 15, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = base.size),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), units = "cm"),
          legend.text = ggplot2::element_text(size = 16),
          legend.title = ggplot2::element_text(size = 16)
        )

      if (any(unique(combined_alpha$group) == "All")) {
        plot <- plot +
          theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none"
          )
      } else if (is.null(strata.var)) {

      }
      else {
        plot <- plot +
          facet_wrap(as.formula(facet_formula), scales = "fixed")
      }

      return(plot)
    })

    change_func_label <- if (is.function(alpha.change.func)) {
      "custom_function"
    } else {
      alpha.change.func
    }

    # Save the plots as a PDF file
    if (pdf) {
      lapply(seq_along(plot_list), function(i) {
        plot <- plot_list[[i]]
        alpha_index <- alpha.name[i]

        pdf_name <- paste0(
          "alpha_change_boxplot_pair_",
          alpha_index,
          "_",
          "subject_",
          subject.var,
          "_",
          "time_",
          time.var,
          "_",
          "change_base_",
          change.base,
          "_",
          "change_func_",
          change_func_label
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

        ggsave(
          filename = pdf_name,
          plot = plot,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      })
    }

    names(plot_list) <- alpha.name

    return(plot_list)
  }
