#' Generate Beta Diversity Ordination Plot
#'
#' This function generates ordination plots (Principal Component Analysis) based on
#' beta-diversity distances. It also allows for stratification and grouping of samples,
#' and calculation of distances at a specific time point.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param dist.obj Distance matrix between samples, usually calculated using
#' \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} function.
#' If NULL, beta diversity will be automatically computed from \code{data.obj}
#' using \code{mStat_calculate_beta_diversity}.
#' @param pc.obj A list containing the results of dimension reduction/Principal Component Analysis.
#' This should be the output from functions like \code{\link[MicrobiomeStat]{mStat_calculate_PC}},
#' containing the PC coordinates and other metadata. If NULL (default), dimension reduction
#' will be automatically performed using metric multidimensional scaling (MDS) via
#' \code{\link[MicrobiomeStat]{mStat_calculate_PC}}. The pc.obj list structure should contain:
#' \describe{
#'   \item{points}{A matrix with samples as rows and PCs as columns containing the coordinates.}
#'   \item{eig}{Eigenvalues for each PC dimension.}
#'   \item{vectors}{Loadings vectors for features onto each PC.}
#'   \item{Other metadata}{like method, dist.name, etc.}
#' }
#' See \code{\link[MicrobiomeStat]{mStat_calculate_PC}} function for details on output format.
#' @param time.var String. Variable to be used for time. Default is NULL.
#' @param t.level Character string specifying the time level/value to subset data to,
#' if a time variable is provided. Default NULL does not subset data.
#' @param group.var String. Variable to be used for grouping. Default is NULL.
#' @param strata.var String. Variable to be used for stratification. Default is NULL.
#' @param adj.vars A character vector of variable names to be used for adjustment.
#' @param dist.name A character vector specifying which beta diversity indices to calculate. Supported indices are "BC" (Bray-Curtis), "Jaccard", "UniFrac" (unweighted UniFrac), "GUniFrac" (generalized UniFrac), "WUniFrac" (weighted UniFrac), and "JS" (Jensen-Shannon divergence). If a name is provided but the corresponding object does not exist within dist.obj, it will be computed internally. If the specific index is not supported, an error message will be returned. Default is c('BC', 'Jaccard').
#' @param base.size Numeric. Base size for plot elements. Default is 16.
#' @param point.size Numeric. The size of the points in the scatter plot. Default is 4.
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
#' @param pdf Logical. If TRUE, the plots are saved as PDF files. Default is TRUE.
#' @param file.ann File annotation. Default is NULL.
#' @param pdf.wid Width of the PDF. Default is 11.
#' @param pdf.hei Height of the PDF. Default is 8.5.
#' @param ... Additional arguments to be passed.
#'
#' @details The function is flexible and allows for various modifications, including the choice of distance measure and stratification factor, providing a comprehensive tool for microbiome beta diversity exploration. It integrates well with other MicrobiomeStat functions and takes their output as input.
#'
#' @return A PCoA plot displaying the beta diversity ordination, stratified by the specified grouping and/or strata variables (if provided). The plot will be saved as a PDF if `pdf` is set to `TRUE`.
#'
#' @seealso \code{\link[MicrobiomeStat]{mStat_calculate_beta_diversity}} for creating the distance object, \code{\link[MicrobiomeStat]{mStat_calculate_PC}} for computing the principal coordinates, and \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{geom_boxplot}} for the underlying plot functions used, and \code{\link[MicrobiomeStat]{mStat_convert_DGEList_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_DESeqDataSet_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_phyloseq_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_convert_SummarizedExperiment_to_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_qiime2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_mothur_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_dada2_as_data_obj}}, \code{\link[MicrobiomeStat]{mStat_import_biom_as_data_obj}} for data conversion.
#'
#' @author Caffery Yang \email{cafferychen7850@@gmail.com}
#'
#' @examples
#' \dontrun{
#' library(aplot)
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC','Jaccard'))
#' generate_beta_ordination_single(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   pc.obj = NULL,
#'   time.var = "time",
#'   t.level = "2",
#'   group.var = "group",
#'   strata.var = "sex",
#'   adj.vars = "sex",
#'   dist.name = c("BC"),
#'   base.size = 20,
#'   point.size = 4,
#'   theme.choice = "bw",
#'   custom.theme = NULL,
#'   palette = NULL,
#'   pdf = TRUE,
#'   file.ann = NULL,
#'   pdf.wid = 11,
#'   pdf.hei = 8.5
#' )
#'
#' data("subset_T2D.obj")
#' dist.obj <- mStat_calculate_beta_diversity(subset_T2D.obj, dist.name = c('BC', 'Jaccard'))
#' pc.obj <- mStat_calculate_PC(dist.obj, method = c('mds'), k = 2, dist.name = c('BC','Jaccard'))
#' generate_beta_ordination_single(
#'   data.obj = subset_T2D.obj,
#'   dist.obj = dist.obj,
#'   pc.obj = pc.obj,
#'   time.var = "visit_number_num",
#'   t.level = NULL,
#'   group.var = "subject_race",
#'   strata.var = "subject_gender",
#'   adj.vars = "sample_body_site",
#'   dist.name = c("BC", 'Jaccard'),
#'   base.size = 20,
#'   point.size = 4,
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
generate_beta_ordination_single <-
  function(data.obj,
           time.var = NULL,
           t.level = NULL,
           group.var = NULL,
           adj.vars = NULL,
           strata.var = NULL,
           dist.obj = NULL,
           dist.name = c('BC', 'Jaccard'),
           pc.obj = NULL,
           base.size = 16,
           point.size = 4,
           theme.choice = "bw",
           custom.theme = NULL,
           palette = NULL,
           pdf = TRUE,
           file.ann = NULL,
           pdf.wid = 11,
           pdf.hei = 8.5,
           ...) {

    # Check if distance metrics are provided. If not, exit the function.
    if (is.null(dist.name)){
      return()
    }

    # This block handles the calculation or retrieval of distance matrices and metadata
    if (is.null(dist.obj)) {
      # Process data based on time variable if it exists
      if (!is.null(time.var)){
        if (!is.null(t.level)){
          # Subset data to specific time point if t.level is provided
          condition <- paste(time.var, "== '", t.level, "'", sep = "")
          data.obj <- mStat_subset_data(data.obj, condition = condition)
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
          dist.obj <-
            mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        } else {
          # If no specific time point is provided, use all time points but warn if multiple exist
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
          if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
            message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
          }
          dist.obj <-
            mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
        }
      } else {
        # If no time variable is provided, proceed with all data
        meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
        dist.obj <-
          mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      }
      # Adjust distances if adjustment variables are provided
      if (!is.null(adj.vars)){
        dist.obj <- mStat_calculate_adjusted_distance(data.obj = data.obj, dist.obj = dist.obj, adj.vars = adj.vars, dist.name = dist.name)
      }
    } else {
      # If distance object is provided, extract metadata accordingly
      if (!is.null(data.obj)){
        if (!is.null(time.var)){
          if (!is.null(t.level)){
            condition <- paste(time.var, "== '", t.level, "'", sep = "")
            data.obj <- mStat_subset_data(data.obj, condition = condition)
            meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
          } else {
            meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
            if (length(levels(as.factor(meta_tab[,time.var]))) != 1){
              message("Multiple time points detected in your dataset. It is recommended to either set t.level or utilize functions for longitudinal data analysis.")
            }
          }
        } else {
          meta_tab <- data.obj$meta.dat %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
        }
      }
      if (!is.null(attr(dist.obj[[dist.name[1]]], "labels"))){
        meta_tab <- attr(dist.obj[[dist.name[1]]], "labels")  %>% dplyr::select(all_of(c(group.var,strata.var,time.var)))
      }
    }

    # Perform dimension reduction if not already done
    # This step reduces the high-dimensional distance data to 2D for visualization
    if (is.null(pc.obj)) {
      pc.obj <-
        mStat_calculate_PC(
          dist.obj = dist.obj,
          method = "mds",
          k = 2,
          dist.name = dist.name
        )
    }

    # Get color palette for plotting
    col <- mStat_get_palette(palette)

    # If no group variable is provided, create a dummy "ALL" group
    if (is.null(group.var)){
      group.var = "ALL"
      meta_tab$ALL <- "ALL"
    }

    # Define aesthetic mapping based on presence of strata variable
    aes_function <- if (!is.null(strata.var)) {
      aes(color = !!sym(group.var),
          shape = !!sym(strata.var))
    } else {
      aes(color = !!sym(group.var))
    }

    # Get appropriate theme for plotting
    theme_to_use <- mStat_get_theme(theme.choice, custom.theme)

    # Generate plots for each distance metric
    plot_list <- lapply(dist.name, function(dist.name) {

      # Extract the first two principal coordinates
      pc.mat <- pc.obj[[dist.name]]$points[, 1:2]

      # Prepare data frame for plotting
      df <-
        pc.mat %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        dplyr::left_join(meta_tab %>%
                           dplyr::select(all_of(c(time.var, group.var, strata.var))) %>%
                           rownames_to_column("sample"), by = "sample")

      # Filter out NA values in time variable if it exists
      if (!is.null(time.var)){
        df <- df %>% dplyr::filter(!is.na(!!sym(time.var)))
      }

      df <- df %>% column_to_rownames("sample")

      colnames(df)[1:2] <- c("PC1", "PC2")

      # Create the main scatter plot
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2)) +
        ggplot2::geom_point(size = point.size, aes_function, show.legend = T) +
        ggplot2::labs(
          x = ifelse(!is.null(pc.obj[[dist.name]]$eig),paste0("Axis 1 (", round(pc.obj[[dist.name]]$eig[1]/sum(pc.obj[[dist.name]]$eig)*100,2),"%)"),"Axis 1"),
          y = ifelse(!is.null(pc.obj[[dist.name]]$eig),paste0("Axis 2 (", round(pc.obj[[dist.name]]$eig[2]/sum(pc.obj[[dist.name]]$eig)*100,2),"%)"),"Axis 2")
        ) +
        # Add 95% confidence ellipses for each group
        ggplot2::stat_ellipse(ggplot2::aes(color = !!sym(group.var)),fill="white",geom = "polygon",
                              level=0.95,alpha = 0.01,show.legend = F) +
        ggplot2::geom_vline(
          xintercept = 0,
          linetype = "dashed",
          color = "black"
        ) +
        ggplot2::geom_hline(
          yintercept = 0,
          linetype = "dashed",
          color = "black"
        ) +
        theme_to_use  +
        ggplot2::theme(
          axis.line.x = ggplot2::element_line(size = 1, colour = "black"),
          axis.line.y = ggplot2::element_line(size = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black"),
          axis.title = ggplot2::element_text(color = "black", size = 20),
          axis.text.x = element_text(color = "black", size = base.size),
          axis.text.y = element_text(color = "black", size = base.size),
          axis.title.x = element_text(size = base.size),
          axis.title.y = element_text(size = base.size),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(
            color = "black",
            size = base.size
          ),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = ggplot2::element_text(size = base.size),
          legend.title = ggplot2::element_text(size = base.size)
        )

      # Set color scale based on grouping
      if (group.var == "ALL") {
        p <- p + scale_color_manual(values = col, guide = "none")
      } else {
        p <- p + scale_color_manual(values = col)
      }

      # Create a boxplot for PC1 values
      Fig1a.taxa.pc1.boxplot <-
        ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x=!!sym(group.var), y=PC1, fill=!!sym(group.var)), color="black", alpha=0.5, show.legend = F) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand = c(0,0.001)) +
        ggplot2::labs(x=NULL, y=NULL) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.text.y =ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())+
        ggplot2::coord_flip()

      # Create a boxplot for PC2 values
      Fig1a.taxa.pc2.boxplot <-
        ggplot2::ggplot(df) +
        ggplot2::geom_boxplot(ggplot2::aes(x=!!sym(group.var), y=PC2, fill=!!sym(group.var)), color="black", alpha=0.5, show.legend = F) +
        ggplot2::scale_fill_manual(values=col) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(expand = c(0,0.001)) +
        ggplot2::labs(x=NULL, y=NULL) +
        ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       axis.text.y =ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())

      # Combine the main plot with the boxplots
      # This provides a comprehensive view of the data distribution along both principal coordinates
      p <- p %>%
        aplot::insert_top(Fig1a.taxa.pc1.boxplot, height = 0.2) %>%
        aplot::insert_right(Fig1a.taxa.pc2.boxplot, width=0.2) %>%
        as.ggplot()

      # Save the plot as a PDF file if requested
      if (pdf) {
        pdf_name <- paste0("beta_ordination_single_", "dist.name_", dist.name)
        if (!is.null(time.var)) {
          pdf_name <- paste0(pdf_name, "_", "time_", time.var)
        }
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
          plot = p,
          width = pdf.wid,
          height = pdf.hei,
          dpi = 300
        )
      }
      return(p)
    })

    # Assign names to the elements of the plot list
    names(plot_list) <- dist.name
    return(plot_list)
  }