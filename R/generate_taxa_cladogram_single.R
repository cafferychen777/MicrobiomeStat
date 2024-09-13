#' Generate a Circular Cladogram with Heatmap for Taxa
#'
#' This function generates a circular cladogram with an integrated heatmap for taxonomic data.
#' It visualizes the phylogenetic relationships between different taxa and their abundances or other
#' coefficients across different taxonomic levels using a tree-like structure (cladogram).
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which includes components like feature.tab, feature.ann, meta.dat, etc.
#' @param test.list A list of test results. If NULL, it will be generated using generate_taxa_test_single.
#' @param group.var The name of the grouping variable in meta.dat.
#' @param feature.level A character vector specifying taxonomic levels to be analyzed.
#' @param feature.mt.method Multiple testing method for features, "none" (default), "fdr", or other methods supported by p.adjust.
#' @param cutoff The p-value cutoff for significance.
#' @param color.group.level The taxonomic level used to color-code the branches of the cladogram.
#' @param pdf Boolean indicating whether to save the plot as a PDF.
#' @param pdf.width The width of the PDF file if saved.
#' @param pdf.height The height of the PDF file if saved.
#'
#' @return A ggplot object representing the circular heatmap with phylogenetic tree.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(subset_T2D.obj)
#'
#' test.list <- generate_taxa_test_single(
#'     data.obj = subset_T2D.obj,
#'     time.var = "visit_number",
#'     t.level = NULL,
#'     group.var = "subject_race",
#'     adj.vars = "subject_gender",
#'     feature.level = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
#'     feature.dat.type = "count",
#'     prev.filter = 0.1,
#'     abund.filter = 0.0001,
#' )
#'
#' plot.list <- generate_taxa_cladogram_single(
#'   data.obj = subset_T2D.obj,
#'   test.list = test.list,
#'   group.var = "subject_gender",
#'   feature.level = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
#'   feature.mt.method = "none",
#'   cutoff = 0.9,
#'   color.group.level = "Order"
#' )
#'
#' data(peerj32.obj)
#'
#' test.list <- generate_taxa_test_single(
#'     data.obj = peerj32.obj,
#'     time.var = "time",
#'     t.level = NULL,
#'     group.var = "group",
#'     adj.vars = "sex",
#'     feature.level = c("Phylum","Family","Genus"),
#'     feature.dat.type = "count",
#'     prev.filter = 0.1,
#'     abund.filter = 0.0001,
#' )
#'
#' plot.list <- generate_taxa_cladogram_single(
#'   data.obj = peerj32.obj,
#'   test.list = test.list,
#'   group.var = "group",
#'   feature.level = c("Phylum", "Family", "Genus"),
#'   cutoff = 0.3,
#'   color.group.level = "Family"
#' )
#' }
generate_taxa_cladogram_single <- function(
    data.obj,
    test.list,
    group.var = NULL,
    feature.level,
    feature.mt.method = "none",
    cutoff = 1,
    color.group.level = NULL,
    pdf = FALSE,
    pdf.width = 10,
    pdf.height = 10
) {

  # Process the data
  level_seq <- feature.level
  min_label <- level_seq[[length(level_seq)]]

  link_frame <- data.obj$feature.ann %>% as.data.frame()

  fix_frame <- function(inputframe, level_i) {
    colnames(inputframe)[[1]] <- level_i
    inputframe$Sites_layr <- level_i
    inputframe
  }

  join_frames <- function(test.list, level_seq, link_frame) {
    join_order <- rev(level_seq)
    link_frame <- link_frame[join_order]

    result <- data.frame()

    for (level in join_order) {
      for (comparison in names(test.list[[level]])) {
        current_frame <- fix_frame(test.list[[level]][[comparison]], level)

        needed_columns <- c(level, "Coefficient", "P.Value", "Adjusted.P.Value", "Sites_layr")
        current_frame <- current_frame[, needed_columns]

        current_frame <- dplyr::left_join(current_frame, link_frame, by = level)
        current_frame$Comparison <- comparison

        result <- dplyr::bind_rows(result, current_frame)
      }
    }

    return(result)
  }

  build_tree <- function(link_frame, level_seq) {
    link_frame <- as.data.frame(link_frame)
    link_frame <- link_frame[level_seq]
    frm <- as.formula(paste0("~", paste0(level_seq, collapse = "/")))
    for (i in 1:ncol(link_frame)) {
      link_frame[[i]] <- as.factor(link_frame[[i]])
    }
    treex <- ape::as.phylo(frm, data = link_frame, collapse = FALSE)
    treex
  }

  filter_h <- function(inputframe_linked, level_seq, feature.mt.method) {
    del_cutoff <- setNames(rep(cutoff, length(level_seq)), level_seq)
    for (i in 1:length(level_seq)) {
      level_i <- rev(level_seq)[[i]]
      if (level_i %in% names(del_cutoff)) {
        tmp_cut_off <- del_cutoff[[level_i]]
        if (feature.mt.method == "none") {
          inputframe_linked$Coefficient[inputframe_linked$Sites_layr == level_i & inputframe_linked$P.Value > tmp_cut_off] <- 0
        } else if (feature.mt.method == "fdr") {
          inputframe_linked$Coefficient[inputframe_linked$Sites_layr == level_i & inputframe_linked$Adjusted.P.Value > tmp_cut_off] <- 0
        }
      }
      if (level_i != min_label) {
        lower_i <- rev(level_seq)[[i-1]]
        keep_level_i <- inputframe_linked %>%
          dplyr::filter(Sites_layr == {{lower_i}} & Coefficient != 0) %>%
          pull(level_i) %>%
          unique()
        inputframe_linked$Coefficient[inputframe_linked$Sites_layr == level_i & (!inputframe_linked[[level_i]] %in% keep_level_i)] <- 0
      }
    }
    inputframe_linked
  }

  # Main processing
  fix_link_frame <- link_frame
  fix_link_frame[[min_label]] <- stringr::str_replace_all(fix_link_frame[[min_label]], pattern = " |\\(|\\)", replacement = "_")
  fix_link_frame[[min_label]] <- stringr::str_replace_all(fix_link_frame[[min_label]], pattern = "\\.", replacement = "")
  treex <- build_tree(fix_link_frame, level_seq)

  inputframe_linked <- join_frames(test.list, level_seq, link_frame)
  inputframe_linked$Variable <- stringr::str_replace_all(inputframe_linked[[min_label]], pattern = "\\.", replacement = "")
  inputframe_linked$Variable <- stringr::str_replace_all(inputframe_linked$Variable, pattern = " |\\(|\\)", replacement = "_")
  inputframe_linked$Sites_layr <- factor(inputframe_linked$Sites_layr, levels = level_seq)
  inputframe_linked <- filter_h(inputframe_linked, level_seq, feature.mt.method)

  sub_inputframe <- inputframe_linked %>% dplyr::filter(Sites_layr == {{color.group.level}})
  split_group <- split(sub_inputframe$Variable, f = sub_inputframe[[color.group.level]])
  treex <- ape::keep.tip(treex, sub_inputframe$Variable)
  treexx <- ggtree::groupOTU(treex, .node = split_group, group_name = color.group.level)

  # Generate plot
  # Custom color palette
  palette <- c(
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#FB8072", "#80B1D3", "#FDB462", "#BC80BD",
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#A6D854",
    "#1E88E5", "#D81B60", "#004D40", "#FFC107",
    "#5E35B1", "#00ACC1", "#3949AB", "#8E24AA",
    "#00897B", "#7CB342", "#C0CA33", "#FB8C00",
    "#6D4C41", "#546E7A", "#B71C1C", "#880E4F",
    "#4A148C", "#311B92", "#0D47A1", "#006064"
  )

  calculate_offset <- function(feature_level_length) {
    base_points <- list(
      c(1, 0.5),
      c(2, 0.8),
      c(3, 1.0),
      c(6, 1.7),
      c(9, 2.3)
    )

    if (feature_level_length <= base_points[[1]][1]) {
      return(base_points[[1]][2])
    }

    if (feature_level_length >= base_points[[length(base_points)]][1]) {
      return(base_points[[length(base_points)]][2])
    }

    for (i in 1:(length(base_points)-1)) {
      if (feature_level_length > base_points[[i]][1] && feature_level_length <= base_points[[i+1]][1]) {
        x1 <- base_points[[i]][1]
        y1 <- base_points[[i]][2]
        x2 <- base_points[[i+1]][1]
        y2 <- base_points[[i+1]][2]

        return(y1 + (feature_level_length - x1) * (y2 - y1) / (x2 - x1))
      }
    }
  }

  calculated_offset <- calculate_offset(length(feature.level))

  # Create a list to store individual plots
  plot.list <- list()

  for (comparison in unique(inputframe_linked$Comparison)) {
    comparison_data <- inputframe_linked %>% dplyr::filter(Comparison == comparison)

    p <- ggtree::ggtree(treexx, layout = "circular", open.angle = 5) +
      ggtreeExtra::geom_fruit(data = comparison_data, geom = geom_tile,
                              mapping = aes(y = Variable, x = Sites_layr,
                                            fill = Coefficient),
                              offset = 0.03, size = 0.02, color = "black") +
      ggtree::geom_tiplab(aes(color = .data[[color.group.level]]), offset = calculated_offset, align = TRUE,
                          linetype = "blank", size = 2, show.legend = FALSE) +
      scale_fill_gradient2(low = "#0571b0", high = "#ca0020") +
      geom_text(mapping = aes(label = "",
                              color = .data[[color.group.level]]),
                key_glyph = draw_key_rect) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8)
      ) +
      scale_color_manual(values = unique(palette)[1:length(unique(comparison_data[[color.group.level]]))]) +
      guides(
        fill = guide_colorbar(title = "Coefficient", barwidth = 10, barheight = 0.5),
        color = guide_legend(
          title = color.group.level,
          override.aes = list(size = 2),
          byrow = TRUE,
          keywidth = unit(0.5, "lines"),
          keyheight = unit(0.5, "lines")
        )
      )

    # Save as PDF if requested
    if (pdf) {
      filename <- paste0("taxa_cladogram_single_", gsub(" ", "_", comparison), ".pdf")
      ggsave(filename = filename, plot = p, width = pdf.width, height = pdf.height)
    }

    plot.list[[comparison]] <- p
  }

  return(plot.list)
}
