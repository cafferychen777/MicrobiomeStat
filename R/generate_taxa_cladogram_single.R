#' Generate a Circular Cladogram with Heatmap for Taxa
#'
#' This function generates a circular cladogram with an integrated heatmap for taxonomic data.
#' 
#' @importFrom rlang expr eval_tidy
#' @importFrom ggplot2 aes scale_fill_gradient2 geom_text element_text guide_colorbar unit scale_color_manual guides
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
#' @param palette An optional vector of colors to be used for the plot. If NULL, a default color palette will be used.
#' @param pdf Boolean indicating whether to save the plot as a PDF.
#' @param pdf.width The width of the PDF file if saved.
#' @param pdf.height The height of the PDF file if saved.
#' @param time.var Character string specifying the column name in metadata containing time variable. Used when test.list is NULL.
#' @param t.level Character string specifying the time level/value to subset data to. Used when test.list is NULL.
#' @param adj.vars Character vector specifying column names in metadata containing covariates. Used when test.list is NULL.
#' @param prev.filter Numeric value specifying the minimum prevalence threshold for filtering taxa. Used when test.list is NULL.
#' @param abund.filter Numeric value specifying the minimum abundance threshold for filtering taxa. Used when test.list is NULL.
#' @param feature.dat.type The type of the feature data, which determines data handling. Should be one of "count", "proportion", or "other". For CLR-transformed data, use "other". Used when test.list is NULL.
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
#' test.list <- generate_taxa_test_single(
#'     data.obj = subset_T2D.obj,
#'     time.var = "visit_number",
#'     t.level = NULL,
#'     group.var = "subject_race",
#'     adj.vars = "subject_gender",
#'     feature.level = c("Order"),
#'     feature.dat.type = "count",
#'     prev.filter = 0.1,
#'     abund.filter = 0.0001,
#' )
#'
#' plot.list <- generate_taxa_cladogram_single(
#'   data.obj = subset_T2D.obj,
#'   test.list = test.list,
#'   group.var = "subject_gender",
#'   feature.level = c("Order"),
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
    test.list = NULL,
    group.var = NULL,
    feature.level,
    feature.mt.method = "none",
    cutoff = 1,
    color.group.level = NULL,
    palette = NULL,
    pdf = FALSE,
    pdf.width = 10,
    pdf.height = 10,
    # Add other parameters that might be needed for generate_taxa_test_single
    time.var = NULL,
    t.level = NULL,
    adj.vars = NULL,
    prev.filter = 0.1,
    abund.filter = 0.0001,
    feature.dat.type = "count"
) {
  # Ensure necessary packages are available
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    stop("Package 'ggtree' is required but not installed.")
  }
  if (!requireNamespace("ggtreeExtra", quietly = TRUE)) {
    stop("Package 'ggtreeExtra' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  # Ensure required namespaces are available
  requireNamespace("ggplot2", quietly = TRUE)

  # Check group.var parameter
  if (is.null(group.var)) {
    stop("'group.var' must be provided")
  }

  # If test.list is NULL, generate it using generate_taxa_test_single
  if (is.null(test.list)) {
    message("Generating test results using generate_taxa_test_single...")
    test.list <- generate_taxa_test_single(
      data.obj = data.obj,
      time.var = time.var,
      t.level = t.level,
      group.var = group.var,
      adj.vars = adj.vars,
      feature.level = feature.level,
      feature.dat.type = feature.dat.type,
      prev.filter = prev.filter,
      abund.filter = abund.filter
    )
    
    # Check if test results were successfully generated
    if (length(test.list) == 0 || all(sapply(test.list, length) == 0)) {
      stop("Failed to generate test results. Please check your input parameters.")
    }
  }

  # If color.group.level is not specified, use the first feature.level
  if (length(color.group.level) == 0){
    color.group.level <- feature.level[1]  # Use the first feature.level as the default value
  }

  # Process the data
  # Extract the sequence of taxonomic levels to be analyzed
  level_seq <- feature.level
  # Identify the most specific taxonomic level (e.g., Species)
  min_label <- level_seq[[length(level_seq)]]

  # Extract the feature annotation data from the input object
  link_frame <- data.obj$feature.ann %>% as.data.frame()

  # Function to standardize the format of input frames
  fix_frame <- function(inputframe, level_i) {
    # Rename the first column to the current taxonomic level
    colnames(inputframe)[[1]] <- level_i
    # Add a column to indicate the taxonomic level
    inputframe$Sites_layr <- level_i
    inputframe
  }

  # Function to join test results with feature annotations
  join_frames <- function(test.list, level_seq, link_frame) {
    # Reverse the order of taxonomic levels for processing
    join_order <- rev(level_seq)
    link_frame <- link_frame[join_order]

    result <- data.frame()

    # Iterate through each taxonomic level and comparison
    for (level in join_order) {
      for (comparison in names(test.list[[level]])) {
        current_frame <- fix_frame(test.list[[level]][[comparison]], level)

        # Select relevant columns
        needed_columns <- c(level, "Coefficient", "P.Value", "Adjusted.P.Value", "Sites_layr")
        current_frame <- current_frame[, needed_columns]

        # Join with feature annotations
        current_frame <- dplyr::left_join(current_frame, link_frame, by = level)
        current_frame$Comparison <- comparison

        # Combine results
        result <- dplyr::bind_rows(result, current_frame)
      }
    }

    return(result)
  }

  #' Get appropriate phylogenetic tree for cladogram
  #'
  #' This function determines the appropriate phylogenetic tree to use for the cladogram.
  #' It first checks if a valid tree exists in the data object, and if so, ensures the
  #' tree tips match the feature names. If no valid tree is found or no matching tips
  #' are found, it builds a taxonomy-based tree.
  #'
  #' @param data.obj A MicrobiomeStat data object
  #' @param fix_link_frame Processed feature annotation data frame
  #' @param min_label The most specific taxonomic level
  #' @param level_seq Vector of taxonomic levels in order
  #' @param verbose Logical, whether to print detailed messages
  #'
  #' @return A phylogenetic tree object
  #' Process Unclassified taxonomic labels
  #'
  #' This function processes "Unclassified" labels in a data frame, ensuring they contain
  #' taxonomic level information by appending the taxonomic level and a unique identifier.
  #'
  #' @param data_frame A data frame containing taxonomic data
  #' @param group_col Column name containing taxonomic level information
  #' @param target_col Column name containing labels to process (default: "Variable")
  #' @param use_grouping Logical, whether to group by group_col before processing (default: TRUE)
  #'
  #' @return A data frame with processed "Unclassified" labels
  process_unclassified_labels <- function(data_frame, group_col, target_col = "Variable", use_grouping = TRUE) {
    # Convert column names to symbols for non-standard evaluation
    group_sym <- rlang::sym(group_col)
    target_sym <- rlang::sym(target_col)
    
    # Define the processing logic
    process_fn <- function(df) {
      df %>%
        dplyr::mutate(!!target_sym := dplyr::case_when(
          !!target_sym == "Unclassified" ~ paste0("Unclassified_", !!group_sym, "_", 
                                                dplyr::row_number()),
          TRUE ~ !!target_sym
        ))
    }
    
    # Apply the processing with or without grouping
    if (use_grouping) {
      data_frame %>%
        dplyr::group_by(!!group_sym) %>%
        process_fn() %>%
        dplyr::ungroup()
    } else {
      process_fn(data_frame)
    }
  }
  
  get_phylogenetic_tree <- function(data.obj, fix_link_frame, min_label, level_seq, verbose = TRUE) {
    # Check if tree exists in data object
    if (is.null(data.obj$tree)) {
      if (verbose) message("No tree found in data object. Building taxonomy-based tree.")
      return(build_tree(fix_link_frame, level_seq))
    }
    
    # Tree exists, check if it's valid
    tree <- data.obj$tree
    if (!inherits(tree, "phylo")) {
      if (verbose) message("Tree in data object is not a valid phylogenetic tree. Building taxonomy-based tree.")
      return(build_tree(fix_link_frame, level_seq))
    }
    
    # Check for matching tips
    common_tips <- intersect(tree$tip.label, fix_link_frame[[min_label]])
    total_tips <- length(tree$tip.label)
    total_features <- length(unique(fix_link_frame[[min_label]]))
    match_percent_tree <- round(length(common_tips) / total_tips * 100, 1)
    match_percent_features <- round(length(common_tips) / total_features * 100, 1)
    
    if (length(common_tips) > 0) {
      if (verbose) {
        message(sprintf("Using phylogenetic tree from data object:"))
        message(sprintf("- Matched %d of %d tree tips (%.1f%%)", 
                       length(common_tips), total_tips, match_percent_tree))
        message(sprintf("- Matched %d of %d features (%.1f%%)", 
                       length(common_tips), total_features, match_percent_features))
      }
      return(ape::keep.tip(tree, common_tips))
    } else {
      if (verbose) {
        message("No matching tips found between phylogenetic tree and features.")
        message(sprintf("- Tree has %d tips, features has %d unique values", 
                       total_tips, total_features))
        message("Using taxonomy-based tree instead.")
      }
      return(build_tree(fix_link_frame, level_seq))
    }
  }
  
  # Function to build phylogenetic tree from feature annotations
  build_tree <- function(link_frame, level_seq) {
    link_frame <- as.data.frame(link_frame)
    link_frame <- link_frame[level_seq]
    # Create a formula for tree construction
    frm <- as.formula(paste0("~", paste0(level_seq, collapse = "/")))
    # Convert all columns to factors
    for (i in 1:ncol(link_frame)) {
      link_frame[[i]] <- as.factor(link_frame[[i]])
    }
    # Generate phylogenetic tree
    treex <- ape::as.phylo(frm, data = link_frame, collapse = FALSE)
    treex
  }

  # Function to filter results based on statistical significance
  filter_h <- function(inputframe_linked, level_seq, feature.mt.method) {
    # Set cutoff values for each taxonomic level
    del_cutoff <- setNames(rep(cutoff, length(level_seq)), level_seq)
    for (i in 1:length(level_seq)) {
      level_i <- rev(level_seq)[[i]]
      if (level_i %in% names(del_cutoff)) {
        tmp_cut_off <- del_cutoff[[level_i]]
        # Apply cutoff based on the multiple testing method
        if (feature.mt.method == "none") {
          # Use raw p-values if no multiple testing correction
          inputframe_linked$Coefficient[inputframe_linked$Sites_layr == level_i & inputframe_linked$P.Value > tmp_cut_off] <- 0
        } else if (feature.mt.method == "fdr") {
          # Use FDR-adjusted p-values for multiple testing correction
          inputframe_linked$Coefficient[inputframe_linked$Sites_layr == level_i & inputframe_linked$Adjusted.P.Value > tmp_cut_off] <- 0
        }
      }
      # Propagate filtering to higher taxonomic levels
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
  # Prepare feature annotations for tree construction
  fix_link_frame <- link_frame %>%
    dplyr::mutate(
      !!rlang::sym(min_label) := stringr::str_replace_all(
        !!rlang::sym(min_label), 
        pattern = " |\\(|\\)", 
        replacement = "_"
      ) %>%
      stringr::str_replace_all(
        pattern = "\\." , 
        replacement = ""
      )
    )

  # Get the appropriate phylogenetic tree using the helper function
  treex <- get_phylogenetic_tree(
    data.obj = data.obj,
    fix_link_frame = fix_link_frame,
    min_label = min_label,
    level_seq = level_seq,
    verbose = TRUE
  )

  # Join and process test results
  inputframe_linked <- join_frames(test.list, level_seq, link_frame) %>%
    dplyr::mutate(
      Variable = stringr::str_replace_all(!!rlang::sym(min_label), pattern = "\\.", replacement = ""),
      Variable = stringr::str_replace_all(Variable, pattern = " |\\(|\\)", replacement = "_"),
      Sites_layr = factor(Sites_layr, levels = level_seq)
    )
  # Apply statistical filtering
  inputframe_linked <- filter_h(inputframe_linked, level_seq, feature.mt.method)

  # Subset data for the chosen color grouping level
  sub_inputframe <- inputframe_linked %>% dplyr::filter(Sites_layr == {{color.group.level}})

  # Process "Unclassified" labels in both data frames
  inputframe_linked <- process_unclassified_labels(
    inputframe_linked, 
    group_col = "Sites_layr", 
    use_grouping = TRUE
  )
  
  sub_inputframe <- process_unclassified_labels(
    sub_inputframe, 
    group_col = "Sites_layr", 
    use_grouping = FALSE
  )

  # Ensure tree labels match data labels
  common_labels <- intersect(treex$tip.label, sub_inputframe$Variable)
  treex <- ape::keep.tip(treex, common_labels)
  sub_inputframe <- sub_inputframe[sub_inputframe$Variable %in% common_labels, ]

  # Group tree nodes by taxonomic level for coloring
  split_group <- split(sub_inputframe$Variable, f = sub_inputframe[[color.group.level]])
  treexx <- ggtree::groupOTU(treex, .node = split_group, group_name = color.group.level)

  # Generate plot
  # Define custom color palette if not provided
  if (is.null(palette)){
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
  }

  # Function to calculate offset for tip labels based on the number of taxonomic levels
  calculate_offset <- function(feature_level_length) {
    if (feature_level_length == 1) {
      return(0.2)
    }

    # Base points for interpolation: (num_levels, offset)
    base_points <- list(
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

    # Linear interpolation between base points
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

  # Generate a plot for each comparison
  for (comparison in unique(inputframe_linked$Comparison)) {
    comparison_data <- inputframe_linked %>% dplyr::filter(Comparison == comparison)

    # Ensure ggplot2 namespace is available in this scope
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required but not installed.")
    }
    
    # Ensure geom_tile function is available in the global environment
    # ggtreeExtra::geom_fruit needs to find this function in the global environment
    if (!exists("geom_tile", envir = .GlobalEnv)) {
      assign("geom_tile", ggplot2::geom_tile, envir = .GlobalEnv)
    }
    
    # Create the circular cladogram plot
    p <- ggtree::ggtree(treexx, layout = "circular", open.angle = 5) +
      # Use the string "geom_tile" as the value for the geom parameter
      # This complies with the requirements of the ggtreeExtra::geom_fruit function
      ggtreeExtra::geom_fruit(
        data = comparison_data,
        geom = "geom_tile",  # Use a string, not a function object or unevaluated symbol
        mapping = aes(y = Variable, x = Sites_layr, fill = Coefficient),
        offset = 0.03,
        size = 0.02,
        color = "black"
      ) +
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

    # Add the plot to the list
    plot.list[[comparison]] <- p
  }

  # Return the list of plots
  return(plot.list)
}