#' Generate a Circular Cladogram with Heatmap for Taxa
#'
#' This function generates a circular cladogram with an integrated heatmap for taxonomic data.
#' It visualizes the phylogenetic relationships between different taxa and their abundances or other
#' coefficients across different taxonomic levels using a tree-like structure (cladogram).
#' 
#' @importFrom rlang expr eval_tidy
#' @importFrom ggplot2 aes scale_fill_gradient2 geom_text element_text guide_colorbar unit scale_color_manual guides
#' @importFrom tidytree as_tibble as.treedata
#' @importFrom dplyr filter mutate group_by ungroup
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
  taxonomic_hierarchy <- feature.level
  # Identify the most specific taxonomic level (e.g., Species)
  finest_taxonomic_level <- taxonomic_hierarchy[[length(taxonomic_hierarchy)]]

  # Extract the feature annotation data from the input object
  taxonomy_annotation <- data.obj$feature.ann %>% as.data.frame()

  # Function to standardize the format of input frames
  standardize_test_frame <- function(inputframe, level_i) {
    # Rename the first column to the current taxonomic level
    colnames(inputframe)[[1]] <- level_i
    # Add a column to indicate the taxonomic level
    inputframe$taxonomic_level <- level_i
    inputframe
  }

  # Function to join test results with feature annotations
  merge_test_with_taxonomy <- function(test.list, taxonomic_hierarchy, taxonomy_annotation) {
    # Reverse the order of taxonomic levels for processing
    join_order <- rev(taxonomic_hierarchy)
    processed_taxonomy <- taxonomy_annotation[join_order]

    result <- data.frame()

    # Iterate through each taxonomic level and comparison
    for (level in join_order) {
      for (comparison in names(test.list[[level]])) {
        current_frame <- standardize_test_frame(test.list[[level]][[comparison]], level)

        # Select relevant columns
        needed_columns <- c(level, "Coefficient", "P.Value", "Adjusted.P.Value", "taxonomic_level")
        current_frame <- current_frame[, needed_columns]

        # Join with feature annotations
        current_frame <- dplyr::left_join(current_frame, processed_taxonomy, by = level)
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
  # Internal function parameters:
  # @param data.obj A MicrobiomeStat data object
  # @param fix_link_frame Processed feature annotation data frame
  # @param min_label The most specific taxonomic level
  # @param level_seq Vector of taxonomic levels in order
  # @param verbose Logical, whether to print detailed messages
  #'
  #' @return A phylogenetic tree object
  # Process Unclassified taxonomic labels
  #
  # This function processes "Unclassified" labels in a data frame, ensuring they contain
  # taxonomic level information by appending the taxonomic level and a unique identifier.
  #
  # @param data_frame A data frame containing taxonomic data
  # @param group_col Column name containing taxonomic level information
  # @param target_col Column name containing labels to process (default: "Variable")
  # @param use_grouping Logical, whether to group by group_col before processing (default: TRUE)
  #
  # @return A data frame with processed "Unclassified" labels
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

  # Function to get the appropriate phylogenetic tree
  get_phylogenetic_tree <- function(data.obj, fix_link_frame, min_label, level_seq, verbose = TRUE) {
    # Check if a tree is available in the data object
    if (!is.null(data.obj$tree)) {
      if (verbose) message("Found tree in data object. Checking compatibility...")
      tree <- data.obj$tree
      
      # Check if the tree is already in phylo format
      if (inherits(tree, "phylo")) {
        # Get the leaf node names from the tree
        tree_tips <- tree$tip.label
        # Get the feature names from the processed taxonomy
        feature_names <- fix_link_frame[[min_label]]
        
        # Find matching nodes between tree and features
        matching_nodes <- intersect(tree_tips, feature_names)
        matching_count <- length(matching_nodes)
        
        # Calculate match percentages
        tree_match_percent <- round(matching_count / length(tree_tips) * 100, 2)
        feature_match_percent <- round(matching_count / length(feature_names) * 100, 2)
        
        if (verbose) {
          message(sprintf("Tree has %d leaf nodes, %d match with features (%.2f%%)", 
                         length(tree_tips), matching_count, tree_match_percent))
          message(sprintf("Data has %d features, %d match with tree (%.2f%%)", 
                         length(feature_names), matching_count, feature_match_percent))
        }
        
        # If there are matching nodes, use the tree from the data object
        if (matching_count > 0) {
          if (verbose) message("Using phylogenetic tree from data object.")
          # Keep only the tips that match with features
          tree <- ape::keep.tip(tree, matching_nodes)
          return(tree)
        } else {
          if (verbose) message("No matching nodes found between tree and features.")
        }
      } else {
        if (verbose) message("Tree in data object is not in phylo format.")
      }
      
      if (verbose) {
        message("Using taxonomy-based tree instead.")
      }
      return(build_taxonomy_tree(taxonomy_annotation, taxonomic_hierarchy))
    }
    return(build_taxonomy_tree(taxonomy_annotation, taxonomic_hierarchy))
  }
  
  # Function to build phylogenetic tree from feature annotations
  build_taxonomy_tree <- function(taxonomy_annotation, taxonomic_hierarchy) {
    taxonomy_annotation <- as.data.frame(taxonomy_annotation)
    taxonomy_annotation <- taxonomy_annotation[taxonomic_hierarchy]
    # Create a formula for tree construction
    frm <- as.formula(paste0("~", paste0(taxonomic_hierarchy, collapse = "/")))
    # Convert all columns to factors
    for (i in seq_len(ncol(taxonomy_annotation))) {
      taxonomy_annotation[[i]] <- as.factor(taxonomy_annotation[[i]])
    }
    # Generate phylogenetic tree
    phylogenetic_tree <- ape::as.phylo(frm, data = taxonomy_annotation, collapse = FALSE)
    
    # Check if tidytree is available
    if (requireNamespace("tidytree", quietly = TRUE)) {
      # Convert phylo object to tibble
      tree_tibble <- tidytree::as_tibble(phylogenetic_tree)
      
      # Check if the tibble contains the necessary columns
      if (!all(c("parent", "node") %in% colnames(tree_tibble))) {
        # If parent column is missing, add it
        if (!("parent" %in% colnames(tree_tibble))) {
          # Create parent column based on tree structure
          # For root node, parent is NA
          # For other nodes, parent is the node number of their parent
          tree_tibble$parent <- NA
          # Fill parent column using edge matrix
          if (!is.null(phylogenetic_tree$edge)) {
            edge_df <- as.data.frame(phylogenetic_tree$edge)
            colnames(edge_df) <- c("parent", "node")
            
            # Find parent for each node
            for (i in seq_len(nrow(edge_df))) {
              node_idx <- which(tree_tibble$node == edge_df$node[i])
              if (length(node_idx) > 0) {
                tree_tibble$parent[node_idx] <- edge_df$parent[i]
              }
            }
          }
        }
        
        # If node column is missing, add it
        if (!("node" %in% colnames(tree_tibble))) {
          # Create node numbers
          tree_tibble$node <- seq_len(nrow(tree_tibble))
        }
      }
      
      # Convert processed tibble back to treedata object
      tree_data <- tidytree::as.treedata(tree_tibble)
      
      # Return phylo object
      return(tree_data@phylo)
    } else {
      # If tidytree is not available, return the original phylo object
      return(phylogenetic_tree)
    }
  }

  # Function to filter results based on statistical significance
  filter_by_significance <- function(inputframe_linked, taxonomic_hierarchy, feature.mt.method, finest_taxonomic_level, cutoff_value = cutoff) {
    # Define significance thresholds
    significance_thresholds <- setNames(rep(cutoff_value, length(taxonomic_hierarchy)), taxonomic_hierarchy)
    for (i in seq_along(taxonomic_hierarchy)) {
      level_i <- rev(taxonomic_hierarchy)[[i]]
      if (level_i %in% names(significance_thresholds)) {
        tmp_cut_off <- significance_thresholds[[level_i]]
        # Apply cutoff based on the multiple testing method
        if (feature.mt.method == "none") {
          # Use raw p-values if no multiple testing correction
          inputframe_linked$Coefficient[inputframe_linked$.data$taxonomic_level == level_i & inputframe_linked$P.Value > tmp_cut_off] <- 0
        } else if (feature.mt.method == "fdr") {
          # Use FDR-adjusted p-values for multiple testing correction
          inputframe_linked$Coefficient[inputframe_linked$.data$taxonomic_level == level_i & inputframe_linked$Adjusted.P.Value > tmp_cut_off] <- 0
        }
      }
      # Propagate filtering to higher taxonomic levels
      if (level_i != finest_taxonomic_level) {
        lower_i <- rev(taxonomic_hierarchy)[[i-1]]
        keep_level_i <- inputframe_linked %>%
          dplyr::filter(.data$taxonomic_level == {{lower_i}} & .data$Coefficient != 0) %>%
          dplyr::pull(level_i) %>%
          unique()
        inputframe_linked$Coefficient[inputframe_linked$.data$taxonomic_level == level_i & (!inputframe_linked[[level_i]] %in% keep_level_i)] <- 0
      }
    }
    inputframe_linked
  }

  # Main processing
  # Prepare feature annotations for tree construction
  processed_taxonomy <- taxonomy_annotation
  processed_taxonomy[[finest_taxonomic_level]] <- stringr::str_replace_all(processed_taxonomy[[finest_taxonomic_level]], pattern = " |\\(|\\)", replacement = "_") %>%
    stringr::str_replace_all(pattern = "\\.", replacement = "")

  # Get the appropriate phylogenetic tree using the helper function
  phylogenetic_tree <- get_phylogenetic_tree(
    data.obj = data.obj,
    fix_link_frame = processed_taxonomy,
    min_label = finest_taxonomic_level,
    level_seq = taxonomic_hierarchy,
    verbose = TRUE
  )

  # Join and process test results
  merged_test_data <- merge_test_with_taxonomy(test.list, taxonomic_hierarchy, taxonomy_annotation) %>%
    dplyr::mutate(
      Variable = stringr::str_replace_all(!!rlang::sym(finest_taxonomic_level), pattern = "\\.", replacement = ""),
      Variable = stringr::str_replace_all(Variable, pattern = " |\\(|\\)", replacement = "_"),
      taxonomic_level = factor(taxonomic_level, levels = taxonomic_hierarchy)
    )
  # Apply statistical filtering
  merged_test_data <- filter_by_significance(merged_test_data, taxonomic_hierarchy, feature.mt.method, finest_taxonomic_level, cutoff)

  # Subset data for the chosen color grouping level
  level_specific_data <- dplyr::filter(merged_test_data, .data$taxonomic_level == color.group.level)

  # Process "Unclassified" labels in the main dataframe
  merged_test_data <- process_unclassified_labels(merged_test_data, "taxonomic_level", use_grouping = TRUE)
  
  # Process "Unclassified" labels in the sub-dataframe
  level_specific_data <- process_unclassified_labels(level_specific_data, "taxonomic_level", use_grouping = FALSE)

  # Ensure tree labels match data labels
  matching_labels <- intersect(phylogenetic_tree$tip.label, level_specific_data$Variable)
  phylogenetic_tree <- ape::keep.tip(phylogenetic_tree, matching_labels)
  level_specific_data <- level_specific_data[level_specific_data$Variable %in% matching_labels, ]

  # Group tree nodes by taxonomic level for coloring
  split_group <- split(level_specific_data$Variable, f = level_specific_data[[color.group.level]])
  annotated_tree <- ggtree::groupOTU(phylogenetic_tree, .node = split_group, group_name = color.group.level)
  
  # Ensure tree is in the correct format for ggtree
  if (requireNamespace("tidytree", quietly = TRUE)) {
    # Convert to treedata format if not already
    if (!inherits(annotated_tree, "treedata")) {
      annotated_tree <- tidytree::as.treedata(annotated_tree)
    }
  }

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
  calculate_level_based_offset <- function(feature_level_length) {
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

  calculated_offset <- calculate_level_based_offset(length(feature.level))

  # Create a list to store individual plots
  plot.list <- list()

  # Generate a plot for each comparison
  for (comparison in unique(merged_test_data$Comparison)) {
    comparison_data <- merged_test_data %>% dplyr::filter(Comparison == comparison)

    # Ensure ggplot2 namespace is available in this scope
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required but not installed.")
    }
    
    # ggtreeExtra::geom_fruit需要在全局环境中查找这个函数或者使用字符串形式的函数名
    # 这是一个已知的问题，在ggtreeExtra包中的非标准求值机制导致的
    if (!exists("geom_tile", envir = .GlobalEnv)) {
      assign("geom_tile", ggplot2::geom_tile, envir = .GlobalEnv)
    }
    
    # Create the circular cladogram plot
    p <- ggtree::ggtree(annotated_tree, layout = "circular", open.angle = 5) +
      # 使用字符串"geom_tile"作为geom参数的值
      # 这符合ggtreeExtra::geom_fruit函数的要求
      ggtreeExtra::geom_fruit(
        data = comparison_data,
        geom = "geom_tile",  # 使用字符串，而非函数对象或未求值的符号
        mapping = ggplot2::aes(y = Variable, x = taxonomic_level, fill = Coefficient),
        offset = 0.03,
        size = 0.02,
        color = "black"
      ) +
      ggtree::geom_tiplab(ggplot2::aes(color = .data[[color.group.level]]), offset = calculated_offset, align = TRUE,
                          linetype = "blank", size = 2, show.legend = FALSE) +
      ggplot2::scale_fill_gradient2(low = "#0571b0", high = "#ca0020") +
      ggplot2::geom_text(mapping = ggplot2::aes(label = "",
                              color = .data[[color.group.level]]),
                key_glyph = ggplot2::draw_key_rect) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = ggplot2::element_text(size = 6),
        legend.title = ggplot2::element_text(size = 8)
      ) +
      ggplot2::scale_color_manual(values = unique(palette)[seq_along(unique(comparison_data[[color.group.level]]))]) +
      ggplot2::guides(
        fill = ggplot2::guide_colorbar(title = "Coefficient", barwidth = 10, barheight = 0.5),
        color = ggplot2::guide_legend(
          title = color.group.level,
          override.aes = list(size = 2),
          byrow = TRUE,
          keywidth = ggplot2::unit(0.5, "lines"),
          keyheight = ggplot2::unit(0.5, "lines")
        )
      )

    # Save as PDF if requested
    if (pdf) {
      filename <- paste0("taxa_cladogram_single_", gsub(" ", "_", comparison), ".pdf")
      ggplot2::ggsave(filename = filename, plot = p, width = pdf.width, height = pdf.height)
    }

    # Add the plot to the list
    plot.list[[comparison]] <- p
  }

  # Return the list of plots
  return(plot.list)
}