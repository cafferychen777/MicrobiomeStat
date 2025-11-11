#' Visualize the relative abundance of selected taxa or functions
#'
#' This function generates various types of plots to visualize the relative abundance
#' of selected taxa or functions in microbiome data. It supports different visualization
#' schemes and plot types for both cross-sectional and longitudinal data.
#'
#' @param data.obj A list object in a format specific to MicrobiomeStat, which can include components such as feature.tab (matrix), feature.ann (matrix), meta.dat (data.frame), tree, and feature.agg.list (list). The data.obj can be converted from other formats using several functions from the MicrobiomeStat package, including: 'mStat_convert_DGEList_to_data_obj', 'mStat_convert_DESeqDataSet_to_data_obj', 'mStat_convert_phyloseq_to_data_obj', 'mStat_convert_SummarizedExperiment_to_data_obj', 'mStat_import_qiime2_as_data_obj', 'mStat_import_mothur_as_data_obj', 'mStat_import_dada2_as_data_obj', and 'mStat_import_biom_as_data_obj'. Alternatively, users can construct their own data.obj. Note that not all components of data.obj may be required for all functions in the MicrobiomeStat package.
#' @param group.var A string indicating the variable for group identifiers. Default is NULL.
#' @param strata.var A string indicating the variable for strata identifiers. Default is NULL.
#' @param subject.var Character string specifying the subject variable in meta.dat
#' @param time.var Character string specifying the time variable in meta.dat
#' @param time.point.plot Character vector specifying time points to plot
#' @param is.plot.change Logical, whether to plot change from baseline.
#'   This parameter is only effective when there are multiple time points
#'   (i.e., length(time.point.plot) > 1). If there's only one time point,
#'   this parameter has no effect.
#' @param feature.change.func The method or function used to calculate the change in feature abundance between time points.
#' The following options are supported:
#'
#' - "relative change": Computes the relative change as (time_2 - time_1) / (time_2 + time_1). If both values are zero, the result is zero.
#' - "log fold change": Computes the log2 fold change between time points. Zero values are imputed as half the minimum nonzero value of the respective feature across BOTH time points combined. The same pseudocount is used at both time points to ensure unbiased log fold change calculations.
#' - "absolute change": Computes the absolute difference between time points.
#' - A custom function: The provided function should take two numeric vectors as input (values at time 1 and time 2) and return a numeric vector of differences. Users should ensure that their function handles zero values appropriately.
#'
#' If an unrecognized value or no value is provided for `feature.change.func`, the default behavior will be to compute the absolute difference between time points.
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
#' @param features.plot A character vector specifying which feature IDs (e.g. OTU IDs) to plot.
#' Default is NULL.
#' @param prop.to.lump Numeric, features with mean proportion less than this will be lumped into "other"
#' @param top.k.plot Integer, plot top k features based on some criterion
#' @param top.k.func Function to order the features
#' @param plot.other Logical, whether to plot the "other" category
#' @param renormalize Logical, whether to renormalize data after removing "other" category
#' @param plot.scheme Character string specifying the plot scheme, one of c("combined", "individual")
#' @param plot.type Character string specifying the plot type, one of c("barplot", "dotplot", "areaplot", "heatmap", "spaghettiplot", "boxplot", "scatterplot")
#'
#' @return A list of ggplot objects
#'
#' @details This function provides a flexible framework for visualizing microbiome data.
#' It can handle both cross-sectional and longitudinal data, and offers various plot types
#' including bar plots, dot plots, area plots, heatmaps, spaghetti plots, box plots, and scatter plots.
#' The function can also visualize changes over time and supports different grouping and stratification options.
#' @examples
#' \donttest{
#' data(ecam.obj)
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = NULL,
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' feature.dat.type = "proportion",
#' features.plot = NULL,
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' feature.dat.type = "proportion",
#' features.plot = NULL,
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum"),
#' feature.dat.type = "proportion",
#' features.plot = ecam.obj$feature.ann[,"Phylum"][1:3],
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' feature.dat.type = "proportion",
#' features.plot = NULL,
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "dotplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' feature.dat.type = "proportion",
#' features.plot = NULL,
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "heatmap"
#' )
#'
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' feature.dat.type = "proportion",
#' features.plot = NULL,
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "boxplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' feature.dat.type = "proportion",
#' features.plot = NULL,
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "individual",
#' plot.type = "boxplot"
#' )
#'
#' data(peerj32.obj)
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = TRUE,
#' feature.change.func = "relative change",
#' feature.level = c("Phylum", "Family", "Genus"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "boxplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = TRUE,
#' feature.change.func = "relative change",
#' feature.level = c("Phylum", "Family", "Genus"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "individual",
#' plot.type = "boxplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = TRUE,
#' feature.change.func = "relative change",
#' feature.level = c("Phylum", "Family", "Genus"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = "combined",
#' plot.type = "heatmap"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = FALSE,
#' feature.change.func = "relative change",
#' feature.level = c("Phylum", "Family", "Genus"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "heatmap"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = FALSE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'individual',
#' plot.type = "boxplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = FALSE,
#' feature.level = c("Phylum", "Family", "Genus"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' features.plot = unique(peerj32.obj$feature.ann[, "Family"])[11:20],
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = peerj32.obj,
#' group.var = "group",
#' strata.var = "sex",
#' subject.var = "subject",
#' time.var = "time",
#' time.point.plot = c("1", "2"),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "dotplot"
#' )
#'
#' data(subset_T2D.obj)
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' features.plot = unique(subset_T2D.obj$feature.ann[, "Family"])[1:10],
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "barplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' features.plot = unique(subset_T2D.obj$feature.ann[, "Family"])[1:10],
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "areaplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "spaghettiplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'individual',
#' plot.type = "spaghettiplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "heatmap"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = TRUE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "heatmap"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'combined',
#' plot.type = "boxplot"
#' )
#'
#' plot_feature_diversity(
#' data.obj = subset_T2D.obj,
#' group.var = "sample_body_site",
#' strata.var = "subject_race",
#' subject.var = "subject_id",
#' time.var = "visit_number_num",
#' time.point.plot = unique(subset_T2D.obj$meta.dat$visit_number_num),
#' is.plot.change = FALSE,
#' feature.level = c("Family"),
#' prop.to.lump = 0.0001,
#' top.k.plot = NULL,
#' top.k.func = NULL,
#' renormalize = FALSE,
#' plot.other = TRUE,
#' plot.scheme = 'individual',
#' plot.type = "boxplot"
#' )
#' }
#' @export
plot_feature_diversity <- function (data.obj,
                                    group.var = NULL,
                                    strata.var = NULL,
                                    subject.var = NULL,
                                    time.var = NULL,
                                    time.point.plot,
                                    # The first pt will be the reference.

                                    is.plot.change = TRUE,

                                    feature.change.func = "relative change",
                                    feature.level = NULL,
                                    # The taxonomic/functional level to be plotted
                                    feature.dat.type = c("count", "proportion", "other"),
                                    features.plot = NULL,
                                    # Features to be plotted if specified

                                    # The following are for taxa selection strategy when features.plot is not specified
                                    prop.to.lump = 0.0001,
                                    # Mean proportion less than that will be lumped into "other"
                                    top.k.plot = NULL,
                                    # Or top k based on some criterion, other will be lumped into "other"
                                    top.k.func = NULL,
                                    # function to order the features
                                    plot.other = TRUE,
                                    # whether "other" group will be plotted
                                    renormalize = FALSE,
                                    # if plot.other = FALSE, whether renormalization to sum 1 will be performed

                                    plot.scheme = c('combined', 'individual'),
                                    plot.type = c(
                                      "barplot",
                                      "dotplot",
                                      "areaplot",
                                      "heatmap",
                                      "spaghettiplot",
                                      "boxplot",
                                      "scatterplot"
                                    )) {
  # Input validation
  feature.dat.type <- match.arg(feature.dat.type)
  plot.scheme <- match.arg(plot.scheme)
  plot.type <- match.arg(plot.type)

  plot.list <- lapply(feature.level, function(feature.level) {
    # Process the data to lump taxa with low mean proportions into "other"
    lump_low_abundance_taxa <- function(data.obj,
                                        prop.to.lump,
                                        feature.level) {
      feature.tab <- data.obj$feature.tab
      original_feature.tab <- feature.tab  # Store the original count data

      # Convert to proportions for identifying low abundance taxa
      feature.tab_prop <- sweep(feature.tab, 2, colSums(feature.tab), '/')

      # Calculate mean proportions for each taxon
      mean_props <- rowMeans(feature.tab_prop)

      # Identify taxa to lump into "other"
      taxa_to_lump <- names(mean_props[mean_props < prop.to.lump])

      # Create a new row for "other" and remove the lumped taxa
      if (length(taxa_to_lump) > 0) {
        other_row <- colSums(original_feature.tab[taxa_to_lump, , drop = FALSE])
        feature.tab <- original_feature.tab[!rownames(original_feature.tab) %in% taxa_to_lump, ]
        feature.tab <- rbind(feature.tab, other = other_row)

        # Update feature annotation
        if (!is.null(data.obj$feature.ann)) {
          other_ann <- rep("Other", ncol(data.obj$feature.ann))
          data.obj$feature.ann <- rbind(data.obj$feature.ann[!rownames(data.obj$feature.ann) %in% taxa_to_lump, ], other = other_ann)
        }
      }

      # Update the data object
      data.obj$feature.tab <- feature.tab
      return(data.obj)
    }

    if (!is.null(features.plot)){
      prop.to.lump <- 0
    }

    # Process the data
    data.obj <- lump_low_abundance_taxa(data.obj, prop.to.lump, feature.level)

    if (!plot.other & renormalize) {
      data.obj <- mStat_remove_feature(data.obj, "Other", feature.level)
      data.obj <- mStat_normalize_data(data.obj, "TSS")$data.obj.norm
      feature.dat.type <- "proportion"
    } else if (!plot.other & !renormalize) {
      data.obj <- mStat_remove_feature(data.obj, "Other", feature.level)
    }

    if (is.null(time.var) | length(time.point.plot) == 1) {
      if (plot.scheme == "combined") {
        if (plot.type == "barplot") {
          p <- generate_taxa_barplot_single(
            data.obj = data.obj,
            subject.var = subject.var,
            time.var = time.var,
            t.level = time.point.plot[1],
            group.var = group.var,
            strata.var = strata.var,
            feature.level = feature.level,
            feature.dat.type = feature.dat.type,
            features.plot = features.plot
          )
        } else if (plot.type == "dotplot") {
          p <- generate_taxa_dotplot_single(
            data.obj = data.obj,
            subject.var = subject.var,
            time.var = time.var,
            t.level = time.point.plot[1],
            group.var = group.var,
            strata.var = strata.var,
            feature.level = feature.level,
            feature.dat.type = feature.dat.type,
            features.plot = features.plot,
            top.k.plot = top.k.plot,
            top.k.func = top.k.func,
          )
        } else if (plot.type == "heatmap") {
          p <- generate_taxa_heatmap_single(
            data.obj = data.obj,
            subject.var = subject.var,
            time.var = time.var,
            t.level = time.point.plot[1],
            group.var = group.var,
            strata.var = strata.var,
            feature.level = feature.level,
            feature.dat.type = feature.dat.type,
            features.plot = features.plot,
            top.k.func = top.k.func,
            top.k.plot = top.k.func
          )
        } else if (plot.type == "boxplot") {
          p <- generate_taxa_boxplot_single(
            data.obj = data.obj,
            subject.var = subject.var,
            time.var = time.var,
            t.level = time.point.plot[1],
            group.var = group.var,
            strata.var = strata.var,
            feature.level = feature.level,
            feature.dat.type = feature.dat.type,
            features.plot = features.plot,
            top.k.func = top.k.func,
            top.k.plot = top.k.plot
          )
        } else {
          message(paste(
            "Currently, we do not support",
            plot.type,
            "output for this scenario."
          ))
          return()
        }
      } else {
        if (plot.type == "boxplot") {
          p <- generate_taxa_indiv_boxplot_single(
            data.obj = data.obj,
            subject.var = subject.var,
            time.var = time.var,
            t.level = time.point.plot[1],
            group.var = group.var,
            strata.var = strata.var,
            feature.level = feature.level,
            features.plot = features.plot,
            feature.dat.type = feature.dat.type,
            top.k.func = top.k.func,
            top.k.plot = top.k.plot
          )
        } else {
          message(paste(
            "Currently, we do not support",
            plot.type,
            "output for this scenario."
          ))
          return()
        }
      }
    } else if (!is.null(time.var) &
               length(time.point.plot) == 2) {
      # Plot two time points and sample pair
      if (is.null(subject.var)) {
        message("Subject variable not specified!")
        return()
      }

      if (is.plot.change) {
        if (plot.scheme == "combined") {
          if (plot.type == "boxplot") {
            p <- generate_taxa_change_boxplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              change.base = time.point.plot[1],
              feature.change.func = feature.change.func,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "heatmap") {
            p <- generate_taxa_change_heatmap_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              change.base = time.point.plot[1],
              feature.change.func = feature.change.func,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "dotplot") {
            p <- generate_taxa_change_dotplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              change.base = time.point.plot[1],
              feature.change.func = feature.change.func,
              fetaure.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "scatterplot") {
            p <- generate_taxa_change_scatterplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              change.base = time.point.plot[1],
              feature.change.func = feature.change.func,
              fetaure.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else {
            message(
              paste(
                "Currently, we do not support",
                plot.type,
                "output for this scenario."
              )
            )
            return()
          }
        } else {
          if (plot.type == "boxplot") {
            p <- generate_taxa_indiv_change_boxplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              change.base = time.point.plot[1],
              feature.change.func = feature.change.func,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "scatterplot") {
            p <- generate_taxa_indiv_change_scatterplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              change.base = time.point.plot[1],
              feature.change.func = feature.change.func,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.func = top.k.func,
              top.k.plot = top.k.plot
            )
          } else {
            message(
              paste(
                "Currently, we do not support",
                plot.type,
                "output for this scenario."
              )
            )
            return()
          }
        }
      } else {
        if (plot.scheme == "combined") {
          if (plot.type == "barplot") {
            p <- generate_taxa_barplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot
            )
          } else if (plot.type == "heatmap") {
            p <- generate_taxa_heatmap_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.func = top.k.func,
              top.k.plot = top.k.plot
            )
          } else if (plot.type == "dotplot") {
            p <- generate_taxa_dotplot_pair(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.func = top.k.func,
              top.k.plot = top.k.plot
            )
          } else {
            message(
              paste(
                "Currently, we do not support",
                plot.type,
                "output for this scenario."
              )
            )
            return()
          }
        } else {
          if (plot.type == "boxplot"){
            p <- generate_taxa_indiv_boxplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[2],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.func = top.k.func,
              top.k.plot = top.k.plot
            )
          } else {
            message(
              paste(
                "Currently, we do not support",
                plot.type,
                "output for this scenario."
              )
            )
            return()
          }
        }
      }
    } else if (!is.null(time.var) & length(time.point.plot) > 2) {
      # Plot more than two time points, which are truly longitudinal
      if (is.null(subject.var)) {
        message("Subject variable not specified!")
        return()
      }

      if (is.plot.change) {
        if (plot.type == "heatmap") {
          p <- generate_taxa_change_heatmap_long(
            data.obj = data.obj,
            subject.var = subject.var,
            time.var = time.var,
            t0.level = time.point.plot[1],
            ts.levels = time.point.plot[-1],
            group.var = group.var,
            strata.var = strata.var,
            feature.level = feature.level,
            feature.change.func = feature.change.func,
            fetaure.dat.type = feature.dat.type,
            features.plot = features.plot,
            top.k.plot = top.k.plot,
            top.k.func = top.k.func
          )
        } else {
          message(paste(
            "Currently, we do not support",
            plot.type,
            "output for this scenario."
          ))
          return()
        }
      } else {
        if (plot.scheme == "combined") {
          if (plot.type == "areaplot") {
            p <- generate_taxa_areaplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot
            )
          } else if (plot.type == "barplot") {
            p <- generate_taxa_barplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot
            )
          } else if (plot.type == "boxplot") {
            p <- generate_taxa_boxplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "heatmap") {
            p <- generate_taxa_heatmap_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "spaghettiplot") {
            p <- generate_taxa_spaghettiplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else {
            message(
              paste(
                "Currently, we do not support",
                plot.type,
                "output for this scenario."
              )
            )
            return()
          }
        } else {
          if (plot.type == "boxplot") {
            p <- generate_taxa_indiv_boxplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.plot = top.k.plot,
              top.k.func = top.k.func
            )
          } else if (plot.type == "spaghettiplot") {
            p <- generate_taxa_indiv_spaghettiplot_long(
              data.obj = data.obj,
              subject.var = subject.var,
              time.var = time.var,
              t0.level = time.point.plot[1],
              ts.levels = time.point.plot[-1],
              group.var = group.var,
              strata.var = strata.var,
              feature.level = feature.level,
              feature.dat.type = feature.dat.type,
              features.plot = features.plot,
              top.k.func = top.k.func,
              top.k.plot = top.k.plot
            )
          } else {
            message(
              paste(
                "Currently, we do not support",
                plot.type,
                "output for this scenario."
              )
            )
            return()
          }
        }
      }
    }
    else {
      return()
    }

    return(p)
  })

  names(plot.list) <- feature.level
  return(plot.list)
}
