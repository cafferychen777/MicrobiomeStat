#' Visualize Feature Diversity
#'
#' Generates various plots to visualize relative abundance of taxa or functions
#' for both cross-sectional and longitudinal microbiome data.
#'
#' @inheritParams mStat_data_obj_doc
#' @param time.point.plot Character vector of time points to plot.
#' @param is.plot.change Logical, plot change from baseline (only for multiple time points).
#' @param feature.change.func Method for calculating change: "relative change", "log fold change", or "absolute change".
#' @param features.plot Character vector of specific feature IDs to plot (default NULL plots all).
#' @param prop.to.lump Numeric, features below this mean proportion are lumped into "other".
#' @param top.k.plot Integer, number of top features to plot.
#' @param top.k.func Function to rank features for top.k selection.
#' @param plot.other Logical, whether to include "other" category in plots.
#' @param renormalize Logical, renormalize after removing "other" category.
#' @param plot.scheme Character, "combined" or "individual" plot layout.
#' @param plot.type Character, plot type: "barplot", "dotplot", "areaplot", "heatmap", "spaghettiplot", "boxplot", or "scatterplot".
#'
#' @return A list of ggplot objects.
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
#' # Get valid Phylum values (excluding empty strings and NAs)
#' valid_phyla <- unique(ecam.obj$feature.ann[,"Phylum"])
#' valid_phyla <- valid_phyla[!valid_phyla %in% c("__", "") & !is.na(valid_phyla)]
#' plot_feature_diversity(
#' data.obj = ecam.obj,
#' group.var = "antiexposedall",
#' strata.var = "delivery",
#' time.var = "month",
#' time.point.plot = unique(ecam.obj$meta.dat$month)[1],
#' is.plot.change = TRUE,
#' feature.level = c("Phylum"),
#' feature.dat.type = "proportion",
#' features.plot = valid_phyla[1:3],
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
            top.k.plot = top.k.plot
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
              feature.level = feature.level,
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
            feature.dat.type = feature.dat.type,
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
