#' Compute and Analyze Pairwise Beta Diversity Measures
#'
#' This function calculates beta diversity measures, performs linear modeling on them, and generates a report.
#' If no precalculated distance object is provided, it calculates the distance object from the data.
#'
#' @param data.obj A list object in MicrobiomeStat format.
#' @param dist.obj A pre-calculated distance object (optional). If not provided, it is computed using the data.obj.
#' @param time.var The name of the time variable column in the metadata.
#' @param subject.var The name of the subject variable column in the metadata.
#' @param group.var The name of the grouping variable column for linear modeling in the metadata.
#' @param adj.vars Names of additional variables to be used as covariates in the analysis.
#' @param change.base The baseline time point for detecting changes in beta diversity.
#' @param dist.name Character vector specifying which beta diversity indices to compute. Defaults to c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS').
#'
#' @examples
#'
#' library(vegan)
#' library(GUniFrac)
#' library(ape)
#' library(philentropy)
#' library(MicrobiomeStat)
#'
#' # Load example data
#' data(peerj32.obj)
#' dist.obj <- mStat_calculate_beta_diversity(peerj32.obj, dist.name = c('BC', 'Jaccard'))
#'
#' # Generate beta diversity report
#' beta_diversity_report <- generate_beta_change_test_pair(
#'   data.obj = peerj32.obj,
#'   dist.obj = NULL,
#'   time.var = "time",
#'   subject.var = "subject",
#'   group.var = "group",
#'   adj.vars = c("sex"),
#'   change.base = "1",
#'   dist.name = c('BC', 'Jaccard')
#' )
#'
#'
#' @return A named list where each element is a coefficient table (with p-values) from a linear model of the corresponding beta diversity measure against the grouping variable and covariates.
#' @export
#' @name generate_beta_change_test_pair
generate_beta_change_test_pair <-
  function(data.obj,
           dist.obj = NULL,
           time.var = NULL,
           subject.var,
           group.var,
           adj.vars,
           change.base = NULL,
           dist.name = c('BC', 'Jaccard', 'UniFrac', 'GUniFrac', 'WUniFrac', 'JS')) {

    if (is.null(dist.obj)&!is.null(data.obj)) {
      dist.obj <-
        mStat_calculate_beta_diversity(data.obj = data.obj, dist.name = dist.name)
      metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,time.var, adj.vars))) %>% rownames_to_column("sample")
    } else {
      if (!is.null(data.obj)){
        metadata <- load_data_obj_metadata(data.obj) %>% select(all_of(c(subject.var,group.var,time.var, adj.vars))) %>% rownames_to_column("sample")
      } else {
        metadata <- attr(dist.obj[[dist.name[1]]], "labels")  %>% select(all_of(c(subject.var,group.var,time.var,adj.vars))) %>% rownames_to_column("sample")
      }
    }

    if (is.null(change.base)){
      change.base <- unique(metadata %>% select(all_of(c(time.var))))[1,]
      message("The 'change.base' variable was NULL. It has been set to the first unique value in the 'time.var' column of the 'meta.dat' data frame: ", change.base)
    }

    change.after <-
      unique(metadata %>% select(all_of(c(time.var))))[unique(metadata %>% select(all_of(c(time.var)))) != change.base]

    test.list <- lapply(dist.name, function(dist.name){
      dist.df <- as.matrix(dist.obj[[dist.name]]) %>%
        as.data.frame() %>%
        rownames_to_column("sample")

      long.df <- dist.df %>%
        gather(key = "sample2", value = "distance", -sample) %>%
        left_join(metadata, by = "sample") %>%
        left_join(metadata, by = c("sample2" = "sample"), suffix = c(".subject", ".sample")) %>%
        filter(!!sym(paste0(subject.var, ".subject")) == !!sym(paste0(subject.var, ".sample"))) %>%
        group_by(!!sym(paste0(subject.var, ".subject"))) %>%
        filter(!!sym(paste0(time.var,".sample")) == change.base) %>%
        filter(!!sym(paste0(time.var,".subject")) != !!sym(paste0(time.var,".sample"))) %>%
        ungroup() %>%
        select(!!sym(paste0(subject.var, ".subject")), !!sym(paste0(time.var, ".subject")), distance) %>%
        dplyr::rename(!!sym(subject.var) := !!sym(paste0(subject.var, ".subject")), !!sym(time.var) := !!sym(paste0(time.var, ".subject")))

      long.df <- long.df %>% left_join(metadata %>% select(-any_of(c(time.var))) %>% distinct(), by = subject.var)

      # Create a formula for lm
      formula <-
        as.formula(paste0("distance", "~", paste(c(
          adj.vars, group.var
        ), collapse = "+")))

      # Run lm and create a coefficient table
      lm.model <- lm(formula, data = long.df)
      coef.tab <- broom::tidy(summary(lm.model))

      # Rearrange the table
      coef.tab <-
        coef.tab %>% select(
          Term = term,
          Estimate = estimate,
          Std.Error = std.error,
          Statistic = statistic,
          P.Value = p.value
        )

      # Run ANOVA on the model if group.var is multi-categorical
      if (length(unique(metadata[[group.var]])) > 2) {
        anova.tab <- broom::tidy(anova(lm.model))

        # Rearrange the table and add missing columns
        anova.tab <- anova.tab %>%
          select(
            Term = term,
            Statistic = statistic,
            df = df,
            P.Value = p.value
          ) %>%
          mutate(Estimate = NA, Std.Error = NA)

        # Reorder the columns to match coef.tab
        anova.tab <- anova.tab %>%
          select(
            Term = term,
            Estimate = Estimate,
            Std.Error = Std.Error,
            Statistic = Statistic,
            P.Value = P.Value
          )

        coef.tab <-
          rbind(coef.tab, anova.tab) # Append the anova.tab to the coef.tab
      }

      return(coef.tab)
    })

    # Assign names to the elements of test.list
    names(test.list) <- dist.name

    return(test.list)

  }
