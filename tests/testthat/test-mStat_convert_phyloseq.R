# Test for mStat_convert_phyloseq_to_data_obj function
# Specifically testing the taxa_are_rows handling (Issue #83)

test_that("mStat_convert_phyloseq_to_data_obj handles taxa_are_rows correctly", {
  # Skip if phyloseq is not installed
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create test data
  set.seed(123)
  n_taxa <- 10
  n_samples <- 5

  # Create OTU matrix (taxa as rows, samples as columns)
  otu_mat <- matrix(
    sample(1:100, n_taxa * n_samples, replace = TRUE),
    nrow = n_taxa,
    ncol = n_samples
  )
  rownames(otu_mat) <- paste0("OTU", 1:n_taxa)
  colnames(otu_mat) <- paste0("Sample", 1:n_samples)

  # Create sample data
  sample_data_df <- data.frame(
    Group = rep(c("A", "B"), length.out = n_samples),
    row.names = paste0("Sample", 1:n_samples)
  )

  # Create taxonomy table
  tax_mat <- matrix(
    paste0("Tax", 1:(n_taxa * 7)),
    nrow = n_taxa,
    ncol = 7
  )
  rownames(tax_mat) <- paste0("OTU", 1:n_taxa)
  colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Test 1: taxa_are_rows = TRUE (standard case)
  otu_table_rows <- otu_table(otu_mat, taxa_are_rows = TRUE)
  sample_data_obj <- sample_data(sample_data_df)
  tax_table_obj <- tax_table(tax_mat)
  phylo_obj_rows <- phyloseq(otu_table_rows, sample_data_obj, tax_table_obj)

  data_obj_rows <- mStat_convert_phyloseq_to_data_obj(phylo_obj_rows)

  # Check dimensions
  expect_equal(nrow(data_obj_rows$feature.tab), n_taxa)
  expect_equal(ncol(data_obj_rows$feature.tab), n_samples)

  # Check row names are taxa
  expect_true(all(grepl("^OTU", rownames(data_obj_rows$feature.tab))))

  # Check column names are samples
  expect_true(all(grepl("^Sample", colnames(data_obj_rows$feature.tab))))

  # Test 2: taxa_are_rows = FALSE (the case that was broken)
  otu_mat_transposed <- t(otu_mat)  # samples as rows, taxa as columns
  otu_table_cols <- otu_table(otu_mat_transposed, taxa_are_rows = FALSE)
  phylo_obj_cols <- phyloseq(otu_table_cols, sample_data_obj, tax_table_obj)

  # Suppress the message about transposing
  suppressMessages({
    data_obj_cols <- mStat_convert_phyloseq_to_data_obj(phylo_obj_cols)
  })

  # Check dimensions (should be same as taxa_are_rows = TRUE)
  expect_equal(nrow(data_obj_cols$feature.tab), n_taxa)
  expect_equal(ncol(data_obj_cols$feature.tab), n_samples)

  # Check row names are taxa (not samples!)
  expect_true(all(grepl("^OTU", rownames(data_obj_cols$feature.tab))))

  # Check column names are samples (not taxa!)
  expect_true(all(grepl("^Sample", colnames(data_obj_cols$feature.tab))))

  # Test 3: Both methods should produce identical data
  expect_equal(data_obj_rows$feature.tab, data_obj_cols$feature.tab)

  # Test 4: Check feature.ann alignment
  expect_true(all(rownames(data_obj_rows$feature.ann) %in% rownames(data_obj_rows$feature.tab)))
  expect_true(all(rownames(data_obj_cols$feature.ann) %in% rownames(data_obj_cols$feature.tab)))

  # Test 5: meta.dat should be identical
  expect_equal(data_obj_rows$meta.dat, data_obj_cols$meta.dat)
})

test_that("mStat_convert_phyloseq_to_data_obj handles square matrices correctly", {
  # Skip if phyloseq is not installed
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create square matrix test (same number of taxa and samples)
  set.seed(456)
  n <- 8

  otu_mat_sq <- matrix(
    sample(1:100, n * n, replace = TRUE),
    nrow = n,
    ncol = n
  )
  rownames(otu_mat_sq) <- paste0("OTU", 1:n)
  colnames(otu_mat_sq) <- paste0("Sample", 1:n)

  sample_data_sq <- data.frame(
    Group = rep(c("A", "B"), length.out = n),
    row.names = paste0("Sample", 1:n)
  )

  tax_mat_sq <- matrix(
    paste0("Tax", 1:(n * 7)),
    nrow = n,
    ncol = 7
  )
  rownames(tax_mat_sq) <- paste0("OTU", 1:n)
  colnames(tax_mat_sq) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Test with taxa_are_rows = TRUE
  otu_table_sq_rows <- otu_table(otu_mat_sq, taxa_are_rows = TRUE)
  phylo_obj_sq_rows <- phyloseq(otu_table_sq_rows, sample_data(sample_data_sq), tax_table(tax_mat_sq))
  data_obj_sq_rows <- mStat_convert_phyloseq_to_data_obj(phylo_obj_sq_rows)

  # Test with taxa_are_rows = FALSE (need to transpose for phyloseq object)
  otu_mat_sq_t <- t(otu_mat_sq)  # samples as rows, taxa as columns
  otu_table_sq_cols <- otu_table(otu_mat_sq_t, taxa_are_rows = FALSE)
  phylo_obj_sq_cols <- phyloseq(otu_table_sq_cols, sample_data(sample_data_sq), tax_table(tax_mat_sq))

  suppressMessages({
    data_obj_sq_cols <- mStat_convert_phyloseq_to_data_obj(phylo_obj_sq_cols)
  })

  # Both should have taxa as rows
  expect_true(all(grepl("^OTU", rownames(data_obj_sq_rows$feature.tab))))
  expect_true(all(grepl("^OTU", rownames(data_obj_sq_cols$feature.tab))))

  # Both should have samples as columns
  expect_true(all(grepl("^Sample", colnames(data_obj_sq_rows$feature.tab))))
  expect_true(all(grepl("^Sample", colnames(data_obj_sq_cols$feature.tab))))
})

test_that("mStat_convert_SummarizedExperiment_to_data_obj keeps aligned feature annotations", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      c(1, 2,
        0, 0,
        3, 4),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2"))
    )),
    rowData = S4Vectors::DataFrame(phylum = c("p1", "p2", "p3"), row.names = c("f1", "f2", "f3")),
    colData = S4Vectors::DataFrame(group = c("A", "B"), row.names = c("s1", "s2"))
  )

  data.obj <- mStat_convert_SummarizedExperiment_to_data_obj(se)

  expect_identical(rownames(data.obj$feature.tab), c("f1", "f3"))
  expect_identical(rownames(data.obj$feature.ann), c("f1", "f3"))
  expect_equal(dim(data.obj$feature.ann), c(2L, 1L))
  expect_identical(rownames(data.obj$meta.dat), c("s1", "s2"))
})

test_that("mStat_convert_DESeqDataSet_to_data_obj keeps aligned feature annotations", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("S4Vectors")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = matrix(
      c(1, 2,
        0, 0,
        3, 4),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2"))
    ),
    colData = S4Vectors::DataFrame(group = c("A", "B"), row.names = c("s1", "s2")),
    design = ~ 1
  )
  SummarizedExperiment::rowData(dds)$phylum <- c("p1", "p2", "p3")

  data.obj <- mStat_convert_DESeqDataSet_to_data_obj(dds)

  expect_identical(rownames(data.obj$feature.tab), c("f1", "f3"))
  expect_identical(rownames(data.obj$feature.ann), c("f1", "f3"))
  expect_equal(dim(data.obj$feature.ann), c(2L, 1L))
  expect_identical(rownames(data.obj$meta.dat), c("s1", "s2"))
})

test_that("mStat_convert_MRExperiment_to_data_obj keeps aligned feature annotations", {
  skip_if_not_installed("metagenomeSeq")
  skip_if_not_installed("Biobase")

  mr <- metagenomeSeq::newMRexperiment(
    counts = matrix(
      c(1, 2,
        0, 0,
        3, 4),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2"))
    ),
    phenoData = Biobase::AnnotatedDataFrame(
      data.frame(group = c("A", "B"), row.names = c("s1", "s2"))
    ),
    featureData = Biobase::AnnotatedDataFrame(
      data.frame(phylum = c("p1", "p2", "p3"), row.names = c("f1", "f2", "f3"))
    )
  )

  data.obj <- mStat_convert_MRExperiment_to_data_obj(mr)

  expect_identical(rownames(data.obj$feature.tab), c("f1", "f3"))
  expect_identical(rownames(data.obj$feature.ann), c("f1", "f3"))
  expect_equal(dim(data.obj$feature.ann), c(2L, 1L))
  expect_identical(rownames(data.obj$meta.dat), c("s1", "s2"))
})

test_that("mStat_convert_MultiAssayExperiment_to_data_obj extracts named experiment with namespace-safe access", {
  skip_if_not_installed("MultiAssayExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  se1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      c(1, 2,
        0, 0,
        3, 4),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("f1", "f2", "f3"), c("s1", "s2"))
    ))
  )
  se2 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(
      c(5, 6,
        7, 8),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("g1", "g2"), c("s1", "s2"))
    ))
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rna = se1, taxa = se2),
    colData = S4Vectors::DataFrame(group = c("A", "B"), row.names = c("s1", "s2"))
  )

  data.obj <- mStat_convert_MultiAssayExperiment_to_data_obj(mae, experiment_name = "taxa")

  expect_identical(rownames(data.obj$feature.tab), c("g1", "g2"))
  expect_identical(colnames(data.obj$feature.tab), c("s1", "s2"))
  expect_identical(rownames(data.obj$meta.dat), c("s1", "s2"))
})

test_that("mStat_convert_phyloseq_to_data_obj removes zero-sum features correctly", {
  # Skip if phyloseq is not installed
  skip_if_not_installed("phyloseq")

  library(phyloseq)

  # Create OTU matrix with some zero-sum rows
  set.seed(789)
  n_taxa <- 12
  n_samples <- 5

  otu_mat <- matrix(
    sample(0:50, n_taxa * n_samples, replace = TRUE),
    nrow = n_taxa,
    ncol = n_samples
  )
  rownames(otu_mat) <- paste0("OTU", 1:n_taxa)
  colnames(otu_mat) <- paste0("Sample", 1:n_samples)

  # Force some rows to be all zeros
  otu_mat[c(3, 7, 10), ] <- 0

  sample_data_df <- data.frame(
    Group = rep(c("A", "B"), length.out = n_samples),
    row.names = paste0("Sample", 1:n_samples)
  )

  # Test with taxa_are_rows = TRUE
  otu_table_obj <- otu_table(otu_mat, taxa_are_rows = TRUE)
  phylo_obj <- phyloseq(otu_table_obj, sample_data(sample_data_df))
  data_obj <- mStat_convert_phyloseq_to_data_obj(phylo_obj)

  # Should have removed 3 zero-sum features
  expect_equal(nrow(data_obj$feature.tab), n_taxa - 3)

  # All remaining features should have sum > 0
  expect_true(all(rowSums(data_obj$feature.tab) > 0))

  # Zero-sum OTUs should not be present
  expect_false("OTU3" %in% rownames(data_obj$feature.tab))
  expect_false("OTU7" %in% rownames(data_obj$feature.tab))
  expect_false("OTU10" %in% rownames(data_obj$feature.tab))

  # Test with taxa_are_rows = FALSE
  otu_mat_t <- t(otu_mat)
  otu_table_obj_t <- otu_table(otu_mat_t, taxa_are_rows = FALSE)
  phylo_obj_t <- phyloseq(otu_table_obj_t, sample_data(sample_data_df))

  suppressMessages({
    data_obj_t <- mStat_convert_phyloseq_to_data_obj(phylo_obj_t)
  })

  # Should also have removed 3 zero-sum features
  expect_equal(nrow(data_obj_t$feature.tab), n_taxa - 3)

  # Results should be identical
  expect_equal(data_obj$feature.tab, data_obj_t$feature.tab)
})
