test_that("generate_taxa_cladogram_single works with and without phylogenetic tree", {
  # Create mock data objects with more realistic taxonomic structure
  features <- paste0("ASV", 1:5)  # Use simple ASV identifiers

  # Create feature annotation data frame - avoid using column names that might conflict
  feature.ann <- data.frame(
    tax_Kingdom = rep("Bacteria", 5),
    tax_Phylum = c("Firmicutes", "Firmicutes", "Bacteroidetes", "Bacteroidetes", "Proteobacteria"),
    tax_Class = c("Bacilli", "Clostridia", "Bacteroidia", "Flavobacteriia", "Gammaproteobacteria"),
    row.names = features,
    stringsAsFactors = FALSE
  )

  # Create abundance matrix
  feature.tab <- matrix(
    rpois(50, lambda = 10),
    nrow = 5,
    ncol = 10,
    dimnames = list(features, paste0("Sample", 1:10))
  )

  # Create metadata
  meta.dat <- data.frame(
    Group = rep(c("A", "B"), each = 5),
    row.names = paste0("Sample", 1:10)
  )

  # Create base data object
  data.obj.no.tree <- list(
    feature.tab = feature.tab,
    feature.ann = feature.ann,
    meta.dat = meta.dat
  )

  # Create mock test results
  test.list <- list(
    tax_Phylum = list(
      "A_vs_B" = data.frame(
        Feature = unique(feature.ann$tax_Phylum),
        Coefficient = rnorm(3),
        P.Value = runif(3),
        Adjusted.P.Value = runif(3)
      )
    ),
    tax_Class = list(
      "A_vs_B" = data.frame(
        Feature = unique(feature.ann$tax_Class),
        Coefficient = rnorm(5),
        P.Value = runif(5),
        Adjusted.P.Value = runif(5)
      )
    )
  )

  # Create data object with tree - ensure tree labels completely match features
  library(ape)
  tree <- ape::read.tree(text = paste0("(", paste(features, collapse = ","), ");"))
  tree$edge.length <- rep(1, nrow(tree$edge))
  # Add necessary attributes
  tree$node.label <- paste0("Node", 1:(tree$Nnode))

  data.obj.with.tree <- data.obj.no.tree
  data.obj.with.tree$tree <- tree

  # Test cases
  suppressWarnings({
    # Test case 1: Without tree
    result_no_tree <- capture_messages({
      plot.list.no.tree <- generate_taxa_cladogram_single(
        data.obj = data.obj.no.tree,
        test.list = test.list,
        group.var = "Group",
        feature.level = c("tax_Phylum", "tax_Class"),
        color.group.level = "tax_Phylum"
      )
    })
    expect_true(any(grepl("Building taxonomy-based tree", result_no_tree)))

    # Test case 2: With tree
    result_with_tree <- capture_messages({
      plot.list.with.tree <- generate_taxa_cladogram_single(
        data.obj = data.obj.with.tree,
        test.list = test.list,
        group.var = "Group",
        feature.level = c("tax_Phylum", "tax_Class"),
        color.group.level = "tax_Phylum"
      )
    })
    cat("\nActual messages with tree:\n")
    print(result_with_tree)

    # Test case 3: With mismatched tree
    mismatched_tree <- ape::read.tree(text = "(Wrong1,Wrong2,Wrong3,Wrong4,Wrong5);")
    mismatched_tree$edge.length <- rep(1, nrow(mismatched_tree$edge))
    data.obj.mismatched <- data.obj.no.tree
    data.obj.mismatched$tree <- mismatched_tree

    result_mismatched <- capture_messages({
      plot.list.mismatched <- generate_taxa_cladogram_single(
        data.obj = data.obj.mismatched,
        test.list = test.list,
        group.var = "Group",
        feature.level = c("tax_Phylum", "tax_Class"),
        color.group.level = "tax_Phylum"
      )
    })
  })

  # Basic checks
  expect_type(plot.list.no.tree, "list")
  expect_type(plot.list.with.tree, "list")
  expect_type(plot.list.mismatched, "list")

  # Check list contents
  expect_true(all(sapply(plot.list.no.tree, function(x) inherits(x, "ggplot"))))
  expect_true(all(sapply(plot.list.with.tree, function(x) inherits(x, "ggplot"))))
  expect_true(all(sapply(plot.list.mismatched, function(x) inherits(x, "ggplot"))))

  # Check list length
  expect_length(plot.list.no.tree, 1)
  expect_length(plot.list.with.tree, 1)
  expect_length(plot.list.mismatched, 1)
})
