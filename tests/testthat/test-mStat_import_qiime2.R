test_that("mStat_import_qiime2_as_data_obj works correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("biomformat")
  skip_if_not_installed("yaml")
  skip_if_not_installed("Biostrings")

  # Load required package
  require(biomformat)

  # Download test files if they don't exist
  test_dir <- file.path(tempdir(), "qiime2_test_data")
  dir.create(test_dir, showWarnings = FALSE)

  # Define URLs
  urls <- c(
    table = "https://docs.qiime2.org/2023.9/data/tutorials/moving-pictures/table.qza",
    taxonomy = "https://docs.qiime2.org/2023.9/data/tutorials/moving-pictures/taxonomy.qza",
    tree = "https://docs.qiime2.org/2023.9/data/tutorials/moving-pictures/rooted-tree.qza",
    metadata = "https://data.qiime2.org/2024.10/tutorials/moving-pictures/sample_metadata.tsv"
  )

  # Download files
  files <- list()
  for(name in names(urls)) {
    files[[name]] <- file.path(test_dir, basename(urls[[name]]))
    if(!file.exists(files[[name]])) {
      tryCatch({
        download.file(urls[[name]], files[[name]], mode = "wb")
      }, error = function(e) {
        skip(paste("Could not download test file:", name))
      })
    }
  }

  # Test full import
  test_that("Full import works with real qza files", {
    skip_if_not_installed("rhdf5")
    skip_if_not_installed("biomformat")

    message("\n=== Test Environment Info ===")
    message("R version: ", R.version.string)
    message("Platform: ", Sys.info()["sysname"])
    message("Package versions:")
    message("- biomformat: ", packageVersion("biomformat"))
    message("- rhdf5: ", packageVersion("rhdf5"))

    message("\nTest files:")
    message("OTU: ", files$table)
    message("Taxonomy: ", files$taxonomy)
    message("Metadata: ", files$metadata)
    message("Tree: ", files$tree)

    message("\nFile existence checks:")
    message("OTU file exists: ", file.exists(files$table))
    message("Taxonomy file exists: ", file.exists(files$taxonomy))
    message("Metadata file exists: ", file.exists(files$metadata))
    message("Tree file exists: ", file.exists(files$tree))

    message("\nFile sizes:")
    message("OTU file: ", file.size(files$table), " bytes")
    message("Taxonomy file: ", file.size(files$taxonomy), " bytes")
    message("Metadata file: ", file.size(files$metadata), " bytes")
    message("Tree file: ", file.size(files$tree), " bytes")

    # Import data with error handling
    data_obj <- tryCatch({
      mStat_import_qiime2_as_data_obj(
        otu_qza = files$table,
        taxa_qza = files$taxonomy,
        sam_tab = files$metadata,
        tree_qza = files$tree
      )
    }, error = function(e) {
      message("\nError during import:")
      message("Error message: ", e$message)
      message("Error class: ", class(e))
      message("Full error:")
      print(e)
      stop(e)
    })

    message("\nValidation results:")
    message("data_obj is null: ", is.null(data_obj))
    if (!is.null(data_obj)) {
      message("feature.tab dimensions: ",
              nrow(data_obj$feature.tab), " x ", ncol(data_obj$feature.tab))
      message("meta.dat dimensions: ",
              nrow(data_obj$meta.dat), " x ", ncol(data_obj$meta.dat))
      message("feature.ann dimensions: ",
              nrow(data_obj$feature.ann), " x ", ncol(data_obj$feature.ann))
    }

    message("=== End of Test Environment Info ===\n")

    expect_true(!is.null(data_obj))
    expect_true(nrow(data_obj$feature.tab) > 0)
    expect_true(ncol(data_obj$feature.tab) > 0)
  })

  # Clean up
  unlink(test_dir, recursive = TRUE)
})
