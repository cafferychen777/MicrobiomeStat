# MicrobiomeStat 1.5.0

## New Features

* Added a unified analysis interface and expanded data handling options, including explicit support for `"other"` feature data flows and consistent feature-level requirements across analysis paths.
* Added `linda2()` with detection-depth weighting and tree-guided smoothing support, and exported it for direct package use.
* Added `ref.level` support to taxa differential testing helpers so reference groups can be chosen explicitly instead of relying on factor ordering.

## Statistical and Workflow Fixes

* Fixed time semantics, metadata alignment, and interface contract regressions across single, paired, and longitudinal workflows.
* Fixed multiple plotting and report-generation regressions so examples, reports, and exported figures remain runnable and consistent with package defaults.
* Fixed remaining R CMD check and dependency handling issues, including a cleanup of GitHub Actions dependency resolution for package checks.

## Internal Improvements

* Refactored taxa helpers, filter handling, and report-section generators to reduce duplication and make plotting/report code paths more consistent.
* Centralized documentation templates and regenerated Rd files to keep exported interfaces aligned with implementation changes.

# MicrobiomeStat 1.4.4

## New Features

* Added `ref.level` parameter to differential abundance testing functions for specifying reference group:
  - `generate_taxa_test_single()`: Single time point differential abundance testing
  - `generate_taxa_test_pair()`: Paired/longitudinal differential abundance testing
  - `generate_taxa_trend_test_long()`: Longitudinal trend testing
  - `generate_taxa_change_test_pair()`: Change score analysis between time points

* The `ref.level` parameter allows users to specify which group level should be used as the reference for comparisons, instead of relying on alphabetical ordering (first level alphabetically). This provides more flexibility and control over the statistical comparisons.

* Example usage:
  ```r
  # Set "Control" as reference instead of alphabetically first group

  test.list <- generate_taxa_test_single(
    data.obj = data.obj,
    group.var = "treatment",
    ref.level = "Control",  # New parameter
    feature.level = c("Genus"),
    ...
  )
  ```

---

# MicrobiomeStat 1.4.3

## Bug Fixes

* Fixed variable type handling in `generate_taxa_test_single()`: The function now correctly distinguishes between categorical (factor/character) and continuous (numeric/integer) variables. Previously, all variables were unconditionally treated as categorical.

## Documentation Improvements

* Updated roxygen documentation for `generate_taxa_test_single()`:
  - Fixed `@return` description to accurately reflect actual output columns
  - Added comprehensive explanation of variable type handling (categorical vs continuous)
  - Added detailed statistical methods documentation
  - Clarified difference between LinDA (for count/proportion data) and linear models (for other data)

## Code Quality

* Removed legacy code from `generate_taxa_volcano_single()`:
  - Removed unused `meta_tab` variable extraction
  - Removed unused `group_level` and `reference_level` variable definitions
  - Code is now cleaner and more maintainable

## Testing

* Verified continuous variable support in volcano plots
* All tests pass with both categorical and continuous variables

---

# MicrobiomeStat 1.4.2

## Bug Fixes

* Fixed pseudocount imputation for taxa change analysis
* Fixed log fold change pseudocount bias in taxa change functions
* Fixed critical statistical issues in normalization and CLR transformation
* Fixed barplot/areaplot mean calculation

## New Features

* Added linda2() function with experimental sample weighting support
* Exported linda2 function for external use

---

# MicrobiomeStat 1.4.1 and earlier

See git commit history for details.
