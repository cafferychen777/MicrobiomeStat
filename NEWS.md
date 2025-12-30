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
