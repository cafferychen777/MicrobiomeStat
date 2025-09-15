# MicrobiomeStat 1.4.1

## Bug Fixes

### Fixed Jensen-Shannon divergence calculation with TSS normalization
* **Resolved Issue #45: NA values in Jensen-Shannon divergence after TSS normalization**:
  - Fixed Kullback-Leibler divergence function to properly handle zero values in probability distributions
  - Added validation to ensure both distributions have positive values before log calculation
  - Implemented safety checks to prevent NA/infinite values in distance matrices
  - Added comprehensive test coverage for Jensen-Shannon divergence with various edge cases
  - Enhanced error handling and user warnings for invalid distance calculations

### Fixed vignette execution errors
* **Resolved vignette build failures in R CMD check**:
  - Fixed duplicate factor level error in `mStat_process_time_variable()` function
  - Added proper handling of duplicate values in `t0.level` and `ts.levels` parameters
  - Improved type conversion for numeric time variables with character level specifications
  - Fixed CRAN mirror issues in vignette installation examples by adding `eval=FALSE`
  - Updated longitudinal analysis vignette examples to use correct parameter combinations

## Additional Improvements

### Complete Normalization Optimization for Barplot Functions
* **Extended normalization skip to proportion data type**:
  - All barplot functions (single, pair, long) now skip normalization for proportion data
  - ggplot2's `position="fill"` handles proportion conversion automatically
  - Eliminates redundant TSS normalization for proportion data
  - Preserves pre-computed `feature.agg.list` for all data types
  - Improved user messaging for proportion data type

### Consistent Behavior Across All Barplot Functions
* **Unified normalization strategy**:
  - `count` data: Skip normalization, let ggplot2 handle conversion
  - `proportion` data: Skip additional normalization, already in proportion format
  - `other` data: Skip normalization, user-processed data
  - All functions now have identical normalization handling logic

### Performance and Compatibility Enhancements
* **Optimal performance**: Eliminated all unnecessary normalization steps in barplot functions
* **Enhanced fake feature.ann compatibility**: All barplot functions work seamlessly with placeholder annotations
* **Data integrity**: Consistent preservation of pre-computed aggregations across all functions
* **Better user experience**: Clear, informative messages for all data types

---

# MicrobiomeStat 1.4.0

## Major Improvements

### Enhanced Normalization Handling in Visualization and Analysis Functions
* **Optimized normalization strategy for barplot visualization**:
  - Removed unnecessary TSS normalization in `generate_taxa_barplot_single` for count data
  - ggplot2's `position="fill"` automatically handles proportion conversion for barplots
  - Preserves pre-computed `feature.agg.list` to avoid re-aggregation with fake `feature.ann`
  - Improved performance by eliminating redundant normalization steps
  - Added informative messages about normalization handling

* **Fixed normalization issues in differential abundance testing**:
  - Removed redundant pre-normalization in `generate_taxa_test_single` for count data
  - LinDA function handles normalization internally based on `feature.dat.type` parameter
  - Fixed hardcoded `feature.dat.type = "proportion"` in LinDA calls - now passes user-specified type
  - Preserves pre-computed `feature.agg.list` when using fake `feature.ann` scenarios
  - Ensures appropriate zero-value handling based on actual data type

### Technical Enhancements
* **Improved parameter validation**:
  - Added `match.arg()` validation for `feature.dat.type` parameter
  - Better error handling and user feedback

* **Enhanced compatibility with fake feature annotation scenarios**:
  - Functions now work correctly when using placeholder `feature.ann` with pre-computed `feature.agg.list`
  - Prevents incorrect re-aggregation that could lead to wrong analysis results
  - Maintains data integrity across different usage patterns

### Bug Fixes
* **Fixed color mapping issues in barplot visualization**:
  - Improved color palette assignment for better ggplot2 compatibility
  - Reduced warning messages about color scale mismatches

---

# MicrobiomeStat 1.3.9

## Bug Fixes

### Fixed distance filtering in `generate_beta_test_single`
* **Fixed issue where all distances in `dist.obj` were tested instead of only those specified in `dist.name`**:
  - Added filtering to ensure only requested distance metrics are processed when `dist.obj` is provided
  - Added validation to check if all requested distances are available in `dist.obj`
  - Clear error messages when requested distances are not available
  - Improved code clarity with better variable naming (snake_case)
  - Fixed logical operators from `&` to `&&` for proper scalar evaluation

---

# MicrobiomeStat 1.3.6

## Major Enhancements

### Enhanced Data Aggregation Function
* **Completely redesigned `mStat_aggregate_data` function**:
  - Added `meta.handle.conflict` parameter with three strategies:
    - `"first"` (default): Use first record's metadata, issue warnings for conflicts
    - `"stop"`: Stop execution immediately when metadata conflicts are detected
    - `"summarise"`: Calculate mean for numeric variables, check consistency for non-numeric
  - Fixed hardcoded `subject.var` limitation - now accepts any column name as subject variable
  - Preserved original feature names during aggregation (no longer modified)
  - Enhanced conflict detection with detailed reporting of which groups and variables have conflicts
  - Improved data integrity validation and error handling

### Technical Improvements
* **Robust conflict detection**: Comprehensive metadata consistency checking across all variables
* **Detailed error reporting**: Clear identification of conflicting groups and variables
* **Backward compatibility**: All existing code continues to work with default settings
* **Enhanced validation**: Stricter data integrity checks prevent silent errors

### User Experience
* **Flexible subject variable**: No longer restricted to columns named "subject"
* **Informative warnings**: Detailed conflict reports help users identify data quality issues
* **Professional error handling**: Clear, actionable error messages for data problems
* **Strategy-based processing**: Choose the most appropriate conflict handling for your analysis

## Notes
This major update transforms `mStat_aggregate_data` from a basic aggregation tool into a sophisticated, configurable data processing function. The new conflict handling strategies ensure data quality and provide transparency in how metadata inconsistencies are managed, making the entire analysis pipeline more robust and trustworthy.

---

# MicrobiomeStat 1.3.3

## Bug Fixes and Improvements

### Negative Data Handling
* **Fixed negative value handling in `generate_taxa_heatmap_single`**:
  - Added automatic detection of negative values in "other" data type
  - Disabled abundance filtering when negative values are detected
  - Updated prevalence calculation to use `!= 0` instead of `> 0` for proper negative value handling
  - Enhanced `mStat_filter` function to support `-Inf` abundance filter for negative data

* **Fixed negative value handling in `generate_taxa_test_single`**:
  - Added automatic detection of negative values in "other" data type
  - Implemented linear model (lm) analysis for "other" data type instead of LinDA
  - Updated prevalence calculation to properly handle negative values
  - Added `perform_lm_analysis` helper function for robust linear model analysis

### Enhanced Features
* **Improved data type handling**: Functions now automatically detect and appropriately handle log-transformed or other preprocessed data containing negative values
* **User-friendly messaging**: Added informative messages when negative values are detected and filtering is adjusted
* **Backward compatibility**: All changes maintain full compatibility with existing positive data workflows

### Technical Improvements
* Enhanced error handling and validation for edge cases
* Improved code documentation and inline comments
* Added comprehensive test coverage for negative data scenarios

## Notes
These improvements specifically address issues encountered when analyzing log-transformed microbiome data or other preprocessed data containing negative values. The package now provides robust analysis capabilities for both traditional count data and transformed abundance data.

---

# MicrobiomeStat 1.3.2

Previous version features and bug fixes...

---

# MicrobiomeStat 1.3.1

Previous version features and bug fixes...

---

# MicrobiomeStat 1.3.0

Previous version features and bug fixes...
