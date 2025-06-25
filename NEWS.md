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
