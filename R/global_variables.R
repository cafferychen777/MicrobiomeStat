# Global variables used in the package
# This file is used to declare global variables to avoid R CMD check notes

utils::globalVariables(c(
  "taxonomic_level", "Variable", "Coefficient", "Significant",
  "Adjusted.P.Value", "subject_prevalence", "group_key",
  "group_prevalence", "reject", ".mean_value", ".prevalence",
  "Q1", "Q3", "IQR", "lower_fence", "upper_fence"
))
