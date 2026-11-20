# Packages the app needs at run time that dependency scanners
# (renv::dependencies(), rsconnect::writeManifest()) cannot detect, because
# MsCoreUtils only loads them when the matching imputation method is chosen.
# Listing them here keeps them in manifest.json. This file is not sourced.
library(impute)  # knn imputation
library(norm)    # MLE imputation
