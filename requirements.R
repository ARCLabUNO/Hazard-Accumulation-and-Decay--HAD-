## =============================================================================
## requirements.R
##
## Installs the R packages this repository's scripts need. Run this once
## before any of them: Rscript requirements.R (or source("requirements.R")
## from an R console / RStudio).
##
## `haven` is only needed if you run 00_prepare_public_data.R to convert a
## raw SPSS (.sav) extract into this repo's public CSV schema -- it is a
## free, open-source R package and reading a .sav file with it does NOT
## require an SPSS license. 01_had_decay_pipeline.R and check_data_schema.R
## never read .sav and don't need haven at all; installed here anyway so
## this one script covers every script in the repo without you having to
## remember which one needs what.
## =============================================================================

required_packages <- c("data.table", "lme4", "haven")
optional_packages <- c("rstudioapi")   # only improves repo-root auto-detection when run from RStudio -- everything degrades gracefully without it

installed <- rownames(installed.packages())

to_install <- setdiff(required_packages, installed)
if (length(to_install) > 0) {
  cat("Installing required packages:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else {
  cat("All required packages already installed:", paste(required_packages, collapse = ", "), "\n")
}

opt_to_install <- setdiff(optional_packages, installed)
if (length(opt_to_install) > 0) {
  cat("Installing optional package(s):", paste(opt_to_install, collapse = ", "), "\n")
  tryCatch(install.packages(opt_to_install), error = function(e) cat("Optional package install failed (non-fatal):", conditionMessage(e), "\n"))
}

cat("\nDone. You can now run check_data_schema.R, then 01_had_decay_pipeline.R.\n")
