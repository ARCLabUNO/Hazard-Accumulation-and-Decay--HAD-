## =============================================================================
## check_data_schema.R
##
## Quick, aggregate-only sanity check for data/conviction_data.csv BEFORE
## running 01_had_decay_pipeline.R -- confirms the required columns are
## present and look like what the pipeline expects. Prints only counts,
## ranges, and a 5-row PREVIEW (not representative, not individual-level in
## any meaningful sense -- just enough to eyeball column formatting). Does
## not fit anything and does not write any cache files.
##
## If your raw records start out in SPSS, run 00_prepare_public_data.R
## first to produce data/conviction_data.csv -- this script only ever
## reads the plain CSV, never a .sav file, and needs no SPSS-related
## package at all.
##
## Run this first if you're setting the pipeline up against a new file or
## a different cohort -- it will tell you exactly what's missing or looks
## wrong before you sink time into the long Stage 1 fit.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

get_repo_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    return(dirname(script_path))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(ctx_path)) return(dirname(normalizePath(ctx_path, mustWork = FALSE)))
  }
  d <- getwd()
  for (i in 1:6) {
    if (dir.exists(file.path(d, "data")) && file.exists(file.path(d, "01_had_decay_pipeline.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  warning("Could not auto-detect the repo root -- falling back to the current working directory.")
  getwd()
}
REPO_ROOT <- get_repo_root()

RAW_DATA_FILENAME <- "conviction_data.csv"   # must match CONFIG$RAW_DATA_FILENAME in 01_had_decay_pipeline.R
RAW_FILE_PATH <- file.path(REPO_ROOT, "data", RAW_DATA_FILENAME)

REQUIRED_COLS <- c("PersonID", "ConvictionYear", "ConvictionType", "YearsSinceConviction",
                    "ReconvictedThisYear", "PriorIncarcerationYears", "AgeAtAssessment")
KNOWN_CATEGORIES <- c("Property", "Alcohol Offense", "Violent", "Other", "Drug", "Sex Offense/FTR", "Weapons", "Escape")

if (!file.exists(RAW_FILE_PATH)) {
  stop(sprintf("Can't find a data file at:\n  %s\nIf your records start out in SPSS, run 00_prepare_public_data.R first. See README.md's \"Getting the data\" section.", RAW_FILE_PATH))
}

cat(sprintf("Reading %s...\n", RAW_FILE_PATH))
dt <- fread(RAW_FILE_PATH)
setnames(dt, names(dt), trimws(names(dt)))
cat(sprintf("Read %s rows, %s columns.\n\n", format(nrow(dt), big.mark = ","), ncol(dt)))

cat("---- Required columns present? ----\n")
missing <- setdiff(REQUIRED_COLS, names(dt))
if (length(missing) == 0) cat("All of:", paste(REQUIRED_COLS, collapse = ", "), "-- PRESENT.\n") else cat("MISSING:", paste(missing, collapse = ", "), "\n")

if ("ConvictionType" %in% names(dt)) {
  cat("\n---- ConvictionType category values ----\n")
  print(table(dt$ConvictionType))
  bad_cats <- setdiff(unique(dt$ConvictionType), KNOWN_CATEGORIES)
  if (length(bad_cats) > 0) cat("UNEXPECTED category value(s) found (not one of the 8 known categories):", paste(bad_cats, collapse = ", "), "\n") else cat("All values match the expected 8-category scheme.\n")
}

if ("ReconvictedThisYear" %in% names(dt)) {
  cat("\n---- ReconvictedThisYear values ----\n")
  print(table(dt$ReconvictedThisYear, useNA = "ifany"))
  cat("Expect strictly 0/1 (no other values, minimal/no NA).\n")
}

if ("YearsSinceConviction" %in% names(dt)) {
  cat(sprintf("\n---- YearsSinceConviction range ----\nmin = %s, max = %s, NA count = %s\n", min(dt$YearsSinceConviction, na.rm = TRUE), max(dt$YearsSinceConviction, na.rm = TRUE), sum(is.na(dt$YearsSinceConviction))))
}

if ("PriorIncarcerationYears" %in% names(dt)) {
  cat(sprintf("\n---- PriorIncarcerationYears range ----\nmin = %.2f, max = %.2f, NA count = %s\n", min(dt$PriorIncarcerationYears, na.rm = TRUE), max(dt$PriorIncarcerationYears, na.rm = TRUE), sum(is.na(dt$PriorIncarcerationYears))))
}

if ("AgeAtAssessment" %in% names(dt)) {
  cat(sprintf("\n---- AgeAtAssessment range ----\nmin = %.1f, max = %.1f, NA count = %s\n", min(dt$AgeAtAssessment, na.rm = TRUE), max(dt$AgeAtAssessment, na.rm = TRUE), sum(is.na(dt$AgeAtAssessment))))
}

## Informational only -- this cannot VERIFY row order is correct (the
## public schema has no full date, only a bare year, so there's nothing
## finer to check against), but it can at least tell you how many
## same-year ties exist, i.e. how many places 01_had_decay_pipeline.R's
## row-order-dependent tie-breaking (see README.md's "Row order matters"
## note) actually comes into play in this file.
if (all(c("PersonID", "ConvictionYear", "YearsSinceConviction") %in% names(dt))) {
  cat("\n---- Same-year conviction row-order note ----\n")
  event_level <- dt[YearsSinceConviction == 1, .N, by = .(PersonID, ConvictionYear)]
  n_tied_groups <- sum(event_level$N > 1)
  if (n_tied_groups > 0) {
    cat(sprintf("%s (PersonID, ConvictionYear) combination(s) have 2+ distinct convictions in the exact same calendar year. 01_had_decay_pipeline.R relies on THIS FILE'S OWN row order to break ties between them correctly. This script can only report how many such ties exist, not verify the order is correct -- if you built or re-sorted this CSV yourself (rather than using 00_prepare_public_data.R's output untouched), double-check README.md's row-order note before trusting results.\n", format(n_tied_groups, big.mark = ",")))
  } else {
    cat("No person has 2+ distinct convictions in the same calendar year in this file -- the row-order tie-breaking rule documented in README.md doesn't come into play here.\n")
  }
}

cat("\n---- Preview (first 5 rows -- not representative) ----\n")
preview_cols <- intersect(c(REQUIRED_COLS, "Gender", "Race", "Ethnicity"), names(dt))
print(dt[1:min(5, .N), ..preview_cols])

cat("\nDone. If everything above looks right, you're ready to run 01_had_decay_pipeline.R.\n")
