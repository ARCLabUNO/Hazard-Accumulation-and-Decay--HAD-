## =============================================================================
## 00_prepare_public_data.R
##
## RUN THIS ONCE, LOCALLY, BEFORE check_data_schema.R / 01_had_decay_
## pipeline.R -- it is NOT part of the pipeline those two scripts run, and
## it is the only script in this repository that reads SPSS (.sav) or
## touches a raw date-of-birth field.
##
## WHAT IT DOES: reads your own raw conviction-event file (typically an
## SPSS export) and writes a single, clean, public-schema CSV --
## data/conviction_data.csv by default -- that the rest of this repository
## reads from then on. Along the way it:
##   1. Keeps ONLY the columns the pipeline actually uses (drops everything
##      else your raw extract may carry -- case numbers, court codes,
##      addresses, free-text fields, anything not listed below).
##   2. Renames those columns to a plain-English public schema (see the
##      table in README.md) instead of the internal names the manuscript's
##      own working pipeline uses.
##   3. Converts a raw offense date into a bare conviction YEAR (no month/
##      day), and a birthdate (if that's what your file has, instead of an
##      already-computed age) into AGE AT ASSESSMENT (age at conviction
##      PLUS the number of years elapsed as of this row, i.e. the person's
##      age during the specific year this row represents) -- so no literal
##      birthdate ever appears in the file this repository or anyone
##      downstream of it reads. Age at conviction on its own is never
##      written out, since the pipeline and manuscript never use it by
##      itself -- only age at the time of each assessment.
##
## ON SPSS: this script uses the `haven` R package to read a .sav file.
## haven is a free, open-source R package -- reading a .sav file with it
## does NOT require an SPSS license or installation. The reason this
## conversion step exists isn't to work around SPSS licensing; it's so
## that (a) the file everyone else interacts with is a plain CSV anyone
## can open in Excel/a text editor without any statistical software at
## all, (b) it carries only the columns actually needed, under names that
## are self-explanatory without this repo's own internal jargon, and (c)
## a birthdate never has to leave your machine in any form.
##
## INPUT: edit CONFIG$RAW_SAV_PATH below to point at your real raw file
## (an absolute path is simplest -- this script does not assume it lives
## inside this repository, since real criminal-justice microdata should
## stay wherever your own data-handling rules require it to live, which is
## almost never a git repository's own folder). Accepts either a .sav file
## or a .csv/.tsv if your raw extract is already tabular text.
##
## REQUIRED RAW COLUMNS (internal/working names, i.e. what your raw file's
## own column headers are expected to be called -- see README.md for how
## these map to the public names this script writes out):
##   ResearchID, OffDate, CharType2, LookBack, Recid, PriorIncarYearsCount,
##   and EITHER AgeAtOffense OR Dob (a birthdate), plus Gender, Race,
##   Ethnicity.
##
## OUTPUT: data/conviction_data.csv (relative to this repo's own root,
## auto-detected the same way every other script here does) with columns:
##   PersonID, ConvictionYear, ConvictionType, YearsSinceConviction,
##   ReconvictedThisYear, PriorIncarcerationYears, AgeAtAssessment, Gender,
##   Race, Ethnicity.
## Rows are written in each person's TRUE chronological offense order
## (sorted by the real date, before it's reduced to a bare year) -- don't
## reshuffle this file afterward; 01_had_decay_pipeline.R's handling of
## people with 2+ convictions in the same calendar year relies on this
## row order to break that tie correctly. See the note in that script's
## own build_event_level_full().
##
## STATUS: synthetic-tested end-to-end (fake data mimicking the real raw
## schema) before delivery. NOT YET run against the real raw file -- run
## it, look at the printed summary/preview, and confirm it looks right
## before trusting data/conviction_data.csv.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
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

CONFIG <- list(
  RAW_SAV_PATH    = "/path/to/your/raw_conviction_extract.sav",   # <-- EDIT THIS to your real file's absolute path
  OUTPUT_CSV_PATH = file.path(REPO_ROOT, "data", "conviction_data.csv")
)

KNOWN_CATEGORIES <- c("Property", "Alcohol Offense", "Violent", "Other", "Drug", "Sex Offense/FTR", "Weapons", "Escape")

log_line <- function(fmt, ...) cat(sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt, "\n"), ...))
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE)

if (!file.exists(CONFIG$RAW_SAV_PATH)) {
  stop(sprintf("Can't find your raw file at:\n  %s\nEdit CONFIG$RAW_SAV_PATH at the top of this script to point at your real file.", CONFIG$RAW_SAV_PATH))
}

## ---------------------------------------------------------------------------
## STEP 1 -- read the raw file (.sav or plain .csv/.tsv).
## ---------------------------------------------------------------------------
ext <- tolower(tools::file_ext(CONFIG$RAW_SAV_PATH))
log_line("Reading raw file: %s", CONFIG$RAW_SAV_PATH)
if (ext == "sav") {
  d <- as.data.table(haven::read_sav(CONFIG$RAW_SAV_PATH))
} else if (ext %in% c("csv", "txt", "tsv")) {
  d <- fread(CONFIG$RAW_SAV_PATH)
} else {
  stop(sprintf("Unrecognized file extension '%s' -- expected .sav, .csv, or .tsv.", ext))
}
setnames(d, names(d), trimws(names(d)))
log_line("Raw file read: %s rows, %s columns.", fmt_n(nrow(d)), ncol(d))

## ---------------------------------------------------------------------------
## STEP 2 -- validate presence of required raw columns. Age-source presence
## is checked per-COLUMN here (at least one of AgeAtOffense/Dob has to
## exist at all) -- the per-ROW fallback (what happens when AgeAtOffense
## exists as a column but is NA for some rows while Dob is also present)
## is handled in STEP 5 below, not here.
## ---------------------------------------------------------------------------
has_age_col <- "AgeAtOffense" %in% names(d)
has_dob_col <- "Dob" %in% names(d)
required_base <- c("ResearchID", "OffDate", "CharType2", "LookBack", "Recid", "PriorIncarYearsCount", "Gender", "Race", "Ethnicity")
if (!has_age_col && !has_dob_col) {
  stop("Raw file is missing required column(s): AgeAtOffense or Dob (need at least one).")
}
missing_cols <- setdiff(required_base, names(d))
if (length(missing_cols) > 0) {
  stop(sprintf("Raw file is missing required column(s): %s.", paste(missing_cols, collapse = ", ")))
}
log_line("All required raw columns present. Age column(s) found: %s.",
         paste(c(if (has_age_col) "AgeAtOffense", if (has_dob_col) "Dob"), collapse = " and "))

## ---------------------------------------------------------------------------
## STEP 3 -- normalize CharType2 (handle SPSS labelled values) and validate
## its category values.
## ---------------------------------------------------------------------------
if (haven::is.labelled(d$CharType2)) d$CharType2 <- as.character(haven::as_factor(d$CharType2))
d[, CharType2 := as.character(CharType2)]
bad_categories <- setdiff(unique(d$CharType2), KNOWN_CATEGORIES)
if (length(bad_categories) > 0) {
  stop(sprintf("Found CharType2 value(s) not in the expected 8-category scheme (%s): %s", paste(KNOWN_CATEGORIES, collapse = "/"), paste(bad_categories, collapse = ", ")))
}

## ---------------------------------------------------------------------------
## STEP 4 -- derive a TRUE sortable offense date (used only to order rows
## correctly before we reduce OffDate down to a bare year below -- see the
## header note and 01_had_decay_pipeline.R's build_event_level_full()).
## Handles the same raw-date encodings the manuscript's own working
## pipeline does: an R Date, an SPSS labelled numeric, or a bare year.
## ---------------------------------------------------------------------------
if (inherits(d$OffDate, "Date")) {
  d[, OffDateSortKey := OffDate]
  d[, OffYear := as.integer(format(OffDate, "%Y"))]
} else if (haven::is.labelled(d$OffDate)) {
  d[, OffDateSortKey := as.numeric(haven::zap_labels(OffDate))]
  d[, OffYear := as.integer(haven::zap_labels(OffDate))]
} else {
  d[, OffDateSortKey := as.numeric(OffDate)]
  d[, OffYear := as.integer(OffDate)]
}

## ---------------------------------------------------------------------------
## STEP 5 -- make sure AgeAtOffense is populated, deriving it from Dob
## PER-ROW wherever it's actually missing -- not only when the AgeAtOffense
## column is absent entirely. A raw extract merged or extended across
## data-entry eras can plausibly have AgeAtOffense directly for most rows
## but only a birthdate for a subset; the previous file-level check never
## attempted the Dob fallback for those rows once the AgeAtOffense column
## existed at all, so any row-level gaps stayed NA even when a usable Dob
## was sitting right there. Same branching as the manuscript's own
## clean_and_derive() (Date / POSIXct-with-1900-sentinel / labelled / bare
## birth year), just scoped to the rows that need it.
## ---------------------------------------------------------------------------
if (!has_age_col) {
  d[, AgeAtOffense := NA_real_]
} else {
  d[, AgeAtOffense := as.numeric(AgeAtOffense)]
}
n_age_missing_before_fallback <- sum(is.na(d$AgeAtOffense))
if (has_dob_col && n_age_missing_before_fallback > 0) {
  log_line("AgeAtOffense is missing for %s row(s) -- deriving it from Dob as OffYear - year(Dob) (calendar-year approximation, not exact-birthday age) for just those rows.",
           fmt_n(n_age_missing_before_fallback))
  missing_idx <- which(is.na(d$AgeAtOffense))
  dob_missing <- d$Dob[missing_idx]
  off_year_missing <- d$OffYear[missing_idx]
  derived_age <- if (inherits(d$Dob, "Date")) {
    off_year_missing - as.integer(format(dob_missing, "%Y"))
  } else if (inherits(d$Dob, "POSIXct") || inherits(d$Dob, "POSIXt")) {
    dob_secs <- as.numeric(dob_missing)
    sentinel_secs <- as.numeric(as.POSIXct("1900-01-01 00:00:00", tz = "UTC"))
    is_sentinel <- !is.na(dob_secs) & abs(dob_secs - sentinel_secs) < 86400
    if (sum(is_sentinel) > 0) dob_missing[is_sentinel] <- as.POSIXct(NA)
    off_year_missing - as.integer(format(dob_missing, "%Y"))
  } else if (haven::is.labelled(d$Dob)) {
    off_year_missing - as.integer(haven::zap_labels(dob_missing))
  } else if (is.numeric(d$Dob) && median(as.numeric(d$Dob), na.rm = TRUE) > 1900 && median(as.numeric(d$Dob), na.rm = TRUE) < 2100) {
    off_year_missing - as.integer(dob_missing)
  } else {
    stop("Dob is numeric but its typical value is not a plausible birth year. Confirm Dob's real encoding before proceeding.")
  }
  d[missing_idx, AgeAtOffense := derived_age]
}
n_age_missing_final <- sum(is.na(d$AgeAtOffense))
if (n_age_missing_final > 0) {
  log_line("NOTE: %s row(s) still have no AgeAtOffense (no usable Dob available for them either). These will carry an NA AgeAtAssessment and will later be dropped by 01_had_decay_pipeline.R's under-18-at-assessment exclusion -- NOT because they are confirmed minors, but because their age is genuinely unknown.",
           fmt_n(n_age_missing_final))
}

## ---------------------------------------------------------------------------
## STEP 6 -- coerce plain types, clean PriorIncarYearsCount. The pre-
## imputation NA count is captured here (not re-derived after the fact in
## STEP 9) since PriorIncarYearsCount is overwritten with 0 below -- a
## diagnostic computed afterward would always read 0 by construction,
## regardless of how many rows were actually missing.
## ---------------------------------------------------------------------------
d[, ResearchID := as.character(ResearchID)]
d[, LookBack := as.numeric(LookBack)]
d[, Recid := as.integer(Recid)]
if (haven::is.labelled(d$PriorIncarYearsCount)) d$PriorIncarYearsCount <- as.numeric(haven::zap_labels(d$PriorIncarYearsCount))
attributes(d$PriorIncarYearsCount) <- NULL
n_incar_na_before_impute <- sum(is.na(d$PriorIncarYearsCount))
if (n_incar_na_before_impute > 0) d[is.na(PriorIncarYearsCount), PriorIncarYearsCount := 0]

## ---------------------------------------------------------------------------
## STEP 6b -- validate Recid is a strict 0/1 indicator, the same way
## CharType2 is validated above. A miscoded Recid (e.g. a reconviction
## COUNT instead of an indicator) would otherwise sail through this script
## silently and only surface later as an opaque crash deep inside
## 01_had_decay_pipeline.R's glm()/glmer() calls.
## ---------------------------------------------------------------------------
bad_recid_vals <- setdiff(unique(d$Recid), c(0L, 1L, NA_integer_))
if (length(bad_recid_vals) > 0) {
  stop(sprintf("Recid contains value(s) other than 0/1: %s. This column must be a strict binary reconviction indicator -- check for a miscoded count or categorical value before proceeding.", paste(bad_recid_vals, collapse = ", ")))
}
n_recid_na <- sum(is.na(d$Recid))
if (n_recid_na > 0) {
  stop(sprintf("Recid has %s NA/missing value(s) -- this column must be a complete 0/1 indicator with no missing values.", fmt_n(n_recid_na)))
}

## ---------------------------------------------------------------------------
## STEP 7 -- sort by each person's TRUE chronological offense order (full
## date resolution), THEN drop down to a bare year. This is what lets
## 01_had_decay_pipeline.R correctly break ties for people with 2+
## convictions in the same calendar year, even though the public CSV only
## ever shows the year -- see the header note above.
## ---------------------------------------------------------------------------
setorder(d, ResearchID, OffDateSortKey)

## ---------------------------------------------------------------------------
## STEP 8 -- select, rename to the public schema, and write. AgeAtAssessment
## (age at conviction PLUS years elapsed as of this row) is written out
## directly rather than AgeAtConviction on its own -- that's the only age
## value 01_had_decay_pipeline.R actually uses (the under-18-at-assessment
## exclusion, the age generalizability table, and Table 1's age
## descriptives all key off age AT THE ASSESSMENT, not age at conviction).
## ---------------------------------------------------------------------------
out <- d[, .(
  PersonID                 = ResearchID,
  ConvictionYear           = OffYear,
  ConvictionType           = CharType2,
  YearsSinceConviction     = LookBack,
  ReconvictedThisYear      = Recid,
  PriorIncarcerationYears  = PriorIncarYearsCount,
  AgeAtAssessment          = AgeAtOffense + LookBack,
  Gender                   = Gender,
  Race                     = Race,
  Ethnicity                = Ethnicity
)]

dir.create(dirname(CONFIG$OUTPUT_CSV_PATH), showWarnings = FALSE, recursive = TRUE)
fwrite(out, CONFIG$OUTPUT_CSV_PATH)
log_line("Wrote %s rows, %s columns, %s unique people to %s.", fmt_n(nrow(out)), ncol(out), fmt_n(uniqueN(out$PersonID)), CONFIG$OUTPUT_CSV_PATH)

## ---------------------------------------------------------------------------
## STEP 9 -- print a summary so you can sanity-check before trusting it.
## ---------------------------------------------------------------------------
cat("\n---- ConvictionType value counts ----\n")
print(table(out$ConvictionType))

cat("\n---- ReconvictedThisYear value counts (expect strictly 0/1) ----\n")
print(table(out$ReconvictedThisYear, useNA = "ifany"))

cat(sprintf("\n---- YearsSinceConviction range ----\nmin = %s, max = %s, NA count = %s\n",
            min(out$YearsSinceConviction, na.rm = TRUE), max(out$YearsSinceConviction, na.rm = TRUE), sum(is.na(out$YearsSinceConviction))))

cat(sprintf("\n---- AgeAtAssessment range ----\nmin = %.1f, max = %.1f, NA count = %s\n",
            min(out$AgeAtAssessment, na.rm = TRUE), max(out$AgeAtAssessment, na.rm = TRUE), sum(is.na(out$AgeAtAssessment))))

cat(sprintf("\n---- PriorIncarcerationYears range ----\nmin = %.2f, max = %.2f, %s row(s) were missing and imputed to 0 (this column cannot have NA by this point in the script)\n",
            min(out$PriorIncarcerationYears, na.rm = TRUE), max(out$PriorIncarcerationYears, na.rm = TRUE), fmt_n(n_incar_na_before_impute)))

cat("\n---- Preview (first 5 rows -- not representative) ----\n")
print(out[1:min(5, .N)])

cat(sprintf("\nDone. %s is ready for check_data_schema.R and 01_had_decay_pipeline.R.\n", CONFIG$OUTPUT_CSV_PATH))
cat("No birthdate, case number, or any column beyond the 10 above was written to this file.\n")
