## =============================================================================
## 01_had_decay_pipeline.R
##
## Public replication script for "Half-Life of a Conviction: A Decay-
## Weighted Approach to Criminal History Risk Assessment" (Hamilton et al.,
## Criminology), using the manuscript's 8-category offense classification
## (CharType2: Property, Alcohol Offense, Violent, Other, Drug, Sex
## Offense/FTR, Weapons, Escape).
##
## WHAT THIS SCRIPT DOES, START TO FINISH:
##   0. Config and shared helper functions.
##   1. Loads the raw conviction-event file and derives the analytic columns
##      (offense category, calendar year, etc.). Age arrives already
##      precomputed as age at THIS row's assessment (not age at
##      conviction, which this pipeline never uses on its own).
##   2. Applies each offense category's own "active window" (how many years
##      after a conviction it remains part of a person's active criminal
##      history) and splits people into a 75% development / 25% held-out
##      sample, with no person appearing on both sides.
##   3. Fits the single-event decay curve (Stage 1): a mixed-effects
##      logistic regression of one-year reconviction on time-since-
##      conviction, interacted with offense category, with a person-level
##      random intercept.
##   4. Builds the full decay-weighted criminal history score (HAD) by
##      summing each active conviction's own decay-adjusted hazard, with an
##      empirical correction for the redundancy introduced when multiple
##      convictions are active at once.
##   5. Calibrates HAD's predicted probabilities (Platt scaling) and reports
##      its held-out predictive performance (AUC, RMSE, calibration) against
##      two comparators: a static logistic-regression score and a simple,
##      unit-weighted ("Burgess-style") point count.
##   6. Reports HAD's generalizability across age, gender, race/ethnicity,
##      and prior-incarceration exposure (held-out sample).
##   7. Reports full descriptive statistics for the analytic sample (Table 1).
##   8. Reports a three-year supplementary predictive-performance comparison
##      (Appendix Table A2 analog).
##
## OUTPUT: every table this script produces is an AGGREGATE table (sample
## sizes, fixed-effects estimates, AUC/calibration summaries, etc.) written
## to the output/ folder as CSV. No individual-level rows, and no p-values,
## are written or printed anywhere in this file -- by design, not omission.
## This matches how the manuscript itself reports results and keeps this
## script safe to run against sensitive criminal-history microdata.
##
## DATA: this script expects a single, plain conviction-event-year CSV (see
## README.md for the exact required columns and their meaning) at
## data/<RAW_DATA_FILENAME>. If your own records live in SPSS (.sav) first,
## run 00_prepare_public_data.R once, locally, to convert and rename them
## into this CSV schema -- see that script and README.md's "Data" section.
## This script itself never reads .sav and never needs SPSS or the haven
## package. The real, fully de-identified data behind the manuscript is NOT
## included directly in this git repository -- not because of any access
## restriction, but because it's several gigabytes at real-data scale, too
## large for a normal git repo. See README.md's "Data" section for the
## external download link. A small synthetic file with the right shape is
## NOT included either -- this script will tell you exactly what's missing
## if the data file isn't found or doesn't have the right columns.
##
## RUNTIME: Sections 2-4 (load, clean, split) are fast. Section 5 (the
## Stage 1 mixed-effects fit) is the slow step -- on the manuscript's real,
## ~70-million-row file this took roughly 8 hours on an ordinary desktop.
## Every expensive step is cached under .pipeline_cache/ (gitignored) and
## keyed by PIPELINE_VERSION/ROW_CAP, so re-running this script after a
## first successful run reloads cached results in seconds rather than
## recomputing them -- change PIPELINE_VERSION below if you want to force a
## full rebuild (e.g., after changing the input data).
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

## =============================================================================
## SECTION 0 -- CONFIG AND REPO-ROOT DETECTION
## =============================================================================

## Locates this repo's own root directory regardless of how the script is
## run -- from the command line (Rscript), from RStudio's Source button, or
## via source() from a console already positioned somewhere inside the
## repo. Means you never need to setwd() by hand before running this.
get_repo_root <- function() {
  ## Strategy 1: command-line invocation (Rscript path/to/this_file.R)
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)
    return(dirname(script_path))
  }
  ## Strategy 2: RStudio's Source button / Ctrl+Shift+Enter
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(ctx_path)) return(dirname(normalizePath(ctx_path, mustWork = FALSE)))
  }
  ## Strategy 3: walk up from the current working directory looking for
  ## this repo's own marker (a "data" folder sitting next to this script).
  d <- getwd()
  for (i in 1:6) {
    if (dir.exists(file.path(d, "data")) && file.exists(file.path(d, "01_had_decay_pipeline.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  warning("Could not auto-detect the repo root -- falling back to the current working directory. If paths below fail, run this script from the repo's own folder, or open it via RStudio's Source button.")
  getwd()
}
REPO_ROOT <- get_repo_root()

CONFIG <- list(
  RAW_DATA_FILENAME  = "conviction_data.csv",   # the file 00_prepare_public_data.R produces -- change this if your file has a different name, see README.md
  CACHE_DIR          = file.path(REPO_ROOT, ".pipeline_cache"),
  OUTPUT_DIR         = file.path(REPO_ROOT, "output"),
  PIPELINE_VERSION   = "v1_chartype2",   # bump this to force a full rebuild (invalidates all caches below)
  STUDY_YEAR_MIN     = 1980L,
  STUDY_YEAR_MAX     = 2021L,
  MIN_ASSESSMENT_AGE = 18,
  TRAIN_FRACTION     = 0.75,
  SPLIT_SEED         = 20260924,
  ROW_CAP_SEED       = 20260926,
  CORRECTION_SEED    = 20260927,
  ROW_CAP            = Inf,   # set to a finite number ONLY for a quick sanity-check run on a subset -- uses its own, separately-tagged cache files
  N_AVERAGING_DRAWS  = 1000,
  FOLLOW_YEARS       = 3L,
  BURGESS_HEADLINE_CAP = 3   # prior-conviction cap for the Burgess-style comparator, matching LSI-R's own cap
)
CONFIG$RAW_FILE_PATH <- file.path(REPO_ROOT, "data", CONFIG$RAW_DATA_FILENAME)
dir.create(CONFIG$CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

## The manuscript's 8 offense categories and their empirically-derived
## active windows (how many years after a conviction it stays part of a
## person's active criminal history). See the manuscript's Appendix Table
## A1 and README.md's "Methodology notes" section for how these were
## derived.
ACTIVE_WINDOWS <- c(Property = 13, `Alcohol Offense` = 12, Violent = 13, Other = 12,
                     Drug = 14, `Sex Offense/FTR` = 14, Weapons = 13, Escape = 15)
KNOWN_CATEGORIES <- names(ACTIVE_WINDOWS)
CAT_COL   <- "CharType2"
INCAR_COL <- "PriorIncarYearsCount"

## =============================================================================
## SECTION 1 -- SHARED HELPERS
## =============================================================================

log_line <- function(fmt, ...) cat(sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt, "\n"), ...))
fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE)

## Reads the Stage 1 fit's fixed effects and returns a function computing one
## conviction's own linear predictor at a given category and elapsed-time
## value. age_z/incar_z are accepted for interface compatibility but are
## INERT here, since Stage 1's formula (below) has no CurrentAge_z/
## IncarYears_z term -- written defensively (checks whether a term exists
## before using it) so the same function would still work if a future refit
## added one.
make_linear_predictor <- function(fixed_effects, lb_mean, lb_sd) {
  function(category, lookback_year, age_z = 0, incar_z = 0) {
    lb_z <- (lookback_year - lb_mean) / lb_sd
    x_names <- names(fixed_effects)
    val <- fixed_effects["(Intercept)"]
    val <- val + fixed_effects["LookBack_z"] * lb_z
    if ("CurrentAge_z" %in% x_names) val <- val + fixed_effects["CurrentAge_z"] * age_z
    if ("IncarYears_z" %in% x_names) val <- val + fixed_effects["IncarYears_z"] * incar_z
    cat_term <- paste0(CAT_COL, category)
    if (cat_term %in% x_names) val <- val + fixed_effects[cat_term]
    inter_term <- paste0("LookBack_z:", cat_term)
    if (inter_term %in% x_names) val <- val + fixed_effects[inter_term] * lb_z
    unname(val)
  }
}

## Manual rank-based AUC (Mann-Whitney form). Counts are cast to double
## before multiplying so n_pos*n_neg can't silently integer-overflow at
## real-data scale.
compute_auc <- function(prob, outcome) {
  keep <- !is.na(prob) & !is.na(outcome)
  prob <- prob[keep]; outcome <- outcome[keep]
  n_pos <- as.double(sum(outcome == 1)); n_neg <- as.double(sum(outcome == 0))
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(prob, ties.method = "average")
  (sum(r[outcome == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

safe_logit <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  log(p / (1 - p))
}

## Maps the PUBLIC column names (what the CSV on disk actually calls its
## columns -- see README.md's "Required columns" table) onto the INTERNAL
## working names the rest of this file uses. Kept as one named translation
## table in one place rather than renaming every reference throughout this
## file, for two reasons: (1) it's far lower-risk than touching ~1000 lines
## of already-validated logic, and (2) CAT_COL's internal value ("CharType2")
## is baked into the Stage 1 model's own fitted term names the moment it's
## fit in Section 4 (e.g. "CharType2Property") -- as long as the column this
## script hands glmer() is consistently named whatever CAT_COL says, and
## every downstream lookup also goes through the same CAT_COL variable
## (which it does, throughout this file), the public-facing rename below is
## purely cosmetic for the input file and never touches model correctness.
PUBLIC_TO_INTERNAL_NAMES <- c(
  PersonID                = "ResearchID",
  ConvictionYear          = "OffDate",   # a plain year, not a real date -- see the note in clean_and_derive()
  ConvictionType          = "CharType2",
  YearsSinceConviction    = "LookBack",
  ReconvictedThisYear     = "Recid",
  PriorIncarcerationYears = "PriorIncarYearsCount",
  AgeAtAssessment         = "CurrentAge"
  ## Gender, Race, Ethnicity, ResearchID (if a file already uses the
  ## internal name) pass through unchanged -- intelligible as-is.
)

## Reads the raw conviction-event CSV/TSV produced by 00_prepare_public_
## data.R (or your own file matching that same public schema -- see
## README.md), then translates its public column names to this script's
## internal working names via PUBLIC_TO_INTERNAL_NAMES above. Used both for
## the main pipeline (Section 2/3) and for the separate demographics-only
## extract (Section 11), so there's exactly one place this project's
## raw-file format/naming logic lives. Does not read .sav -- see
## 00_prepare_public_data.R if your own records start out in SPSS.
load_raw_file <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf(
      "Can't find the data file at:\n  %s\nSee README.md's \"Getting the data\" section -- this repository does not include the underlying criminal-history microdata. If your records start out in SPSS, run 00_prepare_public_data.R once to produce this CSV; otherwise place your own CSV/TSV there (matching the required columns listed in README.md), or change CONFIG$RAW_DATA_FILENAME at the top of this script.",
      path
    ))
  }
  ext <- tolower(tools::file_ext(path))
  if (!ext %in% c("csv", "txt", "tsv")) {
    stop(sprintf("Unrecognized data file extension '%s' -- expected .csv or .tsv. This script no longer reads .sav directly -- run 00_prepare_public_data.R once to convert an SPSS file into the CSV schema documented in README.md.", ext))
  }
  log_line("Loading data file: %s", path)
  dt <- fread(path)
  setnames(dt, names(dt), trimws(names(dt)))
  setnames(dt, names(PUBLIC_TO_INTERNAL_NAMES), unname(PUBLIC_TO_INTERNAL_NAMES), skip_absent = TRUE)
  log_line("Data file loaded: %s rows, %s columns.", fmt_n(nrow(dt)), ncol(dt))
  dt
}

## =============================================================================
## SECTION 2 -- CLEAN AND DERIVE
## =============================================================================

## NOTE: with the standard public input (00_prepare_public_data.R's output
## CSV, after load_raw_file()'s translation above), every column here is
## already a plain CSV type -- ConvictionYear/OffDate is already a bare
## integer year (not a real date) and AgeAtAssessment/CurrentAge is
## already precomputed (age at conviction plus years elapsed, i.e. the
## person's age during THIS row's assessment year -- age at conviction on
## its own is never used anywhere downstream, so the public CSV carries
## the already-combined value directly instead) -- so there is no SPSS
## labelled-value handling, no Date/POSIXct parsing, and no Dob fallback
## needed here; all of that lives in 00_prepare_public_data.R instead, the
## one place it's actually needed (and the only script in this repo that
## still uses haven). This function just validates presence/category
## values and coerces plain types.
clean_and_derive <- function(dt_real, CONFIG) {
  d <- copy(dt_real)

  required <- c("ResearchID", "OffDate", CAT_COL, "LookBack", "Recid", INCAR_COL, "CurrentAge")
  missing_cols <- setdiff(required, names(d))
  if (length(missing_cols) > 0) {
    stop(sprintf("clean_and_derive: data file is missing required column(s): %s. See README.md for the full list of required columns (in their public CSV names).", paste(missing_cols, collapse = ", ")))
  }

  d[[CAT_COL]] <- as.character(d[[CAT_COL]])
  bad_categories <- setdiff(unique(d[[CAT_COL]]), KNOWN_CATEGORIES)
  if (length(bad_categories) > 0) {
    stop(sprintf("clean_and_derive: found %s value(s) not in the expected 8-category scheme (%s): %s",
                 CAT_COL, paste(KNOWN_CATEGORIES, collapse = "/"), paste(bad_categories, collapse = ", ")))
  }

  d[, OffYear := as.integer(OffDate)]
  d[, ResearchID := as.character(ResearchID)]
  d[, LookBack := as.numeric(LookBack)]
  d[, CurrentAge := as.numeric(CurrentAge)]
  d[, Recid := as.integer(Recid)]

  ## Validate Recid is a strict 0/1 indicator, the same way CAT_COL is
  ## validated above. Without this, a miscoded Recid (e.g. a reconviction
  ## COUNT instead of an indicator) sails through silently and only
  ## surfaces later as an opaque "y values must be 0 <= y <= 1" error deep
  ## inside one of this script's several glm()/glmer() calls.
  bad_recid_vals <- setdiff(unique(d$Recid), c(0L, 1L, NA_integer_))
  if (length(bad_recid_vals) > 0) {
    stop(sprintf("clean_and_derive: Recid/ReconvictedThisYear contains value(s) other than 0/1: %s. This column must be a strict binary reconviction indicator.", paste(bad_recid_vals, collapse = ", ")))
  }
  if (sum(is.na(d$Recid)) > 0) {
    stop(sprintf("clean_and_derive: Recid/ReconvictedThisYear has %s NA/missing value(s) -- this column must be a complete 0/1 indicator with no missing values.", fmt_n(sum(is.na(d$Recid)))))
  }

  if (sum(is.na(d[[INCAR_COL]])) > 0) d[is.na(get(INCAR_COL)), (INCAR_COL) := 0]

  keep_cols <- unique(c(required, "OffYear"))
  out <- d[, ..keep_cols]
  log_line("Cleaned data: %s rows, %s people, columns kept: %s", fmt_n(nrow(out)), fmt_n(uniqueN(out$ResearchID)), paste(names(out), collapse = ", "))
  out
}

## =============================================================================
## SECTION 3 -- ACTIVE WINDOW, TRAIN/TEST SPLIT (PERSON-LEVEL, LEAKAGE-CHECKED)
## =============================================================================

rowcap_tag <- if (is.finite(CONFIG$ROW_CAP)) sprintf("cap%s", format(CONFIG$ROW_CAP, scientific = FALSE)) else "nocap"
split_cache_path <- file.path(CONFIG$CACHE_DIR, sprintf("section3_split_%s_%s.rds", CONFIG$PIPELINE_VERSION, rowcap_tag))
if (file.exists(split_cache_path)) {
  log_line("Reloading the cached Section 3 output (raw load + clean + train/test split) from %s.", split_cache_path)
  split_cache <- readRDS(split_cache_path)
  dt_train <- split_cache$dt_train
  dt_test  <- split_cache$dt_test
  rm(split_cache)
} else {
  dt_real  <- load_raw_file(CONFIG$RAW_FILE_PATH)
  clean_dt <- clean_and_derive(dt_real, CONFIG)

  dt_full <- copy(clean_dt)

  n_before_period <- nrow(dt_full)
  dt_full <- dt_full[OffYear >= CONFIG$STUDY_YEAR_MIN & OffYear <= CONFIG$STUDY_YEAR_MAX]
  log_line("Applied study-period restriction (OffYear %d-%d): %s rows kept of %s (%.4f%% dropped).",
           CONFIG$STUDY_YEAR_MIN, CONFIG$STUDY_YEAR_MAX, fmt_n(nrow(dt_full)), fmt_n(n_before_period),
           100 * (n_before_period - nrow(dt_full)) / n_before_period)

  dt_full[, ActiveWindow := ACTIVE_WINDOWS[get(CAT_COL)]]
  if (sum(is.na(dt_full$ActiveWindow)) > 0) stop(sprintf("Found row(s) with a %s value not covered by ACTIVE_WINDOWS.", CAT_COL))
  n_before_window <- nrow(dt_full)
  dt_full <- dt_full[LookBack >= 1 & LookBack <= ActiveWindow]
  log_line("Applied each category's own active window: %s rows kept of %s (%.1f%%).",
           fmt_n(nrow(dt_full)), fmt_n(n_before_window), 100 * nrow(dt_full) / n_before_window)

  all_ids <- unique(dt_full$ResearchID)
  set.seed(CONFIG$SPLIT_SEED)
  split_dt <- data.table(
    ResearchID = all_ids,
    Split = sample(c("train", "test"), length(all_ids), replace = TRUE,
                   prob = c(CONFIG$TRAIN_FRACTION, 1 - CONFIG$TRAIN_FRACTION))
  )
  dt_full <- merge(dt_full, split_dt, by = "ResearchID", all.x = TRUE)
  if (sum(is.na(dt_full$Split)) > 0) stop("Internal error: unassigned rows after merging the train/test split.")
  if (length(intersect(dt_full[Split == "train", unique(ResearchID)], dt_full[Split == "test", unique(ResearchID)])) > 0) {
    stop("LEAKAGE DETECTED: a person appears on both sides of the train/test split.")
  }
  log_line("Leakage check passed: no person appears on both sides of the split.")

  dt_train_full <- dt_full[Split == "train"]
  dt_test       <- dt_full[Split == "test"]

  ## CurrentAge (age at THIS row's assessment) already arrived directly
  ## from the input file via clean_and_derive() above -- no derivation
  ## needed here. (Previously this was computed as AgeAtOffense + LookBack;
  ## the public CSV now ships the already-combined value directly, since
  ## age at conviction on its own is never used anywhere in this pipeline.)

  n_train_full <- nrow(dt_train_full)
  if (is.finite(CONFIG$ROW_CAP) && n_train_full > CONFIG$ROW_CAP) {
    train_ids_all <- unique(dt_train_full$ResearchID)
    set.seed(CONFIG$ROW_CAP_SEED)
    keep_ids <- sample(train_ids_all, size = min(length(train_ids_all),
                                                  ceiling(CONFIG$ROW_CAP / (n_train_full / length(train_ids_all)))))
    dt_train <- dt_train_full[ResearchID %in% keep_ids]
  } else {
    dt_train <- dt_train_full
  }

  for (nm in c("dt_train", "dt_test")) {
    dd <- get(nm)
    dd[, CalendarYear := OffYear + LookBack]
    dd[, NActive := .N, by = .(ResearchID, CalendarYear)]
    assign(nm, dd)
  }

  ## Guard: every category present in TEST must also have at least one row
  ## in TRAIN. Section 9's Static comparator fits a factor/character term
  ## on CAT_COL in TRAIN and predict()s it on TEST -- if a category is thin
  ## enough that this split/seed happens to assign it entirely to TEST,
  ## that predict() call fails deep in the pipeline with a generic
  ## "factor ... has new levels ..." error with no indication of which
  ## category or why. Fail fast here instead, with an actionable message.
  train_cats_check <- unique(dt_train[[CAT_COL]])
  test_cats_check  <- unique(dt_test[[CAT_COL]])
  cats_missing_from_train <- setdiff(test_cats_check, train_cats_check)
  if (length(cats_missing_from_train) > 0) {
    stop(sprintf(
      "Categor(ies) present in TEST but with zero rows in TRAIN after the train/test split: %s. This usually means a category is too sparse, relative to CONFIG$TRAIN_FRACTION and this person-level random split, to reliably land on both sides. Try a different CONFIG$SPLIT_SEED, or investigate whether this category has enough people in this dataset at all.",
      paste(cats_missing_from_train, collapse = ", ")))
  }

  saveRDS(list(dt_train = dt_train, dt_test = dt_test), split_cache_path)
  log_line("Section 3 output complete and cached to %s.", split_cache_path)
}
log_line("TRAIN (for fitting): %s rows, %s people. TEST (held out): %s rows, %s people.",
         fmt_n(nrow(dt_train)), fmt_n(uniqueN(dt_train$ResearchID)), fmt_n(nrow(dt_test)), fmt_n(uniqueN(dt_test$ResearchID)))

## =============================================================================
## SECTION 4 -- STANDARDIZATION CONSTANT + STAGE 1 (SINGLE) DECAY-CURVE FIT.
## Fit ONLY on convictions with no other simultaneously-active conviction
## (NActive == 1), TRAIN only. ONE glmer() fit on the full purified
## population at once -- no partitioning, no pooling.
## =============================================================================

purified_train_raw <- dt_train[NActive == 1]

n_purified_before_age <- nrow(purified_train_raw)
n_purified_na_age <- sum(is.na(purified_train_raw$CurrentAge))
purified_train_raw <- purified_train_raw[CurrentAge >= CONFIG$MIN_ASSESSMENT_AGE]
log_line("Applied under-18-at-assessment exclusion to purified_train_raw: %s rows kept of %s (%.4f%% dropped; of those dropped, %s had a missing/NA age rather than a confirmed under-18 age).",
         fmt_n(nrow(purified_train_raw)), fmt_n(n_purified_before_age),
         100 * (n_purified_before_age - nrow(purified_train_raw)) / n_purified_before_age,
         fmt_n(n_purified_na_age))

lb_mean <- mean(purified_train_raw$LookBack, na.rm = TRUE)
lb_sd   <- sd(purified_train_raw$LookBack, na.rm = TRUE)
log_line("Standardization constant -- LookBack: mean=%.4f sd=%.4f", lb_mean, lb_sd)
if (anyNA(c(lb_mean, lb_sd))) stop("LookBack standardization constants include NA -- investigate before fitting.")

purified_train_raw[, LookBack_z := (LookBack - lb_mean) / lb_sd]

cat("\n---- Stage 1 sufficiency check (purified population, before fitting) ----\n")
suff_tbl <- purified_train_raw[, .(Rows = .N, People = uniqueN(ResearchID)), by = c(CAT_COL)]
print(suff_tbl)
cat("\nIf any category above looks thin (low hundreds or fewer), stop here and investigate before trusting the fit below.\n\n")

single_fit_cache_path <- file.path(CONFIG$CACHE_DIR, sprintf("stage1_fit_%s_%s.rds", CONFIG$PIPELINE_VERSION, rowcap_tag))
stage1_formula <- as.formula(sprintf("Recid ~ LookBack_z * %s + (1 | ResearchID)", CAT_COL))

if (file.exists(single_fit_cache_path)) {
  log_line("Reloading the cached Stage 1 fit from %s.", single_fit_cache_path)
  stage1_fit <- readRDS(single_fit_cache_path)
} else {
  log_line("Fitting Stage 1 on the full purified population (%s rows, %s people). This is the slow step -- see README.md for expected runtime.",
           fmt_n(nrow(purified_train_raw)), fmt_n(uniqueN(purified_train_raw$ResearchID)))
  t0 <- Sys.time()
  stage1_prefit <- glmer(stage1_formula,
                          data = purified_train_raw, family = binomial, nAGQ = 0,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  stage1_fit <- update(stage1_prefit, nAGQ = 1, start = lme4::getME(stage1_prefit, c("theta", "fixef")))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  conv_msgs <- stage1_fit@optinfo$conv$lme4$messages
  conv_status <- if (is.null(conv_msgs) || length(conv_msgs) == 0) "CONVERGED cleanly" else paste(conv_msgs, collapse = " | ")
  log_line("Stage 1 fit done in %.1f minutes -- %s.", elapsed / 60, conv_status)
  saveRDS(stage1_fit, single_fit_cache_path)
  log_line("Stage 1 fit cached to %s.", single_fit_cache_path)
}

stage1_coef_tbl <- summary(stage1_fit)$coefficients
stage1_fixef <- stage1_coef_tbl[, "Estimate"]
stage1_se    <- stage1_coef_tbl[, "Std. Error"]

stage1_tau00  <- as.numeric(lme4::VarCorr(stage1_fit)$ResearchID[1, 1])
stage1_tau_sd <- sqrt(stage1_tau00)

stage1_linear_predictor <- make_linear_predictor(stage1_fixef, lb_mean, lb_sd)
set.seed(CONFIG$CORRECTION_SEED)
stage1_u_draws <- rnorm(CONFIG$N_AVERAGING_DRAWS, mean = 0, sd = stage1_tau_sd)

## =============================================================================
## SECTION 5 -- STAGE 1 RESULTS: SAMPLE SIZES, FIXED EFFECTS, RANDOM EFFECT,
## PREDICTED SINGLE-EVENT DECAY
## =============================================================================

cat("\n\n==================== STAGE 1 RESULTS: SINGLE-EVENT PREDICTION DECAY ====================\n")

purified_train_raw[, ConvictionID := paste(ResearchID, OffYear, sep = "||")]

s1_sample_tbl <- data.table(
  Individuals      = uniqueN(purified_train_raw$ResearchID),
  ConvictionEvents = uniqueN(purified_train_raw$ConvictionID),
  PersonYears      = nrow(purified_train_raw)
)
cat("\n---- 1. Single-event estimation sample ----\n")
print(s1_sample_tbl)

fixef_names <- names(stage1_fixef)
category_dummy_names <- grep(paste0("^", CAT_COL), fixef_names, value = TRUE)
category_dummy_names <- category_dummy_names[!grepl(":", category_dummy_names)]
categories_with_dummy <- gsub(paste0("^", CAT_COL), "", category_dummy_names)
reference_category <- setdiff(names(ACTIVE_WINDOWS), categories_with_dummy)
cat(sprintf("\n---- 2. Stage 1 fixed effects: Recid ~ LookBack_z * %s + (1 | ResearchID) ----\n", CAT_COL))
cat(sprintf("Reference offense category (omitted %s dummy, absorbed into the Intercept): %s\n\n",
            CAT_COL, paste(reference_category, collapse = ", ")))

s1_fixed_tbl <- data.table(
  Term = names(stage1_fixef),
  B    = as.numeric(stage1_fixef),
  SE   = as.numeric(stage1_se[names(stage1_fixef)])
)
s1_fixed_tbl[, OR         := exp(B)]
s1_fixed_tbl[, OR_CI_low  := exp(B - 1.96 * SE)]
s1_fixed_tbl[, OR_CI_high := exp(B + 1.96 * SE)]
print(copy(s1_fixed_tbl)[, `:=`(B = round(B, 4), SE = round(SE, 4),
                                 OR = round(OR, 4), OR_CI_low = round(OR_CI_low, 4), OR_CI_high = round(OR_CI_high, 4))])
cat("No p-values printed anywhere in this file by design. 95% OR CI is the standard Wald interval, exp(B +/- 1.96*SE).\n")

resid_var_logistic <- (pi^2) / 3
icc <- stage1_tau00 / (stage1_tau00 + resid_var_logistic)
s1_random_tbl <- data.table(
  RandomIntercept_Variance        = round(stage1_tau00, 4),
  RandomIntercept_SD              = round(stage1_tau_sd, 4),
  ResidualVariance_LogisticLatent = round(resid_var_logistic, 4),
  ICC                             = round(icc, 4)
)
cat("\n---- 3. Random effect (person-level intercept) ----\n")
print(s1_random_tbl)

## Only report a decay row for a category that actually has rows in the
## Stage 1 fitting population (purified_train_raw). Without this guard, a
## category with ZERO fitting rows has no fitted dummy coefficient, so
## make_linear_predictor()'s `if (cat_term %in% x_names)` check silently
## falls through to the bare intercept -- scoring that category identically
## to the true omitted reference category, with no flag anywhere that this
## is a fabricated fallback rather than a real, data-supported estimate.
categories_in_fit <- unique(purified_train_raw[[CAT_COL]])
categories_absent_from_fit <- setdiff(names(ACTIVE_WINDOWS), categories_in_fit)
if (length(categories_absent_from_fit) > 0) {
  warning(sprintf("Categor(ies) with ZERO rows in the Stage 1 fitting population (purified_train_raw, after the under-18 exclusion): %s. These have no fitted coefficient and are OMITTED from s1_predicted_decay.csv -- there is no data to estimate a decay curve for them.", paste(categories_absent_from_fit, collapse = ", ")))
}

lookback_years_report <- c(1, 5, 10, 15)
decay_grid <- CJ(CategoryTmp = setdiff(names(ACTIVE_WINDOWS), categories_absent_from_fit), LookBack = lookback_years_report)
setnames(decay_grid, "CategoryTmp", CAT_COL)
decay_grid[, LinearPredictor := mapply(stage1_linear_predictor, get(CAT_COL), LookBack, MoreArgs = list(age_z = 0, incar_z = 0))]
decay_grid[, Probability_pct := round(100 * sapply(LinearPredictor, function(lp) mean(plogis(lp + stage1_u_draws))), 2)]

s1_decay_tbl <- dcast(decay_grid, as.formula(paste(CAT_COL, "~ LookBack")), value.var = "Probability_pct")
setnames(s1_decay_tbl, as.character(lookback_years_report), paste0("Year", lookback_years_report, "_pct"))
s1_decay_tbl[, AbsoluteDeclinePP      := round(Year1_pct - Year15_pct, 2)]
s1_decay_tbl[, ProportionalDeclinePct := round(100 * (Year1_pct - Year15_pct) / Year1_pct, 1)]
cat("\n---- 4. Predicted 1-year reconviction probability by category and LookBack year (population-averaged) ----\n")
print(s1_decay_tbl)

fwrite(s1_sample_tbl, file.path(CONFIG$OUTPUT_DIR, "s1_sample.csv"))
fwrite(s1_fixed_tbl,  file.path(CONFIG$OUTPUT_DIR, "s1_fixed_effects.csv"))
fwrite(s1_random_tbl, file.path(CONFIG$OUTPUT_DIR, "s1_random_effect.csv"))
fwrite(s1_decay_tbl,  file.path(CONFIG$OUTPUT_DIR, "s1_predicted_decay.csv"))
cat("\n===========================================================================================\n\n")

## =============================================================================
## SECTION 6 -- HAZARD ACCUMULATION, REDUNDANCY CORRECTION
## =============================================================================

dt_train[, IncarRounded := round(get(INCAR_COL), 1)]
dt_test[,  IncarRounded := round(get(INCAR_COL), 1)]

build_event_level_minimal <- function(d) {
  ev <- d[LookBack == 1, .(ResearchID, OffYear, CurrentAge, Recid)]
  setorder(ev, ResearchID, OffYear)
  ev <- ev[ev[, .I[.N], by = .(ResearchID, OffYear)]$V1]
  ev
}
event_train <- build_event_level_minimal(dt_train); event_train[, EventID := .I]
event_test  <- build_event_level_minimal(dt_test);  event_test[, EventID := .I]

py_outcome_year <- function(d) d[, .(Recid = max(Recid)), by = .(ResearchID, CalendarYear)]
py_train_year <- py_outcome_year(dt_train)
py_test_year  <- py_outcome_year(dt_test)

## The concurrency correction is fit against a 3-year "did this conviction
## reoffend within 3 years" outcome, on the once-per-conviction population --
## a deliberately separate design choice from the 1-year staircase outcome
## the score is ultimately evaluated against below.
check_full_followup_observed <- function(events, py_year, follow_years) {
  lookup <- rbindlist(lapply(seq_len(follow_years), function(k) events[, .(EventID, ResearchID, CalendarYear = OffYear + k)]))
  lookup <- merge(lookup, py_year, by = c("ResearchID", "CalendarYear"), all.x = TRUE)
  lookup[, .(FullyObserved = all(!is.na(Recid))), by = EventID]
}
event_train <- merge(event_train, check_full_followup_observed(event_train, py_train_year, CONFIG$FOLLOW_YEARS), by = "EventID")
event_test  <- merge(event_test,  check_full_followup_observed(event_test,  py_test_year,  CONFIG$FOLLOW_YEARS), by = "EventID")
event_train <- event_train[FullyObserved == TRUE]
event_test  <- event_test[FullyObserved == TRUE]

build_recid_follow <- function(events, py_year, follow_years) {
  lookup <- rbindlist(lapply(seq_len(follow_years), function(k) events[, .(EventID, ResearchID, CalendarYear = OffYear + k)]))
  lookup <- merge(lookup, py_year, by = c("ResearchID", "CalendarYear"), all.x = TRUE)
  if (sum(is.na(lookup$Recid)) > 0) stop("Missing lookup(s) building the follow-up outcome -- should be impossible.")
  lookup[, .(Recid_Follow = max(Recid)), by = EventID]
}
event_train <- merge(event_train, build_recid_follow(event_train, py_train_year, CONFIG$FOLLOW_YEARS), by = "EventID")
event_test  <- merge(event_test,  build_recid_follow(event_test,  py_test_year,  CONFIG$FOLLOW_YEARS), by = "EventID")

## Every distinct (category, elapsed-time) cell gets one hazard value,
## averaged over the random-intercept distribution (population-averaged,
## not a single "typical person" conditional probability) -- then every
## currently-active conviction's hazard is summed per person-year (the
## union/noisy-OR rule: summing on the hazard scale equals the probability
## at least one of several independent risks materializes). Grouped by
## (category, LookBack) ONLY -- mathematically equivalent to also grouping
## by age/incarceration time, since Stage 1's fitted formula has no age or
## incarceration term (confirmed directly from the fitted coefficients'
## own term names above), and far cheaper at real-data scale.
score_side <- function(d) {
  distinct_cells <- unique(d[, c(CAT_COL, "LookBack"), with = FALSE])
  distinct_cells[, LinearPredictor := mapply(stage1_linear_predictor, get(CAT_COL), LookBack, MoreArgs = list(age_z = 0, incar_z = 0))]
  distinct_cells[, Probability := sapply(LinearPredictor, function(lp) mean(plogis(lp + stage1_u_draws)))]
  distinct_cells[, Hazard := -log(pmax(1 - Probability, 1e-12))]
  merge(d, distinct_cells[, c(CAT_COL, "LookBack", "Hazard", "Probability"), with = FALSE],
        by = c(CAT_COL, "LookBack"), all.x = TRUE)
}
build_person_year <- function(d) d[, .(NActive = .N, SumHazard = sum(Hazard), SumProb = sum(Probability)), by = .(ResearchID, CalendarYear)]

log_line("Scoring every active conviction (hazard per cell, summed per person-year) -- TRAIN...")
py_train_sh <- build_person_year(score_side(dt_train))
log_line("Scoring every active conviction -- TEST...")
py_test_sh  <- build_person_year(score_side(dt_test))

assess_train <- merge(event_train[, .(EventID, ResearchID, CalendarYear = OffYear + 1L, Recid_Follow)], py_train_sh, by = c("ResearchID", "CalendarYear"))
assess_train[, UnionP_Pure := 1 - exp(-SumHazard)]

## The redundancy correction: when 2+ convictions are simultaneously active,
## a plain sum of their hazards double-counts overlapping risk information.
## Fit as a continuous function of both the summed hazard and how many
## convictions are overlapping.
concurrent_train <- assess_train[NActive >= 2]
log_line("Fitting the concurrency (redundancy) correction on %s TRAIN assess-once person-years with 2+ active convictions...", fmt_n(nrow(concurrent_train)))
continuous_fit <- glm(Recid_Follow ~ log(SumHazard) * log(NActive), data = concurrent_train, family = binomial(link = "cloglog"))

apply_continuous_correction <- function(dt_scored, outcome_col) {
  d <- copy(dt_scored)
  d[, UnionP_Corrected := UnionP_Pure]
  concurrent_idx <- which(d$NActive >= 2)
  if (length(concurrent_idx) > 0) d[concurrent_idx, UnionP_Corrected := predict(continuous_fit, newdata = d[concurrent_idx], type = "response")]
  d[, UnionP_Final := fifelse(NActive == 1, UnionP_Pure, UnionP_Corrected)]
  list(auc = compute_auc(d$UnionP_Final, d[[outcome_col]]), scored = d)
}

## =============================================================================
## SECTION 7 -- UNDER-18 EXCLUSION ON FINAL SCORING TABLES, CALIBRATION
## (PLATT SCALING), AND HAD'S OWN PERFORMANCE REPORT
## =============================================================================

attach_currentage <- function(py_tbl, d, label) {
  age_lookup <- unique(d[, .(ResearchID, CalendarYear, CurrentAge)])
  chk <- age_lookup[, .(MinAge = min(CurrentAge), MaxAge = max(CurrentAge)), by = .(ResearchID, CalendarYear)]
  n_inconsistent <- sum(abs(chk$MaxAge - chk$MinAge) > 1, na.rm = TRUE)
  log_line("CurrentAge self-consistency check (%s): %s of %s person-years disagree by more than 1 year (should be 0 or near-0).",
           label, fmt_n(n_inconsistent), fmt_n(nrow(chk)))
  age_lookup_unique <- unique(age_lookup, by = c("ResearchID", "CalendarYear"))
  merge(py_tbl, age_lookup_unique, by = c("ResearchID", "CalendarYear"), all.x = TRUE)
}

py_train_scored <- attach_currentage(merge(py_train_sh, py_train_year, by = c("ResearchID", "CalendarYear")), dt_train, "TRAIN")
staircase_test  <- attach_currentage(merge(py_test_sh,  py_test_year,  by = c("ResearchID", "CalendarYear")), dt_test,  "TEST")

n_py_train_before_age <- nrow(py_train_scored)
n_test_before_age     <- nrow(staircase_test)
n_py_train_na_age     <- sum(is.na(py_train_scored$CurrentAge))
n_test_na_age         <- sum(is.na(staircase_test$CurrentAge))
py_train_scored <- py_train_scored[CurrentAge >= CONFIG$MIN_ASSESSMENT_AGE]
staircase_test  <- staircase_test[CurrentAge >= CONFIG$MIN_ASSESSMENT_AGE]
log_line("Applied under-18-at-assessment exclusion to final scoring tables -- TRAIN: %s rows kept of %s (%s of the dropped had missing/NA age, not a confirmed under-18 age). TEST: %s rows kept of %s (%s of the dropped had missing/NA age).",
         fmt_n(nrow(py_train_scored)), fmt_n(n_py_train_before_age), fmt_n(n_py_train_na_age),
         fmt_n(nrow(staircase_test)), fmt_n(n_test_before_age), fmt_n(n_test_na_age))

py_train_scored[, UnionP_Pure := 1 - exp(-SumHazard)]
py_train_scored <- apply_continuous_correction(py_train_scored, "Recid")$scored

staircase_test[, UnionP_Pure := 1 - exp(-SumHazard)]
decay_staircase_1yr <- apply_continuous_correction(staircase_test, "Recid")
staircase_test <- decay_staircase_1yr$scored

cat("\n\n==================== HAZARD ACCUMULATION / REDUNDANCY CORRECTION ====================\n")
h1_concur_tbl <- staircase_test[, .(N = .N, Observed_pct = round(100 * mean(Recid), 1),
                                     MeanPred_UnionFinal_pct = round(100 * mean(UnionP_Final), 1)),
                                 by = .(ActiveBin = fifelse(NActive == 1, "1", fifelse(NActive == 2, "2", fifelse(NActive == 3, "3", "4+"))))]
setorder(h1_concur_tbl, ActiveBin)
print(h1_concur_tbl)
cat("=============================================================================================\n\n")

cat("\n\n==================== CALIBRATION CORRECTION FOR HAD (PLATT SCALING) ====================\n")

py_train_scored[, LogitP_Decay := safe_logit(UnionP_Final)]
staircase_test[, LogitP_Decay := safe_logit(UnionP_Final)]
platt_fit_decay <- glm(Recid ~ LogitP_Decay, data = py_train_scored, family = binomial(link = "logit"))
staircase_test[, UnionP_Platt := predict(platt_fit_decay, newdata = staircase_test, type = "response")]

auc_decay_before <- compute_auc(staircase_test$UnionP_Final, staircase_test$Recid)
auc_decay_after  <- compute_auc(staircase_test$UnionP_Platt, staircase_test$Recid)
cat(sprintf("HAD score AUC (1-year outcome, TEST) -- BEFORE Platt = %.4f, AFTER Platt = %.4f (difference = %.6f -- should be ~0, since Platt scaling cannot change ranking).\n",
            auc_decay_before, auc_decay_after, auc_decay_after - auc_decay_before))

pooled_observed <- mean(staircase_test$Recid)
pooled_before <- mean(staircase_test$UnionP_Final)
pooled_after  <- mean(staircase_test$UnionP_Platt)
cat(sprintf("Pooled calibration: observed = %.4f. Mean predicted BEFORE = %.4f (ratio %.3f). Mean predicted AFTER = %.4f (ratio %.3f).\n",
            pooled_observed, pooled_before, pooled_before / pooled_observed, pooled_after, pooled_after / pooled_observed))

rmse_decay_before <- sqrt(mean((staircase_test$Recid - staircase_test$UnionP_Final)^2))
rmse_decay_after  <- sqrt(mean((staircase_test$Recid - staircase_test$UnionP_Platt)^2))
cat(sprintf("Pooled RMSE: BEFORE = %.4f, AFTER = %.4f (%.1f%% reduction).\n",
            rmse_decay_before, rmse_decay_after, 100 * (rmse_decay_before - rmse_decay_after) / rmse_decay_before))
cat("============================================================================================\n\n")

h1_summary_tbl <- data.table(AUC_before_Platt = round(auc_decay_before, 4), AUC_after_Platt = round(auc_decay_after, 4),
                              RMSE_before_Platt = round(rmse_decay_before, 4), RMSE_after_Platt = round(rmse_decay_after, 4),
                              Observed_pct = round(100 * pooled_observed, 2), MeanPredicted_after_Platt_pct = round(100 * pooled_after, 2))
fwrite(h1_concur_tbl, file.path(CONFIG$OUTPUT_DIR, "concurrency_table.csv"))
fwrite(h1_summary_tbl, file.path(CONFIG$OUTPUT_DIR, "had_performance_summary.csv"))

## =============================================================================
## SECTION 8 -- ASSESS-ONCE CONCURRENCY TABLE (MANUSCRIPT TABLE 3 ANALOG)
## One row per real conviction, scored once at CalendarYear = OffYear + 1,
## evaluated against the 3-year Recid_Follow outcome -- deliberately a
## DIFFERENT population/outcome than Section 7's full-staircase/1-year
## design above. No Platt calibration anywhere in this section.
## =============================================================================

assess_train[, AdditiveP := pmin(1, SumProb)]

assess_test <- merge(event_test[, .(EventID, ResearchID, CalendarYear = OffYear + 1L, Recid_Follow)], py_test_sh, by = c("ResearchID", "CalendarYear"))
assess_test[, UnionP_Pure := 1 - exp(-SumHazard)]
assess_test[, AdditiveP := pmin(1, SumProb)]

apply_assess_correction <- function(dt_scored) {
  d <- copy(dt_scored)
  d[, UnionP_Corrected := UnionP_Pure]
  concurrent_idx <- which(d$NActive >= 2)
  if (length(concurrent_idx) > 0) d[concurrent_idx, UnionP_Corrected := predict(continuous_fit, newdata = d[concurrent_idx], type = "response")]
  d[, UnionP_Final := fifelse(NActive == 1, UnionP_Pure, UnionP_Corrected)]
  d
}
assess_train <- apply_assess_correction(assess_train)
assess_test  <- apply_assess_correction(assess_test)
concurrent_test <- assess_test[NActive >= 2]

cat("\n\n==================== ASSESS-ONCE CONCURRENCY TABLE (TABLE 3 ANALOG) ====================\n")

smart_round_chr <- function(x, digits = 2) {
  r <- round(x, digits)
  needs_more <- (r == 0 & x != 0) | (abs(r) == 1 & abs(x) != 1 & abs(x) <= 1.5)
  out <- sprintf(paste0("%.", digits, "f"), r)
  if (any(needs_more, na.rm = TRUE)) out[needs_more] <- sprintf("%.4f", x[needs_more])
  out
}
build_concurrency_table <- function(d) {
  d[, ConcurrencyLevel := factor(fifelse(NActive == 1, "1", fifelse(NActive == 2, "2", fifelse(NActive == 3, "3", "4+"))), levels = c("1", "2", "3", "4+"))]
  d[, .(N = .N, Observed = mean(Recid_Follow), Predicted_Additive = mean(AdditiveP),
        Predicted_Union = mean(UnionP_Pure), Predicted_Corrected = mean(UnionP_Final)), by = ConcurrencyLevel][order(ConcurrencyLevel)]
}
print_concurrency_table <- function(tbl) {
  p <- copy(tbl)
  for (col in c("Observed", "Predicted_Additive", "Predicted_Union", "Predicted_Corrected")) p[[col]] <- smart_round_chr(p[[col]])
  p[["N"]] <- fmt_n(p[["N"]])
  print(p)
}

item_train_tbl <- build_concurrency_table(assess_train)
cat("---- Observed vs. predicted reconviction risk by concurrency level, TRAIN (the population the correction model was fit on -- a close match here is expected, not independent validation) ----\n")
print_concurrency_table(item_train_tbl)

item_test_tbl <- build_concurrency_table(assess_test)
cat("\n---- Observed vs. predicted reconviction risk by concurrency level, TEST (out-of-sample -- this is Table 3's analog) ----\n")
print_concurrency_table(item_test_tbl)
cat("=======================================================================================================\n\n")

fwrite(item_test_tbl, file.path(CONFIG$OUTPUT_DIR, "table3_concurrency_test.csv"))

## =============================================================================
## SECTION 9 -- COMPARATORS: STATIC (REGRESSION) AND BURGESS (UNIT-WEIGHTED)
## =============================================================================

## NOTE on OffDate: with the public CSV schema, OffDate/ConvictionYear is a
## bare year, not a full date -- so for a person with 2+ distinct
## convictions in the SAME calendar year, this setorder() can only break
## the tie using whatever relative row order the data already has, via
## data.table's stable sort. 00_prepare_public_data.R sorts its OUTPUT rows
## by each person's TRUE original offense date before writing the CSV
## specifically so that order is correct here -- don't reshuffle the CSV's
## row order before running this pipeline, or same-year ties may resolve
## differently than the manuscript's own (full-date) pipeline would.
build_event_level_full <- function(d) {
  cols <- c("ResearchID", CAT_COL, "OffDate", "OffYear", "CurrentAge", "Recid")
  ev <- d[LookBack == 1, ..cols]
  setorder(ev, ResearchID, OffDate)
  ev[, RealPriorCount := seq_len(.N), by = ResearchID]
  ev[, TrueOffYear := OffYear]
  ev <- ev[ev[, .I[.N], by = .(ResearchID, OffYear)]$V1]
  ev
}
event_train_full <- build_event_level_full(dt_train)
event_test_full  <- build_event_level_full(dt_test)

## Marks every conviction from a person's second real Violent-category
## conviction onward (and that conviction itself) as "violence ever" -- the
## running lifetime-violence flag Burgess uses as its second item. With the
## 8-category scheme, "Violent" excludes Weapons (which is its own category).
build_violent_ever_flag <- function(ev) {
  ev <- copy(ev)
  setorder(ev, ResearchID, OffDate)
  ev[, IsViolentConv := as.integer(get(CAT_COL) == "Violent")]
  ev[, ViolentEverFlag := as.integer(cummax(IsViolentConv) > 0L), by = ResearchID]
  ev
}
event_train_anchor <- build_violent_ever_flag(event_train_full)
event_test_anchor  <- build_violent_ever_flag(event_test_full)

build_burgess_score_from_cols <- function(prior_count_raw, violent_flag, cap, include_violence) {
  capped_prior_points <- pmin(prior_count_raw, cap)
  if (include_violence) capped_prior_points + violent_flag else capped_prior_points
}

## The "anchor" conviction for a person-year with several simultaneously
## active convictions is the most recently incurred one -- what Static and
## Burgess score off of, since neither has HAD's own mechanism for combining
## several convictions' information at once.
build_anchor_lookup <- function(d, event_anchor_tbl) {
  anchor <- d[, .SD[which.max(OffYear)], by = .(ResearchID, CalendarYear),
              .SDcols = c(CAT_COL, "OffYear", "LookBack", "CurrentAge", INCAR_COL)]
  anchor[, IsNewConvictionYear := as.integer(LookBack == 1L)]
  anchor <- merge(anchor, event_anchor_tbl[, .(ResearchID, OffYear = TrueOffYear, RealPriorCount, ViolentEverFlag)],
                   by = c("ResearchID", "OffYear"), all.x = TRUE)
  if (sum(is.na(anchor$RealPriorCount)) > 0) stop("Missing anchor-conviction lookup -- every active conviction should trace back to a real event-level row.")
  setnames(anchor, c(CAT_COL, "RealPriorCount", "ViolentEverFlag", "LookBack"),
                   c("RunningCategory", "RunningPriorCount", "RunningViolentEverFlag", "RunningLookBack"))
  setnames(anchor, INCAR_COL, "RunningIncarYears")
  anchor[, `:=`(OffYear = NULL, CurrentAge = NULL)]
  anchor
}
anchor_train <- build_anchor_lookup(dt_train, event_train_anchor)
anchor_test  <- build_anchor_lookup(dt_test, event_test_anchor)

py_train_scored <- merge(py_train_scored, anchor_train, by = c("ResearchID", "CalendarYear"))
staircase_test  <- merge(staircase_test,  anchor_test,  by = c("ResearchID", "CalendarYear"))

running_priorcount_mean <- mean(py_train_scored$RunningPriorCount)
running_priorcount_sd   <- sd(py_train_scored$RunningPriorCount)
py_train_scored[, RunningPriorCount_z := (RunningPriorCount - running_priorcount_mean) / running_priorcount_sd]
staircase_test[,  RunningPriorCount_z := (RunningPriorCount - running_priorcount_mean) / running_priorcount_sd]

## Static comparator: a plain logistic regression (no person-level random
## intercept) on offense category and standardized running prior-conviction
## count.
static_fit_1yr <- glm(Recid ~ RunningCategory + RunningPriorCount_z, data = py_train_scored, family = binomial)
py_train_scored[, StaticP_1yr := predict(static_fit_1yr, newdata = py_train_scored, type = "response")]
staircase_test[,  StaticP_1yr := predict(static_fit_1yr, newdata = staircase_test,  type = "response")]

py_train_scored[, BurgessScore := build_burgess_score_from_cols(RunningPriorCount - 1L, RunningViolentEverFlag, cap = CONFIG$BURGESS_HEADLINE_CAP, include_violence = TRUE)]
staircase_test[,  BurgessScore := build_burgess_score_from_cols(RunningPriorCount - 1L, RunningViolentEverFlag, cap = CONFIG$BURGESS_HEADLINE_CAP, include_violence = TRUE)]

burgess_prob_fit_1yr <- glm(Recid ~ BurgessScore, data = py_train_scored, family = binomial)
py_train_scored[, BurgessP_1yr := predict(burgess_prob_fit_1yr, newdata = py_train_scored, type = "response")]
staircase_test[,  BurgessP_1yr := predict(burgess_prob_fit_1yr, newdata = staircase_test,  type = "response")]

log_line("Static/Burgess comparators built -- TRAIN: %s rows. TEST: %s rows.", fmt_n(nrow(py_train_scored)), fmt_n(nrow(staircase_test)))

## =============================================================================
## SECTION 10 -- H2: HAD VS. STATIC VS. BURGESS (THE HEADLINE COMPARISON)
## =============================================================================

SCORE_COLS_DESC <- c(HAD = "UnionP_Final", HAD_Platt = "UnionP_Platt", Static = "StaticP_1yr", Burgess = "BurgessP_1yr")

cat("\n\n==================== HAD VS. STATIC VS. BURGESS (HELD-OUT, 1-YEAR OUTCOME) ====================\n")
h2_tbl <- rbindlist(lapply(names(SCORE_COLS_DESC), function(nm) {
  col <- SCORE_COLS_DESC[[nm]]
  data.table(Score = nm,
             AUC = round(compute_auc(staircase_test[[col]], staircase_test$Recid), 4),
             MeanPredicted_pct = round(100 * mean(staircase_test[[col]]), 1),
             RMSE = round(sqrt(mean((staircase_test$Recid - staircase_test[[col]])^2)), 4))
}))
h2_tbl[, Observed_pct := round(100 * mean(staircase_test$Recid), 1)]
print(h2_tbl)
fwrite(h2_tbl, file.path(CONFIG$OUTPUT_DIR, "h2_headline.csv"))
cat("=====================================================================================\n\n")

## =============================================================================
## SECTION 11 -- GENERALIZABILITY ACROSS AGE, GENDER, RACE/ETHNICITY, AND
## PRIOR INCARCERATION EXPOSURE
## =============================================================================

build_group_row <- function(sub) {
  data.table(
    N                = nrow(sub),
    Observed_pct     = round(100 * mean(sub$Recid), 1),
    AUC_HAD          = round(compute_auc(sub$UnionP_Platt, sub$Recid), 4),
    AUC_Static       = round(compute_auc(sub$StaticP_1yr,  sub$Recid), 4),
    AUC_Burgess      = round(compute_auc(sub$BurgessP_1yr, sub$Recid), 4),
    MeanPred_HAD     = round(100 * mean(sub$UnionP_Platt), 1),
    MeanPred_Static  = round(100 * mean(sub$StaticP_1yr),  1),
    MeanPred_Burgess = round(100 * mean(sub$BurgessP_1yr), 1),
    RMSE_HAD         = round(sqrt(mean((sub$Recid - sub$UnionP_Platt)^2)),  4),
    RMSE_Static      = round(sqrt(mean((sub$Recid - sub$StaticP_1yr)^2)),  4),
    RMSE_Burgess     = round(sqrt(mean((sub$Recid - sub$BurgessP_1yr)^2)), 4)
  )
}

AGE_BREAKS <- c(-Inf, 24, 34, 44, 54, Inf)
AGE_LABELS <- c("<25", "25-34", "35-44", "45-54", "55+")
staircase_test[, AgeGroup := cut(CurrentAge, breaks = AGE_BREAKS, labels = AGE_LABELS, right = TRUE)]

cat("\n\n==================== GENERALIZABILITY ACROSS AGE ====================\n")
h3_tbl <- rbindlist(lapply(AGE_LABELS, function(g) { row <- build_group_row(staircase_test[AgeGroup == g]); row[, AgeGroup := g]; row }))
setcolorder(h3_tbl, c("AgeGroup", setdiff(names(h3_tbl), "AgeGroup")))
h3_tbl[, AUCdiff_HAD_minus_Static  := round(AUC_HAD - AUC_Static, 4)]
h3_tbl[, AUCdiff_HAD_minus_Burgess := round(AUC_HAD - AUC_Burgess, 4)]
print(h3_tbl)
fwrite(h3_tbl, file.path(CONFIG$OUTPUT_DIR, "generalizability_age.csv"))
cat("=====================================================================================\n\n")

## Gender/Race/Ethnicity are not present on the main cache above, so this
## does a separate, lightweight raw-file read for ResearchID/Gender/Race/
## Ethnicity only -- cached once (reused by this section and Section 13
## below).
demo_cache_path <- file.path(CONFIG$CACHE_DIR, sprintf("demographics_%s.rds", CONFIG$PIPELINE_VERSION))
if (file.exists(demo_cache_path)) {
  log_line("Reloading the cached demographics-only extract from %s.", demo_cache_path)
  full_pop_demo <- readRDS(demo_cache_path)
} else {
  log_line("No cached demographics extract found -- doing a one-time, separate raw-file read for ResearchID/Gender/Race/Ethnicity only. This only has to happen once.")
  raw_demo <- load_raw_file(CONFIG$RAW_FILE_PATH)
  demo_cols_needed <- c("ResearchID", "Gender", "Race", "Ethnicity")
  missing_demo_cols <- setdiff(demo_cols_needed, names(raw_demo))
  if (length(missing_demo_cols) > 0) stop(sprintf("Raw file is missing expected demographic column(s): %s -- cannot build the gender/race generalizability tables without them.", paste(missing_demo_cols, collapse = ", ")))
  full_pop_demo <- raw_demo[, ..demo_cols_needed]
  full_pop_demo[, ResearchID := as.character(ResearchID)]
  full_pop_demo <- unique(full_pop_demo, by = "ResearchID")
  saveRDS(full_pop_demo, demo_cache_path)
  log_line("Demographics-only extract cached to %s (%s unique people).", demo_cache_path, fmt_n(nrow(full_pop_demo)))
  rm(raw_demo)
}

demo_lookup <- full_pop_demo[ResearchID %in% unique(staircase_test$ResearchID)]
demo_lookup <- demo_lookup[is.na(Race) | Race != "U"]
demo_lookup[, RaceEthnicityGroup := fifelse(!is.na(Ethnicity) & Ethnicity == "H", "Hispanic",
                                     fifelse(!is.na(Race) & Race == "W", "White",
                                     fifelse(!is.na(Race) & Race == "B", "Black", "Other")))]

staircase_test_demo <- merge(staircase_test, demo_lookup, by = "ResearchID")

cat("\n\n==================== GENERALIZABILITY ACROSS GENDER ====================\n")
GENDER_LEVELS <- c(Male = "M", Female = "F")
h4_gender_tbl <- rbindlist(lapply(names(GENDER_LEVELS), function(nm) { row <- build_group_row(staircase_test_demo[Gender == GENDER_LEVELS[[nm]]]); row[, Gender := nm]; row }))
setcolorder(h4_gender_tbl, c("Gender", setdiff(names(h4_gender_tbl), "Gender")))
h4_gender_tbl[, AUCdiff_HAD_minus_Static  := round(AUC_HAD - AUC_Static, 4)]
h4_gender_tbl[, AUCdiff_HAD_minus_Burgess := round(AUC_HAD - AUC_Burgess, 4)]
print(h4_gender_tbl)
fwrite(h4_gender_tbl, file.path(CONFIG$OUTPUT_DIR, "generalizability_gender.csv"))
cat("=====================================================================================\n\n")

cat("\n\n==================== GENERALIZABILITY ACROSS RACE/ETHNICITY ====================\n")
RACE_LEVELS <- c("White", "Black", "Hispanic", "Other")
h4_race_tbl <- rbindlist(lapply(RACE_LEVELS, function(g) { row <- build_group_row(staircase_test_demo[RaceEthnicityGroup == g]); row[, RaceEthnicityGroup := g]; row }))
setcolorder(h4_race_tbl, c("RaceEthnicityGroup", setdiff(names(h4_race_tbl), "RaceEthnicityGroup")))
h4_race_tbl[, AUCdiff_HAD_minus_Static  := round(AUC_HAD - AUC_Static, 4)]
h4_race_tbl[, AUCdiff_HAD_minus_Burgess := round(AUC_HAD - AUC_Burgess, 4)]
print(h4_race_tbl)
fwrite(h4_race_tbl, file.path(CONFIG$OUTPUT_DIR, "generalizability_race.csv"))
cat("=====================================================================================\n\n")
rm(staircase_test_demo)

dose_bin_3 <- function(dose) fifelse(dose <= 0, "None", fifelse(dose <= 1, "Low (<=1yr)", "High (>1yr)"))
staircase_test[, DoseBin3 := dose_bin_3(RunningIncarYears)]
DOSE_LEVELS <- c("None", "Low (<=1yr)", "High (>1yr)")

cat("\n\n==================== GENERALIZABILITY ACROSS PRIOR INCARCERATION EXPOSURE (3-BIN) ====================\n")
h5_dose_tbl <- rbindlist(lapply(DOSE_LEVELS, function(g) { row <- build_group_row(staircase_test[DoseBin3 == g]); row[, DoseBin3 := g]; row }))
setcolorder(h5_dose_tbl, c("DoseBin3", setdiff(names(h5_dose_tbl), "DoseBin3")))
h5_dose_tbl[, AUCdiff_HAD_minus_Static  := round(AUC_HAD - AUC_Static, 4)]
h5_dose_tbl[, AUCdiff_HAD_minus_Burgess := round(AUC_HAD - AUC_Burgess, 4)]
print(h5_dose_tbl)
fwrite(h5_dose_tbl, file.path(CONFIG$OUTPUT_DIR, "generalizability_incarceration_3bin.csv"))
cat("=====================================================================================\n\n")

DOSE6_BREAKS <- c(-Inf, 0, 1, 3, 5, 10, Inf)
DOSE6_LABELS <- c("None (0)", "Low (0-1]", "Mid (1-3]", "Mid-high (3-5]", "High (5-10]", "Very high (10+)")
staircase_test[, DoseBin6 := cut(RunningIncarYears, breaks = DOSE6_BREAKS, labels = DOSE6_LABELS, right = TRUE)]

cat("\n\n==================== GENERALIZABILITY ACROSS PRIOR INCARCERATION EXPOSURE (6-BIN) ====================\n")
h5b_dose_tbl <- rbindlist(lapply(DOSE6_LABELS, function(g) { row <- build_group_row(staircase_test[DoseBin6 == g]); row[, DoseBin6 := g]; row }))
setcolorder(h5b_dose_tbl, c("DoseBin6", setdiff(names(h5b_dose_tbl), "DoseBin6")))
h5b_dose_tbl[, AUCdiff_HAD_minus_Static  := round(AUC_HAD - AUC_Static, 4)]
h5b_dose_tbl[, AUCdiff_HAD_minus_Burgess := round(AUC_HAD - AUC_Burgess, 4)]
print(h5b_dose_tbl)
fwrite(h5b_dose_tbl, file.path(CONFIG$OUTPUT_DIR, "generalizability_incarceration_6bin.csv"))
cat("=====================================================================================\n\n")

## =============================================================================
## SECTION 12 -- TABLE 1: FULL DESCRIPTIVE STATISTICS
## =============================================================================

cat("\n\n==================== TABLE 1: FULL DESCRIPTIVE STATISTICS ====================\n")

t1_dev_individuals   <- uniqueN(py_train_scored$ResearchID)
t1_dev_person_years  <- nrow(py_train_scored)
t1_held_individuals  <- uniqueN(staircase_test$ResearchID)
t1_held_person_years <- nrow(staircase_test)
t1_sample_tbl <- data.table(
  Panel = c("Development (TRAIN)", "Held-out (TEST)", "Total"),
  UniqueIndividuals = c(t1_dev_individuals, t1_held_individuals, t1_dev_individuals + t1_held_individuals),
  PersonYearObservations = c(t1_dev_person_years, t1_held_person_years, t1_dev_person_years + t1_held_person_years)
)
cat("\n---- Sample ----\n"); print(t1_sample_tbl)

t1_final_ids <- unique(c(py_train_scored$ResearchID, staircase_test$ResearchID))
t1_demo_lookup <- full_pop_demo[ResearchID %in% t1_final_ids]
t1_demo_lookup <- t1_demo_lookup[is.na(Race) | Race != "U"]
t1_demo_lookup[, RaceEthnicityGroup := fifelse(!is.na(Ethnicity) & Ethnicity == "H", "Hispanic",
                                        fifelse(!is.na(Race) & Race == "W", "White",
                                        fifelse(!is.na(Race) & Race == "B", "Black", "Other")))]
t1_race_tbl <- t1_demo_lookup[, .N, by = RaceEthnicityGroup]
t1_race_tbl[, Pct := round(100 * N / nrow(t1_demo_lookup), 2)]
t1_race_tbl[, RaceEthnicityGroup := factor(RaceEthnicityGroup, levels = c("White", "Hispanic", "Black", "Other"))]
setorder(t1_race_tbl, RaceEthnicityGroup)
cat("\n---- Demographics (race/ethnicity) ----\n"); print(t1_race_tbl)

t1_gender_lookup <- full_pop_demo[ResearchID %in% t1_final_ids, .(ResearchID, Gender)]
t1_n_gender_mf <- t1_gender_lookup[Gender %in% c("M", "F"), .N]
t1_n_male <- t1_gender_lookup[Gender == "M", .N]
cat(sprintf("Male: %.2f%% (N = %s of %s people with Gender recorded as M or F).\n", 100 * t1_n_male / t1_n_gender_mf, fmt_n(t1_n_male), fmt_n(t1_n_gender_mf)))

t1_age_vec <- c(py_train_scored$CurrentAge, staircase_test$CurrentAge)
cat("\n---- Age ----\n"); print(data.table(N = length(t1_age_vec), Mean = round(mean(t1_age_vec), 2), SD = round(sd(t1_age_vec), 2)))

t1_event_all <- rbind(event_train_full, event_test_full)
cat("\n---- Prior convictions (event-level) ----\n")
print(data.table(ConvictionEvents = nrow(t1_event_all), Mean = round(mean(t1_event_all$RealPriorCount), 2), SD = round(sd(t1_event_all$RealPriorCount), 2)))

t1_offense_tbl <- t1_event_all[, .N, by = c(CAT_COL)]
t1_offense_tbl[, Pct := round(100 * N / nrow(t1_event_all), 2)]
setorder(t1_offense_tbl, -N)
cat("\n---- Offense category (event-level, descending by size) ----\n"); print(t1_offense_tbl)

t1_nactive_vec <- c(py_train_scored$NActive, staircase_test$NActive)
cat("\n---- Concurrency ----\n")
print(data.table(N = length(t1_nactive_vec), Mean = round(mean(t1_nactive_vec), 2), SD = round(sd(t1_nactive_vec), 2)))
t1_concur_band_tbl <- data.table(Band = factor(fifelse(t1_nactive_vec == 1, "1", fifelse(t1_nactive_vec == 2, "2", fifelse(t1_nactive_vec == 3, "3", "4+"))), levels = c("1", "2", "3", "4+")))[, .N, by = Band]
t1_concur_band_tbl[, Pct := round(100 * N / length(t1_nactive_vec), 2)]
setorder(t1_concur_band_tbl, Band)
print(t1_concur_band_tbl)

t1_outcome_tbl <- data.table(
  Panel = c("Overall", "Development (TRAIN)", "Held-out (TEST)"),
  N = c(t1_dev_person_years + t1_held_person_years, t1_dev_person_years, t1_held_person_years),
  ObservedOneYearRecidRate_Pct = round(100 * c(mean(c(py_train_scored$Recid, staircase_test$Recid)), mean(py_train_scored$Recid), mean(staircase_test$Recid)), 2)
)
cat("\n---- Outcome ----\n"); print(t1_outcome_tbl)

fwrite(t1_sample_tbl,      file.path(CONFIG$OUTPUT_DIR, "table1_sample.csv"))
fwrite(t1_race_tbl,        file.path(CONFIG$OUTPUT_DIR, "table1_race.csv"))
fwrite(t1_offense_tbl,     file.path(CONFIG$OUTPUT_DIR, "table1_offense_category.csv"))
fwrite(t1_concur_band_tbl, file.path(CONFIG$OUTPUT_DIR, "table1_concurrency.csv"))
fwrite(t1_outcome_tbl,     file.path(CONFIG$OUTPUT_DIR, "table1_outcome.csv"))
cat("\n===========================================================================================\n\n")

## =============================================================================
## SECTION 13 -- APPENDIX TABLE A2 ANALOG: THREE-YEAR PREDICTIVE PERFORMANCE
## Restricts the full-staircase scoring tables to person-years with a fully
## observed 3-year follow-up and re-evaluates HAD/Static/Burgess against
## that longer outcome window, with dedicated 3-year Static/Burgess fits
## and a dedicated 3-year Platt rescaling for HAD.
## =============================================================================

cat("\n\n==================== APPENDIX TABLE A2 ANALOG: 3-YEAR PREDICTIVE PERFORMANCE ====================\n")

py_train_scored[, RowID := .I]
staircase_test[, RowID := .I]

check_full_followup_observed_py <- function(py_tbl, py_year_lookup, follow_years) {
  lookup <- rbindlist(lapply(seq_len(follow_years), function(k) py_tbl[, .(RowID, ResearchID, CalendarYear = CalendarYear + k)]))
  lookup <- merge(lookup, py_year_lookup, by = c("ResearchID", "CalendarYear"), all.x = TRUE)
  lookup[, .(FullyObserved3 = all(!is.na(Recid))), by = RowID]
}
build_recid_follow_py <- function(py_tbl, py_year_lookup, follow_years) {
  lookup <- rbindlist(lapply(seq_len(follow_years), function(k) py_tbl[, .(RowID, ResearchID, CalendarYear = CalendarYear + k)]))
  lookup <- merge(lookup, py_year_lookup, by = c("ResearchID", "CalendarYear"), all.x = TRUE)
  if (sum(is.na(lookup$Recid)) > 0) stop("Missing lookup building the 3-year follow-up outcome -- should be impossible among FullyObserved3 rows.")
  lookup[, .(Recid_Follow3 = max(Recid)), by = RowID]
}

py_train_3yr <- merge(py_train_scored, check_full_followup_observed_py(py_train_scored, py_train_year, CONFIG$FOLLOW_YEARS), by = "RowID")
py_train_3yr <- py_train_3yr[FullyObserved3 == TRUE]
test_3yr     <- merge(staircase_test,  check_full_followup_observed_py(staircase_test,  py_test_year,  CONFIG$FOLLOW_YEARS), by = "RowID")
test_3yr     <- test_3yr[FullyObserved3 == TRUE]

py_train_3yr <- merge(py_train_3yr, build_recid_follow_py(py_train_3yr, py_train_year, CONFIG$FOLLOW_YEARS), by = "RowID")
test_3yr     <- merge(test_3yr,     build_recid_follow_py(test_3yr,     py_test_year,  CONFIG$FOLLOW_YEARS), by = "RowID")

log_line("3-year full-follow-up-observed restriction -- TRAIN: %s rows kept of %s. TEST: %s rows kept of %s.",
         fmt_n(nrow(py_train_3yr)), fmt_n(nrow(py_train_scored)), fmt_n(nrow(test_3yr)), fmt_n(nrow(staircase_test)))

static_fit_3yr <- glm(Recid_Follow3 ~ RunningCategory + RunningPriorCount_z, data = py_train_3yr, family = binomial)
test_3yr[, StaticP_3yr := predict(static_fit_3yr, newdata = test_3yr, type = "response")]

burgess_prob_fit_3yr <- glm(Recid_Follow3 ~ BurgessScore, data = py_train_3yr, family = binomial)
test_3yr[, BurgessP_3yr := predict(burgess_prob_fit_3yr, newdata = test_3yr, type = "response")]

platt_fit_decay_3yr <- glm(Recid_Follow3 ~ LogitP_Decay, data = py_train_3yr, family = binomial(link = "logit"))
test_3yr[, UnionP_Platt_3yr := predict(platt_fit_decay_3yr, newdata = test_3yr, type = "response")]

auc_had_3yr_raw   <- compute_auc(test_3yr$UnionP_Final, test_3yr$Recid_Follow3)
auc_had_3yr_platt <- compute_auc(test_3yr$UnionP_Platt_3yr, test_3yr$Recid_Follow3)
log_line("Confirming the 3-year Platt rescaling cannot change ranking -- raw HAD AUC = %.4f, rescaled HAD AUC = %.4f (difference = %.6f).",
         auc_had_3yr_raw, auc_had_3yr_platt, auc_had_3yr_platt - auc_had_3yr_raw)

SCORE_COLS_3YR <- c(HAD = "UnionP_Platt_3yr", Static = "StaticP_3yr", Burgess = "BurgessP_3yr")
observed_3yr   <- mean(test_3yr$Recid_Follow3)

a2_tbl <- rbindlist(lapply(names(SCORE_COLS_3YR), function(nm) {
  col <- SCORE_COLS_3YR[[nm]]
  mean_pred <- mean(test_3yr[[col]])
  data.table(Score = nm, AUC = round(compute_auc(test_3yr[[col]], test_3yr$Recid_Follow3), 4),
             RMSE = round(sqrt(mean((test_3yr$Recid_Follow3 - test_3yr[[col]])^2)), 4),
             MeanPredicted = round(mean_pred, 4), Observed = round(observed_3yr, 4), CalGap = round(mean_pred - observed_3yr, 4))
}))
cat("\nAppendix Table A2 analog -- 3-year predictive performance, held-out TEST sample:\n")
print(a2_tbl)
cat(sprintf("\nHeld-out person-year N (3-year outcome, full-follow-up-observed TEST sample): %s\n", fmt_n(nrow(test_3yr))))
fwrite(a2_tbl, file.path(CONFIG$OUTPUT_DIR, "appendix_a2_threeyear.csv"))
cat("================================================================================================================\n\n")

py_train_scored[, RowID := NULL]
staircase_test[, RowID := NULL]

log_line("Pipeline complete. All output tables written to %s.", CONFIG$OUTPUT_DIR)
log_line("No individual-level rows and no p-values were written or printed anywhere in this file.")
