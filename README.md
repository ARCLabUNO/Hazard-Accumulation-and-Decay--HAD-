# HAD: a decay-weighted criminal-history risk score — replication code

This repository contains the replication code for the decay-weighted criminal-history risk score ("HAD") described in:

> Hamilton, Z., et al. *Half-Life of a Conviction: A Decay-Weighted Approach to Criminal History Risk Assessment.* Criminology. [Full citation / DOI to be added once available.]

It reproduces every table reported in the manuscript's Results section (the single-event decay curve, the hazard-accumulation/redundancy correction, HAD's predictive performance against two comparators, its generalizability across age/gender/race/incarceration exposure, the sample's full descriptive statistics, and the three-year supplementary performance comparison), using the manuscript's 8-category offense classification.

**This replaces an earlier version of this repository built around a 5-category offense classification.** The manuscript's published analyses now use the 8-category scheme below throughout; the 5-category code is retired and should not be used to reproduce the published results.

## What's in this repository

| File | Purpose |
|---|---|
| `00_prepare_public_data.R` | OPTIONAL. Only needed if you want to regenerate `data/conviction_data.csv` from a different or updated raw extract (e.g. SPSS) — converts, strips, and renames a raw file into this repo's public CSV schema. The only script here that touches SPSS or a raw birthdate. Not required to just run the pipeline below, once you've downloaded `data/conviction_data.csv` — see "Data" below. |
| `01_had_decay_pipeline.R` | The main script. Reads `data/conviction_data.csv`, fits the model, and produces every output table. |
| `check_data_schema.R` | Run this before the main script — confirms the required columns are present and sane before you sink time into the long model fit. |
| `requirements.R` | Installs the R packages the scripts above need. |
| `data/conviction_data.csv` | The processed dataset this repository uses — hosted externally due to file size; see "Data" below for the download link. |
| `output/` | The pipeline's result tables (CSV), included in this repository for reference — see "What it produces" below. |

## Requirements

- R (tested on R 4.x)
- Packages: `data.table`, `lme4` (required by `01_had_decay_pipeline.R`/`check_data_schema.R`); `haven` (only needed if you use `00_prepare_public_data.R` to convert from SPSS — run `requirements.R` to install all of these automatically)

**A note on SPSS:** nothing in this repository requires an SPSS license or installation. `haven`, the R package `00_prepare_public_data.R` uses to read a `.sav` file, is free and open-source — it reads the SPSS file format directly, without SPSS itself. `01_had_decay_pipeline.R` and `check_data_schema.R` never read `.sav` at all; they only ever read the plain CSV `00_prepare_public_data.R` produces (or one you build yourself in that same schema), so anyone can run them without any statistical software beyond R.

## Data

`data/conviction_data.csv` is the processed dataset this pipeline runs against: every column not needed by the pipeline has been dropped, offense category and the other working columns have been renamed to the plain-English public schema below, and no raw birthdate appears anywhere in the file (age is precomputed as age AT THE TIME OF EACH ASSESSMENT — not age at conviction, which the manuscript's analysis never uses on its own). It is fully de-identified — no names, birthdates, or case numbers appear in any column — and there is no restriction on its public provision for replication purposes. If you use this data, please cite the manuscript above.

At real-data scale this file is several gigabytes (~70.5 million rows), too large to host directly in this git repository, so it's hosted externally:

**[https://uofnebraska-my.sharepoint.com/:x:/r/personal/49639061_nebraska_edu/Documents/Washington%20DOC/Decay/Manuscripts/Implementing%20Decay%20in%20RNA/Scripts/decay-github-chartype2/data/conviction_data.csv?d=w8ef71d82f0af4c67bca936bb7c0d5e3b&csf=1&web=1&e=p3fc3d]**

After downloading, place the file at:

```
data/conviction_data.csv
```

### Regenerating it from a different or updated raw extract (optional)

You don't need to do this to run the pipeline — just download `data/conviction_data.csv` from the link above. This is only for producing an updated version (e.g. after a data refresh), and is the only part of this repository that touches SPSS or a raw birthdate.

If your own extract lives in SPSS (or has different column names, extra columns you don't want to carry forward, or a raw birthdate instead of an age), open `00_prepare_public_data.R`, set `CONFIG$RAW_SAV_PATH` to your real file's location, and run it once. It reads your raw file, keeps only the columns below, converts a birthdate to an age at conviction and then adds the years elapsed to get age at THIS row's assessment (never writing the birthdate, or a standalone age-at-conviction value, anywhere), converts a full offense date down to a bare conviction year, renames everything to the public schema, and overwrites `data/conviction_data.csv`. If you re-upload a regenerated file to the external host above, consider whether downstream users need to know the version changed (see the `CONFIG$PIPELINE_VERSION` note under "Running it" below).

**Required RAW columns** (what your own extract's headers need to be called for this conversion step — see the table below for what they become):

| Raw column | Meaning |
|---|---|
| `ResearchID` | A de-identified person identifier. |
| `OffDate` | The conviction's date (or a year). |
| `CharType2` | The offense category. Must be one of: `Property`, `Alcohol Offense`, `Violent`, `Other`, `Drug`, `Sex Offense/FTR`, `Weapons`, `Escape`. |
| `LookBack` | How many years after this conviction the current row represents (1 = the conviction year itself). |
| `Recid` | Whether this person was reconvicted in this calendar year (0/1). |
| `PriorIncarYearsCount` | Running count of this person's prior years incarcerated, as of this row. |
| `AgeAtOffense` **or** `Dob` | The person's age at the time of this conviction, or a birthdate to derive it from. |
| `Gender`, `Race`, `Ethnicity` | Used only for the generalizability tables. |

### The included file's schema

`data/conviction_data.csv` has these columns, one row per (person, conviction, year-since-that-conviction):

| Column | Meaning |
|---|---|
| `PersonID` | A de-identified person identifier. |
| `ConvictionYear` | The calendar year of the conviction (not a full date). |
| `ConvictionType` | The offense category. Must be one of: `Property`, `Alcohol Offense`, `Violent`, `Other`, `Drug`, `Sex Offense/FTR`, `Weapons`, `Escape`. |
| `YearsSinceConviction` | How many years after this conviction the current row represents (1 = the conviction year itself). |
| `ReconvictedThisYear` | Whether this person was reconvicted in this calendar year (0/1). |
| `PriorIncarcerationYears` | Running count of this person's prior years incarcerated, as of this row. |
| `AgeAtAssessment` | The person's age at the time of THIS row's assessment (age at conviction plus `YearsSinceConviction`) — not age at conviction on its own, which the pipeline and manuscript never use by itself. No birthdate appears anywhere in this file. |
| `Gender`, `Race`, `Ethnicity` | Used only for the generalizability tables. |

**Row order matters.** For a person with two or more distinct convictions in the same calendar year, `01_had_decay_pipeline.R` breaks the tie using the file's own row order (it only ever sees a bare year, not a full date). `00_prepare_public_data.R` writes rows in each person's true chronological order for exactly this reason — if you ever rebuild or re-sort this CSV by hand, preserve that same ordering (by person, then by true offense date).

`check_data_schema.R` validates column presence, category values, and value ranges, and prints a small, non-representative preview — run it yourself after downloading the data (see "Running it" below).

## Running it

Once you've downloaded `data/conviction_data.csv` from the link under "Data" above and placed it at that path, go straight to the pipeline — no conversion step needed:

```r
# From an R console, with the working directory set to this repo (or just
# use RStudio's Source button on each file in turn — both scripts
# auto-detect their own location):
source("requirements.R")
source("check_data_schema.R")
source("01_had_decay_pipeline.R")
```

or from a terminal:

```
Rscript requirements.R
Rscript check_data_schema.R
Rscript 01_had_decay_pipeline.R
```

`00_prepare_public_data.R` is only needed if you're regenerating `data/conviction_data.csv` from a different or updated raw extract — see "Data" above.

**Runtime:** loading and cleaning the data is fast. The single-event decay-curve fit (a mixed-effects logistic regression with a person-level random intercept) is the slow step — on the manuscript's real ~70-million-row file this took roughly 8 hours on an ordinary desktop. Everything else in the pipeline runs in minutes. Every expensive step is cached under `.pipeline_cache/` (not tracked by git — see `.gitignore`), so re-running the script after a first successful run reloads cached results in seconds. If you change the input data, bump `CONFIG$PIPELINE_VERSION` at the top of `01_had_decay_pipeline.R` to force a full rebuild rather than reusing a stale cache.

## What it produces

Every table the pipeline reports is written to `output/` as a plain CSV, in addition to being printed to the console as it runs. The copies already committed in `output/` are the exact tables generated from the manuscript's real data and back its reported results directly; re-running the pipeline against the same data file reproduces them.

- `s1_sample.csv`, `s1_fixed_effects.csv`, `s1_random_effect.csv`, `s1_predicted_decay.csv` — the single-event decay-curve model.
- `concurrency_table.csv`, `had_performance_summary.csv` — the hazard-accumulation/redundancy correction and HAD's calibration.
- `table3_concurrency_test.csv` — the assess-once concurrency table.
- `h2_headline.csv` — HAD vs. Static vs. Burgess, held-out, one-year outcome.
- `generalizability_age.csv`, `generalizability_gender.csv`, `generalizability_race.csv`, `generalizability_incarceration_3bin.csv`, `generalizability_incarceration_6bin.csv` — HAD's generalizability across subgroups.
- `table1_sample.csv`, `table1_race.csv`, `table1_offense_category.csv`, `table1_concurrency.csv`, `table1_outcome.csv` — the sample's full descriptive statistics.
- `appendix_a2_threeyear.csv` — the three-year supplementary predictive-performance comparison.

**A note on column names in these output tables:** the tables above are produced by the manuscript's own internal working pipeline, which these scripts fit directly (rather than reload) — a few output columns (for example, the offense-category column in `table1_offense_category.csv`) still use that pipeline's internal short name (`CharType2`) rather than the public `ConvictionType` name used in the input file. Only `data/conviction_data.csv` itself uses the renamed public schema throughout.

**Every output table is aggregate-only** (sample sizes, coefficient estimates, AUC/calibration summaries). No individual-level rows and no p-values are written or printed anywhere in this pipeline, by design — this keeps the script safe to run against sensitive criminal-history microdata and matches how the manuscript itself reports results (95% confidence intervals around estimates, not significance tests).

## Methodology notes

- **Study period.** `01_had_decay_pipeline.R` restricts analysis to convictions with a `ConvictionYear` between `CONFIG$STUDY_YEAR_MIN` and `CONFIG$STUDY_YEAR_MAX` (1980–2021 by default, set near the top of that script). This bounds the conviction year only, not the calendar year a given row's assessment falls in — a long-window category's later assessment years (e.g. `Escape`'s 15-year window) can run past 2021. Adjust these two constants if you're analyzing a different cohort or extending the data forward.
- **Active windows.** Each offense category has its own "active window" — how many years after a conviction it remains part of a person's active criminal history (`Property`: 13, `Alcohol Offense`: 12, `Violent`: 13, `Other`: 12, `Drug`: 14, `Sex Offense/FTR`: 14, `Weapons`: 13, `Escape`: 15 years), derived empirically from the point at which each category's own reconviction hazard converges to the population base rate. See the manuscript's Appendix Table A1.
- **Train/test split.** People (not rows) are randomly split 75%/25% into a development and held-out sample; the script verifies no person appears on both sides.
- **The decay curve (Stage 1).** A single mixed-effects logistic regression — `Recid ~ LookBack_z * CharType2 + (1 | ResearchID)` — fit once on every conviction with no other conviction simultaneously active ("purified" convictions), restricted to the development sample and to assessments at age 18+. (These are the manuscript's own internal variable names, as fit directly by this script — see the note under "What it produces" above.)
- **Combining multiple active convictions (HAD itself).** Each active conviction contributes its own decay-adjusted hazard; a person's active convictions' hazards are summed (the union/noisy-OR rule), then corrected for the redundancy introduced when multiple convictions are simultaneously active (fit empirically against a three-year outcome, on the once-per-conviction population), then Platt-recalibrated against the one-year outcome.
- **Comparators.** *Static*: a plain logistic regression on offense category and a running prior-conviction count, no decay. *Burgess*: a simple unit-weighted point count (prior convictions, capped, plus a lifetime-violence flag — matching the LSI-R's own prior-conviction cap), converted to a probability via its own single-predictor logistic fit.

## Contact

Zachary Hamilton
zhamilton@unomaha.edu

If you use this code or data, please cite the manuscript above.
