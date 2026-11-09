# Aerial Surveillance and the Geography of Selective Violence — replication materials

## Status and scope

**14 saved-panel model checks passed in the author's RStudio session on 9 November 2026.** These comprise H1 (full and anticipation exclusion), six original restricted ACLED H2 appendix estimates, two preferred joint ACLED–LWJ H2 estimates, and four H3 estimates (local comparison plus three global rings). The corrected H3 uses the 241-cell **FINAL30 geography CSV overlaid on the original Afghanistan core monthly panel**; the embedded H3 geography in either panel is not the publication definition.

**This is not a certified from-raw reproduction.** The full final 30-airfield panel construction step, remaining appendix A2–A4 checks, and all manuscript figure rebuilds have **not** been end-to-end validated. Please do not describe the package as a complete from-source JPR replication archive.

## Local quick start (Windows, RStudio)

1. Extract or clone this repository and **open its root folder** in RStudio (use `setwd('C:/Users/<YOUR USERNAME>/Downloads/aerial-surveillance-selective-violence')` if needed).
2. Install `data.table` and `fixest` if absent: `install.packages(c('data.table', 'fixest'))`.
3. `source('PREPARE_LOCAL_DATA.R')` copies required inputs and saved publication panels from `Downloads/JPR_COMPLETE_RESEARCH_ARCHIVE` into Git-ignored local folders. **This uses saved panels and does not rebuild them from raw data.**
4. `source('RUN_VERIFIED_14.R')` re-estimates and compares 14 publication statistics. It writes `outputs/publication_model_replication/ALL_14_PUBLISHED_MODEL_CHECKS.csv` and stops if any do not match.
5. The files `RUN_RAW_UPSTREAM.R`, `RUN_APPENDIX.R` and `code/raw_upstream/` preserve earlier research steps. Their complete integration with the FINAL30 data construction is pending verification; do not run them expecting a single finished paper rebuild.

## Data provenance and sharing

Original local inputs include UCDP-GED v25.1, ACLED, PRIO-GRID static and yearly data, Wikidata airbases, and Afghanistan territorial-control data. The package also depends on two historical LWJ CSV inputs. Raw third-party files, processed event data, and saved panels are excluded from the GitHub upload by `.gitignore` pending a licence review. `reference/` contains only compact geographic and canonical coefficient-reference CSV files; review their provenance before making the repository public. The two private installed data folders are not evidence of reproducibility without a correct from-raw pipeline.

The code files are preserved for provenance; direct Windows paths in archived scripts have been anonymised to `USERNAME`, **not yet made portable**. Scripts must be configured appropriately for other systems.

## Reproduced model targets

H1 full β −2.770786 (SE 1.089019, N 54,312), anticipation exclusion β −2.71921 (SE 1.12358, N 49,956). Joint H2 full β −4.997379 (SE 2.441607, N 19,920), exclusion β −6.309611 (SE 2.512772, N 16,536). Local H3 β −0.886986 (SE 0.158919, N 10,208). Global rings: −0.743051 (SE 0.765509), 0.171800 (SE 0.809994), 0.336591 (SE 0.876770), each N 21,952. The remaining six checks appear in `RUN_VERIFIED_14.R`.

The results were shown by the author's local console, not independently executed in this container; see `RUN_VERIFIED_14.R` for exact programmatic comparison. **Do not replace computational verification with these target numbers.**
