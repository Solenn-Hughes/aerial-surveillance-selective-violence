# ============================================================================
# R28_FIXED.R  — Shock-based selective vs indiscriminate rebel VAC split
#
# Purpose: Re-estimate R27's Test A using the SHOCK variable rather than
#   post_eli × eli_feas. This is the CAUSAL version of the mechanism test.
#
#   R27 runs the descriptive DiD on all 8 theatres. That's correlational
#   because post_eli is theatre-level and endogenous to US deployment
#   decisions. R28 restricts to the 3 shock theatres (Afghanistan, Syria,
#   Iraq) and uses any_isr_removed × eli_feas as the treatment. The shock
#   timing is exogenous (Doha Agreement + S-400 dispute) — so the
#   selective-vs-indiscriminate pattern that emerges here is the version
#   defensible for top-journal submission.
#
# Key design:
#   — Sample: 3 shock theatres (Afghanistan, Syria, Iraq)
#   — Treatment: any_isr_removed × eli_feas (shock variable × feasibility)
#   — Outcomes: log(1 + rebel_vac_selective), log(1 + rebel_vac_indis.)
#   — FE: pg_id + theatre_id × month_int
#   — Inference: Conley 200km spatial HAC
#
#   Prediction: β(shock × eli_feas) should be NEGATIVE for selective
#   (ISR removal contracts selective violence — reverse denunciation
#   incentive drops). β(shock × eli_feas) for indiscriminate should be
#   null (indiscriminate wasn't driving the pattern in the first place).
#
# Notes on R27 bug fix:
#   The n_nonzero column in R27's output was always 0 due to a string
#   mismatch (`gsub` produced `rebel_vac_sel` instead of
#   `rebel_vac_selective`). This is cosmetic; the regressions were fine.
#   R28 computes n_nonzero correctly.
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})
source("R00_config.R")
log_msg("R28_FIXED: shock-based selective/indiscriminate split")

OUT <- file.path("outputs", "paper1", "shock_mechanism_split")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1) LOAD PANEL + RAW ACLED ──────────────────────────────────────────
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(panel_path); setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]

ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT,  "data", "raw", "acled", "acled_full.csv")
)
acled_path <- Filter(file.exists, ACLED_PATHS)[1]
if (is.na(acled_path)) stop("ACLED raw CSV not found.")
ac <- fread(acled_path, showProgress = FALSE, encoding = "UTF-8")
setnames(ac, tolower(names(ac)))
ac[, month := as.IDate(format(as.Date(event_date), "%Y-%m-01"))]
ac[, fatalities := as.integer(fatalities)]
ac[is.na(fatalities), fatalities := 0L]

# Map ACLED events to PRIO-GRID cells
ac[, pg_row := floor((latitude  + 90)  / 0.5) + 1L]
ac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
ac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# ── 2) BUILD FATALITY-CLASS OUTCOMES INLINE (same as R27) ──────────────
rebel_vac <- ac[
  inter1 == 2L &
    grepl("Violence against civilians", event_type, ignore.case = TRUE)
]
rebel_vac[, fat_class := fcase(
  fatalities <= 2, "selective",
  fatalities <= 4, "mid",
  fatalities >= 5, "indiscriminate",
  default         = NA_character_
)]

agg_fat <- rebel_vac[, .(
  rebel_vac_selective      = sum(fat_class == "selective",       na.rm = TRUE),
  rebel_vac_indiscriminate = sum(fat_class == "indiscriminate",  na.rm = TRUE),
  rebel_vac_mid            = sum(fat_class == "mid",             na.rm = TRUE)
), by = .(pg_id, month)]

for (cc in names(agg_fat)[-(1:2)]) if (cc %in% names(dt)) dt[, (cc) := NULL]
dt <- merge(dt, agg_fat, by = c("pg_id", "month"), all.x = TRUE)
for (cc in c("rebel_vac_selective", "rebel_vac_indiscriminate", "rebel_vac_mid")) {
  dt[is.na(get(cc)), (cc) := 0L]
}

# Coordinates for Conley (if missing)
if (!all(c("lat", "lon") %in% names(dt))) {
  if (all(c("ycoord", "xcoord") %in% names(dt))) {
    dt[, `:=`(lat = ycoord, lon = xcoord)]
  } else stop("Cannot find lat/lon columns.")
}

# Log-transforms
dt[, log1p_rv_sel   := log1p(rebel_vac_selective)]
dt[, log1p_rv_indis := log1p(rebel_vac_indiscriminate)]

# ── 3) BUILD any_isr_removed (robust — construct here if missing) ──────
if (!"any_isr_removed" %in% names(dt)) {
  dt[, any_isr_removed := 0L]
  dt[theatre_id == "Afghanistan"         & month >= as.IDate("2021-07-01"),
     any_isr_removed := 1L]
  dt[theatre_id == "Syria"               & month >= as.IDate("2019-10-01"),
     any_isr_removed := 1L]
  dt[theatre_id == "Iraq"                & month >= as.IDate("2019-10-01"),
     any_isr_removed := 1L]
}

# ── 4) RESTRICT TO SHOCK THEATRES (the causal subset) ──────────────────
SHOCK_THEATRES <- c("Afghanistan", "Syria", "Iraq")
dt_shock <- dt[theatre_id %in% SHOCK_THEATRES]
log_msg("Shock-theatre sample: ", nrow(dt_shock), " obs, ",
        uniqueN(dt_shock$pg_id), " cells")

CTRL <- c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(CTRL, names(dt_shock)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── 5) RUN SHOCK-BASED SPEC ────────────────────────────────────────────
# Full specification includes both interactions:
#   β_1 = post_eli × eli_feas         (descriptive DiD, for reference)
#   β_2 = any_isr_removed × eli_feas  (the CAUSAL estimate — shock reversal)
run_shock <- function(yvar, label, col_for_nz) {
  f <- as.formula(paste0(
    yvar, " ~ post_eli:eli_feas + any_isr_removed:eli_feas + ",
    ctrl_str, " | ", FE
  ))
  m <- feols(f, data = dt_shock, notes = FALSE)
  cv <- vcov_conley(m, lat = "lat", lon = "lon", cutoff = 200)
  ct <- as.data.table(coeftable(m, vcov = cv), keep.rownames = "term")
  
  causal_row <- ct[grepl("any_isr_removed", term)]
  panel_row  <- ct[grepl("^post_eli", term) & !grepl("any_isr", term)]
  
  data.table(
    outcome           = label,
    beta_shock_causal = round(causal_row$Estimate[1],   5),
    se_shock_causal   = round(causal_row$`Std. Error`[1], 5),
    p_shock_causal    = round(causal_row$`Pr(>|t|)`[1],   4),
    beta_panel_desc   = round(panel_row$Estimate[1],    5),
    p_panel_desc      = round(panel_row$`Pr(>|t|)`[1],    4),
    n_nonzero         = sum(dt_shock[[col_for_nz]] > 0, na.rm = TRUE),
    n_obs             = nrow(dt_shock)
  )
}

cat("\n── R28: SHOCK-BASED selective vs indiscriminate (3-theatre subset) ──\n")
cat("  Causal coef = β(any_isr_removed × eli_feas) — ISR REMOVAL effect.\n")
cat("  Panel coef  = β(post_eli × eli_feas) — descriptive (for reference).\n\n")

results <- rbindlist(list(
  run_shock("log1p_rv_sel",   "Rebel VAC — SELECTIVE (≤2 fatal.)",      "rebel_vac_selective"),
  run_shock("log1p_rv_indis", "Rebel VAC — INDISCRIMINATE (≥5 fatal.)", "rebel_vac_indiscriminate")
))
print(results)
fwrite(results, file.path(OUT, "R28_shock_selective_indiscriminate.csv"))

# ── 6) INTERPRETATION DIAGNOSTIC ───────────────────────────────────────
b_sel_shock   <- results[grepl("SELECTIVE",      outcome), beta_shock_causal]
b_indis_shock <- results[grepl("INDISCRIMINATE", outcome), beta_shock_causal]
p_sel_shock   <- results[grepl("SELECTIVE",      outcome), p_shock_causal]
p_indis_shock <- results[grepl("INDISCRIMINATE", outcome), p_shock_causal]

cat("\n── CAUSAL DIAGNOSTIC ──\n")
cat(sprintf("β(shock × eli_feas) on SELECTIVE     = %.5f  p = %.4f\n",
            b_sel_shock,   p_sel_shock))
cat(sprintf("β(shock × eli_feas) on INDISCRIMINATE = %.5f  p = %.4f\n",
            b_indis_shock, p_indis_shock))

if (!is.na(b_sel_shock) && b_sel_shock < 0 && p_sel_shock < 0.10) {
  cat("\n✓ Selective rebel VAC CONTRACTS when ISR is exogenously removed.\n")
  cat("  This is the causal version of the R27 finding.\n")
  cat("  Reverse-denunciation mechanism supported under clean identification.\n")
} else {
  cat("\n⚠ Shock-based selective effect not significant — check subset size and\n")
  cat("  consider whether theatre×month FE absorb too much within-theatre\n")
  cat("  variation. Fallback: cell + month FE within each theatre.\n")
}

log_msg("R28_FIXED complete. Outputs in: ", OUT)