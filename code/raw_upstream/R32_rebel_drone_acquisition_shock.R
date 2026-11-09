# ============================================================================
# R32_rebel_drone_acquisition_shock.R
#
#   NEW CAUSAL IDENTIFICATION — DIFFERENT FROM ISR WITHDRAWAL
#
# Purpose: You have two ISR-withdrawal shocks (Bagram, Incirlik). Both
#   share the same structural logic: exogenous removal of state ISR.
#   Reviewers may argue you identify only one TYPE of variation.
#
#   This script adds a STRUCTURALLY DIFFERENT source of causal variation:
#   the documented acquisition of drones BY REBELS. When a rebel group
#   acquires its own ISR capability, the surveillance asymmetry shrinks.
#   Under the displacement (reverse-denunciation) mechanism, the
#   incentive to coerce civilians should weaken. The prediction:
#
#       Post-acquisition, in cells that rebel group controls,
#       rebel VAC (and especially SELECTIVE rebel VAC) should decline.
#
# The acquisition dates are documented by SIPRI, West Point CTC, Airwars:
#   — ISIS: mid-2016 (Iraq, Syria)        → first confirmed Jul 2016
#   — Houthis: early 2017 (Yemen)         → first confirmed Jan 2017
#   — Taliban: 2019+ (Afghanistan)        → less clean; exclude
#
# These acquisitions were driven by capability diffusion, not by cell-level
# violence patterns. They provide an INDEPENDENT causal test of the
# asymmetric-surveillance prediction.
#
# Design:
#   Sample   : ISIS-active cells in Iraq/Syria (2014+) + Houthi cells in Yemen (2015+)
#   Treatment: rebel_has_drones × eli_feas
#   Outcome  : log(1 + rebel_vac) and log(1 + rebel_vac_selective)
#   FE       : pg_id + theatre_id × month_int
#   Inference: Conley 200km
#
# Prediction:
#   β(rebel_has_drones × eli_feas) < 0 for selective rebel VAC.
#   Under displacement: when asymmetry shrinks, coercion pressure drops.
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})
source("R00_config.R")
log_msg("R32: rebel drone acquisition as causal shock")

OUT <- file.path("outputs", "paper1", "rebel_acquisition_shock")
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
ac[, pg_row := floor((latitude  + 90)  / 0.5) + 1L]
ac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
ac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# ── 2) IDENTIFY REBEL-VAC EVENTS BY ACTOR GROUP ────────────────────────
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

# Tag events by which rebel group
# ISIS variants
rebel_vac[, is_isis := grepl("Islamic State|ISIS|ISIL|Daesh|Wilayat", actor1,
                              ignore.case = TRUE)]
# Houthi variants
rebel_vac[, is_houthi := grepl("Houthi|Ansar Allah|Ansarallah", actor1,
                                 ignore.case = TRUE)]

# Acquisition dates (pre-determined from external sources)
ISIS_ACQ   <- as.IDate("2016-07-01")  # West Point CTC, Airwars
HOUTHI_ACQ <- as.IDate("2017-01-01")  # SIPRI, Conflict Armament Research

rebel_vac[, post_rebel_drones := FALSE]
rebel_vac[is_isis   == TRUE & month >= ISIS_ACQ,   post_rebel_drones := TRUE]
rebel_vac[is_houthi == TRUE & month >= HOUTHI_ACQ, post_rebel_drones := TRUE]

cat("\nRebel VAC event counts:\n")
cat("  ISIS (Iraq/Syria):  ", sum(rebel_vac$is_isis),
    "  — post-acq ", sum(rebel_vac$is_isis & rebel_vac$post_rebel_drones), "\n")
cat("  Houthi (Yemen):     ", sum(rebel_vac$is_houthi),
    "  — post-acq ", sum(rebel_vac$is_houthi & rebel_vac$post_rebel_drones), "\n")

# ── 3) AGGREGATE TO PG × MONTH WITH acquisition indicator ───────────────
# Count events where the rebel perpetrator had drones at time of event.
agg_rebel <- rebel_vac[, .(
  rebel_vac_predrone       = sum(!post_rebel_drones,                          na.rm = TRUE),
  rebel_vac_postdrone      = sum( post_rebel_drones,                          na.rm = TRUE),
  rebel_vac_sel_predrone   = sum(!post_rebel_drones & fat_class == "selective",      na.rm = TRUE),
  rebel_vac_sel_postdrone  = sum( post_rebel_drones & fat_class == "selective",      na.rm = TRUE),
  rebel_vac_ind_predrone   = sum(!post_rebel_drones & fat_class == "indiscriminate", na.rm = TRUE),
  rebel_vac_ind_postdrone  = sum( post_rebel_drones & fat_class == "indiscriminate", na.rm = TRUE)
), by = .(pg_id, month)]

for (cc in names(agg_rebel)[-(1:2)]) if (cc %in% names(dt)) dt[, (cc) := NULL]
dt <- merge(dt, agg_rebel, by = c("pg_id", "month"), all.x = TRUE)
for (cc in names(agg_rebel)[-(1:2)]) dt[is.na(get(cc)), (cc) := 0L]

# Coordinates
if (!all(c("lat", "lon") %in% names(dt))) {
  dt[, `:=`(lat = ycoord, lon = xcoord)]
}

# ── 4) BUILD CELL-LEVEL TREATMENT INDICATOR ────────────────────────────
# A cell is in the "rebel-drone regime" if at time t, the dominant rebel
# group active in that cell has acquired drones.
#
# Simplification: we use THEATRE-level timing because ISIS operated
# broadly across Iraq/Syria, and Houthis across Yemen. A cell in Iraq in
# July 2016 is in the ISIS-drone regime if ISIS was ever active there.

dt[, rebel_has_drones := 0L]
# ISIS theatres: Iraq, Syria
dt[theatre_id %in% c("Iraq", "Syria") & month >= ISIS_ACQ,
   rebel_has_drones := 1L]
# Houthi theatre: Yemen (ACLED-matched name varies)
dt[grepl("Yemen", theatre_id) & month >= HOUTHI_ACQ,
   rebel_has_drones := 1L]

log_msg("Cell-months with rebel_has_drones == 1: ",
        sum(dt$rebel_has_drones == 1L))

# ── 5) RESTRICT SAMPLE TO RELEVANT THEATRES ────────────────────────────
# We want cells in theatres where rebel drone acquisition actually happened.
# Compare pre-acquisition vs post-acquisition WITHIN those theatres.
REL_THEATRES <- unique(c(
  dt[theatre_id %in% c("Iraq", "Syria"),              theatre_id],
  dt[grepl("Yemen", theatre_id),                      theatre_id]
))
dt_acq <- dt[theatre_id %in% REL_THEATRES]
log_msg("Acquisition-relevant subset: ", nrow(dt_acq), " obs, ",
        uniqueN(dt_acq$pg_id), " cells, theatres: ",
        paste(REL_THEATRES, collapse = ", "))

# Log-transform outcomes
dt_acq[, log1p_rebel_vac     := log1p(rebel_vac_predrone + rebel_vac_postdrone)]
dt_acq[, log1p_rebel_vac_sel := log1p(rebel_vac_sel_predrone +
                                       rebel_vac_sel_postdrone)]
dt_acq[, log1p_rebel_vac_ind := log1p(rebel_vac_ind_predrone +
                                       rebel_vac_ind_postdrone)]

CTRL <- c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(CTRL, names(dt_acq)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── 6) RUN ACQUISITION-SHOCK REGRESSIONS ───────────────────────────────
run_acq <- function(yvar, label) {
  f <- as.formula(paste0(
    yvar, " ~ rebel_has_drones:eli_feas + ", ctrl_str, " | ", FE
  ))
  m  <- feols(f, data = dt_acq, notes = FALSE)
  cv <- vcov_conley(m, lat = "lat", lon = "lon", cutoff = 200)
  ct <- as.data.table(coeftable(m, vcov = cv), keep.rownames = "term")
  r  <- ct[grepl("rebel_has_drones", term)]
  data.table(
    outcome = label,
    beta    = round(r$Estimate[1],   5),
    se      = round(r$`Std. Error`[1], 5),
    p_value = round(r$`Pr(>|t|)`[1],   4)
  )
}

cat("\n── R32: rebel drone acquisition shock ──\n\n")
results <- rbindlist(list(
  run_acq("log1p_rebel_vac",     "Rebel VAC — all"),
  run_acq("log1p_rebel_vac_sel", "Rebel VAC — SELECTIVE (≤2 fatal.)"),
  run_acq("log1p_rebel_vac_ind", "Rebel VAC — INDISCRIMINATE (≥5 fatal.)")
))
print(results)
fwrite(results, file.path(OUT, "R32_rebel_acquisition_shock.csv"))

# ── 7) INTERPRETATION ──────────────────────────────────────────────────
b_sel <- results[grepl("SELECTIVE",      outcome), beta]
p_sel <- results[grepl("SELECTIVE",      outcome), p_value]
b_ind <- results[grepl("INDISCRIMINATE", outcome), beta]

cat("\n── INTERPRETATION ──\n")
cat(sprintf("β(rebel_has_drones × eli_feas) on SELECTIVE rebel VAC = %.5f  p = %.4f\n",
            b_sel, p_sel))
cat(sprintf("β(rebel_has_drones × eli_feas) on INDISCRIMINATE VAC  = %.5f\n", b_ind))

cat("\nDisplacement prediction: β_selective < 0 (coercion pressure drops when\n")
cat("asymmetry shrinks). β_indiscriminate ≈ 0 (wasn't driving the pattern).\n\n")

if (!is.na(b_sel) && b_sel < 0 && p_sel < 0.10) {
  cat("✓ Supported. Rebels reduce selective civilian violence after acquiring\n")
  cat("  their own drones — consistent with the reverse-denunciation story.\n")
  cat("  This is INDEPENDENT causal evidence: a structurally different shock\n")
  cat("  (rebel capability acquisition) from the ISR-withdrawal shocks, and\n")
  cat("  it points the same way.\n")
} else if (!is.na(b_sel) && b_sel > 0) {
  cat("⚠ Opposite sign. Rebel drone capability associated with MORE selective\n")
  cat("  VAC. Possible explanations:\n")
  cat("    (i) drone acquisition reflects organisational capacity that enables\n")
  cat("        more sophisticated selective targeting — substitution within\n")
  cat("        the rebel camp.\n")
  cat("    (ii) geographic range expansion post-acquisition.\n")
  cat("  Either way, this deserves discussion in the paper.\n")
} else {
  cat("→ Null or imprecisely estimated. Sample may be too small or acquisition\n")
  cat("  timing coincides with confounds (e.g. ISIS territorial decline).\n")
}

log_msg("R32 complete. Outputs in: ", OUT)
