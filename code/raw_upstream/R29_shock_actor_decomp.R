# ============================================================================
# R29_shock_actor_decomp.R — CAUSAL version of the actor decomposition
#
# Purpose: The 6:1 rebel-vs-state ratio is your main substantive finding.
#   Currently it's estimated on the descriptive DiD (post_eli × eli_feas)
#   across all 8 theatres. This script re-estimates the decomposition
#   using the shock variable (any_isr_removed × eli_feas) across the
#   3 shock theatres only. If the rebel effect remains negative and
#   significant while the state effect remains null, the perpetrator
#   signature is causally identified — not just correlational.
#
# Outcomes tested:
#   1. Rebel civilian harm        (main: drives displacement mechanism)
#   2. State civilian harm (corr) (null under displacement)
#   3. State drone strikes        (null under displacement)
#   4. Non-state civilian harm    (rebel + militia broadly)
#
# Design:
#   Sample   : 3 shock theatres (Afghanistan, Syria, Iraq)
#   Treatment: any_isr_removed × eli_feas
#   FE       : pg_id + theatre_id × month_int
#   Inference: Conley 200km spatial HAC
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})
source("R00_config.R")
log_msg("R29: shock-based actor decomposition")

OUT <- file.path("outputs", "paper1", "shock_actor_decomp")
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

# ── 2) BUILD ACTOR-DISAGGREGATED OUTCOMES (with ACLED correction) ──────
# The correction: state drone-strike casualties are coded as
# Explosions/Remote violence in ACLED, not Violence against civilians.
# For the state civilian harm measure, we include both.
ac[, is_vac      := grepl("Violence against civilians", event_type, ignore.case = TRUE)]
ac[, is_remote   := grepl("Explosions|Remote violence",  event_type, ignore.case = TRUE)]
ac[, state_inv   := (inter1 == 1L | inter2 == 1L)]
ac[, rebel_inv   := (inter1 == 2L | inter2 == 2L)]
ac[, militia_inv := (inter1 %in% c(3L,4L) | inter2 %in% c(3L,4L))]
ac[, civ_target  := (inter1 == 7L | inter2 == 7L)]
if ("sub_event_type" %in% names(ac)) {
  ac[, is_airdrone := grepl("Air/drone strike", sub_event_type, ignore.case = TRUE)]
} else ac[, is_airdrone := FALSE]

# Corrected state civilian harm = state VAC + state remote violence against civilians
ac[, state_civ_harm := (state_inv & (is_vac | (is_remote & civ_target)))]
ac[, rebel_civ_harm := (rebel_inv & (is_vac | (is_remote & civ_target)))]
ac[, nonstate_civ_harm := ((rebel_inv | militia_inv) & (is_vac | (is_remote & civ_target)))]
ac[, state_drone := (state_inv & is_airdrone)]

agg_actor <- ac[, .(
  state_civ_harm_n    = sum(state_civ_harm,    na.rm = TRUE),
  rebel_civ_harm_n    = sum(rebel_civ_harm,    na.rm = TRUE),
  nonstate_civ_harm_n = sum(nonstate_civ_harm, na.rm = TRUE),
  state_drone_n       = sum(state_drone,       na.rm = TRUE)
), by = .(pg_id, month)]

# Merge into panel
for (cc in names(agg_actor)[-(1:2)]) if (cc %in% names(dt)) dt[, (cc) := NULL]
dt <- merge(dt, agg_actor, by = c("pg_id", "month"), all.x = TRUE)
for (cc in c("state_civ_harm_n", "rebel_civ_harm_n",
             "nonstate_civ_harm_n", "state_drone_n")) {
  dt[is.na(get(cc)), (cc) := 0L]
}

# Coordinates
if (!all(c("lat", "lon") %in% names(dt))) {
  if (all(c("ycoord", "xcoord") %in% names(dt))) {
    dt[, `:=`(lat = ycoord, lon = xcoord)]
  } else stop("Cannot find lat/lon columns.")
}

# Log outcomes
dt[, log1p_state := log1p(state_civ_harm_n)]
dt[, log1p_rebel := log1p(rebel_civ_harm_n)]
dt[, log1p_ns    := log1p(nonstate_civ_harm_n)]
dt[, log1p_drone := log1p(state_drone_n)]

# ── 3) BUILD any_isr_removed IF NEEDED ─────────────────────────────────
if (!"any_isr_removed" %in% names(dt)) {
  dt[, any_isr_removed := 0L]
  dt[theatre_id == "Afghanistan" & month >= as.IDate("2021-07-01"), any_isr_removed := 1L]
  dt[theatre_id == "Syria"       & month >= as.IDate("2019-10-01"), any_isr_removed := 1L]
  dt[theatre_id == "Iraq"        & month >= as.IDate("2019-10-01"), any_isr_removed := 1L]
}

# ── 4) RESTRICT TO SHOCK THEATRES ──────────────────────────────────────
SHOCK_THEATRES <- c("Afghanistan", "Syria", "Iraq")
dt_shock <- dt[theatre_id %in% SHOCK_THEATRES]
log_msg("Shock subset: ", nrow(dt_shock), " obs, ",
        uniqueN(dt_shock$pg_id), " cells")

CTRL <- c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(CTRL, names(dt_shock)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── 5) RUN EACH OUTCOME ────────────────────────────────────────────────
run_actor <- function(yvar, label, count_col) {
  f <- as.formula(paste0(
    yvar, " ~ post_eli:eli_feas + any_isr_removed:eli_feas + ",
    ctrl_str, " | ", FE
  ))
  m  <- feols(f, data = dt_shock, notes = FALSE)
  cv <- vcov_conley(m, lat = "lat", lon = "lon", cutoff = 200)
  ct <- as.data.table(coeftable(m, vcov = cv), keep.rownames = "term")

  r_causal <- ct[grepl("any_isr_removed", term)]
  r_panel  <- ct[grepl("^post_eli", term) & !grepl("any_isr", term)]

  data.table(
    outcome           = label,
    beta_shock_causal = round(r_causal$Estimate[1],   5),
    se_shock_causal   = round(r_causal$`Std. Error`[1], 5),
    p_shock_causal    = round(r_causal$`Pr(>|t|)`[1],   4),
    beta_panel        = round(r_panel$Estimate[1],   5),
    p_panel           = round(r_panel$`Pr(>|t|)`[1],   4),
    n_nonzero         = sum(dt_shock[[count_col]] > 0, na.rm = TRUE)
  )
}

cat("\n── R29: SHOCK-BASED ACTOR DECOMPOSITION (3-theatre subset) ──\n\n")
results <- rbindlist(list(
  run_actor("log1p_rebel", "Rebel civilian harm",            "rebel_civ_harm_n"),
  run_actor("log1p_state", "State civilian harm (corrected)","state_civ_harm_n"),
  run_actor("log1p_ns",    "Non-state civilian harm",        "nonstate_civ_harm_n"),
  run_actor("log1p_drone", "State drone strikes",            "state_drone_n")
))
print(results)
fwrite(results, file.path(OUT, "R29_shock_actor_decomposition.csv"))

# ── 6) CAUSAL RATIO DIAGNOSTIC ────────────────────────────────────────
b_rebel <- results[outcome == "Rebel civilian harm",             beta_shock_causal]
b_state <- results[outcome == "State civilian harm (corrected)", beta_shock_causal]
p_rebel <- results[outcome == "Rebel civilian harm",             p_shock_causal]
p_state <- results[outcome == "State civilian harm (corrected)", p_shock_causal]

cat("\n── CAUSAL RATIO DIAGNOSTIC ──\n")
cat(sprintf("β(shock × eli_feas) on REBEL civilian harm  = %.5f  p = %.4f\n",
            b_rebel, p_rebel))
cat(sprintf("β(shock × eli_feas) on STATE civilian harm  = %.5f  p = %.4f\n",
            b_state, p_state))

if (!is.na(b_rebel) && !is.na(b_state) && abs(b_state) > 0) {
  ratio <- b_rebel / b_state
  cat(sprintf("\nCausal ratio |rebel/state| = %.2f\n", abs(ratio)))
  cat("\nUnder displacement: rebel effect should be LARGE and significant;\n")
  cat("                    state effect should be SMALL and insignificant.\n")
  cat("Under substitution: the reverse.\n")
}

log_msg("R29 complete. Outputs in: ", OUT)
