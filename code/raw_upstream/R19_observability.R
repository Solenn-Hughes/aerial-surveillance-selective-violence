# ============================================================
# R19_observability.R — Observability / Reporting-Bias Checks
#
# PURPOSE
#   Address the concern that drone ELI improves event *reporting*
#   rather than (or in addition to) changing actual violence.
#   Three falsification tests, ordered by analytic weight:
#
#   A) Actor decomposition
#      If drones raise counts via improved observation, the effect
#      should be symmetric across ALL actor types (state-VAC,
#      rebel-VAC, and even air/drone strikes themselves).
#      A real displacement effect predicts rebel-VAC >> state-VAC
#      and no effect on air/drone strikes (which are not "hidden").
#
#   B) Stable-coverage subsample
#      Restrict to cells with ACLED activity in ≥ 20% of their
#      pre-treatment months (confirmed observable before ELI onset).
#      If results hold in this subsample, coverage expansion cannot
#      account for the main finding.
#
#   C) Biden pause observability falsification (Somalia + Yemen)
#      During the Biden pause (Jan 2021 – Aug 2022), drone strikes
#      fell sharply.  Observability predicts: acled_airdrone_strike
#      and acled_rebel_vac should BOTH fall (less observation →
#      fewer recorded events).  A real deterrence-rebound effect
#      predicts: acled_rebel_vac RISES (rebels expand) while
#      acled_airdrone_strike falls.  Compare the two β coefficients.
#
# INPUTS
#   data/processed/panels/panel_pg_month_ALLCW_paper1_shocks.rds
#   data/processed/base/acled_pg_month.rds
#
# OUTPUTS (to outputs/paper1/ and outputs/paper1/audit/)
#   R19_observability_actor_decomp_<TS>.csv
#   R19_observability_stable_coverage_<TS>.csv
#   R19_observability_biden_pause_obs_<TS>.csv
#   R19_observability_audit_<TS>.csv
#
# KNOWN ISSUE (inherited from R17/R18)
#   Single-theatre FE spec "theatre_id^month_int" fails when there
#   is only one theatre in the data.  Section C detects this and
#   falls back to "pg_id + month_int" automatically.
#
# Run AFTER R16 (shocks panel) and R12 (ACLED pg-month).
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found (R00_config.R or R/00_config.R).")
source(cfg_path)

if (!requireNamespace("fixest", quietly = TRUE))
  stop("Package 'fixest' required. Install via install.packages('fixest').")
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_msg("Step 19: Observability / reporting-bias checks — ", TS)

# ── PATHS ────────────────────────────────────────────────────
P_SHOCKS <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ACLED  <- file.path(DIRS$base,   "acled_pg_month.rds")

if (!file.exists(P_SHOCKS)) stop("Missing: ", P_SHOCKS, " — run R16 first.")
if (!file.exists(P_ACLED))  stop("Missing: ", P_ACLED,  " — run R12 first.")

OUT_ACTOR   <- file.path(DIRS$out_paper1, paste0("R19_observability_actor_decomp_",     TS, ".csv"))
OUT_STABLE  <- file.path(DIRS$out_paper1, paste0("R19_observability_stable_coverage_",  TS, ".csv"))
OUT_BIDEN   <- file.path(DIRS$out_paper1, paste0("R19_observability_biden_pause_obs_",  TS, ".csv"))
OUT_AUDIT   <- file.path(DIRS$audit,      paste0("R19_observability_audit_",            TS, ".csv"))

# ── LOAD PANELS ──────────────────────────────────────────────
log_msg("Loading shocks panel...")
dt <- readRDS(P_SHOCKS); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) {
  dt[, month := as_month_idate(month)]
} else {
  dt[, month := as.IDate(month)]
}
dt[, month_int := year(month) * 12L + month(month)]
log_msg("  Shocks panel: ", nrow(dt), " rows, ", uniqueN(dt$theatre_id), " theatres")

log_msg("Loading ACLED pg-month...")
ac <- readRDS(P_ACLED); setDT(ac)
if (exists("as_month_idate")) {
  ac[, month := as_month_idate(month)]
} else {
  ac[, month := as.IDate(month)]
}
log_msg("  ACLED panel: ", nrow(ac), " rows, ", uniqueN(ac$pg_id), " cells")

# Verify ACLED actor columns present
req_acled <- c("acled_vac", "acled_state_vac", "acled_rebel_vac",
               "acled_airdrone_strike", "acled_all_events")
miss_acled <- setdiff(req_acled, names(ac))
if (length(miss_acled)) {
  stop("ACLED panel missing columns: ", paste(miss_acled, collapse = ", "),
       "\nRe-run R12_build_acled_controls.R.")
}

# ── MERGE ACLED INTO MAIN PANEL ──────────────────────────────
log_msg("Merging ACLED into shocks panel...")
ac_sub <- ac[, c("pg_id", "month", req_acled), with = FALSE]
dt <- merge(dt, ac_sub, by = c("pg_id", "month"), all.x = TRUE)

# Zero-fill ACLED NAs (cells with no ACLED events = 0, not missing)
for (cc in req_acled) {
  n_na <- sum(is.na(dt[[cc]]))
  if (n_na > 0) {
    dt[is.na(get(cc)), (cc) := 0L]
    log_msg("  Zero-filled ", n_na, " NAs in ", cc)
  }
}

log_msg("Merge complete: ", nrow(dt), " rows")

# ── CONTROLS ─────────────────────────────────────────────────
# Standard controls; intersect with what's actually in the panel
CTRL_CAND <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
               "prec_gpcp", "temp", "water_ih",
               "pop_gpw_sum", "bdist3", "capdist",
               "ttime_mean", "excluded", "agri_ih")
ctrl_cols  <- intersect(CTRL_CAND, names(dt))
ctrl_str   <- paste(ctrl_cols, collapse = " + ")
log_msg("Controls used (", length(ctrl_cols), "): ", ctrl_str)

# ── COORDINATES (needed for Conley) ──────────────────────────
if ("xcoord" %in% names(dt) && "ycoord" %in% names(dt)) {
  dt[, longitude := xcoord]
  dt[, latitude  := ycoord]
} else {
  stop("xcoord/ycoord not found. Check panel build (R04/R07).")
}

# ── HELPERS ──────────────────────────────────────────────────
# Extract a single coefficient row from a fitted model
extract_coef_row <- function(m, term_pattern, section, label, outcome,
                             vcov_type = "Conley_200km") {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  hit <- ct[grepl(term_pattern, term, perl = TRUE)]
  if (nrow(hit) == 0L) {
    return(data.table(
      section = section, label = label, outcome = outcome,
      vcov_type = vcov_type,
      coef = NA_real_, se = NA_real_, t_val = NA_real_, p = NA_real_,
      n_obs = nobs(m), note = "term not found"
    ))
  }
  data.table(
    section   = section,
    label     = label,
    outcome   = outcome,
    vcov_type = vcov_type,
    coef  = round(hit$Estimate[1],       5),
    se    = round(hit$`Std. Error`[1],   5),
    t_val = round(hit$`t value`[1],      3),
    p     = round(hit$`Pr(>|t|)`[1],     4),
    n_obs = nobs(m),
    note  = "ok"
  )
}

# Safe Conley summary (falls back to two-way clustering if Conley fails)
safe_conley <- function(m, cutoff_km = 200) {
  res <- tryCatch(summary(m, vcov = conley(cutoff = cutoff_km)), error = function(e) NULL)
  if (!is.null(res)) return(list(m = res, vcov_type = paste0("Conley_", cutoff_km, "km")))
  res <- tryCatch(summary(m, vcov = vcov_conley(cutoff = cutoff_km)), error = function(e) NULL)
  if (!is.null(res)) return(list(m = res, vcov_type = paste0("Conley_", cutoff_km, "km")))
  # fallback
  res <- tryCatch(summary(m, cluster = ~pg_id + month_int), error = function(e) NULL)
  if (!is.null(res)) return(list(m = res, vcov_type = "TwoWay_cluster_fallback"))
  list(m = summary(m), vcov_type = "OLS_fallback")
}

# Detect if a data.table has only one theatre (FE collapse issue)
single_theatre <- function(d) uniqueN(d$theatre_id) == 1L

# Choose FE spec based on theatre count
pick_fe <- function(d) {
  if (single_theatre(d)) {
    log_msg("  Single-theatre subset detected — using pg_id + month_int FE")
    "pg_id + month_int"
  } else {
    "pg_id + theatre_id^month_int"
  }
}

# Collector
all_results <- list()
audit_rows  <- list()

# ============================================================
# A) ACTOR DECOMPOSITION
# ============================================================
log_msg("\n── SECTION A: Actor decomposition ──────────────────────")

ACTOR_OUTCOMES <- c(
  "acled_vac",           # All VAC (benchmark)
  "acled_rebel_vac",     # Rebel/non-state VAC — displacement mechanism
  "acled_state_vac",     # State VAC — should NOT rise if displacement
  "acled_airdrone_strike"# Air/drone strikes — should NOT track ISR observability
)

fe_actor  <- pick_fe(dt)   # Should be theatre_month for full panel
rhs_actor <- paste0("post_eli:eli_feas + ", ctrl_str)

actor_rows <- list()

for (yy in ACTOR_OUTCOMES) {
  log_msg("  Running actor decomp: outcome = ", yy)
  
  ylog <- paste0("log1p_", yy)
  dt[, (ylog) := log1p(get(yy))]
  
  fml <- as.formula(paste0(ylog, " ~ ", rhs_actor, " | ", fe_actor))
  
  m_raw <- tryCatch(
    feols(fml, data = dt, notes = FALSE),
    error = function(e) { log_msg("  ERROR: ", e$message); NULL }
  )
  if (is.null(m_raw)) next
  
  # Primary: Conley 200km
  res_con  <- safe_conley(m_raw, 200)
  hit_con  <- extract_coef_row(res_con$m, "post_eli.*eli_feas|eli_feas.*post_eli",
                               "A_actor_decomp", yy, yy, res_con$vcov_type)
  actor_rows[[length(actor_rows) + 1]] <- hit_con
  
  # Secondary: two-way cluster
  m_tw <- tryCatch(summary(m_raw, cluster = ~pg_id + month_int), error = function(e) NULL)
  if (!is.null(m_tw)) {
    hit_tw <- extract_coef_row(m_tw, "post_eli.*eli_feas|eli_feas.*post_eli",
                               "A_actor_decomp", yy, yy, "TwoWay_cluster")
    actor_rows[[length(actor_rows) + 1]] <- hit_tw
  }
  
  log_msg("  Done: ", yy, " β = ", hit_con$coef, "  p = ", hit_con$p)
}

actor_dt <- rbindlist(actor_rows, fill = TRUE)
fwrite(actor_dt, OUT_ACTOR)
log_msg("Wrote: ", OUT_ACTOR)
all_results$actor_decomp <- actor_dt

# Summary print
cat("\n============================================================\n")
cat("SECTION A — ACTOR DECOMPOSITION\n")
cat("============================================================\n")
cat("Observability prediction: β ≈ equal across all actors\n")
cat("Displacement prediction : β(rebel_vac) > β(state_vac); β(airdrone) ≈ 0\n\n")
conley_only <- actor_dt[grepl("Conley", vcov_type)]
print(conley_only[, .(outcome, coef, se, p, n_obs)])

# Audit: ratio of rebel to state VAC β (key diagnostic)
b_rebel <- conley_only[outcome == "acled_rebel_vac", coef]
b_state <- conley_only[outcome == "acled_state_vac",  coef]
b_drone <- conley_only[outcome == "acled_airdrone_strike", coef]
if (length(b_rebel) && length(b_state) && !is.na(b_rebel) && !is.na(b_state) && b_state != 0) {
  log_msg("  β(rebel_vac) / β(state_vac) = ", round(b_rebel / b_state, 2),
          "  [>1 supports displacement; ~1 supports observability]")
}
audit_rows$actor_ratio <- data.table(
  section = "A_actor_decomp",
  metric  = c("beta_rebel_vac", "beta_state_vac", "beta_airdrone",
              "ratio_rebel_to_state"),
  value   = c(b_rebel[1], b_state[1], b_drone[1],
              if (!is.na(b_rebel[1]) && !is.na(b_state[1]) && b_state[1] != 0)
                round(b_rebel[1] / b_state[1], 3) else NA_real_)
)

# ============================================================
# B) STABLE-COVERAGE SUBSAMPLE
# ============================================================
log_msg("\n── SECTION B: Stable-coverage subsample ────────────────")
#
# Define "stable ACLED coverage" per cell:
# For treated cells: pre-treatment months where post_eli == 0
# For control cells: all months (never treated → always pre-treatment)
# A cell is "stable" if ≥ 20% of its pre-treatment months have
# acled_all_events > 0 (i.e., ACLED was actively reporting there).
#
# Threshold choice: 20% is deliberately low to be conservative
# (even one event every 5 months suggests ACLED was monitoring).

STABLE_THRESHOLD <- 0.20  # share of pre-treatment months with ≥ 1 event

# Compute pre-treatment ACLED activity per cell
pre_dt <- dt[post_eli == 0L, .(
  pre_months      = .N,
  pre_active      = sum(acled_all_events > 0),
  pre_share_active = mean(acled_all_events > 0)
), by = pg_id]

log_msg("  Pre-treatment cell-months: ", dt[post_eli == 0L, .N])
log_msg("  Cells with any pre-treatment months: ", nrow(pre_dt))
log_msg("  Share of cells ≥ ", STABLE_THRESHOLD * 100, "% active: ",
        round(mean(pre_dt$pre_share_active >= STABLE_THRESHOLD) * 100, 1), "%")

stable_cells <- pre_dt[pre_share_active >= STABLE_THRESHOLD, pg_id]
log_msg("  Stable cells (≥ ", STABLE_THRESHOLD * 100, "% active): ",
        length(stable_cells), " / ", uniqueN(dt$pg_id))

dt_stable <- dt[pg_id %in% stable_cells]
log_msg("  Stable subsample: ", nrow(dt_stable), " cell-months, ",
        uniqueN(dt_stable$theatre_id), " theatres")

if (nrow(dt_stable) < 1000L) {
  warning("Stable subsample very small (", nrow(dt_stable), " rows). ",
          "Consider reducing STABLE_THRESHOLD.")
}

fe_stable  <- pick_fe(dt_stable)
rhs_stable <- paste0("post_eli:eli_feas + ", ctrl_str)
ylog_main  <- "log1p_events_type3"
if (!("events_type3" %in% names(dt_stable))) {
  stop("events_type3 not found. Check R03/R04 pipeline.")
}
dt_stable[, (ylog_main) := log1p(events_type3)]

fml_stable <- as.formula(paste0(ylog_main, " ~ ", rhs_stable, " | ", fe_stable))

m_stable_raw <- tryCatch(
  feols(fml_stable, data = dt_stable, notes = FALSE),
  error = function(e) { log_msg("  ERROR: ", e$message); NULL }
)

stable_rows <- list()

if (!is.null(m_stable_raw)) {
  res_con_s <- safe_conley(m_stable_raw, 200)
  hit_con_s <- extract_coef_row(res_con_s$m, "post_eli.*eli_feas|eli_feas.*post_eli",
                                "B_stable_coverage", "stable_subsample",
                                "events_type3", res_con_s$vcov_type)
  hit_con_s[, n_cells_stable := length(stable_cells)]
  hit_con_s[, threshold_pct  := STABLE_THRESHOLD * 100]
  stable_rows[[1]] <- hit_con_s
  
  m_tw_s <- tryCatch(summary(m_stable_raw, cluster = ~pg_id + month_int), error = function(e) NULL)
  if (!is.null(m_tw_s)) {
    hit_tw_s <- extract_coef_row(m_tw_s, "post_eli.*eli_feas|eli_feas.*post_eli",
                                 "B_stable_coverage", "stable_subsample",
                                 "events_type3", "TwoWay_cluster")
    hit_tw_s[, n_cells_stable := length(stable_cells)]
    hit_tw_s[, threshold_pct  := STABLE_THRESHOLD * 100]
    stable_rows[[2]] <- hit_tw_s
  }
  
  log_msg("  Stable-coverage DiD: β = ", hit_con_s$coef,
          "  p = ", hit_con_s$p, "  (Conley 200km)")
}

# Also run the full-sample reference to print side-by-side
dt[, (ylog_main) := log1p(events_type3)]
fml_full <- as.formula(paste0(ylog_main, " ~ post_eli:eli_feas + ", ctrl_str,
                              " | pg_id + theatre_id^month_int"))
m_full_raw <- tryCatch(
  feols(fml_full, data = dt, notes = FALSE),
  error = function(e) { log_msg("  ERROR full-sample: ", e$message); NULL }
)
if (!is.null(m_full_raw)) {
  res_full <- safe_conley(m_full_raw, 200)
  hit_full <- extract_coef_row(res_full$m, "post_eli.*eli_feas|eli_feas.*post_eli",
                               "B_full_sample_ref", "full_sample",
                               "events_type3", res_full$vcov_type)
  hit_full[, n_cells_stable := uniqueN(dt$pg_id)]
  hit_full[, threshold_pct  := 0]
  stable_rows[[length(stable_rows) + 1]] <- hit_full
}

stable_dt <- rbindlist(stable_rows, fill = TRUE)
fwrite(stable_dt, OUT_STABLE)
log_msg("Wrote: ", OUT_STABLE)
all_results$stable_coverage <- stable_dt

cat("\n============================================================\n")
cat("SECTION B — STABLE-COVERAGE SUBSAMPLE\n")
cat("============================================================\n")
cat("Interpretation: if β(stable) ≈ β(full), coverage expansion\n")
cat("does not explain the main finding.\n\n")
print(stable_dt[grepl("Conley", vcov_type), .(label, n_cells_stable, coef, se, p)])

audit_rows$stable_coverage <- data.table(
  section = "B_stable_coverage",
  metric  = c("n_stable_cells", "n_full_cells", "threshold_pct",
              "beta_stable_Conley200", "beta_full_Conley200"),
  value   = c(length(stable_cells), uniqueN(dt$pg_id),
              STABLE_THRESHOLD * 100,
              stable_dt[label == "stable_subsample" & grepl("Conley", vcov_type), coef][1],
              stable_dt[label == "full_sample"       & grepl("Conley", vcov_type), coef][1])
)

# ============================================================
# C) BIDEN PAUSE OBSERVABILITY FALSIFICATION
# ============================================================
log_msg("\n── SECTION C: Biden pause observability falsification ──")
#
# Logic:
#   Observability prediction: Drone pause → less surveillance →
#     BOTH acled_airdrone_strike AND acled_rebel_vac FALL (β < 0).
#   Real-effect prediction: Drone pause → deterrence removed →
#     acled_rebel_vac RISES (β > 0) while acled_airdrone_strike
#     falls (β < 0, just documenting fewer strikes).
#
#   We estimate: y ~ shock_biden_pause:eli_feas + controls
#   where "treated" = high eli_feas cells in Somalia + Yemen,
#   and the shock period is post Jan 2021.
#
# NOTE: shock_biden_pause is 1 for Somalia + Yemen post Jan 2021.
# Filter to only Somalia + Yemen so the comparison group is
# control cells within those theatres only.

BIDEN_THEATRES <- c("Somalia", "Yemen")

if (!("shock_biden_pause" %in% names(dt))) {
  log_msg("  shock_biden_pause not found — skipping Section C.")
  log_msg("  Run R16_build_shocks.R first.")
} else {
  dt_biden <- dt[theatre_id %in% BIDEN_THEATRES]
  log_msg("  Biden pause subset: ", nrow(dt_biden), " rows, ",
          uniqueN(dt_biden$theatre_id), " theatres")
  log_msg("  Treated cell-months (shock=1): ",
          dt_biden[shock_biden_pause == 1L, .N])
  
  # Automatic FE fallback for single-theatre subset
  fe_biden  <- pick_fe(dt_biden)
  rhs_biden <- paste0("shock_biden_pause:eli_feas + ", ctrl_str)
  
  BIDEN_OUTCOMES <- c(
    "acled_rebel_vac",      # Real effect: should RISE; Observability: should FALL
    "acled_airdrone_strike",# Both: should FALL (direct causal — fewer drones fly)
    "acled_state_vac",      # Additional check: state VAC
    "events_type3"          # UCDP core outcome for reference
  )
  
  biden_rows <- list()
  
  for (yy in BIDEN_OUTCOMES) {
    if (!(yy %in% names(dt_biden))) {
      log_msg("  SKIP: ", yy, " not found")
      next
    }
    log_msg("  Biden pause falsification: outcome = ", yy)
    
    ylog <- paste0("log1p_", yy)
    dt_biden[, (ylog) := log1p(get(yy))]
    
    fml <- as.formula(paste0(ylog, " ~ ", rhs_biden, " | ", fe_biden))
    
    m_raw <- tryCatch(
      feols(fml, data = dt_biden, notes = FALSE),
      error = function(e) { log_msg("  ERROR: ", e$message); NULL }
    )
    if (is.null(m_raw)) next
    
    # Conley primary
    res_con <- safe_conley(m_raw, 200)
    hit_con <- extract_coef_row(
      res_con$m, "shock_biden_pause.*eli_feas|eli_feas.*shock_biden_pause",
      "C_biden_obs", yy, yy, res_con$vcov_type
    )
    hit_con[, fe_spec  := fe_biden]
    hit_con[, theatres := paste(BIDEN_THEATRES, collapse = "+")]
    biden_rows[[length(biden_rows) + 1]] <- hit_con
    
    # Two-way cluster secondary
    m_tw <- tryCatch(summary(m_raw, cluster = ~pg_id + month_int), error = function(e) NULL)
    if (!is.null(m_tw)) {
      hit_tw <- extract_coef_row(
        m_tw, "shock_biden_pause.*eli_feas|eli_feas.*shock_biden_pause",
        "C_biden_obs", yy, yy, "TwoWay_cluster"
      )
      hit_tw[, fe_spec  := fe_biden]
      hit_tw[, theatres := paste(BIDEN_THEATRES, collapse = "+")]
      biden_rows[[length(biden_rows) + 1]] <- hit_tw
    }
    
    log_msg("  Done: ", yy, " β = ", hit_con$coef, "  p = ", hit_con$p)
  }
  
  if (length(biden_rows)) {
    biden_dt <- rbindlist(biden_rows, fill = TRUE)
    fwrite(biden_dt, OUT_BIDEN)
    log_msg("Wrote: ", OUT_BIDEN)
    all_results$biden_obs <- biden_dt
    
    cat("\n============================================================\n")
    cat("SECTION C — BIDEN PAUSE OBSERVABILITY FALSIFICATION\n")
    cat("============================================================\n")
    cat("Observability prediction:\n")
    cat("  β(rebel_vac) < 0  AND  β(airdrone) < 0  (both fall during pause)\n")
    cat("Real-effect prediction:\n")
    cat("  β(rebel_vac) > 0  AND  β(airdrone) < 0  (rebels expand; fewer strikes)\n\n")
    print(biden_dt[grepl("Conley", vcov_type), .(outcome, coef, se, p, fe_spec)])
    
    # Key diagnostic
    b_rv  <- biden_dt[outcome == "acled_rebel_vac"      & grepl("Conley", vcov_type), coef][1]
    b_ad  <- biden_dt[outcome == "acled_airdrone_strike" & grepl("Conley", vcov_type), coef][1]
    cat("\nKey diagnostic:\n")
    cat("  β(rebel_vac):   ", round(b_rv, 4), "\n")
    cat("  β(airdrone):    ", round(b_ad, 4), "\n")
    if (!is.na(b_rv) && !is.na(b_ad)) {
      if (b_rv > 0 && b_ad < 0) {
        cat("  → Signs DIVERGE: consistent with REAL deterrence-rebound effect.\n")
        cat("    Observability alone cannot explain opposite-sign divergence.\n")
      } else if (b_rv < 0 && b_ad < 0) {
        cat("  → Signs both negative: consistent with observability mechanism.\n")
        cat("    Cannot rule out reporting-bias explanation for Biden pause.\n")
        cat("    NOTE: main result (drone onset → violence UP) still requires\n")
        cat("    replication in stable-coverage subsample (Section B).\n")
      } else {
        cat("  → Mixed/null: inconclusive.\n")
      }
    }
    
    audit_rows$biden_obs <- data.table(
      section = "C_biden_obs",
      metric  = c("beta_rebel_vac_Conley200", "beta_airdrone_Conley200",
                  "diverge_sign"),
      value   = c(b_rv, b_ad,
                  if (!is.na(b_rv) && !is.na(b_ad))
                    as.numeric(sign(b_rv) != sign(b_ad)) else NA_real_)
    )
  } else {
    log_msg("  No Biden pause models succeeded.")
  }
}  # end Biden pause block

# ── COMBINED AUDIT ───────────────────────────────────────────
audit_dt <- rbindlist(audit_rows, fill = TRUE)
fwrite(audit_dt, OUT_AUDIT)
log_msg("Wrote audit: ", OUT_AUDIT)

# ── FINAL SUMMARY ────────────────────────────────────────────
cat("\n============================================================\n")
cat("R19 OBSERVABILITY CHECKS — SUMMARY INTERPRETATION\n")
cat("============================================================\n")
cat("A) Actor decomp: Look at β(rebel_vac) vs β(state_vac) vs β(airdrone).\n")
cat("   If β(rebel) >> β(state) and β(airdrone) ≈ 0, displacement favoured.\n")
cat("B) Stable coverage: If β(stable) ≈ β(full), coverage bias is unlikely.\n")
cat("C) Biden pause: If β(rebel_vac) > 0 and β(airdrone) < 0, real effect.\n")
cat("All three tests jointly address the observability critique.\n")
cat("============================================================\n\n")

log_msg("Step 19 complete.")