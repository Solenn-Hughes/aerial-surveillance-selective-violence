# ============================================================================
# R36_design3_geographic_rd.R
#
# DESIGN 3 — LOCATION SIGNATURE, GEOGRAPHIC RD
#
# Purpose: Identify the causal effect of drone surveillance exposure
#   using spatial discontinuity at the 400 km MQ-9 operational radius.
#
#   At that threshold, cells on either side are observably similar in
#   terrain, population, conflict intensity — different only in being
#   within documented drone operational range. This gives a causal test
#   independent of any shock timing.
#
# METHOD:
#   Running variable: d_i = distance to nearest active military airbase
#                            (in km, during active ISR period)
#   Treatment:         T_i = 1(d_i ≤ 400)
#   Cutoff:            c = 400 km (documented MQ-9 Reaper radius)
#   Outcome:           log(1 + type3)_it, aggregated to cell-level
#                       averages during the ISR active period
#
# SPECIFICATION:
#   Local-linear RD with Calonico-Cattaneo-Titiunik (2014) bias
#   correction and Calonico-Cattaneo-Farrell (2020) coverage-optimal
#   bandwidth selection.
#
# ROBUSTNESS:
#   (1) Bandwidth sensitivity: 200, 300, 400, 500, 600 km
#   (2) Donut-hole RD: drop cells within 50 km of cutoff (McCrary 2008)
#   (3) Placebo cutoffs at 300 km and 500 km (should be null)
#   (4) McCrary density test (manipulation of running variable)
#   (5) Covariate balance tests
#
# IDENTIFYING ASSUMPTION:
#   Potential outcomes are continuous at the 400 km threshold. That is,
#   absent drone surveillance, cells at 399 km and 401 km have the same
#   expected violence level.
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
})

# Install rdrobust if not available
if (!requireNamespace("rdrobust", quietly = TRUE)) {
  install.packages("rdrobust")
}
library(rdrobust)
# Density test for McCrary
if (!requireNamespace("rddensity", quietly = TRUE)) {
  install.packages("rddensity")
}
library(rddensity)

source("R00_config.R")
log_msg("R36: Design 3 — Geographic Regression Discontinuity")

OUT <- file.path("outputs", "paper1", "design3_geographic_rd")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1) LOAD PANEL ──────────────────────────────────────────────────────
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(panel_path); setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_type3 := log1p(events_type3)]

# Coordinates fallback
if (!all(c("lat", "lon") %in% names(dt))) {
  if (all(c("ycoord", "xcoord") %in% names(dt))) {
    dt[, `:=`(lat = ycoord, lon = xcoord)]
  }
}

# ── 2) CONSTRUCT RUNNING VARIABLE: DISTANCE TO NEAREST AIRBASE ────────
# We need distance_to_airbase as a continuous measure.
# If not already in the panel, reconstruct from eli_feas (reverse-engineer).
# eli_feas = exp(-d/400) × observability
# log(eli_feas / observability) = -d/400
# d = -400 × log(eli_feas / observability)
#
# But the cleaner route: reload the airbase data and compute distance directly.

AIRBASE_PATH <- file.path(DIRS$data_raw, "airbases", "airbases.csv")
if (!file.exists(AIRBASE_PATH)) {
  AIRBASE_PATH <- file.path(PROJECT_ROOT, "data", "raw", "airbases",
                             "airbases.csv")
}
if (!file.exists(AIRBASE_PATH)) {
  # Fallback: reverse-engineer from eli_feas
  log_msg("Airbase CSV not found — reverse-engineering from eli_feas.")
  # Need terrain observability component. Approximate using
  # ruggedness-based variables if we have them.
  # Simplest: assume observability varies little and extract -400 * log(eli_feas)
  # as an approximate distance proxy.
  dt[, dist_to_airbase := pmax(1e-6, -400 * log(pmax(eli_feas, 1e-6)))]
  log_msg("WARNING: using approximate distance. Load airbases.csv for precise RD.")
} else {
  airbases <- fread(AIRBASE_PATH)
  setnames(airbases, tolower(names(airbases)))
  stopifnot(all(c("lat", "lon") %in% names(airbases)))

  # Compute minimum distance from each cell to any airbase
  haversine <- function(lat1, lon1, lat2, lon2) {
    R  <- 6371
    to_rad <- pi / 180
    dlat <- (lat2 - lat1) * to_rad
    dlon <- (lon2 - lon1) * to_rad
    a <- sin(dlat/2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) *
         sin(dlon/2)^2
    2 * R * asin(pmin(1, sqrt(a)))
  }

  cells_uniq <- unique(dt[, .(pg_id, lat, lon)])
  log_msg("Computing cell-to-airbase distances for ", nrow(cells_uniq), " cells")
  cells_uniq[, dist_to_airbase := vapply(seq_len(.N), function(i) {
    d <- haversine(lat[i], lon[i], airbases$lat, airbases$lon)
    min(d, na.rm = TRUE)
  }, numeric(1))]
  dt <- merge(dt, cells_uniq[, .(pg_id, dist_to_airbase)],
              by = "pg_id", all.x = TRUE)
}

cat("Distance-to-airbase distribution:\n")
print(summary(dt$dist_to_airbase))

# ── 3) COLLAPSE TO CELL-LEVEL AVERAGES (DURING ACTIVE ISR PERIOD) ─────
# RD at the cell level — one observation per cell, outcome = mean VAC
# during the ISR-active period.
dt_isr <- dt[post_eli == 1L]
cell_avg <- dt_isr[, .(
  mean_type3      = mean(events_type3, na.rm = TRUE),
  mean_log1p      = mean(log1p_type3,  na.rm = TRUE),
  n_obs           = .N,
  dist_to_airbase = dist_to_airbase[1],
  lat             = lat[1], lon = lon[1],
  theatre_id      = theatre_id[1],
  nlights_mean    = mean(nlights_calib_mean, na.rm = TRUE),
  forest_ih       = mean(forest_ih, na.rm = TRUE),
  water_ih        = mean(water_ih, na.rm = TRUE)
), by = pg_id]
cell_avg <- cell_avg[!is.na(dist_to_airbase) & n_obs >= 12L]
log_msg("RD cell sample: ", nrow(cell_avg), " cells")

# ── 4) MAIN RD — 400 KM THRESHOLD ─────────────────────────────────────
cat("\n== MAIN RD — 400 km threshold ==\n")

rd_main <- rdrobust(
  y      = cell_avg$mean_log1p,
  x      = cell_avg$dist_to_airbase,
  c      = 400,
  p      = 1,         # local-linear
  bwselect = "mserd"  # MSE-optimal bandwidth
)
summary(rd_main)

rd_main_out <- data.table(
  cutoff       = 400,
  bandwidth_h  = rd_main$bws[1, 1],
  n_left       = rd_main$N_h[1],
  n_right      = rd_main$N_h[2],
  coef_conv    = rd_main$coef[1, 1],
  coef_bc      = rd_main$coef[2, 1],
  coef_rb      = rd_main$coef[3, 1],
  se_rb        = rd_main$se[3, 1],
  p_rb         = rd_main$pv[3, 1],
  ci_lower_rb  = rd_main$ci[3, 1],
  ci_upper_rb  = rd_main$ci[3, 2]
)
cat("Main RD results (robust bias-corrected):\n")
print(rd_main_out)
fwrite(rd_main_out, file.path(OUT, "R36_rd_main_400km.csv"))

# ── 5) BANDWIDTH SENSITIVITY ──────────────────────────────────────────
cat("\n== BANDWIDTH SENSITIVITY ==\n")
bw_tests <- rbindlist(lapply(c(100, 150, 200, 250, 300, 400, 500), function(bw) {
  res <- tryCatch({
    rdrobust(
      y = cell_avg$mean_log1p, x = cell_avg$dist_to_airbase,
      c = 400, p = 1, h = bw
    )
  }, error = function(e) NULL)
  if (is.null(res)) return(data.table(bandwidth = bw))
  data.table(
    bandwidth   = bw,
    coef_bc     = res$coef[2, 1],
    coef_rb     = res$coef[3, 1],
    se_rb       = res$se[3, 1],
    p_rb        = res$pv[3, 1],
    n_left      = res$N_h[1],
    n_right     = res$N_h[2]
  )
}), fill = TRUE)
print(bw_tests)
fwrite(bw_tests, file.path(OUT, "R36_rd_bandwidth_sensitivity.csv"))

# ── 6) DONUT-HOLE RD (DROP CELLS WITHIN 50 KM OF CUTOFF) ──────────────
cat("\n== DONUT-HOLE RD (exclude ±50 km around 400 km) ==\n")
cell_avg_donut <- cell_avg[abs(dist_to_airbase - 400) > 50]
log_msg("Donut-hole cells: ", nrow(cell_avg_donut))

rd_donut <- tryCatch({
  rdrobust(
    y = cell_avg_donut$mean_log1p,
    x = cell_avg_donut$dist_to_airbase,
    c = 400, p = 1, bwselect = "mserd"
  )
}, error = function(e) { message("Donut RD failed: ", e$message); NULL })

if (!is.null(rd_donut)) {
  donut_out <- data.table(
    design       = "donut-hole (±50 km)",
    coef_rb      = rd_donut$coef[3, 1],
    se_rb        = rd_donut$se[3, 1],
    p_rb         = rd_donut$pv[3, 1],
    bandwidth_h  = rd_donut$bws[1, 1]
  )
  print(donut_out)
  fwrite(donut_out, file.path(OUT, "R36_rd_donut.csv"))
}

# ── 7) PLACEBO CUTOFFS ────────────────────────────────────────────────
cat("\n== PLACEBO CUTOFFS (should be ~0) ==\n")
placebos <- rbindlist(lapply(c(200, 250, 300, 350, 500, 600, 700), function(cc) {
  res <- tryCatch({
    rdrobust(
      y = cell_avg$mean_log1p,
      x = cell_avg$dist_to_airbase,
      c = cc, p = 1, bwselect = "mserd"
    )
  }, error = function(e) NULL)
  if (is.null(res)) return(data.table(cutoff = cc))
  data.table(
    cutoff   = cc,
    coef_rb  = res$coef[3, 1],
    se_rb    = res$se[3, 1],
    p_rb     = res$pv[3, 1]
  )
}), fill = TRUE)
print(placebos)
fwrite(placebos, file.path(OUT, "R36_rd_placebo_cutoffs.csv"))

# ── 8) McCRARY DENSITY TEST ───────────────────────────────────────────
cat("\n== McCRARY DENSITY TEST ==\n")
density_test <- tryCatch({
  rddensity(X = cell_avg$dist_to_airbase, c = 400)
}, error = function(e) { message("Density test failed: ", e$message); NULL })
if (!is.null(density_test)) {
  summary(density_test)
  dens_out <- data.table(
    test       = "McCrary manipulation test",
    t_stat     = density_test$test$t_jk,
    p_value    = density_test$test$p_jk,
    interpretation = ifelse(density_test$test$p_jk > 0.05,
                            "NO manipulation detected (GOOD — running var is continuous)",
                            "MANIPULATION detected (running var has jump — RD invalid)")
  )
  print(dens_out)
  fwrite(dens_out, file.path(OUT, "R36_density_test.csv"))
}

# ── 9) COVARIATE BALANCE AT THRESHOLD ─────────────────────────────────
cat("\n== COVARIATE BALANCE AT CUTOFF ==\n")
balance_out <- rbindlist(lapply(
  c("nlights_mean", "forest_ih", "water_ih"), function(cc) {
    xx <- cell_avg[[cc]]
    if (all(is.na(xx))) return(data.table(covariate = cc))
    res <- tryCatch({
      rdrobust(
        y = xx,
        x = cell_avg$dist_to_airbase,
        c = 400, p = 1, bwselect = "mserd"
      )
    }, error = function(e) NULL)
    if (is.null(res)) return(data.table(covariate = cc))
    data.table(
      covariate = cc,
      coef_rb   = res$coef[3, 1],
      se_rb     = res$se[3, 1],
      p_rb      = res$pv[3, 1]
    )
  }), fill = TRUE)
cat("Expected: all covariates should have p > 0.10 at the threshold.\n")
print(balance_out)
fwrite(balance_out, file.path(OUT, "R36_covariate_balance.csv"))

# ── 10) INTERPRETATION ────────────────────────────────────────────────
cat("\n== INTERPRETATION ==\n")
main_coef  <- rd_main$coef[3, 1]
main_p     <- rd_main$pv[3, 1]
cat(sprintf("Main RD (400km): β = %.5f  p = %.4f\n", main_coef, main_p))
cat(sprintf("Donut-hole: β = %.5f  p = %.4f\n",
            ifelse(!is.null(rd_donut), rd_donut$coef[3, 1], NA),
            ifelse(!is.null(rd_donut), rd_donut$pv[3, 1],   NA)))

if (!is.na(main_p) && main_p < 0.10) {
  cat("\n✓ Causal effect identified at 400 km threshold.\n")
  cat("  Cells within drone operational radius show significantly\n")
  cat("  different violence patterns from cells just outside.\n")
  cat("  This is an identification channel INDEPENDENT of shock timing.\n")
} else {
  cat("\n⚠ No discontinuity at 400 km. Consider:\n")
  cat("  (i) Drone radius varies by platform — may need airbase-specific cutoffs\n")
  cat("  (ii) Coverage gradient is smooth, not sharp — use fuzzy RD with\n")
  cat("       documented coverage probabilities as instrument\n")
}

log_msg("R36 complete. Outputs in: ", OUT)
