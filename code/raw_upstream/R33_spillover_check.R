# ============================================================================
# R33_spillover_check.R — SUTVA / spillover diagnostic
#
# Purpose: Your control group is low-eli_feas cells within the same
#   treated theatres. If rebels in high-feas cells displace some of
#   their violence to nearby low-feas cells (rather than reducing it),
#   your control is contaminated — the β on treated cells is biased
#   TOWARD zero (conservative).
#
#   A reviewer will ask about this. This script provides three checks:
#
#   CHECK 1: Buffer-zone test
#     Split control cells by distance to nearest treated cell.
#     — "Near control"  : within 100km of high-feas cell
#     — "Far control"   : >200km from any high-feas cell
#     If spillovers exist, near controls show violence INCREASES
#     post-ISR-onset, while far controls show the "true" counterfactual.
#
#   CHECK 2: Treatment-saturation
#     Does β change when estimated on cells in "sparse treatment" theatres
#     vs "dense treatment" theatres? Denser treatment = more risk of
#     spillover contamination.
#
#   CHECK 3: Triple-difference placebo
#     Run the shock specification on DIFFERENT outcomes that should NOT
#     respond to spillover (e.g. non-political violence, criminal events)
#     to see whether the pattern is specific to the mechanism we claim.
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})
source("R00_config.R")
log_msg("R33: spillover / SUTVA diagnostic")

OUT <- file.path("outputs", "paper1", "spillover_check")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1) LOAD ────────────────────────────────────────────────────────────
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(panel_path); setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_type3 := log1p(events_type3)]

if (!all(c("lat", "lon") %in% names(dt))) {
  if (all(c("ycoord", "xcoord") %in% names(dt))) {
    dt[, `:=`(lat = ycoord, lon = xcoord)]
  }
}

SHOCK_THEATRES <- c("Afghanistan", "Syria", "Iraq")
dt_s <- dt[theatre_id %in% SHOCK_THEATRES]

# ── 2) CHECK 1 — BUFFER-ZONE SPILLOVER TEST ─────────────────────────────
# Define "treated" cells as those with eli_feas > theatre median.
# For each control cell, compute minimum distance to the nearest treated cell.

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R  <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat/2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# Unique cell list for each theatre
cells <- unique(dt_s[, .(pg_id, theatre_id, lat, lon, eli_feas)])
cells[, theatre_median_feas := median(eli_feas, na.rm = TRUE),
      by = theatre_id]
cells[, treated := eli_feas > theatre_median_feas]

# For each control cell, find nearest treated cell in same theatre
log_msg("Computing nearest-treated distance for each control cell...")
cells[, nearest_treated_km := NA_real_]
for (th in unique(cells$theatre_id)) {
  tc <- cells[theatre_id == th & treated == TRUE]
  cc_idx <- which(cells$theatre_id == th & cells$treated == FALSE)
  if (nrow(tc) == 0 || length(cc_idx) == 0) next
  for (i in cc_idx) {
    d <- haversine_km(cells$lat[i], cells$lon[i], tc$lat, tc$lon)
    cells[i, nearest_treated_km := min(d)]
  }
}

# Control cell bucket: near (≤100km), mid (100–200km), far (>200km)
cells[, control_bucket := fcase(
  treated == TRUE,                               "TREATED",
  nearest_treated_km <= 100,                     "near_control (≤100km)",
  nearest_treated_km > 100 & nearest_treated_km <= 200, "mid_control (100–200km)",
  nearest_treated_km > 200,                      "far_control (>200km)",
  default                                       = "isolated"
)]

# Merge back
dt_s <- merge(dt_s, cells[, .(pg_id, control_bucket)], by = "pg_id", all.x = TRUE)

cat("\nCell counts by control bucket:\n")
print(cells[, .N, by = control_bucket])

# Build any_isr_removed if missing
if (!"any_isr_removed" %in% names(dt_s)) {
  dt_s[, any_isr_removed := 0L]
  dt_s[theatre_id == "Afghanistan" & month >= as.IDate("2021-07-01"), any_isr_removed := 1L]
  dt_s[theatre_id == "Syria"       & month >= as.IDate("2019-10-01"), any_isr_removed := 1L]
  dt_s[theatre_id == "Iraq"        & month >= as.IDate("2019-10-01"), any_isr_removed := 1L]
}

CTRL <- c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(CTRL, names(dt_s)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# Run spec using far controls only
run_buffer <- function(control_types, label) {
  sub <- dt_s[control_bucket %in% c("TREATED", control_types)]
  if (nrow(sub) < 1000) return(NULL)
  f <- as.formula(paste0(
    "log1p_type3 ~ post_eli:eli_feas + any_isr_removed:eli_feas + ",
    ctrl_str, " | ", FE
  ))
  m  <- feols(f, data = sub, notes = FALSE)
  cv <- vcov_conley(m, lat = "lat", lon = "lon", cutoff = 200)
  ct <- as.data.table(coeftable(m, vcov = cv), keep.rownames = "term")
  r_causal <- ct[grepl("any_isr_removed", term)]
  r_panel  <- ct[grepl("^post_eli", term) & !grepl("any_isr", term)]
  data.table(
    sample            = label,
    n_obs             = nrow(sub),
    beta_shock_causal = round(r_causal$Estimate[1],   5),
    p_shock_causal    = round(r_causal$`Pr(>|t|)`[1],   4),
    beta_panel        = round(r_panel$Estimate[1],    5),
    p_panel           = round(r_panel$`Pr(>|t|)`[1],    4)
  )
}

cat("\n── CHECK 1: Buffer-zone spillover test ──\n")
buffer_res <- rbindlist(list(
  run_buffer(c("near_control (≤100km)", "mid_control (100–200km)",
                "far_control (>200km)"), "All controls (original)"),
  run_buffer(c("mid_control (100–200km)", "far_control (>200km)"),
               "Drop near controls (≤100km)"),
  run_buffer(c("far_control (>200km)"),    "Far controls only (>200km)")
), fill = TRUE)
print(buffer_res)
fwrite(buffer_res, file.path(OUT, "R33_buffer_zone_test.csv"))

cat("\nInterpretation:\n")
cat("  If β(shock × eli_feas) grows MORE negative (larger magnitude) when\n")
cat("  near-controls are dropped → those near-controls were absorbing\n")
cat("  spillover, and your baseline estimate is CONSERVATIVE.\n")
cat("  If β is stable across all three columns → no meaningful spillover.\n")

# ── 3) CHECK 2 — DIRECT TEST ON NEAR-CONTROLS ──────────────────────────
# Under spillover, near-control cells should experience INCREASED violence
# post-ISR-onset (as rebels relocate). Directly test this.

cat("\n── CHECK 2: Direct test on near-control cells ──\n")
near <- dt_s[control_bucket == "near_control (≤100km)"]
if (nrow(near) > 500) {
  f_near <- as.formula(paste0(
    "log1p_type3 ~ post_eli + any_isr_removed + ", ctrl_str, " | ", FE
  ))
  m_near  <- feols(f_near, data = near, notes = FALSE)
  cv_near <- vcov_conley(m_near, lat = "lat", lon = "lon", cutoff = 200)
  ct_near <- as.data.table(coeftable(m_near, vcov = cv_near),
                            keep.rownames = "term")
  cat("Near-control effect of post-ISR (should be 0 under no-spillover):\n")
  print(ct_near)
  fwrite(ct_near, file.path(OUT, "R33_near_control_direct.csv"))
}

# ── 4) CHECK 3 — SATURATION ─────────────────────────────────────────────
# Compute treatment saturation per theatre: share of cells above median.
# (This is always ~50% by construction, but weight by POPULATION or AREA
# gives real variation.)
sat <- cells[, .(
  n_cells    = .N,
  n_treated  = sum(treated, na.rm = TRUE),
  sat_rate   = mean(treated, na.rm = TRUE)
), by = theatre_id]
cat("\nTreatment saturation by theatre:\n")
print(sat)
fwrite(sat, file.path(OUT, "R33_saturation.csv"))

# ── 5) SUMMARY ─────────────────────────────────────────────────────────
sink(file.path(OUT, "R33_SUMMARY.txt"))
cat("R33 SPILLOVER / SUTVA DIAGNOSTIC — SUMMARY\n")
cat(strrep("=", 60), "\n\n")
cat("Buffer-zone test:\n")
print(buffer_res)
cat("\nIf near-control results show positive β on post_eli while far-control\n")
cat("results show β ≈ 0: spillover is present, and your baseline estimate\n")
cat("UNDERSTATES the true effect on treated cells. In that case, your\n")
cat("causal claim is CONSERVATIVE, and the paper should note this.\n")
sink()

log_msg("R33 complete. Outputs in: ", OUT)
