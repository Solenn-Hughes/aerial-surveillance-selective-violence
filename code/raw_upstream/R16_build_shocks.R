# ============================================================
# R16_build_shocks.R — Code quasi-experiment shock variables
#
# Adds binary (0/1) shock indicators to the panel for:
#   CORE (presented in full):
#     shock_bagram          — Bagram closure, post July 2021 (Afghanistan)
#     shock_somalia_wd      — Somalia troop withdrawal, post Jan 2021 (Somalia)
#     shock_biden_pause     — Biden strike pause, post Jan 2021 (Somalia + Yemen)
#   ROBUSTNESS (Section 3.4):
#     shock_incirlik        — Incirlik restrictions, post Oct 2019 (Syria + Iraq)
#     shock_starlink        — Starlink geofencing proxy (Ukraine, 2022+)
#     shock_trump2017       — Trump active hostilities expansion, post Mar 2017 (Somalia + Yemen)
#     shock_solar           — May 2024 solar storm (Ukraine, narrow window)
#
# Also adds eli_feas interaction terms for each shock.
#
# Input:  panel_pg_month_ALLCW_paper1_posteli_v2.rds (from R08v2 with FOB distances)
# Output: panel_pg_month_ALLCW_paper1_shocks.rds
#
# Run AFTER R08v2 (with FOB distance addition).
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 16: Build shock variables")

# ── LOAD PANEL ───────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL, " — run R08v2 first.")

dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

# Verify required columns exist
req_cols <- c("pg_id", "month", "theatre_id", "post_eli", "eli_feas")
assert_cols(dt, req_cols, "panel")

# Check FOB distances are present (from R08v2 addition)
dist_cols <- c("dist_bagram_km", "dist_somalia_fob_km", "dist_djibouti_km",
               "dist_incirlik_km", "dist_muwaffaq_km")
miss_dist <- setdiff(dist_cols, names(dt))
if (length(miss_dist)) {
  stop("Missing distance columns: ", paste(miss_dist, collapse = ", "),
       "\nRe-run R08v2 with the FOB distance addition.")
}

log_msg("Panel loaded: ", nrow(dt), " rows, ", uniqueN(dt$theatre_id), " theatres")

# ══════════════════════════════════════════════════════════════
# SHOCK 1: BAGRAM CLOSURE (Afghanistan, July 2021)
# ══════════════════════════════════════════════════════════════
# Exogeneity: Doha Agreement (Feb 2020) + Biden withdrawal (Apr 2021)
# Effect: MQ-9 base removed. Cells within ~1000km lose ISR.
# Coding: binary post July 2021, Afghanistan only.
# Spatial: distance to Bagram already computed.

BAGRAM_DATE <- as.IDate("2021-07-01")
MQ9_RANGE_KM <- 1000  # operational combat radius

dt[, shock_bagram := 0L]
dt[theatre_id == "Afghanistan" & month >= BAGRAM_DATE, shock_bagram := 1L]

# Continuous distance-based intensity (for diff-in-RD)
dt[, bagram_in_range := as.integer(dist_bagram_km <= MQ9_RANGE_KM)]
# Normalised distance (0 at Bagram, 1 at range limit)
dt[theatre_id == "Afghanistan", bagram_range_norm := pmin(dist_bagram_km / MQ9_RANGE_KM, 2)]

log_msg("Bagram shock: ", dt[shock_bagram == 1L, .N], " treated cell-months")

# ══════════════════════════════════════════════════════════════
# SHOCK 2: SOMALIA TROOP WITHDRAWAL (Dec 2020 – Jan 2021)
# ══════════════════════════════════════════════════════════════
# Exogeneity: Trump "end forever wars" + presidential transition
# Effect: Forward operators removed from Mogadishu/Baledogle.
#         Strikes shift to Djibouti (~1500km away).
# Coding: binary post Jan 2021, Somalia only.
# Spatial: distance to nearest former FOB (interior-to-interior per Thomas).

SOMALIA_WD_DATE <- as.IDate("2021-01-01")

dt[, shock_somalia_wd := 0L]
dt[theatre_id == "Somalia" & month >= SOMALIA_WD_DATE, shock_somalia_wd := 1L]

# Interior-to-interior: cells near former FOBs vs cells not near FOBs
# Both groups are "interior" — this avoids coast/interior confounding
# Median split on dist_somalia_fob_km within Somalia
somalia_med_fob <- dt[theatre_id == "Somalia", median(dist_somalia_fob_km, na.rm = TRUE)]
dt[, somalia_near_fob := as.integer(dist_somalia_fob_km <= somalia_med_fob)]

log_msg("Somalia withdrawal shock: ", dt[shock_somalia_wd == 1L, .N], " treated cell-months")
log_msg("Somalia FOB median distance: ", round(somalia_med_fob), " km")

# ══════════════════════════════════════════════════════════════
# SHOCK 3: BIDEN DRONE PAUSE (Somalia + Yemen, Jan 2021)
# ══════════════════════════════════════════════════════════════
# Exogeneity: Sullivan order on inauguration day (Jan 20, 2021)
# Effect: Strikes halted ~6 months; ISR continued from Djibouti.
# Tests: information-only vs information+strike mechanism.
# Coding: binary for Jan–Jul 2021 in Somalia + Yemen.

BIDEN_PAUSE_START <- as.IDate("2021-01-01")
BIDEN_PAUSE_END   <- as.IDate("2021-07-01")  # first Biden strikes in Somalia ~July 2021

dt[, shock_biden_pause := 0L]
dt[theatre_id %in% c("Somalia", "Yemen (North Yemen)") &
     month >= BIDEN_PAUSE_START & month < BIDEN_PAUSE_END,
   shock_biden_pause := 1L]

log_msg("Biden pause shock: ", dt[shock_biden_pause == 1L, .N], " treated cell-months")

# ══════════════════════════════════════════════════════════════
# SHOCK 4 (ROBUSTNESS): INCIRLIK RESTRICTIONS (Syria + Iraq, Oct 2019)
# ══════════════════════════════════════════════════════════════
INCIRLIK_DATE <- as.IDate("2019-10-01")

dt[, shock_incirlik := 0L]
dt[theatre_id %in% c("Syria", "Iraq") & month >= INCIRLIK_DATE, shock_incirlik := 1L]

# Spatial: distance shift from Incirlik to Muwaffaq Salti
# Cells where Incirlik was closer than Muwaffaq are "treated"
dt[, incirlik_was_closer := as.integer(dist_incirlik_km < dist_muwaffaq_km)]

log_msg("Incirlik shock: ", dt[shock_incirlik == 1L, .N], " treated cell-months")

# ══════════════════════════════════════════════════════════════
# SHOCK 5 (ROBUSTNESS): STARLINK GEOFENCING (Ukraine, 2022+)
# ══════════════════════════════════════════════════════════════
# Proxy: we don't have exact geofencing dates by area.
# Coding as a theatre-level indicator for awareness that
# connectivity disruptions occurred. Refine with ISW data later.
# For now: treat all Ukraine post-Feb 2022 as "Starlink-dependent"
# and use reported disruption months as the shock.

# Known major disruption periods (approximate from reporting):
# Sep-Oct 2022: geofencing around Crimea (Sevastopol drone attack)
# Oct 2022: Kherson/Zaporizhzhia/Donetsk/Luhansk outages
STARLINK_DISRUPTIONS <- as.IDate(c("2022-09-01", "2022-10-01"))

dt[, shock_starlink := 0L]
dt[theatre_id == "Ukraine" & month %in% STARLINK_DISRUPTIONS, shock_starlink := 1L]

log_msg("Starlink shock: ", dt[shock_starlink == 1L, .N], " treated cell-months")

# ══════════════════════════════════════════════════════════════
# SHOCK 6 (ROBUSTNESS): TRUMP ACTIVE HOSTILITIES (Somalia + Yemen, Mar 2017)
# ══════════════════════════════════════════════════════════════
# Tests the "ON" direction: strike authority expanded.
TRUMP_AH_DATE <- as.IDate("2017-03-01")

dt[, shock_trump2017 := 0L]
dt[theatre_id %in% c("Somalia", "Yemen (North Yemen)") &
     month >= TRUMP_AH_DATE,
   shock_trump2017 := 1L]

log_msg("Trump 2017 shock: ", dt[shock_trump2017 == 1L, .N], " treated cell-months")

# ══════════════════════════════════════════════════════════════
# SHOCK 7 (ROBUSTNESS): SOLAR STORM (Ukraine, May 2024)
# ══════════════════════════════════════════════════════════════
# Requires ACLED daily data for Ukraine. Monthly proxy here.
SOLAR_STORM_MONTH <- as.IDate("2024-05-01")

dt[, shock_solar := 0L]
dt[theatre_id == "Ukraine" & month == SOLAR_STORM_MONTH, shock_solar := 1L]

log_msg("Solar storm shock: ", dt[shock_solar == 1L, .N], " treated cell-months")

# ══════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════

# Audit: shock activation summary
shock_cols <- c("shock_bagram", "shock_somalia_wd", "shock_biden_pause",
                "shock_incirlik", "shock_starlink", "shock_trump2017", "shock_solar")

shock_summary <- rbindlist(lapply(shock_cols, function(sc) {
  data.table(
    shock = sc,
    n_treated = dt[get(sc) == 1L, .N],
    n_cells   = dt[get(sc) == 1L, uniqueN(pg_id)],
    n_months  = dt[get(sc) == 1L, uniqueN(month)],
    theatres  = dt[get(sc) == 1L, paste(unique(theatre_id), collapse = ", ")]
  )
}))

cat("\n============================================================\n")
cat("SHOCK SUMMARY\n")
cat("============================================================\n")
print(shock_summary)
write_audit_dt(shock_summary, "R16_shock_summary")

OUT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
saveRDS(dt, OUT)
log_msg("Wrote: ", OUT)
log_msg("Step 16 complete.")
