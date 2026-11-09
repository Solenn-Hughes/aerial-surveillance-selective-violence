# ============================================================================
# R34_build_dyadic_panel.R
#
# Purpose: Build the cell × group × month dyadic panel required for
#   Designs 2 (perpetrator) and 4 (asymmetry).
#
# Logic:
#   — Identify active rebel groups per theatre from ACLED.
#   — For each (cell, group, month), compute whether the group was
#     "present" in the cell (lagged 6-month window of documented
#     activity — battles OR VAC events).
#   — For each (cell, group, month), compute VAC events attributable
#     to that group in that cell-month.
#   — Merge in the cell-level eli_feas and shock variables.
#   — Output: dyadic panel ready for fixed-effect regressions.
#
# Design choices:
#   (1) "Active group" = any group with >= 20 events in panel, avoids
#       noise from one-off unnamed actors.
#   (2) "Group present in cell" = lagged 6-month indicator, cumulative.
#       Once a group is documented operating in a cell, it remains
#       present until 12 months of silence. This matches the conflict
#       literature's standard operationalization (Raleigh 2016).
#   (3) We create dyads only for (cell, group) combinations where the
#       group was EVER present in the cell's theatre. Otherwise, panel
#       size explodes to tens of millions of zero rows.
#   (4) Full panel expected size: ~5–8M dyad-months.
#
# INPUT:  panel_pg_month_ALLCW_paper1_shocks.rds (R07)
#         data/raw/acled/acled_full.csv
# OUTPUT: panel_dyadic_cell_group_month.rds
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(zoo)
})

source("R00_config.R")
log_msg("R34: building dyadic cell × group × month panel")

OUT_DIR <- DIRS$panels
OUT_PATH <- file.path(OUT_DIR, "panel_dyadic_cell_group_month.rds")

# ── 1) LOAD BASE CELL-MONTH PANEL ─────────────────────────────────────
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
base_panel <- readRDS(panel_path); setDT(base_panel)
base_panel[, month_int := year(month) * 12L + month(month)]
log_msg("Base panel: ", nrow(base_panel), " cell-months, ",
        uniqueN(base_panel$pg_id), " cells, ",
        uniqueN(base_panel$theatre_id), " theatres")

# Identify which variables we want to bring into the dyadic panel
# (cell-level, will be constant across groups within a cell-month)
CELL_VARS <- c(
  "pg_id", "month", "month_int", "theatre_id",
  "eli_feas", "post_eli", "lat", "lon",
  "nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih"
)
# Coordinates fallback
if (!all(c("lat", "lon") %in% names(base_panel))) {
  if (all(c("ycoord", "xcoord") %in% names(base_panel))) {
    base_panel[, `:=`(lat = ycoord, lon = xcoord)]
  }
}
# any_isr_removed (shock variable)
if (!"any_isr_removed" %in% names(base_panel)) {
  base_panel[, any_isr_removed := 0L]
  base_panel[theatre_id == "Afghanistan" & month >= as.IDate("2021-07-01"),
             any_isr_removed := 1L]
  base_panel[theatre_id == "Syria"       & month >= as.IDate("2019-10-01"),
             any_isr_removed := 1L]
  base_panel[theatre_id == "Iraq"        & month >= as.IDate("2019-10-01"),
             any_isr_removed := 1L]
}
CELL_VARS <- c(CELL_VARS, "any_isr_removed")
CELL_VARS <- intersect(CELL_VARS, names(base_panel))
cell_slim <- base_panel[, ..CELL_VARS]

# ── 2) LOAD ACLED AND CLEAN ───────────────────────────────────────────
ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT,  "data", "raw", "acled", "acled_full.csv")
)
acled_path <- Filter(file.exists, ACLED_PATHS)[1]
if (is.na(acled_path)) stop("ACLED raw CSV not found.")

ac <- fread(acled_path, showProgress = FALSE, encoding = "UTF-8")
setnames(ac, tolower(names(ac)))
ac[, month     := as.IDate(format(as.Date(event_date), "%Y-%m-01"))]
ac[, month_int := year(month) * 12L + month(month)]
ac[, fatalities := as.integer(fatalities)]
ac[is.na(fatalities), fatalities := 0L]

# Map to PRIO-GRID
ac[, pg_row := floor((latitude  + 90)  / 0.5) + 1L]
ac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
ac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# Restrict to rebel-initiated events (inter1 == 2) — this is the mechanism actor
rebel_events <- ac[inter1 == 2L & !is.na(actor1) & actor1 != ""]
log_msg("Rebel events in ACLED: ", nrow(rebel_events))

# ── 3) IDENTIFY ACTIVE GROUPS ─────────────────────────────────────────
# Count events per actor; keep groups with >= 20 events total
actor_counts <- rebel_events[, .(n_events = .N), by = actor1]
setorder(actor_counts, -n_events)
active_groups <- actor_counts[n_events >= 20L, actor1]
log_msg("Active groups (>=20 events): ", length(active_groups),
        " / ", nrow(actor_counts))

cat("\nTop-10 active groups by event count:\n")
print(head(actor_counts[actor1 %in% active_groups], 10))

rebel_events <- rebel_events[actor1 %in% active_groups]

# ── 4) MAP GROUPS TO THEATRES ─────────────────────────────────────────
# For each group, identify which theatres they operate in.
# A group "belongs" to a theatre if >= 5 events there (filter noise).
# We need the base panel's theatre mapping for cells.
cell_theatre <- unique(base_panel[, .(pg_id, theatre_id)])
rebel_events_t <- merge(rebel_events, cell_theatre, by = "pg_id", all.x = FALSE)

group_theatre <- rebel_events_t[, .(n_events = .N),
                                  by = .(actor1, theatre_id)]
group_theatre <- group_theatre[n_events >= 5L]
log_msg("Group-theatre pairs: ", nrow(group_theatre))

# ── 5) BUILD "PRESENCE" INDICATOR (6-month lagged window) ─────────────
# A group is "present" in a cell in month t if it had any documented
# event (battle OR VAC) in that cell in the 6 months leading up to t.
# We use both battles and VAC because presence is about operational
# activity, not specifically about civilian targeting.

# Include battles in the presence measure
battle_events <- ac[inter1 == 2L & !is.na(actor1) &
                     grepl("^Battles", event_type, ignore.case = TRUE) &
                     actor1 %in% active_groups]
all_presence_events <- rbindlist(list(
  rebel_events_t[, .(pg_id, actor1, month_int)],
  merge(battle_events[, .(pg_id, actor1, month_int)], cell_theatre,
        by = "pg_id", all.x = FALSE)[, .(pg_id, actor1, month_int)]
), use.names = TRUE)
all_presence_events <- unique(all_presence_events)

# For each (pg_id, actor1) with activity, create a 6-month forward window
# i.e. presence = 1 in month t if any event in [t-5, t-0]
log_msg("Building 6-month lagged presence indicator...")

# Group's first event month in each cell
group_cell_first <- all_presence_events[, .(
  first_month_int = min(month_int),
  last_month_int  = max(month_int)
), by = .(pg_id, actor1)]
log_msg("Group-cell pairs with any activity: ", nrow(group_cell_first))

# ── 6) CONSTRUCT DYADIC PANEL: CARTESIAN OF CELL-MONTHS × GROUPS-IN-THEATRE ──
# For each cell-month, get list of groups that operate in that theatre
cell_month_ids <- base_panel[, .(pg_id, theatre_id, month_int)]
cell_month_ids <- merge(cell_month_ids, group_theatre[, .(actor1, theatre_id)],
                        by = "theatre_id", allow.cartesian = TRUE)
cell_month_ids[, theatre_id := NULL]
log_msg("Candidate dyad-months (before filtering): ", nrow(cell_month_ids))

# Now add presence flag: 1 if group had event in cell in [t-5, t-0]
setkey(all_presence_events, pg_id, actor1, month_int)

# Create a sparse data structure of (pg_id, actor1, month_int) events
# For each candidate (pg_id, actor1, month_int) dyad, check if any event
# in [month_int - 5, month_int]
log_msg("Computing presence flag...")

# Efficient join: for each event, expand to the next 6 months
pres_expand <- all_presence_events[, .(
  pg_id, actor1,
  month_int = seq(month_int, month_int + 5L)
), by = .I][, !"I"]
pres_expand <- unique(pres_expand)
pres_expand[, group_present := 1L]

# Merge presence flag into candidate dyads
dyad <- merge(cell_month_ids, pres_expand,
              by = c("pg_id", "actor1", "month_int"),
              all.x = TRUE)
dyad[is.na(group_present), group_present := 0L]
log_msg("Dyads with group_present == 1: ", sum(dyad$group_present))

# ── 7) ATTACH OUTCOMES (VAC events per group × cell × month) ──────────
vac_events <- rebel_events_t[grepl("Violence against civilians", event_type,
                                    ignore.case = TRUE)]
vac_events[, fat_class := fcase(
  fatalities <= 2, "selective",
  fatalities <= 4, "mid",
  fatalities >= 5, "indiscriminate",
  default         = NA_character_
)]

dyad_vac <- vac_events[, .(
  vac_events      = .N,
  vac_selective   = sum(fat_class == "selective",      na.rm = TRUE),
  vac_mid         = sum(fat_class == "mid",            na.rm = TRUE),
  vac_indisc      = sum(fat_class == "indiscriminate", na.rm = TRUE),
  vac_fatalities  = sum(fatalities, na.rm = TRUE)
), by = .(pg_id, actor1, month_int)]

dyad <- merge(dyad, dyad_vac,
              by = c("pg_id", "actor1", "month_int"), all.x = TRUE)
for (cc in c("vac_events", "vac_selective", "vac_mid",
             "vac_indisc", "vac_fatalities")) {
  dyad[is.na(get(cc)), (cc) := 0L]
}

# ── 8) MERGE CELL-LEVEL VARIABLES (eli_feas, shocks, controls) ────────
# Use month_int as the join key (not month, since that's Date)
cell_slim_join <- cell_slim[, month_int := year(month) * 12L + month(month)][
  , c("pg_id", "month_int", setdiff(names(cell_slim), c("pg_id", "month", "month_int"))),
  with = FALSE
]
dyad <- merge(dyad, cell_slim_join,
              by = c("pg_id", "month_int"),
              all.x = TRUE)

# Some dyads will have no cell data (cell outside panel) — drop them
dyad <- dyad[!is.na(eli_feas)]
log_msg("Final dyadic panel: ", nrow(dyad), " dyad-months, ",
        uniqueN(dyad$pg_id), " cells, ",
        uniqueN(dyad$actor1), " groups, ",
        uniqueN(dyad$month_int), " months")

# ── 9) ADD GROUP_HAS_DRONES VARIABLE ──────────────────────────────────
# Hand-coded acquisition dates (from R23)
drone_acq <- rbindlist(list(
  data.table(actor_pattern = "Taliban",
             theatre = "Afghanistan", acq_date = as.IDate("2019-07-01")),
  data.table(actor_pattern = "Islamic State|ISIS|ISIL|Daesh|Wilayat",
             theatre = "Iraq", acq_date = as.IDate("2016-07-01")),
  data.table(actor_pattern = "Islamic State|ISIS|ISIL|Daesh|Wilayat",
             theatre = "Syria", acq_date = as.IDate("2016-07-01")),
  data.table(actor_pattern = "Kata.ib Hezbollah|Popular Mobilization|PMF|Hashd",
             theatre = "Iraq", acq_date = as.IDate("2018-01-01")),
  data.table(actor_pattern = "Libyan National Army|LNA|Haftar",
             theatre = "Libya", acq_date = as.IDate("2019-06-01")),
  data.table(actor_pattern = "Government of National Accord|GNA",
             theatre = "Libya", acq_date = as.IDate("2019-10-01")),
  data.table(actor_pattern = "Hayat Tahrir|HTS|Jabhat|Nusra",
             theatre = "Syria", acq_date = as.IDate("2018-06-01")),
  data.table(actor_pattern = "DPR|LPR|Donetsk|Luhansk",
             theatre = "Ukraine", acq_date = as.IDate("2014-04-01")),
  data.table(actor_pattern = "Houthi|Ansar Allah|Ansarallah",
             theatre = "Yemen (North Yemen)", acq_date = as.IDate("2017-01-01"))
))
drone_acq[, acq_month_int := year(acq_date) * 12L + month(acq_date)]

# Need theatre for each dyad
dyad_theatre <- merge(dyad[, .(pg_id, actor1, month_int)], cell_theatre,
                       by = "pg_id", all.x = TRUE)
dyad[, theatre_id := dyad_theatre$theatre_id]

# Match each dyad to its acquisition date based on actor1 + theatre_id
dyad[, group_has_drones := 0L]
for (i in seq_len(nrow(drone_acq))) {
  pat <- drone_acq$actor_pattern[i]
  thr <- drone_acq$theatre[i]
  amo <- drone_acq$acq_month_int[i]
  dyad[grepl(pat, actor1, ignore.case = TRUE) &
       theatre_id == thr &
       month_int >= amo,
       group_has_drones := 1L]
}
log_msg("Dyad-months with group_has_drones == 1: ",
        sum(dyad$group_has_drones))

# ── 10) DERIVED OUTCOMES (LOG-TRANSFORMS) ─────────────────────────────
dyad[, log1p_vac        := log1p(vac_events)]
dyad[, log1p_vac_sel    := log1p(vac_selective)]
dyad[, log1p_vac_ind    := log1p(vac_indisc)]

# ── 11) SAVE ──────────────────────────────────────────────────────────
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(dyad, OUT_PATH)
log_msg("Saved dyadic panel to: ", OUT_PATH,
        " (", round(file.size(OUT_PATH) / 1e6, 1), " MB)")

# Summary diagnostics
cat("\n== DYADIC PANEL SUMMARY ==\n")
cat(sprintf("Rows: %d\n", nrow(dyad)))
cat(sprintf("Unique cells: %d\n", uniqueN(dyad$pg_id)))
cat(sprintf("Unique groups: %d\n", uniqueN(dyad$actor1)))
cat(sprintf("Unique months: %d\n", uniqueN(dyad$month_int)))
cat(sprintf("Share with group_present == 1: %.3f\n",
            mean(dyad$group_present)))
cat(sprintf("Share with group_has_drones == 1: %.3f\n",
            mean(dyad$group_has_drones)))
cat(sprintf("Non-zero VAC dyads: %d (%.2f%%)\n",
            sum(dyad$vac_events > 0),
            100 * mean(dyad$vac_events > 0)))
cat(sprintf("Non-zero SELECTIVE VAC dyads: %d\n", sum(dyad$vac_selective > 0)))
cat(sprintf("Non-zero INDISCRIMINATE VAC dyads: %d\n", sum(dyad$vac_indisc > 0)))

log_msg("R34 complete.")
