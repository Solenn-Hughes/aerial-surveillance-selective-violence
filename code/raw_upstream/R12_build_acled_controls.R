# ============================================================
# R12_build_acled_controls.R  (v2 — full codebook awareness)
#
# Ingest ACLED CSV, clean, and build PRIO-GRID × month
# conflict indicators for H2 (contestation) and H3 (rebel cap)
#
# ACLED CODEBOOK REFERENCE (https://acleddata.com/acleddatanew/
#   wp-content/uploads/dlm_uploads/2019/01/ACLED_Codebook_2019FINAL.pdf)
#
# Actor type codes (inter1/inter2):
#   1 = State Forces          5 = Rioters
#   2 = Rebel Groups          6 = Protesters
#   3 = Political Militias    7 = Civilians
#   4 = Identity Militias     8 = External/Other Forces
#
# Interaction codes (two-digit: first=inter1, second=inter2):
#   12/21 = State vs Rebel    17 = State vs Civilian
#   13/31 = State vs Militia  27 = Rebel vs Civilian
#   23/32 = Rebel vs Militia  37 = Militia vs Civilian
#
# Event types:
#   Battles, Violence against civilians, Explosions/Remote violence,
#   Protests, Riots, Strategic developments
#
# Key sub-event types:
#   Battles: Armed clash, Government regains territory,
#            Non-state actor overtakes territory
#   Remote:  Air/drone strike, Shelling/artillery/missile attack,
#            Remote explosive/landmine/IED, Suicide bomb, Grenade
#   VAC:     Sexual violence, Attack, Abduction/forced disappearance
#
# Input:  data/raw/acled/acled_full.csv
# Output: data/processed/base/acled_pg_month.rds
#
# Run AFTER downloading ACLED, BEFORE R13 and R14.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 12: Ingest ACLED and build pg-month conflict indicators (v2)")

# ── LOCATE ACLED FILE ────────────────────────────────────────
ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT, "data", "raw", "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT, "data", "raw", "acled_full.csv")
)

acled_path <- NULL
for (pp in ACLED_PATHS) {
  if (file.exists(pp)) { acled_path <- pp; break }
}

if (is.null(acled_path)) {
  stop("ACLED data not found.\n",
       "Expected at: ", ACLED_PATHS[1], "\n",
       "See DOWNLOAD_AND_RUN_GUIDE.R section D for download instructions.\n")
}

log_msg("Reading ACLED from: ", acled_path)

# ── READ AND CLEAN ───────────────────────────────────────────
ac <- fread(acled_path, showProgress = TRUE, encoding = "UTF-8")
log_msg("Raw ACLED rows: ", nrow(ac))

# Standardise column names (ACLED exports vary in casing)
setnames(ac, tolower(names(ac)))

# ── COLUMN INVENTORY ─────────────────────────────────────────
# Log which columns are present for debugging
all_cols <- names(ac)
log_msg("ACLED columns found: ", paste(all_cols, collapse = ", "))

# Required minimum
req_min <- c("event_date", "event_type", "latitude", "longitude")
miss <- setdiff(req_min, all_cols)
if (length(miss)) stop("ACLED missing critical columns: ", paste(miss, collapse = ", "))

# ── PARSE DATES ──────────────────────────────────────────────
# ACLED date format can be "DD Month YYYY" or "YYYY-MM-DD"
ac[, event_date_raw := event_date]

# Try standard ISO first
ac[, event_date := tryCatch(as.IDate(event_date_raw), error = function(e) NA)]

# If mostly NA, try the "DD Month YYYY" format
if (mean(is.na(ac$event_date)) > 0.5) {
  log_msg("Standard date parse failed. Trying ACLED format (e.g. '15 January 2020')...")
  ac[, event_date := as.IDate(event_date_raw, format = "%d %B %Y")]
}
if (mean(is.na(ac$event_date)) > 0.5) {
  # Also try MM/DD/YYYY
  ac[, event_date := as.IDate(event_date_raw, format = "%m/%d/%Y")]
}

n_bad_date <- sum(is.na(ac$event_date))
if (n_bad_date > 0) log_msg("WARNING: ", n_bad_date, " rows with unparseable dates")
ac <- ac[!is.na(event_date)]

# Year/month
if ("year" %in% names(ac)) {
  ac[, year := as.integer(year)]
} else {
  ac[, year := year(event_date)]
}
ac[, month := as_month_idate(event_date)]

# ── PARSE ACTOR TYPE CODES ───────────────────────────────────
# ACLED's inter1/inter2 can be NUMERIC codes (1-8) if user
# ticked "display interactions as codes", or TEXT labels
# ("State Forces", "Rebel Groups", etc.) if they didn't.
# We handle both cases.

# Mapping from text labels to numeric codes
actor_type_map <- c(
  "State Forces"        = 1L,
  "Rebel Groups"        = 2L,
  "Political Militias"  = 3L,
  "Identity Militias"   = 4L,
  "Rioters"             = 5L,
  "Protesters"          = 6L,
  "Civilians"           = 7L,
  "External/Other Forces" = 8L
)

parse_inter <- function(x) {
  # Try numeric first
  num <- suppressWarnings(as.integer(x))
  if (mean(!is.na(num), na.rm = TRUE) > 0.8) return(num)
  
  # Text labels — map to codes
  log_msg("  inter column appears to be text labels. Mapping to codes...")
  out <- actor_type_map[as.character(x)]
  out <- as.integer(out)
  # Log unmapped values
  unmapped <- unique(x[is.na(out) & !is.na(x)])
  if (length(unmapped) > 0 && length(unmapped) < 20) {
    log_msg("  Unmapped actor types: ", paste(unmapped, collapse = "; "))
  }
  out
}

if ("inter1" %in% names(ac)) {
  ac[, inter1 := parse_inter(inter1)]
  log_msg("inter1 parsed. Distribution:")
  print(ac[, .N, by = inter1][order(inter1)])
} else {
  log_msg("WARNING: inter1 not found. Actor-type classification will be limited.")
  ac[, inter1 := NA_integer_]
}

if ("inter2" %in% names(ac)) {
  ac[, inter2 := parse_inter(inter2)]
} else {
  ac[, inter2 := NA_integer_]
}

# Interaction code (combined two-digit)
if ("interaction" %in% names(ac)) {
  ac[, interaction := as.integer(interaction)]
} else if (!any(is.na(ac$inter1)) && !any(is.na(ac$inter2))) {
  # Reconstruct from inter1 and inter2
  ac[, interaction := inter1 * 10L + inter2]
}

# ── PARSE OTHER COLUMNS ──────────────────────────────────────
ac[, latitude  := as.numeric(latitude)]
ac[, longitude := as.numeric(longitude)]

if ("fatalities" %in% names(ac)) {
  ac[, fatalities := as.integer(fatalities)]
} else {
  ac[, fatalities := 0L]
}

if ("geo_precision" %in% names(ac)) {
  ac[, geo_precision := as.integer(geo_precision)]
} else {
  ac[, geo_precision := NA_integer_]
}

if (!("sub_event_type" %in% names(ac))) {
  ac[, sub_event_type := NA_character_]
}

# ── DROP BAD COORDINATES ─────────────────────────────────────
n_miss_coord <- sum(is.na(ac$latitude) | is.na(ac$longitude))
log_msg("Dropping ", n_miss_coord, " rows with missing coordinates")
ac <- ac[!is.na(latitude) & !is.na(longitude)]

# ── CLASSIFY EVENTS (full codebook) ──────────────────────────
# Event types
ac[, is_battle   := grepl("^Battles", event_type, ignore.case = TRUE)]
ac[, is_vac      := grepl("Violence against civilians", event_type, ignore.case = TRUE)]
ac[, is_remote   := grepl("Explosions|Remote violence", event_type, ignore.case = TRUE)]
ac[, is_riot     := grepl("^Riots", event_type, ignore.case = TRUE)]
ac[, is_protest  := grepl("^Protests", event_type, ignore.case = TRUE)]
ac[, is_strategic := grepl("Strategic developments", event_type, ignore.case = TRUE)]

# Sub-event types (critical for H2 and mechanism)
ac[, is_armed_clash       := grepl("Armed clash", sub_event_type, ignore.case = TRUE)]
ac[, is_territory_change  := grepl("regains territory|overtakes territory",
                                    sub_event_type, ignore.case = TRUE)]
ac[, is_govt_regains      := grepl("Government regains territory",
                                    sub_event_type, ignore.case = TRUE)]
ac[, is_rebel_overtakes   := grepl("Non-state actor overtakes territory",
                                    sub_event_type, ignore.case = TRUE)]
ac[, is_airdrone_strike   := grepl("Air/drone strike", sub_event_type, ignore.case = TRUE)]
ac[, is_shelling          := grepl("Shelling|artillery|missile",
                                    sub_event_type, ignore.case = TRUE)]
ac[, is_remote_explosive  := grepl("Remote explosive|landmine|IED",
                                    sub_event_type, ignore.case = TRUE)]
ac[, is_suicide_bomb      := grepl("Suicide bomb", sub_event_type, ignore.case = TRUE)]
ac[, is_abduction         := grepl("Abduction|forced disappearance",
                                    sub_event_type, ignore.case = TRUE)]

# Actor-type classifications
ac[, state_actor1   := inter1 == 1L]
ac[, rebel_actor1   := inter1 == 2L]
ac[, militia_actor1 := inter1 %in% c(3L, 4L)]

ac[, state_actor2   := inter2 == 1L]
ac[, rebel_actor2   := inter2 == 2L]
ac[, militia_actor2 := inter2 %in% c(3L, 4L)]

# Specific interaction patterns (using interaction code or inter1/inter2)
# State vs Rebel battles (core for H2 contestation)
ac[, is_state_rebel_battle := is_battle & (
  (inter1 == 1L & inter2 == 2L) | (inter1 == 2L & inter2 == 1L)
)]

# State vs Militia battles
ac[, is_state_militia_battle := is_battle & (
  (inter1 == 1L & inter2 %in% c(3L, 4L)) | (inter1 %in% c(3L, 4L) & inter2 == 1L)
)]

# Rebel vs Rebel / inter-rebel fighting (indicates fragmentation)
ac[, is_rebel_rebel := is_battle & inter1 == 2L & inter2 == 2L]

# Rebel-initiated events (any event where rebels are actor 1)
ac[, is_rebel_initiated := rebel_actor1 & (is_battle | is_vac | is_remote)]

# State-initiated violence against civilians (targeting)
ac[, is_state_vac := is_vac & state_actor1]

# Rebel-initiated violence against civilians
ac[, is_rebel_vac := is_vac & rebel_actor1]

# All violence involving non-state armed groups
ac[, is_nonstate_violence := (rebel_actor1 | militia_actor1) & (is_battle | is_vac | is_remote)]

# ── GEO-PRECISION FLAG ───────────────────────────────────────
# geo_precision: 1 = exact location, 2 = near, 3 = admin region centroid
# Flag low-precision events but keep them (dropping would bias)
if (!all(is.na(ac$geo_precision))) {
  geo_tab <- ac[, .N, by = geo_precision][order(geo_precision)]
  log_msg("Geo-precision distribution:")
  print(geo_tab)
  ac[, geo_precise := geo_precision == 1L]
} else {
  ac[, geo_precise := TRUE]  # assume precise if column missing
}

# ── MAP TO PRIO-GRID ────────────────────────────────────────
# PRIO-GRID is 0.5° × 0.5°
# gid = (row - 1) * 720 + col
# row = floor((lat + 90) / 0.5) + 1
# col = floor((lon + 180) / 0.5) + 1

ac[, pg_row := floor((latitude + 90) / 0.5) + 1L]
ac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
ac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# Restrict to cells in our panel
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (file.exists(P_ALLCW)) {
  panel <- readRDS(P_ALLCW); setDT(panel)
  if (!("pg_id" %in% names(panel)) && "gid" %in% names(panel)) panel[, pg_id := gid]
  panel_cells <- unique(as.integer(panel$pg_id))
  panel_months <- sort(unique(panel$month))
  
  n_in <- sum(ac$pg_id %in% panel_cells)
  log_msg("ACLED events in panel cells: ", n_in, " / ", nrow(ac),
          " (", round(100 * n_in / nrow(ac), 1), "%)")
  ac <- ac[pg_id %in% panel_cells]
  log_msg("After spatial filter: ", nrow(ac), " events")
} else {
  log_msg("WARNING: Panel not found. Keeping all ACLED events.")
  panel_cells <- unique(ac$pg_id)
  panel_months <- sort(unique(ac$month))
}

# ── AGGREGATE TO PG × MONTH ─────────────────────────────────
# Count distinct actors using actor1/actor2 columns if available
has_actor1 <- "actor1" %in% names(ac)

agg <- ac[, .(
  # ── Event counts by type ──
  acled_all_events       = .N,
  acled_battles          = sum(is_battle),
  acled_vac              = sum(is_vac),
  acled_remote           = sum(is_remote),
  acled_riots            = sum(is_riot),
  acled_protests         = sum(is_protest),
  acled_strategic        = sum(is_strategic),
  
  # ── Sub-event detail ──
  acled_armed_clash      = sum(is_armed_clash),
  acled_territory_change = sum(is_territory_change),
  acled_govt_regains     = sum(is_govt_regains),
  acled_rebel_overtakes  = sum(is_rebel_overtakes),
  acled_airdrone_strike  = sum(is_airdrone_strike),
  acled_shelling         = sum(is_shelling),
  acled_ied              = sum(is_remote_explosive),
  acled_suicide_bomb     = sum(is_suicide_bomb),
  acled_abduction        = sum(is_abduction),
  
  # ── Actor-specific interactions ──
  acled_state_rebel      = sum(is_state_rebel_battle),
  acled_state_militia    = sum(is_state_militia_battle),
  acled_rebel_rebel      = sum(is_rebel_rebel),
  acled_rebel_initiated  = sum(is_rebel_initiated),
  acled_state_vac        = sum(is_state_vac),
  acled_rebel_vac        = sum(is_rebel_vac),
  acled_nonstate_violence = sum(is_nonstate_violence),
  
  # ── Fatalities ──
  acled_fatalities       = sum(fatalities, na.rm = TRUE),
  acled_battle_fatalities = sum(fatalities[is_battle], na.rm = TRUE),
  acled_vac_fatalities   = sum(fatalities[is_vac], na.rm = TRUE),
  
  # ── Geo-precision share ──
  acled_share_precise    = mean(geo_precise, na.rm = TRUE),
  
  # ── Distinct actors (if actor columns available) ──
  acled_n_distinct_actors = if (has_actor1) {
    uniqueN(c(na.omit(actor1), na.omit(get("actor2"))))
  } else { NA_integer_ }
  
), by = .(pg_id, month)]

# ── FILL ZEROS FOR COMPLETE PANEL ────────────────────────────
grid <- CJ(pg_id = panel_cells, month = panel_months)
agg <- merge(grid, agg, by = c("pg_id", "month"), all.x = TRUE)

# Fill NAs with 0 (no events = 0 events)
num_cols <- setdiff(names(agg), c("pg_id", "month"))
for (cc in num_cols) {
  if (cc == "acled_share_precise") {
    agg[is.na(get(cc)), (cc) := 1]  # assume precise if no events
  } else if (cc == "acled_n_distinct_actors") {
    agg[is.na(get(cc)), (cc) := 0L]
  } else {
    agg[is.na(get(cc)), (cc) := 0L]
  }
}

log_msg("Aggregated ACLED: ", nrow(agg), " cell-months, ",
        uniqueN(agg$pg_id), " cells, ", uniqueN(agg$month), " months")

# ── COMPREHENSIVE AUDIT ──────────────────────────────────────
audit_summary <- data.table(
  metric = c("Total cell-months", "Cells", "Months",
             "Total events", "Total battles", "Total VAC",
             "Total remote/explosive", "Territory changes",
             "Air/drone strikes", "State-rebel battles",
             "Rebel-initiated events", "Total fatalities",
             "Share zero-event months", "Share zero-battle months"),
  value = c(nrow(agg), uniqueN(agg$pg_id), uniqueN(agg$month),
            sum(agg$acled_all_events), sum(agg$acled_battles),
            sum(agg$acled_vac), sum(agg$acled_remote),
            sum(agg$acled_territory_change),
            sum(agg$acled_airdrone_strike),
            sum(agg$acled_state_rebel),
            sum(agg$acled_rebel_initiated),
            sum(agg$acled_fatalities),
            round(mean(agg$acled_all_events == 0), 3),
            round(mean(agg$acled_battles == 0), 3))
)
cat("\n============================================================\n")
cat("ACLED SUMMARY\n")
cat("============================================================\n")
print(audit_summary)
write_audit_dt(audit_summary, "R12_acled_summary")

# Event type breakdown
event_audit <- ac[, .N, by = event_type][order(-N)]
write_audit_dt(event_audit, "R12_acled_event_types")

# Sub-event type breakdown
sub_audit <- ac[, .N, by = .(event_type, sub_event_type)][order(event_type, -N)]
write_audit_dt(sub_audit, "R12_acled_sub_event_types")

# Interaction code breakdown (top 20)
if ("interaction" %in% names(ac)) {
  inter_audit <- ac[, .N, by = interaction][order(-N)][1:min(.N, 20)]
  write_audit_dt(inter_audit, "R12_acled_interaction_codes")
}

# Country breakdown
country_audit <- ac[, .(events = .N, fatalities = sum(fatalities, na.rm = TRUE)),
                    by = country][order(-events)]
write_audit_dt(country_audit, "R12_acled_by_country")

# ── SAVE ─────────────────────────────────────────────────────
out_path <- file.path(DIRS$base, "acled_pg_month.rds")
saveRDS(agg, out_path)
log_msg("Wrote: ", out_path, " (", round(file.size(out_path)/1e6, 1), " MB)")

log_msg("Step 12 complete.")
