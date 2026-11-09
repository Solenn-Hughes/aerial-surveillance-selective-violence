# ============================================================
# R08_define_post_eli_v2.R — theatre-specific onsets (no input CSV)
# Outputs *_posteli_v2.rds + eli_onsets_v2.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found: R00_config.R or R/00_config.R")
source(cfg_path)

P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1.rds")
GED_CLEAN <- file.path(DIRS$base, "ged_clean.rds")

for (pp in c(P_ALLCW, P_STRICT, GED_CLEAN)) if (!file.exists(pp)) stop("Missing file: ", pp)

p1 <- readRDS(P_ALLCW);  setDT(p1)
p2 <- readRDS(P_STRICT); setDT(p2)

# ---- standardize keys
std_keys <- function(dt, label){
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  assert_cols(dt, c("pg_id","month"), label)
  dt[, pg_id := as.integer(pg_id)]
  dt[, month := as_month_idate(month)]
  assert_unique(dt, c("pg_id","month"), paste0(label," keys"))
  dt
}
p1 <- std_keys(p1, "ALLCW_paper1")
p2 <- std_keys(p2, "STRICT_paper1")

# ---- ensure we have a theatre identifier (prefer country name; else country_id)
ensure_country_name <- function(dt){
  if ("country" %in% names(dt)) return(dt)
  if (!("country_id" %in% names(dt))) return(dt)
  
  ged <- readRDS(GED_CLEAN); setDT(ged)
  if (!("country_id" %in% names(ged))) return(dt)
  # GED often has a 'country' string column
  cname_col <- intersect(c("country","country_name"), names(ged))
  if (!length(cname_col)) return(dt)
  
  map <- unique(ged[, .(country_id = as.integer(country_id), country = as.character(get(cname_col[1])))])
  setkey(map, country_id)
  dt[, country_id := as.integer(country_id)]
  setkey(dt, country_id)
  dt <- map[dt]  # left join adds country
  dt
}

p1 <- ensure_country_name(p1)
p2 <- ensure_country_name(p2)

theatre_col <- NULL
if ("country" %in% names(p1)) theatre_col <- "country"
if (is.null(theatre_col) && "country_id" %in% names(p1)) theatre_col <- "country_id"
if (is.null(theatre_col)) {
  stop("Cannot define theatre_id: no 'country' or 'country_id' in panel. Run R98_inspect_theatre_id.R and paste output.")
}

# ---- Define onsets internally (EDIT LATER WITH SOURCED JUSTIFICATION)
# These are INITIAL placeholders to remove global-onset bias.
# You can shrink/expand this list as you validate the design.
eli_onsets_map <- c(
  "Afghanistan" = "2001-10-01",
  "Pakistan"    = "2004-06-01",
  "Iraq"        = "2014-08-01",
  "Syria"       = "2014-09-01",
  "Yemen (North Yemen)" = "2009-12-01",
  "Somalia"     = "2007-01-01",
  "Libya"       = "2011-03-01",
  "Ukraine"     = "2022-02-01"
)

# If theatre_id is numeric country_id, we can only treat by name if we have 'country'
if (theatre_col == "country_id" && !("country" %in% names(p1))) {
  stop("Panel only has country_id (numeric) and no country name. ",
       "Run R98_inspect_theatre_id.R and we will treat by IDs or merge names.")
}

# ---- Apply onsets
apply_onsets <- function(dt, label){
  if (theatre_col == "country") {
    dt[, theatre_id := as.character(country)]
    dt[, eli_onset_month := as.IDate(NA)]
    dt[theatre_id %in% names(eli_onsets_map), eli_onset_month := as.IDate(eli_onsets_map[theatre_id])]
  } else {
    # theatre_col == "country_id" but we also have dt$country from join above
    dt[, theatre_id := as.character(country)]  # keep stable readable label
    dt[, eli_onset_month := as.IDate(NA)]
    dt[theatre_id %in% names(eli_onsets_map), eli_onset_month := as.IDate(eli_onsets_map[theatre_id])]
  }
  
  dt[, treated := as.integer(!is.na(eli_onset_month))]
  dt[, post_eli := as.integer(treated == 1L & month >= eli_onset_month)]
  
  # months since onset (only treated)
  dt[, event_time := NA_integer_]
  dt[treated == 1L, event_time := (year(month) - year(eli_onset_month)) * 12L + (month(month) - month(eli_onset_month))]
  
  # Integrity
  assert_no_na(dt, c("treated","post_eli"), label)
  if (dt[, any(post_eli == 1L & treated == 0L)]) stop(label, ": post_eli==1 for untreated (impossible).")
  dt
}

p1v2 <- apply_onsets(copy(p1), "ALLCW_posteli_v2")
p2v2 <- apply_onsets(copy(p2), "STRICT_posteli_v2")

# ---- Audit onsets
onsets_dt <- unique(p1v2[treated == 1L, .(theatre_id, eli_onset_month)])
setorder(onsets_dt, theatre_id)
onsets_dt[, eli_onset_month := as.character(eli_onset_month)]

audit_summary <- rbindlist(list(
  data.table(panel="ALLCW",  rows=nrow(p1v2), theatres=uniqueN(p1v2$theatre_id), treated_theatres=uniqueN(p1v2[treated==1L]$theatre_id),
             share_post=mean(p1v2$post_eli==1L), et_min=min(p1v2$event_time, na.rm=TRUE), et_max=max(p1v2$event_time, na.rm=TRUE)),
  data.table(panel="STRICT", rows=nrow(p2v2), theatres=uniqueN(p2v2$theatre_id), treated_theatres=uniqueN(p2v2[treated==1L]$theatre_id),
             share_post=mean(p2v2$post_eli==1L), et_min=min(p2v2$event_time, na.rm=TRUE), et_max=max(p2v2$event_time, na.rm=TRUE))
), fill=TRUE)

write_audit_dt(audit_summary, "R08_define_post_eli_v2_summary")
write_audit_dt(onsets_dt, "R08_define_post_eli_v2_onsets")

# ============================================================
# ADDITION TO R08_define_post_eli_v2.R
# Insert this block BEFORE the "# ---- Save outputs" section
# (i.e., before line ~127 "OUT1 <- file.path(...)")
# ============================================================

# ---- Compute distance to nearest Somalia FOB (for withdrawal shock) ----
# These are the forward operating locations US troops were withdrawn FROM
# in Dec 2020 – Jan 2021. Used by R16 for interior-to-interior comparison.

SOMALIA_FOBS <- data.table(
  fob_name = c("Mogadishu_Airport", "Baledogle_Airbase"),
  fob_lat  = c(2.0144, 2.367),
  fob_lon  = c(45.3047, 44.767)
)

haversine_km_fob <- function(lat1, lon1, lat2, lon2) {
  R <- 6371.0; to_rad <- pi / 180
  lat1r <- lat1 * to_rad; lon1r <- lon1 * to_rad
  lat2r <- lat2 * to_rad; lon2r <- lon2 * to_rad
  dlat <- lat2r - lat1r; dlon <- lon2r - lon1r
  a <- sin(dlat/2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

add_fob_distances <- function(dt) {
  # Need cell coordinates — get from eli_feasibility or panel
  eli_path <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
  if (file.exists(eli_path)) {
    eli <- readRDS(eli_path); setDT(eli)
    coords <- eli[, .(pg_id = as.integer(gid), cell_lon = lon, cell_lat = lat)]
  } else {
    # Fallback: derive from PRIO-GRID formula
    dt[, cell_lat := ((pg_id - 1L) %/% 720L) * 0.5 - 90 + 0.25]
    dt[, cell_lon := ((pg_id - 1L) %% 720L) * 0.5 - 180 + 0.25]
    coords <- unique(dt[, .(pg_id, cell_lon, cell_lat)])
  }
  
  # Compute min distance to any Somalia FOB for each cell
  coords[, dist_somalia_fob_km := {
    dists <- vapply(seq_len(nrow(SOMALIA_FOBS)), function(i) {
      haversine_km_fob(cell_lat, cell_lon, SOMALIA_FOBS$fob_lat[i], SOMALIA_FOBS$fob_lon[i])
    }, numeric(.N))
    if (is.matrix(dists)) apply(dists, 1, min) else dists
  }]
  
  # Also compute distance to Djibouti (Camp Lemonnier — fallback base post-withdrawal)
  DJIBOUTI_LAT <- 11.5469; DJIBOUTI_LON <- 43.1556
  coords[, dist_djibouti_km := haversine_km_fob(cell_lat, cell_lon, DJIBOUTI_LAT, DJIBOUTI_LON)]
  
  # Also compute distance to Bagram (for Bagram shock)
  BAGRAM_LAT <- 34.9461; BAGRAM_LON <- 69.2650
  coords[, dist_bagram_km := haversine_km_fob(cell_lat, cell_lon, BAGRAM_LAT, BAGRAM_LON)]
  
  # Also compute distance to Incirlik (for Incirlik shock)
  INCIRLIK_LAT <- 37.0021; INCIRLIK_LON <- 35.4259
  coords[, dist_incirlik_km := haversine_km_fob(cell_lat, cell_lon, INCIRLIK_LAT, INCIRLIK_LON)]
  
  # And distance to Muwaffaq Salti (Jordan — replacement base post-Incirlik)
  MUWAFFAQ_LAT <- 32.3564; MUWAFFAQ_LON <- 36.7822
  coords[, dist_muwaffaq_km := haversine_km_fob(cell_lat, cell_lon, MUWAFFAQ_LAT, MUWAFFAQ_LON)]
  
  # Merge into panel
  merge_cols <- c("pg_id", "dist_somalia_fob_km", "dist_djibouti_km", 
                  "dist_bagram_km", "dist_incirlik_km", "dist_muwaffaq_km")
  coords_merge <- coords[, ..merge_cols]
  
  # Remove if already present (re-run safety)
  for (cc in setdiff(merge_cols, "pg_id")) {
    if (cc %in% names(dt)) dt[, (cc) := NULL]
  }
  
  dt <- merge(dt, coords_merge, by = "pg_id", all.x = TRUE)
  log_msg("Added FOB/base distances. Sample: Somalia FOB dist range = ",
          round(min(dt$dist_somalia_fob_km, na.rm = TRUE)), " - ",
          round(max(dt$dist_somalia_fob_km, na.rm = TRUE)), " km")
  dt
}

p1v2 <- add_fob_distances(p1v2)
p2v2 <- add_fob_distances(p2v2)

# Audit distances
dist_audit <- p1v2[, .(
  n = .N, cells = uniqueN(pg_id),
  bagram_med = median(dist_bagram_km, na.rm = TRUE),
  somalia_fob_med = median(dist_somalia_fob_km, na.rm = TRUE),
  djibouti_med = median(dist_djibouti_km, na.rm = TRUE),
  incirlik_med = median(dist_incirlik_km, na.rm = TRUE),
  muwaffaq_med = median(dist_muwaffaq_km, na.rm = TRUE)
)]
write_audit_dt(dist_audit, "R08_fob_distances_summary")
log_msg("FOB distances added to both panels.")

# ---- Save outputs (new names, do not overwrite v1)
OUT1 <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
OUT2 <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds")
saveRDS(p1v2, OUT1)
saveRDS(p2v2, OUT2)
log_msg("Wrote: ", OUT1)
log_msg("Wrote: ", OUT2)

# also write a convenience CSV (output only; not required as input)
OUTCSV <- file.path(DIRS$eli, "eli_onsets_v2.csv")
fwrite(onsets_dt, OUTCSV)
log_msg("Wrote: ", OUTCSV)

log_msg("Step 8 v2 complete.")