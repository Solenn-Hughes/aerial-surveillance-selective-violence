# ============================================================
# R02_ingest_airbases.R — Ingest & clean airbase/airport points
# FIX: reads quoted CSV correctly (no early stop on commas-in-quotes)
#
# Input:
#   data/raw/wikidata/wikidata_airbases_raw.csv  (auto-detected)
#
# Output:
#   data/processed/base/airbases_points.rds
#
# Audits:
#   outputs/paper1/audit/airbases_inputs_*.csv
#   outputs/paper1/audit/airbases_read_summary_*.csv
#   outputs/paper1/audit/airbases_dropped_missing_coords_*.csv (if any)
#   outputs/paper1/audit/airbases_duplicates_removed_*.csv (if any)
# ============================================================

source("R00_config.R")
log_msg("Step 2: ingest airbases/airports (Wikidata)")

# -------------------------
# 1) Locate input
# -------------------------
wikidata_dir <- file.path(PROJECT_ROOT, "data", "raw", "wikidata")
if (!dir.exists(wikidata_dir)) stop("Wikidata raw dir not found: ", wikidata_dir)

cands <- list.files(
  wikidata_dir,
  pattern = "wikidata_airbases.*\\.csv$",
  ignore.case = TRUE,
  full.names = TRUE
)

if (!length(cands)) stop("No wikidata_airbases*.csv found in: ", wikidata_dir)

pick_latest <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  paths[order(file.info(paths)$mtime, decreasing = TRUE)][1]
}

IN_CSV <- pick_latest(cands)
log_msg("Using: ", IN_CSV)

write_audit_dt(
  data.table(path = IN_CSV, mtime = file.info(IN_CSV)$mtime),
  "airbases_inputs"
)

# -------------------------
# 2) Read (IMPORTANT: allow normal quote parsing)
# -------------------------
dt <- fread(
  IN_CSV,
  showProgress = FALSE,
  fill = TRUE,          # tolerate occasional ragged lines
  encoding = "UTF-8"    # keep Å/é/etc
  # NOTE: do NOT set quote="" here; we WANT quote parsing
)

# quick read audit
write_audit_dt(
  data.table(n_rows_raw = nrow(dt), n_cols_raw = ncol(dt), cols = paste(names(dt), collapse = " | ")),
  "airbases_read_summary"
)

# -------------------------
# 3) Standardize columns
# Expected: item | itemLabel | countryLabel | lat | lon | enwiki
# -------------------------
rename_if_present <- function(old, new) {
  if (old %in% names(dt) && !(new %in% names(dt))) setnames(dt, old, new)
}

rename_if_present("itemLabel", "asset_name")
rename_if_present("enwiki", "enwiki_url")
rename_if_present("item", "qid")

# tolerate alternate column names if present
rename_if_present("latitude", "lat")
rename_if_present("longitude", "lon")

need <- c("qid", "asset_name", "lat", "lon")
miss <- setdiff(need, names(dt))
if (length(miss)) stop("Missing required cols in airbases raw: ", paste(miss, collapse = ", "))

# -------------------------
# 4) Clean types
# -------------------------
dt[, asset_name := trimws(as.character(asset_name))]
dt[, qid := trimws(as.character(qid))]
if ("countryLabel" %in% names(dt)) dt[, countryLabel := trimws(as.character(countryLabel))]
if ("enwiki_url" %in% names(dt)) dt[, enwiki_url := trimws(as.character(enwiki_url))]

dt[, lat := suppressWarnings(as.numeric(lat))]
dt[, lon := suppressWarnings(as.numeric(lon))]

# -------------------------
# 5) Drop missing/invalid coords (audit)
# -------------------------
bad <- dt[is.na(lat) | is.na(lon) | lat < -90 | lat > 90 | lon < -180 | lon > 180]
if (nrow(bad)) {
  write_audit_dt(bad[1:min(200, .N), .(qid, asset_name, lat, lon)], "airbases_dropped_missing_coords")
  dt <- dt[!(is.na(lat) | is.na(lon) | lat < -90 | lat > 90 | lon < -180 | lon > 180)]
}
log_msg("After dropping missing/invalid coords: n=", nrow(dt))

# -------------------------
# 6) Stable IDs + dedupe
# -------------------------
# Extract QID if qid is a full URL; otherwise keep as is.
dt[, qid_clean := fifelse(grepl("Q\\d+", qid), sub(".*(Q\\d+).*", "\\1", qid), qid)]
dt[, base_id := qid_clean]

# if any missing base_id, fallback to name+coords signature
dt[is.na(base_id) | !nzchar(base_id),
   base_id := sprintf("name_%s_%0.5f_%0.5f",
                      gsub("[^A-Za-z0-9]+","_", tolower(asset_name)),
                      lat, lon)]

# primary dedupe: base_id
setorder(dt, base_id)
dups1 <- dt[duplicated(base_id)]
if (nrow(dups1)) {
  write_audit_dt(dups1[1:min(200, .N), .(base_id, qid, asset_name, lat, lon)], "airbases_duplicates_removed")
  dt <- dt[!duplicated(base_id)]
}

# secondary dedupe: same name + same rounded coords (catches repeat entries)
dt[, lat_r := round(lat, 4)]
dt[, lon_r := round(lon, 4)]
setorder(dt, asset_name, lat_r, lon_r)
dups2 <- dt[duplicated(dt[, .(asset_name, lat_r, lon_r)])]
if (nrow(dups2)) {
  write_audit_dt(dups2[1:min(200, .N), .(base_id, qid, asset_name, lat, lon)], "airbases_duplicates_removed_secondary")
  dt <- dt[!duplicated(dt[, .(asset_name, lat_r, lon_r)])]
}
dt[, c("lat_r","lon_r") := NULL]

# -------------------------
# 7) Final schema + save
# -------------------------
dt[, source := "wikidata"]

keep <- c("base_id", "asset_name", "lat", "lon", "source", "qid_clean")
if ("countryLabel" %in% names(dt)) keep <- c(keep, "countryLabel")
if ("enwiki_url" %in% names(dt)) keep <- c(keep, "enwiki_url")

dt_out <- dt[, ..keep]
setnames(dt_out, "qid_clean", "qid")

assert_unique(dt_out, "base_id", "airbases_points")

# summary audit
summ <- data.table(
  n_points = nrow(dt_out),
  n_unique_qid = if ("qid" %in% names(dt_out)) uniqueN(dt_out$qid) else NA_integer_,
  lat_min = min(dt_out$lat), lat_max = max(dt_out$lat),
  lon_min = min(dt_out$lon), lon_max = max(dt_out$lon)
)
write_audit_dt(summ, "airbases_clean_summary")

out_path <- file.path(DIRS$base, "airbases_points.rds")
saveRDS(dt_out, out_path)
log_msg("Wrote: ", out_path)
log_msg("Step 2 complete.")
