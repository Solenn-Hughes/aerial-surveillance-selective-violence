# ============================================================
# R03_ingest_ged.R — Ingest UCDP GED and aggregate to PRIO grid×month
#
# Input:
#   data/raw/ucdp_ged/v25_1/GEDEvent_v25_1.csv
#
# Outputs:
#   data/processed/base/ged_clean.rds
#   data/processed/base/ged_pg_month.rds
#
# Audits:
#   outputs/paper1/audit/ged_inputs_*.csv
#   outputs/paper1/audit/ged_clean_summary_*.csv
#   outputs/paper1/audit/ged_missing_top30_*.csv
#   outputs/paper1/audit/ged_type_counts_*.csv
#   outputs/paper1/audit/ged_pg_month_summary_*.csv
# ============================================================

source("R00_config.R")
log_msg("Step 3: ingest GED and build PRIO grid×month aggregates")

IN_CSV <- file.path(PROJECT_ROOT, "data", "raw", "ucdp_ged", "v25_1", "GEDEvent_v25_1.csv")
if (!file.exists(IN_CSV)) stop("GED file not found: ", IN_CSV)

write_audit_dt(
  data.table(path = IN_CSV, mtime = file.info(IN_CSV)$mtime),
  "ged_inputs"
)

# -------------------------
# 1) Read
# -------------------------
ged <- fread(IN_CSV, showProgress = TRUE, fill = TRUE, encoding = "UTF-8")

# -------------------------
# 2) Require columns we rely on
# -------------------------
need <- c(
  "id","type_of_violence","priogrid_gid","latitude","longitude",
  "date_start","date_end",
  "best","high","low",
  "deaths_civilians","deaths_unknown",
  "country","country_id","region","year"
)
miss <- setdiff(need, names(ged))
if (length(miss)) stop("GED missing required columns: ", paste(miss, collapse = ", "))

# -------------------------
# 3) Type coercion (robust)
# -------------------------
to_int <- function(x) suppressWarnings(as.integer(x))
to_num <- function(x) suppressWarnings(as.numeric(x))

ged[, id := as.character(id)]

ged[, type_of_violence := to_int(type_of_violence)]
ged[, priogrid_gid     := to_int(priogrid_gid)]
ged[, year             := to_int(year)]

ged[, latitude  := to_num(latitude)]
ged[, longitude := to_num(longitude)]

# dates are YYYY-MM-DD strings; parse to IDate
ged[, date_start := as.IDate(date_start)]
ged[, date_end   := as.IDate(date_end)]

ged[, best := to_int(best)]
ged[, high := to_int(high)]
ged[, low  := to_int(low)]

ged[, deaths_civilians := to_int(deaths_civilians)]
ged[, deaths_unknown   := to_int(deaths_unknown)]

# -------------------------
# 4) Basic integrity checks + audits
# -------------------------
# date_start should exist for almost all events
n_bad_date <- ged[is.na(date_start), .N]
if (n_bad_date > 0) {
  write_audit_dt(ged[is.na(date_start)][1:min(200, .N),
                                        .(id, year, type_of_violence, date_start, date_end, country, priogrid_gid)],
                 "ged_missing_date_start_examples")
  stop("GED has missing date_start for ", n_bad_date, " rows. See audit ged_missing_date_start_examples_*.csv")
}

# if priogrid_gid missing, we can compute later, but for now we require it
n_bad_pg <- ged[is.na(priogrid_gid), .N]
if (n_bad_pg > 0) {
  write_audit_dt(ged[is.na(priogrid_gid)][1:min(200, .N),
                                          .(id, year, type_of_violence, latitude, longitude, priogrid_gid, country)],
                 "ged_missing_priogrid_gid_examples")
  stop("GED has missing priogrid_gid for ", n_bad_pg, " rows. See audit ged_missing_priogrid_gid_examples_*.csv")
}

# quick missingness table
miss_top <- audit_missing(ged, "ged")[1:min(30, .N)]
write_audit_dt(miss_top, "ged_missing_top30")

# type counts
type_counts <- ged[, .N, by = .(type_of_violence)][order(type_of_violence)]
write_audit_dt(type_counts, "ged_type_counts")

# -------------------------
# 5) Create month key (YYYY-MM)
# -------------------------
ged[, month := format(date_start, "%Y-%m")]
# also numeric month_id for FE convenience later if desired
ged[, month_id := as.integer(factor(month, levels = sort(unique(month))))]

# -------------------------
# 6) Build clean GED (kept columns)
# -------------------------
ged_clean <- ged[, .(
  id, year, month, month_id,
  type_of_violence,
  priogrid_gid,
  latitude, longitude,
  country, country_id, region,
  date_start, date_end,
  best, high, low,
  deaths_civilians, deaths_unknown
)]

# ensure no duplicates on id (GED id should be unique)
if (anyDuplicated(ged_clean$id)) {
  dup <- ged_clean[duplicated(id) | duplicated(id, fromLast = TRUE)][1:min(200, .N)]
  write_audit_dt(dup, "ged_duplicate_id_examples")
  stop("GED id is not unique. See audit ged_duplicate_id_examples_*.csv")
}

out_clean <- file.path(DIRS$base, "ged_clean.rds")
saveRDS(ged_clean, out_clean)
log_msg("Wrote: ", out_clean)

# summary audit
summ_clean <- ged_clean[, .(
  n_events = .N,
  years_min = min(year, na.rm=TRUE),
  years_max = max(year, na.rm=TRUE),
  months_min = min(month),
  months_max = max(month),
  share_type3 = mean(type_of_violence == 3, na.rm=TRUE),
  share_missing_civ = mean(is.na(deaths_civilians))
)]
write_audit_dt(summ_clean, "ged_clean_summary")

# -------------------------
# 7) Aggregate to PRIO grid×month
# -------------------------
# Create counts + fatality sums, overall and for type 3
ged_pg_month <- ged_clean[, .(
  events_ged = .N,
  fat_best_ged = sum(best, na.rm = TRUE),
  civ_best_ged = sum(fifelse(is.na(deaths_civilians), 0L, deaths_civilians), na.rm = TRUE),
  events_type3 = sum(type_of_violence == 3, na.rm = TRUE),
  fat_best_type3 = sum(fifelse(type_of_violence == 3, best, 0L), na.rm = TRUE),
  civ_best_type3 = sum(fifelse(type_of_violence == 3, fifelse(is.na(deaths_civilians), 0L, deaths_civilians), 0L), na.rm = TRUE)
), by = .(pg_id = priogrid_gid, month)]

# add month_id stable ordering
ged_pg_month[, month_id := as.integer(factor(month, levels = sort(unique(month))))]

# sanity: no missing pg_id or month
if (ged_pg_month[is.na(pg_id) | is.na(month), .N] > 0) stop("Aggregated ged_pg_month has missing pg_id or month.")

# audit summary
summ_pg <- ged_pg_month[, .(
  n_cells_months = .N,
  months = uniqueN(month),
  pg_cells = uniqueN(pg_id),
  tot_events = sum(events_ged),
  tot_events_type3 = sum(events_type3),
  share_type3_events = sum(events_type3) / sum(events_ged)
)]
write_audit_dt(summ_pg, "ged_pg_month_summary")

out_pg <- file.path(DIRS$base, "ged_pg_month.rds")
saveRDS(ged_pg_month, out_pg)
log_msg("Wrote: ", out_pg)

log_msg("Step 3 complete.")
