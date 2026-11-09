# ============================================================
# R05_build_cw_panels.R — Build civil-war restrictions (ALLCW + ACD_STRICT) and filtered panels
# FIXED: month expansion is vector-safe (Map over rows)
# ============================================================

source("R00_config.R")
log_msg("Step 5: build civil-war restrictions (ALLCW + ACD_STRICT) and filtered panels")

# Paths
ACD_CSV <- file.path(PROJECT_ROOT, "data", "raw", "ucdp_acd", "UcdpPrioConflict_v25_1.csv")
GED_RDS <- file.path(DIRS$base, "ged_clean.rds")
PAN_RDS <- file.path(DIRS$base, "panel_pg_month_base.rds")

if (!file.exists(ACD_CSV)) stop("Missing ACD CSV: ", ACD_CSV)
if (!file.exists(GED_RDS)) stop("Missing GED clean: ", GED_RDS)
if (!file.exists(PAN_RDS)) stop("Missing base panel: ", PAN_RDS)

dir_panels <- file.path(PROJECT_ROOT, "data", "processed", "panels")
dir.create(dir_panels, recursive = TRUE, showWarnings = FALSE)

write_audit_dt(data.table(path = ACD_CSV, mtime = file.info(ACD_CSV)$mtime), "acd_inputs")

# -------------------------
# 1) Build pg_id -> country_id map from GED (modal)
# -------------------------
ged <- readRDS(GED_RDS)

ct <- ged[, .N, by = .(pg_id = priogrid_gid, country_id, country)]
setorder(ct, pg_id, -N)

ct_n <- ct[, .(n_countries = uniqueN(country_id), topN = max(N)), by = pg_id]
multi <- ct_n[n_countries > 1]

if (nrow(multi) > 0) {
  ex <- ct[pg_id %in% multi$pg_id][order(pg_id, -N)][, head(.SD, 3), by=pg_id]
  write_audit_dt(ex[1:min(200, .N)], "pg_country_map_multicountry_examples")
}

pg_map <- ct[, .SD[1], by = pg_id]  # modal row

map_summ <- data.table(
  pg_cells_with_events = uniqueN(ged$priogrid_gid),
  pg_cells_mapped = nrow(pg_map),
  share_multicountry = nrow(multi) / max(1, nrow(pg_map))
)
write_audit_dt(map_summ, "pg_country_map_summary")

saveRDS(pg_map[, .(pg_id, country_id, country)], file.path(DIRS$base, "pg_country_map_from_ged.rds"))
log_msg("Wrote: ", file.path(DIRS$base, "pg_country_map_from_ged.rds"))

months_ged <- sort(unique(ged$month))

# -------------------------
# 2) Ingest ACD and build country×month civil-war sets
# -------------------------
acd <- fread(ACD_CSV, showProgress = TRUE, fill = TRUE, encoding = "UTF-8")

need <- c("year","type_of_conflict","intensity_level","gwno_loc","start_date","ep_end_date")
miss <- setdiff(need, names(acd))
if (length(miss)) stop("ACD missing required cols: ", paste(miss, collapse = ", "))

to_int <- function(x) suppressWarnings(as.integer(x))

acd[, year := to_int(year)]
acd[, type_of_conflict := to_int(type_of_conflict)]
acd[, intensity_level := to_int(intensity_level)]
acd[, gwno_loc := to_int(gwno_loc)]
acd[, start_date := as.IDate(start_date)]
acd[, ep_end_date := as.IDate(ep_end_date)]

# Keep only internal + internationalized internal
acd_cw <- acd[type_of_conflict %in% c(3L, 4L) & !is.na(year) & !is.na(gwno_loc)]

y_start <- function(y) as.IDate(paste0(y, "-01-01"))
y_end   <- function(y) as.IDate(paste0(y, "-12-31"))

acd_cw[, start_in_year := fifelse(is.na(start_date), y_start(year), pmax(start_date, y_start(year)))]
acd_cw[, end_in_year   := fifelse(is.na(ep_end_date), y_end(year),   pmin(ep_end_date, y_end(year)))]

bad_span <- acd_cw[end_in_year < start_in_year]
if (nrow(bad_span) > 0) {
  write_audit_dt(bad_span[1:min(200, .N),
                          .(gwno_loc, year, start_date, ep_end_date, start_in_year, end_in_year)],
                 "acd_bad_spans_examples")
  stop("ACD has rows where end_in_year < start_in_year. See audit acd_bad_spans_examples_*.csv")
}

acd_summ <- acd_cw[, .(
  rows_cw = .N,
  countries = uniqueN(gwno_loc),
  years_min = min(year),
  years_max = max(year),
  share_intensity2 = mean(intensity_level == 2, na.rm=TRUE)
)]
write_audit_dt(acd_summ, "acd_clean_summary")

# ---- FIX: vector-safe month expansion ----
expand_months_scalar <- function(s, e) {
  if (is.na(s) || is.na(e)) return(character())
  s0 <- as.Date(format(as.Date(s), "%Y-%m-01"))
  e0 <- as.Date(format(as.Date(e), "%Y-%m-01"))
  if (e0 < s0) return(character())
  format(seq.Date(s0, e0, by = "month"), "%Y-%m")
}

make_cw_months <- function(dt) {
  # dt must have gwno_loc, start_in_year, end_in_year
  out <- dt[, .(
    month = unique(unlist(Map(expand_months_scalar, start_in_year, end_in_year), use.names = FALSE))
  ), by = .(country_id = gwno_loc)]
  
  out <- unique(out[month %in% months_ged])
  out
}

cw_allcw  <- make_cw_months(acd_cw)
cw_strict <- make_cw_months(acd_cw[intensity_level == 2L])

write_audit_dt(cw_allcw[, .(n=.N, countries=uniqueN(country_id),
                            months=uniqueN(month), min_m=min(month), max_m=max(month))],
               "cw_months_ALLCW_summary")
write_audit_dt(cw_strict[, .(n=.N, countries=uniqueN(country_id),
                             months=uniqueN(month), min_m=min(month), max_m=max(month))],
               "cw_months_ACD_STRICT_summary")

saveRDS(cw_allcw,  file.path(DIRS$base, "cw_country_month_ALLCW.rds"))
saveRDS(cw_strict, file.path(DIRS$base, "cw_country_month_ACD_STRICT.rds"))
log_msg("Wrote CW sets: ", file.path(DIRS$base, "cw_country_month_ALLCW.rds"),
        " and ", file.path(DIRS$base, "cw_country_month_ACD_STRICT.rds"))

# -------------------------
# 3) Filter base panel
# -------------------------
panel <- readRDS(PAN_RDS)

# attach country_id (modal from GED) — unmapped cells will be NA
setkey(pg_map, pg_id)
setkey(panel, pg_id)
panel <- pg_map[panel]

map_cov <- panel[, .(
  n_rows = .N,
  n_pg = uniqueN(pg_id),
  share_country_mapped = mean(!is.na(country_id))
)]
write_audit_dt(map_cov, "panel_country_map_coverage")

filter_panel <- function(pan, cw) {
  setkey(cw, country_id, month)
  pan[cw, on = .(country_id, month), nomatch = 0]
}

panel_allcw  <- filter_panel(panel, cw_allcw)
panel_strict <- filter_panel(panel, cw_strict)

audit_panel <- function(pan, name) {
  s <- pan[, .(
    n_rows = .N,
    n_pg = uniqueN(pg_id),
    n_countries = uniqueN(country_id),
    n_months = uniqueN(month),
    tot_events = sum(events_ged, na.rm=TRUE),
    tot_type3 = sum(events_type3, na.rm=TRUE),
    share_zero = mean(events_ged == 0)
  )]
  write_audit_dt(s, name)
}

audit_panel(panel_allcw,  "panel_filter_ALLCW_summary")
audit_panel(panel_strict, "panel_filter_ACD_STRICT_summary")

out1 <- file.path(dir_panels, "panel_pg_month_ALLCW.rds")
out2 <- file.path(dir_panels, "panel_pg_month_ALLCW_ACD_STRICT.rds")
saveRDS(panel_allcw, out1)
saveRDS(panel_strict, out2)

log_msg("Wrote: ", out1)
log_msg("Wrote: ", out2)
log_msg("Step 5 complete.")
