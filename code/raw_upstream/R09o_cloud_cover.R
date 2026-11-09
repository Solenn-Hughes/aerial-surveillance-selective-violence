# ============================================================
# R09o_cloud_cover_robustness.R  (FIXED)
# Adds ERA5 monthly cloud cover as time-varying IMINT component
#
# FIXES:
#   1) ERA5 time dimension can be "time" or "valid_time" —
#      now auto-detects
#   2) Uses xcoord/ycoord from panel (avoids coordinate merge)
#   3) More robust NetCDF dimension handling
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 9o: ERA5 cloud cover robustness")

# ── CHECK FOR ERA5 DATA ──────────────────────────────────────
ERA5_PATHS <- c(
  file.path(DIRS$data_raw, "era5", "era5_tcc_monthly.nc"),
  file.path(PROJECT_ROOT, "data", "raw", "era5", "era5_tcc_monthly.nc"),
  file.path(PROJECT_ROOT, "data", "raw", "era5_tcc_monthly.nc")
)

era5_path <- NULL
for (pp in ERA5_PATHS) {
  if (file.exists(pp)) { era5_path <- pp; break }
}

if (is.null(era5_path)) {
  stop("ERA5 cloud cover data not found.\n",
       "Expected at: ", ERA5_PATHS[1], "\n",
       "See DOWNLOAD_AND_RUN_GUIDE.R for download instructions.\n")
}

if (!requireNamespace("ncdf4", quietly = TRUE)) {
  stop("Package 'ncdf4' required. Install via: install.packages('ncdf4')")
}
library(ncdf4)

# ── READ ERA5 NETCDF ─────────────────────────────────────────
log_msg("Reading ERA5 cloud cover from: ", era5_path)
nc <- nc_open(era5_path)

# Print what's in the file for debugging
cat("── NetCDF contents ──\n")
cat("Variables:", paste(names(nc$var), collapse = ", "), "\n")
cat("Dimensions:", paste(names(nc$dim), collapse = ", "), "\n")

# Find cloud cover variable
tcc_candidates <- c("tcc", "TCC", "total_cloud_cover")
tcc_name <- intersect(names(nc$var), tcc_candidates)
if (!length(tcc_name)) {
  # Try first variable if none of the expected names match
  tcc_name <- names(nc$var)[1]
  log_msg("WARNING: No standard TCC variable found. Using first variable: ", tcc_name)
}
tcc_name <- tcc_name[1]
log_msg("Using variable: ", tcc_name)

# FIX: Auto-detect dimension names (ERA5 uses different names depending on download method)
dim_names <- names(nc$dim)

# Longitude: could be "longitude", "lon", "x"
lon_name <- intersect(dim_names, c("longitude", "lon", "x"))[1]
if (is.na(lon_name)) stop("Cannot find longitude dimension. Available: ", paste(dim_names, collapse = ", "))

# Latitude: could be "latitude", "lat", "y"
lat_name <- intersect(dim_names, c("latitude", "lat", "y"))[1]
if (is.na(lat_name)) stop("Cannot find latitude dimension. Available: ", paste(dim_names, collapse = ", "))

# Time: could be "time", "valid_time", "date", "month"
time_name <- intersect(dim_names, c("valid_time", "time", "date", "forecast_reference_time"))[1]
if (is.na(time_name)) stop("Cannot find time dimension. Available: ", paste(dim_names, collapse = ", "))

log_msg("Dimensions: lon=", lon_name, " lat=", lat_name, " time=", time_name)

lon_era <- ncvar_get(nc, lon_name)
lat_era <- ncvar_get(nc, lat_name)
time_raw <- ncvar_get(nc, time_name)

# Convert time to dates
# ERA5 time encoding varies:
#   - "hours since 1900-01-01" (classic CDS)
#   - Unix timestamps (seconds since 1970-01-01)
#   - ISO date strings
time_units <- ncatt_get(nc, time_name, "units")

if (time_units$hasatt) {
  tu <- time_units$value
  log_msg("Time units: ", tu)
  
  if (grepl("hours since 1900", tu, ignore.case = TRUE)) {
    dates <- as.Date(as.POSIXct("1900-01-01 00:00:00", tz = "UTC") + time_raw * 3600)
  } else if (grepl("seconds since 1970", tu, ignore.case = TRUE)) {
    dates <- as.Date(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC"))
  } else if (grepl("days since", tu, ignore.case = TRUE)) {
    origin <- sub("days since ", "", tu, ignore.case = TRUE)
    dates <- as.Date(time_raw, origin = as.Date(trimws(origin)))
  } else {
    # Try interpreting as numeric dates
    dates <- tryCatch(as.Date(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC")),
                      error = function(e) {
                        log_msg("Could not parse time. Trying as hours since 1900.")
                        as.Date(as.POSIXct("1900-01-01", tz = "UTC") + time_raw * 3600)
                      })
  }
} else {
  # No units attribute — try common encodings
  if (all(time_raw > 1e9)) {
    # Looks like Unix timestamp (seconds since 1970)
    dates <- as.Date(as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC"))
  } else if (all(time_raw > 800000)) {
    # Looks like hours since 1900
    dates <- as.Date(as.POSIXct("1900-01-01", tz = "UTC") + time_raw * 3600)
  } else {
    stop("Cannot determine time encoding. Raw values range: ",
         min(time_raw), " to ", max(time_raw))
  }
}

months_era <- as.IDate(format(dates, "%Y-%m-01"))
log_msg("Date range: ", min(months_era), " to ", max(months_era), " (", length(months_era), " months)")

# Read the cloud cover data
tcc <- ncvar_get(nc, tcc_name)
nc_close(nc)

log_msg("ERA5 array dimensions: ", paste(dim(tcc), collapse = " x "))

# Sanity: ERA5 tcc should be fraction [0, 1]
tcc_range <- range(tcc, na.rm = TRUE)
log_msg("TCC range: [", round(tcc_range[1], 3), ", ", round(tcc_range[2], 3), "]")
if (tcc_range[2] > 1.5) {
  log_msg("Rescaling from percentage to fraction.")
  tcc <- tcc / 100
}

# ── AGGREGATE TO PRIO-GRID ───────────────────────────────────
log_msg("Aggregating ERA5 to PRIO-GRID resolution...")

P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

# Get coordinates from panel (xcoord/ycoord)
if ("xcoord" %in% names(dt) && "ycoord" %in% names(dt)) {
  coords <- unique(dt[, .(pg_id, lon = xcoord, lat = ycoord)])
} else {
  P_ELI <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
  eli <- readRDS(P_ELI); setDT(eli)
  if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := as.integer(gid)]
  coords <- unique(eli[, .(pg_id, lon, lat)])
}

# Only cells actually in our panel
panel_cells <- unique(dt$pg_id)
coords <- coords[pg_id %in% panel_cells]

era5_res <- median(abs(diff(lon_era[1:min(5, length(lon_era))])))
log_msg("ERA5 resolution: ~", round(era5_res, 3), " degrees")

# Vectorised approach: for each month, extract values for all cells at once
cloud_list <- list()

for (t_idx in seq_along(months_era)) {
  m_date <- months_era[t_idx]
  tcc_slice <- tcc[, , t_idx]  # lon x lat matrix
  
  # For each cell, find nearest ERA5 pixel
  vals <- vapply(seq_len(nrow(coords)), function(i) {
    lon_idx <- which.min(abs(lon_era - coords$lon[i]))
    lat_idx <- which.min(abs(lat_era - coords$lat[i]))
    tcc_slice[lon_idx, lat_idx]
  }, numeric(1))
  
  cloud_list[[t_idx]] <- data.table(
    pg_id = coords$pg_id,
    month = m_date,
    tcc_mean = vals
  )
  
  if (t_idx %% 24 == 0) log_msg("  Processed ", t_idx, " / ", length(months_era), " months")
}

cloud_dt <- rbindlist(cloud_list)
log_msg("Cloud cover extracted: ", nrow(cloud_dt), " cell-months")

cloud_dt[, cloud_access := 1 - tcc_mean]

# ── MERGE INTO PANEL ─────────────────────────────────────────
n_before <- nrow(dt)
dt <- merge(dt, cloud_dt[, .(pg_id, month, cloud_access)],
            by = c("pg_id", "month"), all.x = TRUE)
stopifnot(nrow(dt) == n_before)

# Fill NAs
na_count <- sum(is.na(dt$cloud_access))
log_msg("Cloud access NAs: ", na_count, " / ", nrow(dt))
if (na_count > 0) {
  dt[, cloud_access_fill := median(cloud_access, na.rm = TRUE), by = .(theatre_id, month)]
  dt[is.na(cloud_access), cloud_access := cloud_access_fill]
  dt[, cloud_access_fill := NULL]
  
  n_still <- sum(is.na(dt$cloud_access))
  if (n_still > 0) {
    dt[is.na(cloud_access), cloud_access := median(dt$cloud_access, na.rm = TRUE)]
  }
}

# ── BUILD IMINT-AUGMENTED ELI_FEAS ───────────────────────────
if (!("eli_feas_terrain" %in% names(dt)) || !("eli_feas_range" %in% names(dt))) {
  P_ELI <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
  eli <- readRDS(P_ELI); setDT(eli)
  if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := as.integer(gid)]
  eli_merge <- unique(eli[, .(pg_id, eli_feas_terrain, eli_feas_range)])
  n_before <- nrow(dt)
  dt <- merge(dt, eli_merge, by = "pg_id", all.x = TRUE, suffixes = c("", ".eli2"))
  stopifnot(nrow(dt) == n_before)
}

dt[, eli_feas_imint := eli_feas_terrain * cloud_access * eli_feas_range]

log_msg("eli_feas_imint constructed.")
cat("  Mean:", round(mean(dt$eli_feas_imint, na.rm = TRUE), 4), "\n")
cat("  SD:", round(sd(dt$eli_feas_imint, na.rm = TRUE), 4), "\n")
cat("  Cor with eli_feas:", round(cor(dt$eli_feas, dt$eli_feas_imint, use = "complete.obs"), 4), "\n")

# ── RUN MODELS ───────────────────────────────────────────────
CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")
dt[, log1p_y := log1p(events_type3)]

fml_base <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m0 <- feols(fml_base, data = dt, cluster = ~theatre_id, notes = FALSE)

fml_imint <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas_imint + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m1 <- feols(fml_imint, data = dt, cluster = ~theatre_id, notes = FALSE)

fml_both <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + post_eli:cloud_access + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m2 <- feols(fml_both, data = dt, cluster = ~theatre_id, notes = FALSE)

# ── EXTRACT AND COMPARE ──────────────────────────────────────
extract_hit <- function(m, pattern) {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  ct[grepl(pattern, term)]
}

cat("\n============================================================\n")
cat("CLOUD COVER ROBUSTNESS\n")
cat("============================================================\n")

cat("\nM0: Baseline (eli_feas, time-invariant)\n")
h0 <- extract_hit(m0, "post_eli.*eli_feas|eli_feas.*post_eli")
if (nrow(h0)) print(h0[, .(term, Estimate = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

cat("\nM1: IMINT-augmented (eli_feas_terrain * cloud_access * range)\n")
h1 <- extract_hit(m1, "post_eli.*eli_feas_imint|eli_feas_imint.*post_eli")
if (nrow(h1)) print(h1[, .(term, Estimate = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

cat("\nM2: Horse race (eli_feas + cloud_access)\n")
h2 <- extract_hit(m2, "post_eli|eli_feas|cloud")
if (nrow(h2)) print(h2[, .(term, Estimate = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

results <- data.table(
  model = c("M0: eli_feas (baseline)", "M1: eli_feas_imint (+ cloud)", "M2: horse race (eli_feas)"),
  coef = c(
    if (nrow(h0)) round(h0$Estimate[1], 5) else NA,
    if (nrow(h1)) round(h1$Estimate[1], 5) else NA,
    if (nrow(h2)) round(h2$Estimate[1], 5) else NA
  ),
  se = c(
    if (nrow(h0)) round(h0$`Std. Error`[1], 5) else NA,
    if (nrow(h1)) round(h1$`Std. Error`[1], 5) else NA,
    if (nrow(h2)) round(h2$`Std. Error`[1], 5) else NA
  ),
  p = c(
    if (nrow(h0)) round(h0$`Pr(>|t|)`[1], 4) else NA,
    if (nrow(h1)) round(h1$`Pr(>|t|)`[1], 4) else NA,
    if (nrow(h2)) round(h2$`Pr(>|t|)`[1], 4) else NA
  )
)

cat("\n── Summary ──\n")
print(results)

write_audit_dt(results, "R09o_cloud_cover_robustness")
log_msg("Step 9o complete.")