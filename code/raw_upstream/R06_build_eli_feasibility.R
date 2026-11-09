# ============================================================
# R06_build_eli_feasibility.R — Build ELI feasibility (Option C)
# FIXED: robust numeric parsing for terrain vars + fail-fast checks
#
# Inputs:
#   data/processed/base/priogrid_static.rds
#   data/processed/base/airbases_points.rds
#
# Output:
#   data/processed/eli/eli_feasibility_pg.rds
#
# Audits:
#   outputs/paper1/audit/eli_feasibility_summary_*.csv
#   outputs/paper1/audit/eli_feasibility_top_bottom_*.csv
#   outputs/paper1/audit/eli_feasibility_terrain_audit_*.csv
# ============================================================

source("R00_config.R")
log_msg("Step 6: build ELI feasibility (terrain × basing)")

# -------------------------
# 1) Load inputs
# -------------------------
priopath <- file.path(DIRS$base, "priogrid_static.rds")
airpath  <- file.path(DIRS$base, "airbases_points.rds")

if (!file.exists(priopath)) stop("Missing: ", priopath, " (run Step 1)")
if (!file.exists(airpath))  stop("Missing: ", airpath,  " (run Step 2)")

pg <- readRDS(priopath)
ab <- readRDS(airpath)

assert_cols(pg, c("gid","xcoord","ycoord","mountains_mean","forest_gc","water_gc"), "priogrid_static")
assert_unique(pg, "gid", "priogrid_static")

assert_cols(ab, c("base_id","lat","lon"), "airbases_points")
assert_unique(ab, "base_id", "airbases_points")

# -------------------------
# 2) Robust numeric parsing (terrain vars)
# -------------------------
to_num_locale <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u00A0", "", x, fixed = TRUE) # NBSP
  x <- gsub(" ", "", x, fixed = TRUE)
  
  # if comma but no dot, treat comma as decimal separator
  x <- ifelse(grepl(",", x, fixed = TRUE) & !grepl("\\.", x), gsub(",", ".", x, fixed = TRUE), x)
  
  # strip everything except number tokens
  x <- gsub("[^0-9eE\\+\\-\\.]", "", x)
  
  suppressWarnings(as.numeric(x))
}

for (v in c("mountains_mean","forest_gc","water_gc")) {
  pg[, (v) := to_num_locale(get(v))]
}

# If forest/water are percent-like (0..100), bring to 0..1 (not required for z-scores, but nice)
rescale_if_pct <- function(x) {
  ok <- x[!is.na(x)]
  if (!length(ok)) return(x)
  mx <- max(ok); mn <- min(ok)
  if (is.finite(mx) && is.finite(mn) && mn >= 0 && mx > 1.5 && mx <= 100.5) return(x/100)
  x
}
pg[, forest_gc := rescale_if_pct(forest_gc)]
pg[, water_gc  := rescale_if_pct(water_gc)]

terrain_audit <- data.table(
  var = c("mountains_mean","forest_gc","water_gc"),
  class = c(class(pg$mountains_mean)[1], class(pg$forest_gc)[1], class(pg$water_gc)[1]),
  share_na = c(mean(is.na(pg$mountains_mean)), mean(is.na(pg$forest_gc)), mean(is.na(pg$water_gc))),
  min = c(min(pg$mountains_mean, na.rm=TRUE), min(pg$forest_gc, na.rm=TRUE), min(pg$water_gc, na.rm=TRUE)),
  max = c(max(pg$mountains_mean, na.rm=TRUE), max(pg$forest_gc, na.rm=TRUE), max(pg$water_gc, na.rm=TRUE))
)
write_audit_dt(terrain_audit, "eli_feasibility_terrain_audit")

# Fail fast if any terrain var is basically missing everywhere
if (any(terrain_audit$share_na > 0.25)) {
  stop("Terrain vars have too much missingness (>25%). See audit eli_feasibility_terrain_audit_*.csv")
}

# -------------------------
# 3) Ensure PRIO coords are lon/lat
# -------------------------
lon <- pg$xcoord
lat <- pg$ycoord

in_lon <- function(x) quantile(x, 0.999, na.rm=TRUE) <= 180 && quantile(x, 0.001, na.rm=TRUE) >= -180
in_lat <- function(x) quantile(x, 0.999, na.rm=TRUE) <= 90  && quantile(x, 0.001, na.rm=TRUE) >= -90

if (!in_lon(lon) || !in_lat(lat)) {
  if (in_lon(lat) && in_lat(lon)) {
    log_msg("Detected swapped coords in PRIO static. Swapping xcoord/ycoord.")
    tmp <- lon; lon <- lat; lat <- tmp
  } else {
    stop("PRIO coords do not look like lon/lat. Check xcoord/ycoord ranges in priogrid_static.rds.")
  }
}
pg[, `:=`(lon = lon, lat = lat)]

# -------------------------
# 4) Terrain feasibility component
# -------------------------
z_safe <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

scale01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] == rng[1]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

terrain_raw <- (-z_safe(pg$mountains_mean) + -z_safe(pg$forest_gc) + -z_safe(pg$water_gc)) / 3
pg[, eli_feas_terrain := scale01(terrain_raw)]

# Fail fast if terrain is constant (this is what happened to you before)
if (uniqueN(pg$eli_feas_terrain) <= 1) {
  stop("eli_feas_terrain is constant (no variation). Terrain parsing still failed — stop here.")
}

# -------------------------
# 4b) Geographic pre-filter: restrict to ISR-relevant airbases
# -------------------------
# The Wikidata list includes civilian airports worldwide (Heathrow, JFK, etc.).
# These inflate eli_feas_range scores in cells far from any conflict theatre.
# We keep only airbases within AIRBASE_MAX_KM of at least one conflict-theatre centroid.
# This does NOT filter the analysis sample — it only changes which bases drive eli_feas_range.
#
# Theatre centroids (approximate geographic centres):
THEATRE_CENTROIDS <- data.table(
  theatre    = c("Afghanistan","Pakistan","Iraq","Syria","Yemen","Somalia","Libya","Ukraine"),
  lat        = c(33.9,        30.4,      33.2,  35.0,   15.6,   6.1,     26.3,   49.0),
  lon        = c(67.7,        69.3,      43.7,  38.0,   48.5,   46.2,    17.2,   32.0)
)

AIRBASE_MAX_KM <- 2000  # keep bases within 2000 km of any theatre centroid
# covers regional hubs (Ramstein, Diego Garcia, Djibouti, etc.)
# without including wholly irrelevant civilian airports

haversine_km_vec <- function(lat1, lon1, lat2_vec, lon2_vec) {
  R <- 6371.0; to_rad <- pi / 180
  lat1r <- lat1 * to_rad; lon1r <- lon1 * to_rad
  lat2r <- lat2_vec * to_rad; lon2r <- lon2_vec * to_rad
  dlat <- lat2r - lat1r; dlon <- lon2r - lon1r
  a <- sin(dlat/2)^2 + cos(lat1r) * cos(lat2r) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

# For each airbase, compute minimum distance to any theatre centroid
ab[, min_dist_to_theatre := {
  dists <- vapply(seq_len(nrow(THEATRE_CENTROIDS)), function(i) {
    haversine_km_vec(THEATRE_CENTROIDS$lat[i], THEATRE_CENTROIDS$lon[i], lat, lon)
  }, numeric(.N))
  if (is.matrix(dists)) apply(dists, 1, min) else dists
}]

n_before_filter <- nrow(ab)
ab <- ab[min_dist_to_theatre <= AIRBASE_MAX_KM]
n_after_filter  <- nrow(ab)

log_msg(sprintf(
  "Airbase geographic filter: kept %d of %d bases (within %d km of a theatre centroid).",
  n_after_filter, n_before_filter, AIRBASE_MAX_KM
))

if (n_after_filter == 0) stop("No airbases remain after geographic filter. Check THEATRE_CENTROIDS or increase AIRBASE_MAX_KM.")

write_audit_dt(
  data.table(
    n_before = n_before_filter,
    n_after  = n_after_filter,
    n_dropped = n_before_filter - n_after_filter,
    max_km_threshold = AIRBASE_MAX_KM
  ),
  "eli_airbase_geographic_filter"
)

# Also write which bases survived (useful for face-validity check)
write_audit_dt(
  ab[order(min_dist_to_theatre), .(base_id, asset_name, lat, lon, min_dist_to_theatre)][1:min(200, .N)],
  "eli_airbase_geographic_filter_kept_examples"
)

ab[, min_dist_to_theatre := NULL]  # clean up temp column

# -------------------------
# 5) Nearest airbase distance (fast via RANN)
# -------------------------
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371.0
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad; lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad; lon2 <- lon2 * to_rad
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

R_KM <- 400

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("RANN not installed but required for this run (you had it earlier). Install.packages('RANN')")
}

log_msg("Using fast nearest-neighbor via RANN (KD-tree on sphere).")

to_xyz <- function(lat_deg, lon_deg) {
  latr <- lat_deg * pi/180
  lonr <- lon_deg * pi/180
  cbind(cos(latr) * cos(lonr), cos(latr) * sin(lonr), sin(latr))
}

base_xyz <- to_xyz(ab$lat, ab$lon)
cell_xyz <- to_xyz(pg$lat, pg$lon)

nn <- RANN::nn2(data = base_xyz, query = cell_xyz, k = 1)
idx <- nn$nn.idx[,1]

pg[, min_km_airbase := haversine_km(lat, lon, ab$lat[idx], ab$lon[idx])]
pg[, eli_feas_range := exp(-min_km_airbase / R_KM)]
pg[, eli_feas := eli_feas_terrain * eli_feas_range]

# -------------------------
# 6) Output + audits (ROBUST)
# -------------------------

# If any terrain feasibility is NA, replace with median (conservative, avoids dropping cells)
n_na_terrain <- sum(is.na(pg$eli_feas_terrain))
if (n_na_terrain > 0) {
  med_terrain <- median(pg$eli_feas_terrain, na.rm = TRUE)
  pg[is.na(eli_feas_terrain), eli_feas_terrain := med_terrain]
  log_msg("Imputed eli_feas_terrain NA with median. n_na_terrain=", n_na_terrain, " median=", round(med_terrain, 4))
}

# If any min_km_airbase is NA (shouldn’t happen), recompute those rows only
n_na_km <- sum(is.na(pg$min_km_airbase))
if (n_na_km > 0) {
  log_msg("min_km_airbase has NA for ", n_na_km, " rows; recomputing those via RANN.")
  idx_na <- which(is.na(pg$min_km_airbase))
  
  to_xyz <- function(lat_deg, lon_deg) {
    latr <- lat_deg * pi/180
    lonr <- lon_deg * pi/180
    cbind(cos(latr) * cos(lonr), cos(latr) * sin(lonr), sin(latr))
  }
  base_xyz <- to_xyz(ab$lat, ab$lon)
  q_xyz    <- to_xyz(pg$lat[idx_na], pg$lon[idx_na])
  
  nn2 <- RANN::nn2(data = base_xyz, query = q_xyz, k = 1)
  idx2 <- nn2$nn.idx[,1]
  
  pg$min_km_airbase[idx_na] <- haversine_km(
    pg$lat[idx_na], pg$lon[idx_na],
    ab$lat[idx2], ab$lon[idx2]
  )
}

# Recompute range + composite (in case of any imputation)
pg[, eli_feas_range := exp(-min_km_airbase / R_KM)]
pg[, eli_feas := eli_feas_terrain * eli_feas_range]

out <- pg[, .(gid, lon, lat,
              eli_feas_terrain, min_km_airbase, eli_feas_range, eli_feas)]

# Audit missingness in final outputs (should be zero)
miss_out <- out[, .(
  n = .N,
  na_terrain = sum(is.na(eli_feas_terrain)),
  na_km      = sum(is.na(min_km_airbase)),
  na_range   = sum(is.na(eli_feas_range)),
  na_feas    = sum(is.na(eli_feas))
)]
write_audit_dt(miss_out, "eli_feasibility_missingness")

# Final sanity checks (now should never fail)
if (any(!is.finite(out$eli_feas))) {
  # write examples before stopping
  bad <- out[!is.finite(eli_feas)][1:min(200, .N)]
  write_audit_dt(bad, "eli_feasibility_nonfinite_examples")
  stop("eli_feas has non-finite values. See audit eli_feasibility_nonfinite_examples_*.csv")
}
if (any(out$eli_feas < 0 | out$eli_feas > 1, na.rm=TRUE)) stop("eli_feas outside [0,1].")

summ <- out[, .(
  n = .N,
  terrain_p01 = quantile(eli_feas_terrain, 0.01),
  terrain_med = median(eli_feas_terrain),
  terrain_p99 = quantile(eli_feas_terrain, 0.99),
  min_km_p01 = quantile(min_km_airbase, 0.01),
  min_km_med = median(min_km_airbase),
  min_km_p99 = quantile(min_km_airbase, 0.99),
  feas_p01 = quantile(eli_feas, 0.01),
  feas_med = median(eli_feas),
  feas_p99 = quantile(eli_feas, 0.99),
  cor_km_feas = cor(min_km_airbase, eli_feas),
  cor_terrain_feas = cor(eli_feas_terrain, eli_feas)
)]
write_audit_dt(summ, "eli_feasibility_summary")

topbot <- rbind(
  out[order(-eli_feas)][1:20],
  out[order(eli_feas)][1:20]
)
write_audit_dt(topbot, "eli_feasibility_top_bottom")

out_path <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
saveRDS(out, out_path)
log_msg("Wrote: ", out_path)
log_msg("Step 6 complete.")

# ============================================================
# R06_validate_eli_feasibility.R — Validation checks for Step 6
# ============================================================

source("R00_config.R")
log_msg("Validate Step 6: ELI feasibility output")

# Paths
priopath <- file.path(DIRS$base, "priogrid_static.rds")
airpath  <- file.path(DIRS$base, "airbases_points.rds")
efpath   <- file.path(DIRS$eli,  "eli_feasibility_pg.rds")

for (p in c(priopath, airpath, efpath)) if (!file.exists(p)) stop("Missing file: ", p)

pg <- readRDS(priopath)
ab <- readRDS(airpath)
ef <- readRDS(efpath)

# -------------------------
# A) Structure checks
# -------------------------
assert_cols(pg, c("gid","xcoord","ycoord","mountains_mean","forest_gc","water_gc"), "priogrid_static")
assert_unique(pg, "gid", "priogrid_static")

assert_cols(ab, c("base_id","lat","lon"), "airbases_points")
assert_unique(ab, "base_id", "airbases_points")

assert_cols(ef, c("gid","lon","lat","eli_feas_terrain","min_km_airbase","eli_feas_range","eli_feas"), "eli_feasibility_pg")
assert_unique(ef, "gid", "eli_feasibility_pg")

# row counts must match PRIO static
if (nrow(ef) != nrow(pg)) stop("Row mismatch: ef=", nrow(ef), " pg=", nrow(pg))

# gid coverage must match exactly
miss_in_ef <- setdiff(pg$gid, ef$gid)
miss_in_pg <- setdiff(ef$gid, pg$gid)
if (length(miss_in_ef) || length(miss_in_pg)) {
  write_audit_dt(data.table(miss_in_ef = head(miss_in_ef, 200)), "validate_gid_missing_in_ef")
  write_audit_dt(data.table(miss_in_pg = head(miss_in_pg, 200)), "validate_gid_missing_in_pg")
  stop("GID sets differ between PRIO and ELI feasibility outputs.")
}

# -------------------------
# B) Value sanity checks
# -------------------------
# No missingness allowed in final ef columns
miss <- ef[, .(
  na_terrain = sum(is.na(eli_feas_terrain)),
  na_km      = sum(is.na(min_km_airbase)),
  na_range   = sum(is.na(eli_feas_range)),
  na_feas    = sum(is.na(eli_feas))
)]
write_audit_dt(miss, "validate_missingness_final")

if (any(miss != 0)) stop("Final ef has missing values. See audit validate_missingness_final_*.csv")

# Bounds
if (any(ef$eli_feas < 0 | ef$eli_feas > 1)) stop("eli_feas outside [0,1]")
if (any(ef$eli_feas_range <= 0 | ef$eli_feas_range > 1)) stop("eli_feas_range outside (0,1]")
if (any(ef$eli_feas_terrain < 0 | ef$eli_feas_terrain > 1)) stop("eli_feas_terrain outside [0,1]")
if (any(ef$min_km_airbase < 0)) stop("min_km_airbase has negative values")

# Terrain variation must exist (prevents the old 'all 0.5' failure mode)
if (uniqueN(ef$eli_feas_terrain) <= 10) {
  stop("eli_feas_terrain has suspiciously low variation (<=10 unique values).")
}

# Distribution audits
dist <- ef[, .(
  n = .N,
  km_p01 = quantile(min_km_airbase, 0.01),
  km_med = median(min_km_airbase),
  km_p99 = quantile(min_km_airbase, 0.99),
  feas_p01 = quantile(eli_feas, 0.01),
  feas_med = median(eli_feas),
  feas_p99 = quantile(eli_feas, 0.99),
  terrain_p01 = quantile(eli_feas_terrain, 0.01),
  terrain_med = median(eli_feas_terrain),
  terrain_p99 = quantile(eli_feas_terrain, 0.99),
  cor_km_feas = cor(min_km_airbase, eli_feas),
  cor_km_range = cor(min_km_airbase, eli_feas_range),
  cor_terrain_feas = cor(eli_feas_terrain, eli_feas)
)]
write_audit_dt(dist, "validate_distributions")

# -------------------------
# C) Confirm where the 2835 terrain-missing came from (PRIO covariates)
# -------------------------
to_num_locale <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u00A0", "", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)
  x <- ifelse(grepl(",", x, fixed = TRUE) & !grepl("\\.", x), gsub(",", ".", x, fixed = TRUE), x)
  x <- gsub("[^0-9eE\\+\\-\\.]", "", x)
  suppressWarnings(as.numeric(x))
}
rescale_if_pct <- function(x) {
  ok <- x[!is.na(x)]
  if (!length(ok)) return(x)
  mx <- max(ok); mn <- min(ok)
  if (is.finite(mx) && is.finite(mn) && mn >= 0 && mx > 1.5 && mx <= 100.5) return(x/100)
  x
}

mm <- to_num_locale(pg$mountains_mean)
ff <- rescale_if_pct(to_num_locale(pg$forest_gc))
ww <- rescale_if_pct(to_num_locale(pg$water_gc))

na_any <- is.na(mm) | is.na(ff) | is.na(ww)
na_counts <- data.table(
  n = nrow(pg),
  na_any = sum(na_any),
  na_mountains = sum(is.na(mm)),
  na_forest    = sum(is.na(ff)),
  na_water     = sum(is.na(ww))
)
write_audit_dt(na_counts, "validate_prio_terrain_missingness")

# Optional: write a small sample of the NA rows for inspection
if (na_counts$na_any > 0) {
  samp <- pg[na_any == TRUE, .(gid, xcoord, ycoord, mountains_mean, forest_gc, water_gc)][1:min(200, .N)]
  write_audit_dt(samp, "validate_prio_terrain_missing_examples")
}

log_msg("Validation complete: Step 6 output looks consistent.")