# ============================================================
# R_h2_spatial_adm2_diagnostic.R
# Diagnostic branch: improve H2 by spatially matching ACLED to ADM2
# Does NOT overwrite final dissertation results.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(fixest)
})

# -------------------------
# 1) PATHS
# -------------------------
THESIS_DIR <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
BASE_DIR <- file.path(THESIS_DIR, "outputs", "paper1_geo_afghanistan_executable_acled")

PANEL_PATH <- file.path(BASE_DIR, "panels", "panel_core_cell_month.rds")
RESULTS_PATH <- file.path(BASE_DIR, "tables", "all_model_results.csv")
ADM2_PATH <- file.path(THESIS_DIR, "data", "raw", "afg_admin", "afg_admin2.shp")

OUT_DIR <- file.path(BASE_DIR, "diagnostics_h2_spatial_adm2")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 2) HELPERS
# -------------------------
norm_text <- function(x) {
  x <- trimws(as.character(x))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

safe_name <- function(x) {
  x <- trimws(as.character(x))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

pick_col <- function(nms, candidates, required = TRUE) {
  sn <- safe_name(nms)
  sc <- safe_name(candidates)
  m <- match(sc, sn)
  m <- m[!is.na(m)]
  if (!length(m)) {
    if (required) stop("Could not find required column among: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  nms[m[1]]
}

parse_acled_date <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  out <- suppressWarnings(as.IDate(x))
  
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    p <- suppressWarnings(as.Date(x[bad], format = "%d %B %Y"))
    out[bad] <- as.IDate(p)
  }
  
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    p <- suppressWarnings(as.Date(x[bad], format = "%d %b %Y"))
    out[bad] <- as.IDate(p)
  }
  
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    p <- suppressWarnings(as.Date(x[bad], format = "%m/%d/%Y"))
    out[bad] <- as.IDate(p)
  }
  
  as.IDate(format(out, "%Y-%m-01"))
}

extract_term <- function(model, wanted, vcov_spec) {
  ct <- as.data.table(fixest::coeftable(model, vcov = vcov_spec), keep.rownames = "term")
  if ("Estimate" %in% names(ct)) setnames(ct, "Estimate", "estimate")
  if ("Std. Error" %in% names(ct)) setnames(ct, "Std. Error", "std_error")
  if ("z value" %in% names(ct)) setnames(ct, "z value", "stat")
  if ("t value" %in% names(ct)) setnames(ct, "t value", "stat")
  if ("Pr(>|z|)" %in% names(ct)) setnames(ct, "Pr(>|z|)", "p_value")
  if ("Pr(>|t|)" %in% names(ct)) setnames(ct, "Pr(>|t|)", "p_value")
  ct[term == wanted, .(term, estimate, std_error, stat, p_value)]
}

run_h2_model <- function(dt, moderator, label) {
  controls <- intersect(c("prec_gpcp", "temp", "nlights_calib_mean"), names(dt))
  rhs <- c(
    "post_closure:eli_feas",
    paste0("post_closure:eli_feas:", moderator),
    controls
  )
  
  fml <- as.formula(paste0(
    "events_type3 ~ ",
    paste(rhs, collapse = " + "),
    " | cell_id + month"
  ))
  
  m <- fepois(fml, data = dt)
  
  r1 <- extract_term(m, "post_closure:eli_feas", ~ cell_id + month)
  r2 <- extract_term(m, paste0("post_closure:eli_feas:", moderator), ~ cell_id + month)
  
  out <- rbindlist(list(
    data.table(model = label, moderator = moderator, component = "H1-within-H2", r1),
    data.table(model = label, moderator = moderator, component = "H2-triple", r2)
  ), fill = TRUE)
  
  out[, nobs := nobs(m)]
  out
}

# -------------------------
# 3) LOAD PANEL
# -------------------------
if (!file.exists(PANEL_PATH)) stop("Missing panel: ", PANEL_PATH)
core <- as.data.table(readRDS(PANEL_PATH))
core[, month := as.IDate(month)]
core[, cell_id := as.factor(cell_id)]

district_col <- intersect(
  c("district_id", "adm2_id", "adm2_pcode", "district_pcode"),
  names(core)
)[1]

if (is.na(district_col)) {
  stop("Could not find a district ID column in core panel. Available names are:\n",
       paste(names(core), collapse = ", "))
}

core[, district_id_for_h2 := as.character(get(district_col))]

cat("\nUsing core district column:", district_col, "\n")
cat("Core rows:", nrow(core), "\n")
cat("Unique core districts:", uniqueN(core$district_id_for_h2, na.rm = TRUE), "\n")

# -------------------------
# 4) FIND AND LOAD ACLED
# -------------------------
acled_files <- list.files(
  file.path(THESIS_DIR, "data", "raw"),
  pattern = "acled.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (!length(acled_files)) stop("No ACLED CSV found under data/raw.")

# Prefer non-OLD file if available
acled_files_ordered <- acled_files[order(grepl("old", acled_files, ignore.case = TRUE))]
ACLED_PATH <- acled_files_ordered[1]

cat("\nUsing ACLED file:\n", ACLED_PATH, "\n")

acled <- fread(ACLED_PATH, showProgress = FALSE)

country_col <- pick_col(names(acled), c("country"), required = FALSE)
lat_col <- pick_col(names(acled), c("latitude", "lat"), required = TRUE)
lon_col <- pick_col(names(acled), c("longitude", "lon", "long"), required = TRUE)
date_col <- pick_col(names(acled), c("event_date", "eventdate", "date"), required = TRUE)
event_col <- pick_col(names(acled), c("event_type"), required = TRUE)
subevent_col <- pick_col(names(acled), c("sub_event_type"), required = FALSE)

actor_cols <- intersect(
  c("actor1", "actor2", "assoc_actor_1", "assoc_actor_2", "inter1", "inter2"),
  names(acled)
)

if (!is.na(country_col)) {
  acled <- acled[norm_text(get(country_col)) == "afghanistan"]
}

acled[, latitude_num := suppressWarnings(as.numeric(get(lat_col)))]
acled[, longitude_num := suppressWarnings(as.numeric(get(lon_col)))]
acled[, month := parse_acled_date(get(date_col))]

acled <- acled[
  !is.na(latitude_num) &
    !is.na(longitude_num) &
    !is.na(month) &
    latitude_num >= -90 & latitude_num <= 90 &
    longitude_num >= -180 & longitude_num <= 180
]

cat("ACLED Afghanistan rows with usable coordinates/month:", nrow(acled), "\n")

# -------------------------
# 5) LOAD ADM2 AND SPATIALLY JOIN
# -------------------------
if (!file.exists(ADM2_PATH)) stop("Missing ADM2 shapefile: ", ADM2_PATH)

adm2 <- st_read(ADM2_PATH, quiet = TRUE)
adm2 <- st_make_valid(adm2)
adm2 <- st_transform(adm2, 4326)

adm2_id_col <- pick_col(
  names(adm2),
  c("adm2_pcode", "adm2_code", "dist_id", "district_id", "adm2_code_en"),
  required = FALSE
)

adm2_name_col <- pick_col(
  names(adm2),
  c("adm2_en", "adm2_name", "district", "district_name", "name_2"),
  required = TRUE
)

adm1_name_col <- pick_col(
  names(adm2),
  c("adm1_en", "adm1_name", "province", "province_name", "name_1"),
  required = FALSE
)

if (is.na(adm2_id_col)) {
  adm2$adm2_id_temp <- seq_len(nrow(adm2))
  adm2_id_col <- "adm2_id_temp"
}

adm2_small <- adm2[, c(adm2_id_col, adm2_name_col, adm1_name_col)]
names(adm2_small)[names(adm2_small) == adm2_id_col] <- "adm2_id_spatial"
names(adm2_small)[names(adm2_small) == adm2_name_col] <- "adm2_name_spatial"
if (!is.na(adm1_name_col)) {
  names(adm2_small)[names(adm2_small) == adm1_name_col] <- "adm1_name_spatial"
} else {
  adm2_small$adm1_name_spatial <- NA_character_
}

adm2_small$adm2_id_spatial <- as.character(adm2_small$adm2_id_spatial)

pts <- st_as_sf(
  acled,
  coords = c("longitude_num", "latitude_num"),
  crs = 4326,
  remove = FALSE
)

joined <- st_join(pts, adm2_small, join = st_within, left = TRUE)
joined_dt <- as.data.table(st_drop_geometry(joined))

joined_dt[, spatial_within_match := !is.na(adm2_id_spatial)]

# Nearest fallback for points just outside polygons
# Create fallback columns first so the code works even if every point is matched.
joined_dt[, nearest_adm2_id := NA_character_]
joined_dt[, nearest_adm2_name := NA_character_]
joined_dt[, nearest_adm1_name := NA_character_]
joined_dt[, nearest_distance_km := NA_real_]

unmatched_idx <- which(!joined_dt$spatial_within_match)

if (length(unmatched_idx) > 0) {
  nearest <- st_nearest_feature(pts[unmatched_idx, ], adm2_small)
  d_km <- as.numeric(st_distance(
    pts[unmatched_idx, ],
    adm2_small[nearest, ],
    by_element = TRUE
  )) / 1000
  
  joined_dt[unmatched_idx, nearest_adm2_id := adm2_small$adm2_id_spatial[nearest]]
  joined_dt[unmatched_idx, nearest_adm2_name := adm2_small$adm2_name_spatial[nearest]]
  joined_dt[unmatched_idx, nearest_adm1_name := adm2_small$adm1_name_spatial[nearest]]
  joined_dt[unmatched_idx, nearest_distance_km := d_km]
}

# Conservative nearest fallback: accept if within 25 km
joined_dt[, adm2_id_final := adm2_id_spatial]
joined_dt[, adm2_name_final := adm2_name_spatial]
joined_dt[, adm1_name_final := adm1_name_spatial]
joined_dt[, match_method := fifelse(spatial_within_match, "within_polygon", NA_character_)]

joined_dt[
  is.na(adm2_id_final) & !is.na(nearest_distance_km) & nearest_distance_km <= 25,
  `:=`(
    adm2_id_final = nearest_adm2_id,
    adm2_name_final = nearest_adm2_name,
    adm1_name_final = nearest_adm1_name,
    match_method = "nearest_within_25km"
  )
]

cat("Spatial within matches:", sum(joined_dt$spatial_within_match, na.rm = TRUE), "\n")
cat("Unmatched:", sum(!joined_dt$spatial_within_match, na.rm = TRUE), "\n")
cat("Match rate:", mean(joined_dt$spatial_within_match, na.rm = TRUE), "\n")

joined_dt[, final_spatial_match := !is.na(adm2_id_final)]

match_summary <- data.table(
  metric = c(
    "acled_rows_with_coordinates_month",
    "within_polygon_matches",
    "within_polygon_match_rate",
    "nearest_25km_extra_matches",
    "final_spatial_matches",
    "final_spatial_match_rate",
    "unique_adm2_matched"
  ),
  value = c(
    nrow(joined_dt),
    sum(joined_dt$spatial_within_match, na.rm = TRUE),
    round(mean(joined_dt$spatial_within_match, na.rm = TRUE), 4),
    sum(joined_dt$match_method == "nearest_within_25km", na.rm = TRUE),
    sum(joined_dt$final_spatial_match, na.rm = TRUE),
    round(mean(joined_dt$final_spatial_match, na.rm = TRUE), 4),
    uniqueN(joined_dt$adm2_id_final[!is.na(joined_dt$adm2_id_final)])
  )
)

cat("\n================ SPATIAL ACLED-ADM2 MATCH SUMMARY ================\n")
print(match_summary, row.names = FALSE)

fwrite(match_summary, file.path(OUT_DIR, "spatial_acled_adm2_match_summary.csv"))
fwrite(joined_dt, file.path(OUT_DIR, "spatial_acled_adm2_joined_rows.csv"))

# -------------------------
# 6) BUILD DISTRICT-MONTH CONTESTEDNESS FROM SPATIAL MATCH
# -------------------------
matched <- joined_dt[final_spatial_match == TRUE]

if (!length(actor_cols)) stop("No usable actor columns found in ACLED.")

matched[, actor_text := do.call(paste, c(.SD, sep = " ")), .SDcols = actor_cols]
matched[, actor_text_norm := norm_text(actor_text)]
matched[, event_type_norm := norm_text(get(event_col))]

matched[, taliban_involved := grepl("taliban", actor_text_norm)]

matched[, govt_involved := grepl(
  paste(c(
    "military forces of afghanistan",
    "police forces of afghanistan",
    "government of afghanistan",
    "afghan military",
    "afghan national army",
    "afghan national police",
    "afghan security",
    "national directorate of security",
    "nds"
  ), collapse = "|"),
  actor_text_norm
)]

matched[, is_battle := event_type_norm == "battles"]

# ACLED political violence categories used for Taliban activity.
# This is a proxy, not a direct territorial-control measure.
matched[, is_acled_polviol := event_type_norm %in% c(
  "battles",
  "explosions remote violence",
  "violence against civilians"
)]

matched[, battles_GT := as.integer(is_battle & taliban_involved & govt_involved)]
matched[, taliban_polviol := as.integer(taliban_involved & is_acled_polviol)]

dm_raw <- matched[
  ,
  .(
    battles_GT = sum(battles_GT, na.rm = TRUE),
    taliban_polviol = sum(taliban_polviol, na.rm = TRUE),
    acled_events_total = .N
  ),
  by = .(district_id = as.character(adm2_id_final), month)
]

# Complete district-month grid using districts observed in core and panel months
districts_core <- sort(unique(core$district_id_for_h2[!is.na(core$district_id_for_h2)]))
months_core <- sort(unique(core$month))

grid <- CJ(
  district_id = districts_core,
  month = months_core
)

dm <- merge(grid, dm_raw, by = c("district_id", "month"), all.x = TRUE)
for (v in c("battles_GT", "taliban_polviol", "acled_events_total")) {
  dm[is.na(get(v)), (v) := 0L]
}

dm[, contested_any := as.integer(battles_GT > 0)]
dm[, taliban_active_any := as.integer(taliban_polviol > 0)]

SHOCK_DATE <- as.IDate("2021-07-01")

pre_spatial <- dm[
  month < SHOCK_DATE,
  .(
    pre_contested_share_spatial = mean(contested_any, na.rm = TRUE),
    pre_taliban_activity_share_spatial = mean(taliban_active_any, na.rm = TRUE),
    pre_contested_intensity_spatial = log1p(mean(battles_GT, na.rm = TRUE)),
    pre_contested_total_spatial = sum(battles_GT, na.rm = TRUE),
    pre_taliban_polviol_total_spatial = sum(taliban_polviol, na.rm = TRUE),
    pre_months_observed_spatial = .N
  ),
  by = district_id
]

cat("\n================ SPATIAL H2 MODERATOR SUMMARY ================\n")
print(summary(pre_spatial$pre_contested_share_spatial))
print(summary(pre_spatial$pre_taliban_activity_share_spatial))
print(summary(pre_spatial$pre_contested_intensity_spatial))

fwrite(dm, file.path(OUT_DIR, "district_month_contestedness_spatial_adm2.csv"))
fwrite(pre_spatial, file.path(OUT_DIR, "pre_shock_h2_moderators_spatial_adm2.csv"))

# -------------------------
# 7) MERGE SPATIAL H2 MODERATORS TO CORE PANEL
# -------------------------
core_spatial <- merge(
  core,
  pre_spatial,
  by.x = "district_id_for_h2",
  by.y = "district_id",
  all.x = TRUE,
  sort = FALSE
)

merge_audit <- data.table(
  metric = c(
    "core_rows",
    "core_rows_with_existing_district_id",
    "core_rows_with_spatial_h2_moderator",
    "core_row_match_rate_to_spatial_h2",
    "unique_core_districts",
    "unique_core_districts_with_spatial_h2"
  ),
  value = c(
    nrow(core_spatial),
    sum(!is.na(core_spatial$district_id_for_h2)),
    sum(!is.na(core_spatial$pre_contested_share_spatial)),
    round(mean(!is.na(core_spatial$pre_contested_share_spatial)), 4),
    uniqueN(core_spatial$district_id_for_h2, na.rm = TRUE),
    uniqueN(core_spatial$district_id_for_h2[!is.na(core_spatial$pre_contested_share_spatial)], na.rm = TRUE)
  )
)

cat("\n================ CORE MERGE AUDIT ================\n")
print(merge_audit, row.names = FALSE)

fwrite(merge_audit, file.path(OUT_DIR, "core_merge_to_spatial_h2_audit.csv"))
saveRDS(core_spatial, file.path(OUT_DIR, "panel_core_with_spatial_h2_diagnostic.rds"))

# -------------------------
# 8) RERUN H2 DIAGNOSTIC MODELS
# -------------------------
core_spatial[, month := as.IDate(month)]
core_spatial[, cell_id := as.factor(cell_id)]

h2_spatial_results <- rbindlist(list(
  run_h2_model(core_spatial, "pre_contested_share_spatial", "Spatial ADM2 H2: pre-contested share"),
  run_h2_model(core_spatial, "pre_taliban_activity_share_spatial", "Spatial ADM2 H2: Taliban activity share"),
  run_h2_model(core_spatial, "pre_contested_intensity_spatial", "Spatial ADM2 H2: contested intensity")
), fill = TRUE)

cat("\n================ SPATIAL ADM2 H2 RESULTS ================\n")
print(h2_spatial_results, row.names = FALSE)

fwrite(h2_spatial_results, file.path(OUT_DIR, "h2_spatial_adm2_diagnostic_results.csv"))

# -------------------------
# 9) COMPARE WITH CURRENT FINAL H2 RESULTS
# -------------------------
if (file.exists(RESULTS_PATH)) {
  current <- fread(RESULTS_PATH)
  
  current_h2 <- current[
    (result_block == "afghanistan_h2_contested_share_ppml" &
       term == "post_closure:eli_feas:pre_contested_share") |
      (result_block == "afghanistan_h2_taliban_activity_share_ppml" &
         term == "post_closure:eli_feas:pre_taliban_activity_share") |
      (result_block == "afghanistan_h2_contested_intensity_ppml" &
         term == "post_closure:eli_feas:pre_contested_intensity"),
    .(
      old_result_block = result_block,
      old_term = term,
      old_estimate = estimate,
      old_std_error = std_error,
      old_p_value = p_value,
      old_nobs = nobs
    )
  ]
  
  new_h2 <- h2_spatial_results[
    component == "H2-triple",
    .(
      new_model = model,
      new_term = term,
      new_estimate = estimate,
      new_std_error = std_error,
      new_p_value = p_value,
      new_nobs = nobs
    )
  ]
  
  comparison <- data.table(
    old_label = c("pre_contested_share", "pre_taliban_activity_share", "pre_contested_intensity")
  )
  
  comparison <- cbind(comparison, current_h2, new_h2)
  
  cat("\n================ OLD VS SPATIAL-H2 COMPARISON ================\n")
  print(comparison, row.names = FALSE)
  
  fwrite(comparison, file.path(OUT_DIR, "h2_old_vs_spatial_adm2_comparison.csv"))
}

cat("\n================ DONE ================\n")
cat("Diagnostic files written to:\n", OUT_DIR, "\n")
cat("\nPaste back:\n")
cat("1) SPATIAL ACLED-ADM2 MATCH SUMMARY\n")
cat("2) CORE MERGE AUDIT\n")
cat("3) SPATIAL ADM2 H2 RESULTS\n")
cat("4) OLD VS SPATIAL-H2 COMPARISON\n")