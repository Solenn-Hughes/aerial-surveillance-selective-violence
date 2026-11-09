# ============================================================
# 01_build_foundation_geo_afghanistan_only_EXECUTABLE_ACLED.R
# Foundation panel with:
# - ACLED district-contestedness H2 moderators
# - legacy H2 proxy retained for appendix
# ============================================================

source("C:/Users/USERNAME/OneDrive/Desktop/THESIS/Proper code/super proper/00_config_geo_afghanistan_only.R")
suppressPackageStartupMessages(library(data.table))

log_msg("01_build_foundation_geo_afghanistan_only_EXECUTABLE_ACLED: start")

base <- read_any_dt(FOUNDATION_PATH)
setnames(base, names(base), make.unique(safe_name(names(base)), sep = "__dup"))

cell_col    <- pick_col(names(base), c("cell_id", "pg_id", "gid"))
month_col   <- pick_col(names(base), c("month", "month_id", "date"))
theatre_col <- pick_col(names(base), c("theatre_id", "country", "country_name"))
eli_col     <- pick_col(names(base), c("eli_feas", "eli_feasibility", "eli_access", "eli_geo_access"))

base[, cell_id := as.integer(get(cell_col))]
base[, month := coerce_month(get(month_col))]
base[, theatre_id := clean_theatre(get(theatre_col))]
base[, eli_feas := as.numeric(get(eli_col))]

if (all(c("lat", "lon") %in% names(base))) {
  base[, `:=`(lat = as.numeric(lat), lon = as.numeric(lon))]
} else if (all(c("ycoord", "xcoord") %in% names(base))) {
  base[, `:=`(lat = as.numeric(ycoord), lon = as.numeric(xcoord))]
} else {
  base[, `:=`(lat = NA_real_, lon = NA_real_)]
}

keep_cols <- unique(c(
  "cell_id", "month", "theatre_id", "eli_feas", "lat", "lon",
  intersect(c(
    "prec_gpcp", "temp", "nlights_calib_mean",
    "water_ih", "forest_ih", "pop_gpw_sum",
    "control_weak", "state_control",
    "contestation_tv", "contestedness_cts",
    "contest_spatial_tv", "contestedness_spatial",
    "rebel_cap_tv", "rebel_capability",
    "eli_feas_terrain", "min_km_airbase", "eli_feas_range"
  ), names(base))
))

foundation <- copy(base[, ..keep_cols])
foundation <- foundation[!is.na(cell_id) & !is.na(month) & !is.na(theatre_id) & !is.na(eli_feas)]
foundation <- foundation[theatre_id == TARGET_THEATRE]
foundation <- unique(foundation)

if ("control_weak" %in% names(foundation) && "state_control" %notin% names(foundation)) {
  foundation[, state_control := 1 - as.numeric(control_weak)]
}
if ("contestation_tv" %in% names(foundation) && "contestedness_cts" %notin% names(foundation)) {
  foundation[, contestedness_cts := as.numeric(contestation_tv)]
}
assert_unique_keys(foundation, c("cell_id", "month"), "foundation_raw")

xw_path <- file.path(DIRS$base, "cell_district_crosswalk.rds")
ctrl_path <- file.path(DIRS$base, "district_contestation_month.rds")
if (!file.exists(xw_path)) stop("Run 01a_build_cell_district_crosswalk_afg_EXECUTABLE_ACLED.R first.")
if (!file.exists(ctrl_path)) stop("Run 01b_build_district_contestation_month_afg_EXECUTABLE_ACLED.R first.")

xw <- as.data.table(readRDS(xw_path))
ctrl <- as.data.table(readRDS(ctrl_path))
assert_unique_keys(xw, "cell_id", "cell_district_crosswalk")
assert_unique_keys(ctrl, c("district_id", "month"), "district_contestation_month")

foundation <- merge(foundation, xw[, .(cell_id, district_id, district_name, province_name)], by = "cell_id", all.x = TRUE, sort = FALSE)
foundation <- merge(foundation,
                    ctrl[, .(district_id, month, district_battles_gt, district_taliban_polviol, contested_month, taliban_activity_month)],
                    by = c("district_id", "month"), all.x = TRUE, sort = FALSE)
for (cc in c("district_battles_gt", "district_taliban_polviol", "contested_month", "taliban_activity_month")) {
  foundation[is.na(get(cc)), (cc) := 0]
}

foundation[, month_int := month_int_from_idate(month)]
foundation[, shock_date := SHOCK_DATE]
foundation[, shock_var := SHOCK_VAR]
foundation[, shock_label := SHOCK_LABEL]
foundation[, shock_bagram := as.integer(month >= SHOCK_DATE)]
foundation[, shock_active := shock_bagram]
foundation[, post_closure := as.integer(month >= SHOCK_DATE)]

foundation[, eli_high_theatre := as.integer(eli_feas >= median(eli_feas, na.rm = TRUE))]
foundation[, eli_std_theatre := {
  m <- mean(eli_feas, na.rm = TRUE)
  s <- stats::sd(eli_feas, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, .N) else (eli_feas - m) / s
}]
foundation[, eli_q_theatre := as.integer(pmin(4L, pmax(1L, ceiling(4 * frank(eli_feas, ties.method = "average") / .N))))]

pre_base <- foundation[month < SHOCK_DATE]

# Main executable H2 variables
pre_h2 <- pre_base[, .(
  pre_district_obs_n = .N,
  pre_contested_share = mean(contested_month, na.rm = TRUE),
  pre_taliban_activity_share = mean(taliban_activity_month, na.rm = TRUE),
  pre_contested_intensity = mean(log1p(district_battles_gt), na.rm = TRUE),
  pre_taliban_activity_intensity = mean(log1p(district_taliban_polviol), na.rm = TRUE)
), by = cell_id]

# Optional binary version for appendix
cut_ct <- quantile(pre_h2$pre_contested_share, 0.75, na.rm = TRUE)
pre_h2[, pre_high_contested := as.integer(pre_contested_share >= cut_ct)]

# Legacy H2 retained for comparison
pre_legacy <- pre_base[, .(
  pre_obs_n = .N,
  pre_state_control_mean = if ("state_control" %in% names(pre_base)) mean(state_control, na.rm = TRUE) else NA_real_,
  pre_contestedness_mean = if ("contestedness_cts" %in% names(pre_base)) mean(contestedness_cts, na.rm = TRUE) else NA_real_
), by = cell_id]
if (!all(is.na(pre_legacy$pre_state_control_mean))) {
  cut_sc <- quantile(pre_legacy$pre_state_control_mean, 0.75, na.rm = TRUE)
  pre_legacy[, pre_state_control_high := as.integer(pre_state_control_mean >= cut_sc)]
} else pre_legacy[, pre_state_control_high := 0L]
if (!all(is.na(pre_legacy$pre_contestedness_mean))) {
  cut_c2 <- quantile(pre_legacy$pre_contestedness_mean, 0.75, na.rm = TRUE)
  pre_legacy[, pre_contested_high := as.integer(pre_contestedness_mean >= cut_c2)]
} else pre_legacy[, pre_contested_high := 0L]
pre_legacy[, pre_low_state_control := as.integer(1L - pmin(1L, pmax(0L, pre_state_control_high)))]
pre_legacy[, pre_weak_control_or_contested_legacy := as.integer(pre_low_state_control == 1L | pre_contested_high == 1L)]

foundation <- merge(foundation, pre_h2, by = "cell_id", all.x = TRUE, sort = FALSE)
foundation <- merge(foundation, pre_legacy, by = "cell_id", all.x = TRUE, sort = FALSE)
for (cc in c("pre_contested_share", "pre_taliban_activity_share", "pre_contested_intensity", "pre_taliban_activity_intensity")) {
  foundation[is.na(get(cc)), (cc) := 0]
}
foundation[is.na(pre_high_contested), pre_high_contested := 0L]
for (cc in c("pre_state_control_high", "pre_contested_high", "pre_low_state_control", "pre_weak_control_or_contested_legacy")) {
  foundation[is.na(get(cc)), (cc) := 0L]
}

foundation[, low_visibility_q25 := 0L]
foundation[, low_visibility_q10 := 0L]
foundation[, dist_mainshock_km := fifelse(!is.na(lat) & !is.na(lon), haversine_km(lat, lon, SHOCK_SITE_LAT, SHOCK_SITE_LON), NA_real_)]

assert_unique_keys(foundation, c("cell_id", "month"), "panel_foundation_cell_month")
saveRDS(foundation, file.path(DIRS$panels, "panel_foundation_cell_month.rds"))

fwrite(foundation[month < SHOCK_DATE, .(
  cell_id,
  pre_contested_share = unique(pre_contested_share)[1],
  pre_taliban_activity_share = unique(pre_taliban_activity_share)[1],
  pre_contested_intensity = unique(pre_contested_intensity)[1],
  pre_high_contested = unique(pre_high_contested)[1],
  legacy_h2 = unique(pre_weak_control_or_contested_legacy)[1]
)], file.path(DIRS$audit, "foundation_h2_moderators_cell_audit.csv"))

write_audit_text(c(
  "Foundation build summary (executable ACLED branch)",
  paste0("Backbone used: ", FOUNDATION_PATH),
  paste0("District crosswalk used: ", xw_path),
  paste0("District contestation used: ", ctrl_path),
  paste0("Rows: ", nrow(foundation)),
  paste0("Cells: ", uniqueN(foundation$cell_id)),
  paste0("Months: ", uniqueN(foundation$month)),
  "Main H2 variables created: pre_contested_share, pre_taliban_activity_share, pre_contested_intensity, pre_taliban_activity_intensity",
  "Legacy H2 retained as pre_weak_control_or_contested_legacy"
), "foundation_build_summary_executable_acled")

log_msg("01_build_foundation_geo_afghanistan_only_EXECUTABLE_ACLED: done")
