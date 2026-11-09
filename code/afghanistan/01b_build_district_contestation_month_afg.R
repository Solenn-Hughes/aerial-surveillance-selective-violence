# ============================================================
# 01b_build_district_contestation_month_afg_EXECUTABLE_ACLED.R
# Build district-month H2 moderators from ACLED only.
# No hand-built district-control CSV is required.
#
# Main H2 variables built:
#   pre_contested_share       = share of pre-shock months with govt-Taliban battles in district
#   pre_taliban_activity_share= share of pre-shock months with any Taliban political-violence event in district
#   pre_high_contested        = top-quartile pre_contested_share (optional binary appendix)
# ============================================================

source("C:/Users/USERNAME/OneDrive/Desktop/THESIS/Proper code/super proper/00_config_geo_afghanistan_only.R")
suppressPackageStartupMessages(library(data.table))

log_msg("01b_build_district_contestation_month_afg_EXECUTABLE_ACLED: start")

xw_path <- file.path(DIRS$base, "cell_district_crosswalk.rds")
if (!file.exists(xw_path)) stop("Run 01a_build_cell_district_crosswalk_afg_EXECUTABLE_ACLED.R first.")
xw <- as.data.table(readRDS(xw_path))
district_lookup <- unique(xw[, .(district_id, district_name, province_name, district_name_norm, province_name_norm)])
assert_unique_keys(district_lookup, "district_id", "district_lookup")

acled_path <- find_first_recursive(DIRS$data_raw, RAW_ACLED_PATTERNS)
if (is.na(acled_path)) stop("No ACLED CSV found under data/raw.")
ac <- read_any_dt(acled_path)
setnames(ac, names(ac), make.unique(safe_name(names(ac)), sep = "__dup"))

month_col <- pick_col(names(ac), ACLED_DATE_CANDIDATES, required = FALSE)
dist_col  <- pick_col(names(ac), ACLED_DISTRICT_CANDIDATES, required = FALSE)
prov_col  <- pick_col(names(ac), ACLED_PROVINCE_CANDIDATES, required = FALSE)
event_col <- pick_col(names(ac), ACLED_EVENT_CANDIDATES, required = FALSE)
sub_col   <- pick_col(names(ac), ACLED_SUBEVENT_CANDIDATES, required = FALSE)

if (is.na(month_col) || is.na(dist_col) || is.na(event_col)) {
  stop("ACLED file must have event date, district/admin2, and event_type columns.")
}

ac[, month := coerce_month(get(month_col))]
ac[, district_name_raw := normalize_text(as.character(get(dist_col)))]
ac[, province_name_raw := if (!is.na(prov_col)) normalize_text(as.character(get(prov_col))) else ""]
ac[, district_name_norm := district_name_raw]
ac[, province_name_norm := province_name_raw]
ac[, event_type_std := normalize_text(get(event_col))]
ac[, sub_event_type_std := if (!is.na(sub_col)) normalize_text(get(sub_col)) else NA_character_]

# Optional manual district-name cleanup
if (file.exists(AFG_DISTRICT_NAME_CROSSWALK_PATH)) {
  xwc <- fread(AFG_DISTRICT_NAME_CROSSWALK_PATH)
  setnames(xwc, names(xwc), make.unique(safe_name(names(xwc)), sep = "__dup"))
  req <- c("district_name_raw", "district_name_std")
  if (all(req %in% names(xwc))) {
    xwc[, district_name_raw := normalize_text(district_name_raw)]
    xwc[, district_name_std := normalize_text(district_name_std)]
    if ("province_name_raw" %in% names(xwc)) xwc[, province_name_raw := normalize_text(province_name_raw)]
    if ("province_name_std" %in% names(xwc)) xwc[, province_name_std := normalize_text(province_name_std)]
    if (!"province_name_raw" %in% names(xwc)) xwc[, province_name_raw := ""]
    ac <- merge(
      ac,
      xwc[, .(district_name_raw, province_name_raw, district_name_std, province_name_std)],
      by = c("district_name_raw", "province_name_raw"),
      all.x = TRUE,
      sort = FALSE
    )
    ac[!is.na(district_name_std), district_name_norm := district_name_std]
    if ("province_name_std" %in% names(ac)) ac[!is.na(province_name_std), province_name_norm := province_name_std]
  }
}

actor_cols <- intersect(c("actor1", "actor2", "assoc_actor_1", "assoc_actor_2"), names(ac))
if (!length(actor_cols)) stop("ACLED file needs actor columns (actor1/actor2 at minimum).")
for (cc in actor_cols) ac[, (cc) := normalize_text(get(cc))]

combine_actor_text <- function(dt) {
  txt <- rep("", nrow(dt))
  for (cc in actor_cols) {
    v <- dt[[cc]]
    v[is.na(v)] <- ""
    txt <- paste(txt, v)
  }
  normalize_text(txt)
}
ac[, actor_text := combine_actor_text(.SD), .SDcols = actor_cols]

is_taliban <- function(x) grepl("\\btaliban\\b|\\bhaqqani\\b", x)
is_gov_afg <- function(x) grepl("government of afghanistan|military forces of afghanistan|police forces of afghanistan|national directorate of security|pro government militia", x)

ac[, taliban_involved := is_taliban(actor_text)]
ac[, govt_involved := is_gov_afg(actor_text)]
ac[, is_battle_gt := as.integer(event_type_std == "battles" & taliban_involved & govt_involved)]
ac[, is_taliban_polviol := as.integer(taliban_involved & event_type_std %in% c("battles", "violence against civilians", "explosions remote violence"))]

# Restrict to Afghanistan and valid matched districts
country_col <- pick_col(names(ac), c("country", "country_name"), required = FALSE)
if (!is.na(country_col)) ac <- ac[normalize_text(get(country_col)) == normalize_text(TARGET_THEATRE)]

ac <- ac[!is.na(month) & !is.na(district_name_norm) & nzchar(district_name_norm)]
ac[, ac_rowid := .I]
ac[, province_name_norm := fifelse(is.na(province_name_norm), "", province_name_norm)]

# First pass: exact district + province match
ac_m_prov <- merge(
  ac,
  district_lookup,
  by = c("district_name_norm", "province_name_norm"),
  all.x = FALSE,
  sort = FALSE
)
if (nrow(ac_m_prov) > 0L) ac_m_prov[, match_source := "district_plus_province"]

matched_ids_prov <- if (nrow(ac_m_prov) > 0L) unique(ac_m_prov$ac_rowid) else integer()
ac_unmatched <- ac[!.(matched_ids_prov), on = "ac_rowid"]

# Second pass: district-only fallback for still-unmatched rows.
# Keep only unique district-name matches to avoid silently assigning ambiguous names.
dl2 <- unique(district_lookup[, .(district_id, district_name, province_name, district_name_norm)])
ac_m_dist_raw <- merge(
  ac_unmatched,
  dl2,
  by = "district_name_norm",
  all.x = FALSE,
  sort = FALSE,
  allow.cartesian = TRUE
)

ambig_dt <- if (nrow(ac_m_dist_raw) > 0L) ac_m_dist_raw[, .N, by = ac_rowid] else data.table(ac_rowid = integer(), N = integer())
unique_dist_ids <- ambig_dt[N == 1L, ac_rowid]
ambiguous_dist_ids <- ambig_dt[N > 1L, ac_rowid]

ac_m_dist <- ac_m_dist_raw[ac_rowid %in% unique_dist_ids]
if (nrow(ac_m_dist) > 0L) ac_m_dist[, match_source := "district_only_unique"]

ac_m <- rbindlist(list(ac_m_prov, ac_m_dist), use.names = TRUE, fill = TRUE)

if (nrow(ac_m) == 0L) stop("Could not match ACLED admin names to Afghanistan district boundaries. Add a crosswalk CSV.")

# Match audit outputs
matched_ids_all <- unique(ac_m$ac_rowid)
ac_unmatched_final <- ac[!.(matched_ids_all), on = "ac_rowid"]

match_summary <- data.table(
  metric = c(
    "acled_rows_afghanistan_after_basic_cleaning",
    "matched_rows_district_plus_province",
    "matched_rows_district_only_unique",
    "matched_rows_total",
    "unmatched_rows_total",
    "ambiguous_rows_district_only",
    "match_rate_rows_total",
    "unique_admin_pairs_input",
    "unique_admin_pairs_matched",
    "unique_admin_pairs_unmatched"
  ),
  value = c(
    nrow(ac),
    nrow(ac_m_prov),
    nrow(ac_m_dist),
    length(matched_ids_all),
    nrow(ac_unmatched_final),
    length(unique(ambiguous_dist_ids)),
    if (nrow(ac) > 0L) round(length(matched_ids_all) / nrow(ac), 4) else NA_real_,
    uniqueN(ac[, .(district_name_norm, province_name_norm)]),
    uniqueN(ac_m[, .(district_name_norm, province_name_norm)]),
    uniqueN(ac_unmatched_final[, .(district_name_norm, province_name_norm)])
  )
)
fwrite(match_summary, file.path(DIRS$audit, "district_name_match_summary.csv"))

if (nrow(ac_unmatched_final) > 0L) {
  unmatched_names <- ac_unmatched_final[, .(
    n_rows = .N,
    n_months = uniqueN(month)
  ), by = .(district_name_norm, province_name_norm)][order(-n_rows, district_name_norm, province_name_norm)]
  fwrite(unmatched_names, file.path(DIRS$audit, "district_name_match_unmatched.csv"))
}

if (length(ambiguous_dist_ids) > 0L) {
  ambiguous_names <- ac_m_dist_raw[ac_rowid %in% ambiguous_dist_ids, .(
    district_name_norm,
    province_name_norm_input = province_name_norm,
    matched_province_name = province_name,
    district_id
  )][order(district_name_norm, province_name_norm_input, matched_province_name, district_id)]
  fwrite(ambiguous_names, file.path(DIRS$audit, "district_name_match_ambiguous.csv"))
}

# Build district-month panel over the same month support as the backbone
base <- read_any_dt(FOUNDATION_PATH)
setnames(base, names(base), make.unique(safe_name(names(base)), sep = "__dup"))
month_base_col <- pick_col(names(base), c("month", "month_id", "date"))
theatre_col <- pick_col(names(base), c("theatre_id", "country", "country_name"))
base[, month := coerce_month(get(month_base_col))]
base[, theatre_id := clean_theatre(get(theatre_col))]
months_afg <- sort(unique(base[theatre_id == TARGET_THEATRE, month]))

all_district_months <- CJ(district_id = district_lookup$district_id, month = months_afg, unique = TRUE)
agg <- ac_m[, .(
  district_battles_gt = sum(is_battle_gt, na.rm = TRUE),
  district_taliban_polviol = sum(is_taliban_polviol, na.rm = TRUE)
), by = .(district_id, month)]

out <- merge(all_district_months, agg, by = c("district_id", "month"), all.x = TRUE, sort = TRUE)
out[is.na(district_battles_gt), district_battles_gt := 0L]
out[is.na(district_taliban_polviol), district_taliban_polviol := 0L]
out[, contested_month := as.integer(district_battles_gt > 0L)]
out[, taliban_activity_month := as.integer(district_taliban_polviol > 0L)]
out <- merge(out, district_lookup[, .(district_id, district_name, province_name)], by = "district_id", all.x = TRUE, sort = FALSE)
assert_unique_keys(out, c("district_id", "month"), "district_contestation_month")

saveRDS(out, file.path(DIRS$base, "district_contestation_month.rds"))
fwrite(out, file.path(DIRS$audit, "district_contestation_month_preview.csv"))

summ <- out[, .(
  district_months = .N,
  share_contested = mean(contested_month),
  share_taliban_activity = mean(taliban_activity_month),
  mean_battles_gt = mean(district_battles_gt),
  mean_taliban_polviol = mean(district_taliban_polviol)
)]
fwrite(summ, file.path(DIRS$audit, "district_contestation_month_summary.csv"))

write_audit_text(c(
  "District contestation build summary (ACLED-based)",
  paste0("ACLED source: ", acled_path),
  paste0("District-month rows: ", nrow(out)),
  paste0("Districts represented: ", uniqueN(out$district_id)),
  paste0("Matched ACLED rows (all methods): ", length(matched_ids_all)),
  paste0("Unmatched ACLED rows: ", nrow(ac_unmatched_final)),
  paste0("District+province matches: ", nrow(ac_m_prov)),
  paste0("District-only unique fallback matches: ", nrow(ac_m_dist)),
  "Main H2 proxies are now built from ACLED district-month contestedness and Taliban political-violence activity.",
  "District-name matching is audited explicitly in district_name_match_summary.csv.",
  "No hand-built district-control CSV is required in this executable branch."
), "district_contestation_build_summary")

log_msg("01b_build_district_contestation_month_afg_EXECUTABLE_ACLED: done")
