# ============================================================
# 01a_build_cell_district_crosswalk_afg_EXECUTABLE_ACLED.R
# Map each PRIO cell centroid to an Afghanistan ADM2 district.
# ============================================================

source("C:/Users/USERNAME/OneDrive/Desktop/THESIS/Proper code/super proper/00_config_geo_afghanistan_only.R")
suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})

log_msg("01a_build_cell_district_crosswalk_afg_EXECUTABLE_ACLED: start")

if (!file.exists(FOUNDATION_PATH)) stop("Missing FOUNDATION_PATH: ", FOUNDATION_PATH)
if (!file.exists(AFG_DISTRICT_BOUNDARIES_PATH)) stop("Missing district boundaries: ", AFG_DISTRICT_BOUNDARIES_PATH)

base <- read_any_dt(FOUNDATION_PATH)
setnames(base, names(base), make.unique(safe_name(names(base)), sep = "__dup"))

cell_col    <- pick_col(names(base), c("cell_id", "pg_id", "gid"))
theatre_col <- pick_col(names(base), c("theatre_id", "country", "country_name"))
eli_col     <- pick_col(names(base), c("eli_feas", "eli_feasibility", "eli_access", "eli_geo_access"))

base[, cell_id := as.integer(get(cell_col))]
base[, theatre_id := clean_theatre(get(theatre_col))]
base[, eli_feas := as.numeric(get(eli_col))]

if (all(c("lat", "lon") %in% names(base))) {
  base[, `:=`(lat = as.numeric(lat), lon = as.numeric(lon))]
} else if (all(c("ycoord", "xcoord") %in% names(base))) {
  base[, `:=`(lat = as.numeric(ycoord), lon = as.numeric(xcoord))]
} else {
  stop("Need lat/lon or ycoord/xcoord in backbone to build district crosswalk.")
}

cells <- unique(base[theatre_id == TARGET_THEATRE & !is.na(cell_id) & !is.na(lat) & !is.na(lon),
                     .(cell_id, lat, lon)])
if (nrow(cells) == 0L) stop("No Afghanistan cells with coordinates found in backbone.")

cells_sf <- st_as_sf(cells, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

dist_sf <- suppressWarnings(st_read(AFG_DISTRICT_BOUNDARIES_PATH, quiet = TRUE))
if (is.na(st_crs(dist_sf))) st_crs(dist_sf) <- 4326
dist_sf <- st_transform(dist_sf, 4326)

dist_names <- names(dist_sf)
id_col <- pick_col(dist_names, DISTRICT_ID_CANDIDATES, required = FALSE)
dist_col <- pick_col(dist_names, DISTRICT_NAME_CANDIDATES, required = FALSE)
prov_col <- pick_col(dist_names, PROVINCE_NAME_CANDIDATES, required = FALSE)

if (is.na(id_col)) {
  dist_sf$district_id <- as.character(seq_len(nrow(dist_sf)))
} else {
  dist_sf$district_id <- as.character(dist_sf[[id_col]])
}
if (is.na(dist_col)) stop("Could not detect district-name column in boundaries.")
if (is.na(prov_col)) {
  dist_sf$province_name <- NA_character_
} else {
  dist_sf$province_name <- as.character(dist_sf[[prov_col]])
}
dist_sf$district_name <- as.character(dist_sf[[dist_col]])

# Try within, then intersects, then nearest fallback for unmatched cells.
join1 <- suppressWarnings(st_join(cells_sf, dist_sf[, c("district_id", "district_name", "province_name")], join = st_within, left = TRUE))
join1$crosswalk_method <- ifelse(is.na(join1$district_name), NA_character_, "within")

unmatched1 <- is.na(join1$district_name)
if (any(unmatched1)) {
  join_intersects <- suppressWarnings(
    st_drop_geometry(
      st_join(cells_sf[unmatched1, ], dist_sf[, c("district_id", "district_name", "province_name")], join = st_intersects, left = TRUE)
    )[, c("district_id", "district_name", "province_name")]
  )
  if (nrow(join_intersects) > 0L) {
    join1[unmatched1, c("district_id", "district_name", "province_name")] <- join_intersects
    matched_intersects <- unmatched1 & !is.na(join1$district_name)
    join1$crosswalk_method[matched_intersects] <- "intersects"
  }
}

unmatched2 <- is.na(join1$district_name)
if (any(unmatched2)) {
  nearest_idx <- st_nearest_feature(cells_sf[unmatched2, ], dist_sf)
  join1$district_id[unmatched2] <- dist_sf$district_id[nearest_idx]
  join1$district_name[unmatched2] <- dist_sf$district_name[nearest_idx]
  join1$province_name[unmatched2] <- dist_sf$province_name[nearest_idx]
  join1$crosswalk_method[unmatched2] <- "nearest"
}

join1$crosswalk_method[is.na(join1$crosswalk_method) | join1$crosswalk_method == ""] <- "within"

out <- as.data.table(st_drop_geometry(join1))[
  , .(
    cell_id = as.integer(cell_id),
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    district_id = as.character(district_id),
    district_name = as.character(district_name),
    province_name = as.character(province_name),
    district_name_norm = normalize_text(district_name),
    province_name_norm = normalize_text(province_name),
    crosswalk_method = as.character(crosswalk_method)
  )
]
assert_unique_keys(out, "cell_id", "cell_district_crosswalk")

saveRDS(out, file.path(DIRS$base, "cell_district_crosswalk.rds"))
fwrite(out, file.path(DIRS$audit, "cell_district_crosswalk_preview.csv"))

audit <- out[, .(
  n_cells = .N,
  n_districts = uniqueN(district_id),
  n_nearest_fallback = sum(crosswalk_method == "nearest"),
  share_nearest_fallback = mean(crosswalk_method == "nearest")
)]
fwrite(audit, file.path(DIRS$audit, "cell_district_crosswalk_audit.csv"))

write_audit_text(c(
  "Cell-to-district crosswalk summary",
  paste0("Boundary file: ", AFG_DISTRICT_BOUNDARIES_PATH),
  paste0("Cells mapped: ", nrow(out)),
  paste0("Districts represented: ", uniqueN(out$district_id)),
  paste0("Nearest-neighbour fallback cells: ", sum(out$crosswalk_method == 'nearest'))
), "cell_district_crosswalk_summary")

log_msg("01a_build_cell_district_crosswalk_afg_EXECUTABLE_ACLED: done")
