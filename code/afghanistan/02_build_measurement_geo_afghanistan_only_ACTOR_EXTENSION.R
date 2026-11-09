# ============================================================
# 02_build_measurement_geo_afghanistan_only_EXECUTABLE_ACLED.R
# Builds:
# - core GED outcomes
# - ACLED protests / riots negative-control outcomes
# - low-visibility geography indicators
# - H3 local-cluster and global ring spillover variables
# ============================================================

source("C:/Users/USERNAME/OneDrive/Desktop/THESIS/Proper code/super proper/00_config_geo_afghanistan_only.R")
suppressPackageStartupMessages({
  library(data.table)
})

log_msg("02_build_measurement_geo_afghanistan_only_EXECUTABLE_ACLED: start")

foundation_path <- file.path(DIRS$panels, "panel_foundation_cell_month.rds")
if (!file.exists(foundation_path)) stop("Run 01_build_foundation_geo_afghanistan_only_EXECUTABLE_ACLED.R first.")
foundation <- as.data.table(readRDS(foundation_path))
assert_unique_keys(foundation, c("cell_id", "month"), "foundation")

meta_cols <- intersect(c(
  "cell_id", "theatre_id", "month", "month_int",
  "eli_feas", "eli_high_theatre", "eli_std_theatre", "eli_q_theatre",
  "lat", "lon", "eli_feas_terrain", "dist_mainshock_km",
  "pre_contested_share", "pre_taliban_activity_share",
  "pre_contested_intensity", "pre_taliban_activity_intensity",
  "pre_weak_control_or_contested_legacy",
  "prec_gpcp", "temp", "nlights_calib_mean",
  "district_id", "district_name", "province_name",
  "post_closure", "shock_bagram"
), names(foundation))
cell_map <- unique(foundation[, ..meta_cols])
cell_lookup <- unique(cell_map[, .(cell_id, month, month_int, theatre_id)])
assert_unique_keys(cell_map, c("cell_id", "month"), "cell_map")

# ---- GED outcomes ----
ged_path <- find_first_recursive(DIRS$data_raw, RAW_GED_PATTERNS)
if (is.na(ged_path)) stop("No raw GED/UCDP CSV found.")
log_msg("Using GED source: ", ged_path)

ged <- read_any_dt(ged_path)
setnames(ged, names(ged), make.unique(safe_name(names(ged)), sep = "__dup"))

country_col <- pick_col(names(ged), c("country", "country_name", "location"), required = FALSE)
month_col   <- pick_col(names(ged), c("date_start", "date_precise", "date", "month_start", "start_date"), required = FALSE)
lat_col     <- pick_col(names(ged), c("latitude", "lat"), required = FALSE)
lon_col     <- pick_col(names(ged), c("longitude", "lon"), required = FALSE)
pg_col      <- pick_col(names(ged), c("priogrid_gid", "priogrid_id", "pg_id", "gid"), required = FALSE)
viol_col    <- pick_col(names(ged), c("type_of_violence", "type_of_violence_1"), required = TRUE)
fat_col     <- pick_col(names(ged), c("best", "fat_best", "deaths_best", "fatalities"), required = FALSE)

if (!is.na(pg_col)) {
  ged[, cell_id := as.integer(get(pg_col))]
} else if (!is.na(lat_col) && !is.na(lon_col)) {
  ged[, cell_id := compute_pg_id(as.numeric(get(lat_col)), as.numeric(get(lon_col)))]
} else {
  stop("GED needs either priogrid_gid/pg_id or latitude+longitude.")
}

if (!is.na(month_col)) {
  ged[, month := coerce_month(get(month_col))]
} else if (all(c("year", "month") %in% names(ged))) {
  ged[, month := as.IDate(sprintf("%04d-%02d-01", as.integer(year), as.integer(month)))]
} else {
  stop("GED needs a date or month column.")
}

if (!is.na(country_col)) ged[, theatre_id_raw := clean_theatre(get(country_col))]
ged[, type_of_violence := suppressWarnings(as.integer(get(viol_col)))]
ged[, fat_best := if (!is.na(fat_col)) as.numeric(get(fat_col)) else 0]

ged <- ged[cell_id %in% cell_map$cell_id & month %in% cell_map$month]
ged <- merge(ged, cell_lookup, by = c("cell_id", "month"), all.x = FALSE, sort = FALSE)
ged <- ged[theatre_id == TARGET_THEATRE]
saveRDS(ged, file.path(DIRS$base, "ged_events_pg_geo_afghanistan_causal.rds"))

core_agg_type3 <- ged[type_of_violence == 3L, .(
  events_type3 = .N,
  fat_best_type3 = sum(fat_best, na.rm = TRUE),
  any_type3 = as.integer(.N > 0L)
), by = .(cell_id, month)]
core_agg_type1 <- ged[type_of_violence == 1L, .(
  events_type1 = .N,
  fat_best_type1 = sum(fat_best, na.rm = TRUE),
  any_type1 = as.integer(.N > 0L)
), by = .(cell_id, month)]

core_panel <- merge(cell_map, core_agg_type3, by = c("cell_id", "month"), all.x = TRUE, sort = FALSE)
core_panel <- merge(core_panel, core_agg_type1, by = c("cell_id", "month"), all.x = TRUE, sort = FALSE)
for (v in c("events_type3", "fat_best_type3", "any_type3", "events_type1", "fat_best_type1", "any_type1")) {
  core_panel[is.na(get(v)), (v) := 0]
}
core_panel[, theatre_month_events_type3 := sum(events_type3, na.rm = TRUE), by = month]
core_panel[, share_type3_theatre_month := fifelse(theatre_month_events_type3 > 0, events_type3 / theatre_month_events_type3, 0)]

# ---- ACLED negative-control outcomes ----
acled_path <- find_first_recursive(DIRS$data_raw, RAW_ACLED_PATTERNS)
if (is.na(acled_path)) stop("No raw ACLED CSV found using RAW_ACLED_PATTERNS.")
log_msg("Using ACLED source: ", acled_path)

ac <- read_any_dt(acled_path)
setnames(ac, names(ac), make.unique(safe_name(names(ac)), sep = "__dup"))

ac_country_col <- pick_col(names(ac), c("country", "country_name"), required = FALSE)
ac_date_col    <- pick_col(names(ac), c("event_date", "event_date_std", "date"), required = TRUE)
ac_lat_col     <- pick_col(names(ac), c("latitude", "lat"), required = TRUE)
ac_lon_col     <- pick_col(names(ac), c("longitude", "lon"), required = TRUE)
ac_event_col   <- pick_col(names(ac), c("event_type", "event_type_std"), required = TRUE)
ac_sub_col     <- pick_col(names(ac), c("sub_event_type", "sub_event_type_std"), required = FALSE)

ac[, month := coerce_month(get(ac_date_col))]
ac[, cell_id := compute_pg_id(as.numeric(get(ac_lat_col)), as.numeric(get(ac_lon_col)))]
if (!is.na(ac_country_col)) ac[, theatre_id_raw := clean_theatre(get(ac_country_col))]
ac[, event_type_clean := normalize_text(get(ac_event_col))]
ac[, sub_event_type_clean := if (!is.na(ac_sub_col)) normalize_text(get(ac_sub_col)) else NA_character_]

ac <- ac[cell_id %in% cell_map$cell_id & month %in% cell_map$month]
ac <- merge(ac, cell_lookup, by = c("cell_id", "month"), all.x = FALSE, sort = FALSE)
ac <- ac[theatre_id == TARGET_THEATRE]
saveRDS(ac, file.path(DIRS$base, "acled_events_pg_geo_afghanistan_causal.rds"))

# ---- Actor / style extensions for mechanism adjudication ----
parse_interaction_label <- function(x) {
  x <- suppressWarnings(as.integer(x))
  out <- rep(NA_character_, length(x))
  out[x == 1L] <- "state"
  out[x == 2L] <- "rebel"
  out[x == 3L] <- "political_militia"
  out[x == 4L] <- "identity_militia"
  out[x == 5L] <- "rioters"
  out[x == 6L] <- "protesters"
  out[x == 7L] <- "civilians"
  out[x == 8L] <- "other"
  out
}

is_rebel_like <- function(x) {
  x <- normalize_text(x)
  x %chin% c("rebel", "political militia", "identity militia")
}

have_actor_inputs <- all(c("inter1", "actor1") %in% names(ac))
actor_panel <- NULL
style_panel <- NULL
actor_style_panel <- NULL
actor_audit <- data.table(
  component = c("actor_panel", "style_panel", "actor_style_panel"),
  built = 0L,
  reason = "Missing ACLED actor/interaction columns"
)

if (have_actor_inputs) {
  ac[, inter1_norm := parse_interaction_label(get("inter1"))]
  if ("inter2" %in% names(ac)) ac[, inter2_norm := parse_interaction_label(get("inter2"))] else ac[, inter2_norm := NA_character_]
  ac[, actor1_clean := trimws(as.character(get("actor1")))]

  ac[, is_vac := grepl("violence against civilians", event_type_clean, ignore.case = TRUE)]
  ac[, is_remote := grepl("explosions|remote violence", event_type_clean, ignore.case = TRUE)]
  ac[, civilian_target := inter1_norm == "civilians" | inter2_norm == "civilians"]
  ac[, state_involved := inter1_norm == "state" | inter2_norm == "state"]
  ac[, rebel_involved := is_rebel_like(inter1_norm) | is_rebel_like(inter2_norm) | grepl("^unidentified armed group", normalize_text(actor1_clean))]
  ac[, actor_type_primary := fifelse(state_involved & !rebel_involved, "state",
                              fifelse(rebel_involved, "rebel", NA_character_))]

  ac[, state_civ_harm := state_involved & (is_vac | (is_remote & civilian_target))]
  ac[, rebel_civ_harm := rebel_involved & (is_vac | (is_remote & civilian_target))]

  sub_lower <- normalize_text(ac$sub_event_type_clean)
  ac[, selective_flag := (state_civ_harm | rebel_civ_harm) & is_vac & sub_lower %chin% c(
    "attack", "abduction forced disappearance", "sexual violence"
  )]
  ac[, indiscriminate_flag := (state_civ_harm | rebel_civ_harm) & (
    sub_lower %chin% c(
      "air drone strike", "shelling artillery missile attack",
      "remote explosive landmine ied", "grenade", "chemical weapon", "suicide bomb"
    ) |
      (is_vac & sub_lower %chin% c("shelling artillery missile attack", "air drone strike"))
  )]

  actor_long <- rbindlist(list(
    ac[state_civ_harm == TRUE, .(cell_id, month, actor_type = "state", vac_actor = 1L)],
    ac[rebel_civ_harm == TRUE, .(cell_id, month, actor_type = "rebel", vac_actor = 1L)]
  ), use.names = TRUE, fill = TRUE)
  actor_agg <- actor_long[, .(vac_actor = .N), by = .(cell_id, month, actor_type)]

  actor_support <- rbindlist(list(
    copy(cell_map)[, actor_type := "state"],
    copy(cell_map)[, actor_type := "rebel"]
  ), use.names = TRUE, fill = TRUE)
  actor_panel <- merge(actor_support, actor_agg, by = c("cell_id", "month", "actor_type"), all.x = TRUE, sort = FALSE)
  actor_panel[is.na(vac_actor), vac_actor := 0L]
  assert_unique_keys(actor_panel, c("cell_id", "month", "actor_type"), "actor_panel")
  saveRDS(actor_panel, file.path(DIRS$panels, "panel_actor_cell_actor_month.rds"))

  style_long <- rbindlist(list(
    ac[selective_flag == TRUE, .(cell_id, month, style = "selective", vac_style = 1L)],
    ac[indiscriminate_flag == TRUE, .(cell_id, month, style = "indiscriminate", vac_style = 1L)]
  ), use.names = TRUE, fill = TRUE)
  style_agg <- style_long[, .(vac_style = .N), by = .(cell_id, month, style)]

  style_support <- rbindlist(list(
    copy(cell_map)[, style := "selective"],
    copy(cell_map)[, style := "indiscriminate"]
  ), use.names = TRUE, fill = TRUE)
  style_panel <- merge(style_support, style_agg, by = c("cell_id", "month", "style"), all.x = TRUE, sort = FALSE)
  style_panel[is.na(vac_style), vac_style := 0L]
  assert_unique_keys(style_panel, c("cell_id", "month", "style"), "style_panel")
  saveRDS(style_panel, file.path(DIRS$panels, "panel_style_cell_style_month.rds"))

  actor_style_long <- rbindlist(list(
    ac[selective_flag == TRUE & !is.na(actor_type_primary), .(cell_id, month, actor_type = actor_type_primary, style = "selective", vac_actor_style = 1L)],
    ac[indiscriminate_flag == TRUE & !is.na(actor_type_primary), .(cell_id, month, actor_type = actor_type_primary, style = "indiscriminate", vac_actor_style = 1L)]
  ), use.names = TRUE, fill = TRUE)
  actor_style_agg <- actor_style_long[, .(vac_actor_style = .N), by = .(cell_id, month, actor_type, style)]

  actor_style_support <- rbindlist(list(
    copy(cell_map)[, `:=`(actor_type = "state", style = "selective")],
    copy(cell_map)[, `:=`(actor_type = "state", style = "indiscriminate")],
    copy(cell_map)[, `:=`(actor_type = "rebel", style = "selective")],
    copy(cell_map)[, `:=`(actor_type = "rebel", style = "indiscriminate")]
  ), use.names = TRUE, fill = TRUE)
  actor_style_panel <- merge(actor_style_support, actor_style_agg, by = c("cell_id", "month", "actor_type", "style"), all.x = TRUE, sort = FALSE)
  actor_style_panel[is.na(vac_actor_style), vac_actor_style := 0L]
  assert_unique_keys(actor_style_panel, c("cell_id", "month", "actor_type", "style"), "actor_style_panel")
  saveRDS(actor_style_panel, file.path(DIRS$panels, "panel_actor_style_cell_actor_style_month.rds"))

  actor_audit <- data.table(
    component = c("actor_panel", "style_panel", "actor_style_panel"),
    built = 1L,
    rows = c(nrow(actor_panel), nrow(style_panel), nrow(actor_style_panel)),
    nonzero = c(sum(actor_panel$vac_actor > 0L), sum(style_panel$vac_style > 0L), sum(actor_style_panel$vac_actor_style > 0L)),
    reason = NA_character_
  )
}

fwrite(actor_audit, file.path(DIRS$audit, "measurement_actor_extension_audit.csv"))

acled_neg <- ac[, .(
  events_protests = sum(grepl("^protests$", event_type_clean)),
  events_riots = sum(grepl("^riots$", event_type_clean)),
  events_demo = sum(grepl("^protests$|^riots$", event_type_clean))
), by = .(cell_id, month)]

core_panel <- merge(core_panel, acled_neg, by = c("cell_id", "month"), all.x = TRUE, sort = FALSE)
for (v in c("events_protests", "events_riots", "events_demo")) {
  core_panel[is.na(get(v)), (v) := 0L]
}

# ---- Low-visibility geography ----
if ("eli_feas_terrain" %in% names(core_panel)) {
  q25 <- quantile(core_panel$eli_feas_terrain, LOW_VIS_MAIN_Q, na.rm = TRUE)
  q10 <- quantile(core_panel$eli_feas_terrain, LOW_VIS_STRICT_Q, na.rm = TRUE)
  core_panel[, low_visibility_q25 := as.integer(eli_feas_terrain <= q25)]
  core_panel[, low_visibility_q10 := as.integer(eli_feas_terrain <= q10)]
} else {
  core_panel[, `:=`(low_visibility_q25 = 0L, low_visibility_q10 = 0L)]
}

# ---- Old global H3 rings (support only) ----
core_panel[, high_eli_cell := as.integer(eli_high_theatre == 1L)]
core_panel[, low_eli_cell := as.integer(eli_high_theatre == 0L)]

hi_cells <- unique(core_panel[high_eli_cell == 1L & !is.na(lat) & !is.na(lon), .(cell_id_hi = cell_id, lat_hi = lat, lon_hi = lon)])
lo_cells <- unique(core_panel[low_eli_cell == 1L & !is.na(lat) & !is.na(lon), .(cell_id, lat, lon)])

if (nrow(hi_cells) > 0L && nrow(lo_cells) > 0L) {
  lo_cells[, dist_to_high_eli_km := vapply(
    seq_len(.N),
    function(i) min(haversine_km(lat[i], lon[i], hi_cells$lat_hi, hi_cells$lon_hi)),
    numeric(1)
  )]
} else {
  lo_cells[, dist_to_high_eli_km := NA_real_]
}

core_panel <- merge(core_panel, lo_cells[, .(cell_id, dist_to_high_eli_km)], by = "cell_id", all.x = TRUE, sort = FALSE)
core_panel[, low_eli_refuge_primary := as.integer(low_eli_cell == 1L & !is.na(dist_to_high_eli_km) & dist_to_high_eli_km <= H3_LOCAL_MAX_KM)]
core_panel[, refuge_ring_0_50 := as.integer(low_eli_cell == 1L & !is.na(dist_to_high_eli_km) & dist_to_high_eli_km <= 50)]
core_panel[, refuge_ring_50_100 := as.integer(low_eli_cell == 1L & !is.na(dist_to_high_eli_km) & dist_to_high_eli_km > 50 & dist_to_high_eli_km <= 100)]
core_panel[, refuge_ring_100_150 := as.integer(low_eli_cell == 1L & !is.na(dist_to_high_eli_km) & dist_to_high_eli_km > 100 & dist_to_high_eli_km <= 150)]
core_panel[, refuge_ring_gt150 := as.integer(low_eli_cell == 1L & !is.na(dist_to_high_eli_km) & dist_to_high_eli_km > 150)]

# ---- New headline H3: local source clusters ----
# Build contiguous clusters among high-ELI cells using queen adjacency on the PRIO lattice.
hi_nodes <- unique(core_panel[high_eli_cell == 1L, .(cell_id)])
hi_nodes[, `:=`(
  pg_row = ((as.integer(cell_id) - 1L) %/% 720L) + 1L,
  pg_col = ((as.integer(cell_id) - 1L) %% 720L) + 1L
)]
if (nrow(hi_nodes) > 0L) {
  setkey(hi_nodes, pg_row, pg_col)
  offsets <- data.table(dr = rep(-1:1, each = 3), dc = rep(-1:1, 3))
  offsets <- offsets[!(dr == 0 & dc == 0)]

  edges <- rbindlist(lapply(seq_len(nrow(offsets)), function(k) {
    tmp <- copy(hi_nodes)
    tmp[, `:=`(nbr_row = pg_row + offsets$dr[k], nbr_col = pg_col + offsets$dc[k])]
    tmp <- merge(
      tmp[, .(from = cell_id, nbr_row, nbr_col)],
      hi_nodes[, .(to = cell_id, pg_row, pg_col)],
      by.x = c("nbr_row", "nbr_col"), by.y = c("pg_row", "pg_col"),
      all = FALSE, sort = FALSE
    )
    tmp[from != to, .(from, to)]
  }), fill = TRUE)
  if (nrow(edges) == 0L) {
    hi_nodes[, source_cluster_id := paste0("cluster_", seq_len(.N))]
  } else {
    # custom connected components to avoid external igraph dependency
    edges_cc <- unique(rbind(
      data.table(from = as.character(edges$from), to = as.character(edges$to)),
      data.table(from = as.character(edges$to), to = as.character(edges$from))
    ), by = c("from", "to"))
    adj_list <- split(edges_cc$to, edges_cc$from)
    node_ids <- as.character(hi_nodes$cell_id)
    membership <- setNames(rep(NA_integer_, length(node_ids)), node_ids)
    comp_id <- 0L
    for (node in node_ids) {
      if (!is.na(membership[[node]])) next
      comp_id <- comp_id + 1L
      queue <- node
      membership[[node]] <- comp_id
      while (length(queue) > 0L) {
        cur <- queue[[1]]
        if (length(queue) == 1L) {
          queue <- character()
        } else {
          queue <- queue[-1]
        }
        nbrs <- adj_list[[cur]]
        if (length(nbrs) == 0L) next
        new_nbrs <- nbrs[is.na(membership[nbrs])]
        if (length(new_nbrs) > 0L) {
          membership[new_nbrs] <- comp_id
          queue <- c(queue, new_nbrs)
        }
      }
    }
    memb <- data.table(cell_id = as.integer(names(membership)),
                       source_cluster_id = paste0("cluster_", as.integer(membership)))
    hi_nodes <- merge(hi_nodes, memb, by = "cell_id", all.x = TRUE, sort = FALSE)
  }

  hi_cluster_coords <- merge(
    hi_nodes[, .(cell_id, source_cluster_id)],
    unique(core_panel[, .(cell_id, lat, lon)]),
    by = "cell_id", all.x = TRUE, sort = FALSE
  )

  lo_loc <- unique(core_panel[low_eli_cell == 1L & !is.na(lat) & !is.na(lon), .(cell_id, lat, lon)])
  if (nrow(lo_loc) > 0L) {
    cluster_ids <- unique(hi_cluster_coords$source_cluster_id)
    assign_list <- lapply(cluster_ids, function(cl) {
      src <- hi_cluster_coords[source_cluster_id == cl & !is.na(lat) & !is.na(lon)]
      if (nrow(src) == 0L) return(NULL)
      data.table(
        cell_id = lo_loc$cell_id,
        source_cluster_id = cl,
        dist_km = vapply(seq_len(nrow(lo_loc)), function(i) {
          min(haversine_km(lo_loc$lat[i], lo_loc$lon[i], src$lat, src$lon))
        }, numeric(1))
      )
    })
    assign_dt <- rbindlist(assign_list, fill = TRUE)
    setorder(assign_dt, cell_id, dist_km)
    nearest_cluster <- assign_dt[, .SD[1], by = cell_id]
    setnames(nearest_cluster, "dist_km", "dist_to_source_cluster_km")
    core_panel <- merge(core_panel, nearest_cluster, by = "cell_id", all.x = TRUE, sort = FALSE)
  } else {
    core_panel[, `:=`(source_cluster_id = NA_character_, dist_to_source_cluster_km = NA_real_)]
  }
} else {
  core_panel[, `:=`(source_cluster_id = NA_character_, dist_to_source_cluster_km = NA_real_)]
}

core_panel[, adjacent_0_50 := as.integer(low_eli_cell == 1L & !is.na(dist_to_source_cluster_km) & dist_to_source_cluster_km <= H3_ADJ_KM)]
core_panel[, fringe_50_100 := as.integer(low_eli_cell == 1L & !is.na(dist_to_source_cluster_km) & dist_to_source_cluster_km > H3_ADJ_KM & dist_to_source_cluster_km <= H3_LOCAL_MAX_KM)]
core_panel[, h3_local_sample := as.integer(low_eli_cell == 1L & !is.na(dist_to_source_cluster_km) & dist_to_source_cluster_km <= H3_LOCAL_MAX_KM)]

assert_unique_keys(core_panel, c("cell_id", "month"), "core_panel")
saveRDS(core_panel, file.path(DIRS$panels, "panel_core_cell_month.rds"))

# audits
fwrite(data.table(
  metric = c("core_rows", "core_cells", "core_months", "nonzero_type3", "nonzero_protests", "nonzero_riots", "h3_local_cells"),
  value = c(
    nrow(core_panel),
    uniqueN(core_panel$cell_id),
    uniqueN(core_panel$month),
    sum(core_panel$events_type3 > 0L),
    sum(core_panel$events_protests > 0L),
    sum(core_panel$events_riots > 0L),
    uniqueN(core_panel$cell_id[core_panel$h3_local_sample == 1L])
  )
), file.path(DIRS$audit, "measurement_causal_summary.csv"))

write_audit_text(c(
  "Measurement build summary (causal upgrade)",
  paste0("GED source: ", ged_path),
  paste0("ACLED source: ", acled_path),
  "Core outcomes built: events_type3, events_type1, any_type3, share_type3_theatre_month",
  "Negative-control outcomes built: events_protests, events_riots, events_demo",
  "Actor adjudication panels built when ACLED actor/interaction columns are available: actor, style, actor_style",
  "Low-visibility geography built: low_visibility_q25, low_visibility_q10",
  "Headline H3 built from local source clusters: adjacent_0_50, fringe_50_100, h3_local_sample",
  "Global refuge rings retained as appendix support checks."
), "measurement_build_summary_causal_upgrade")

log_msg("02_build_measurement_geo_afghanistan_only_EXECUTABLE_ACLED: done")
