# ============================================================
# R04_build_base_panel_pg_month.R — Build full PRIO grid×month panel
# Inputs:
#   data/processed/base/priogrid_static.rds
#   data/processed/base/ged_pg_month.rds
#   data/processed/eli/eli_feasibility_pg.rds
#
# Output:
#   data/processed/base/panel_pg_month_base.rds
#
# Audits:
#   outputs/paper1/audit/panel_base_dims_*.csv
#   outputs/paper1/audit/panel_base_months_*.csv
#   outputs/paper1/audit/panel_base_eli_missing_*.csv
# ============================================================

source("R00_config.R")
log_msg("Step 4: build full PRIO grid×month base panel (zeros + ELI feasibility)")

pg   <- readRDS(file.path(DIRS$base, "priogrid_static.rds"))
pm   <- readRDS(file.path(DIRS$base, "ged_pg_month.rds"))
eli  <- readRDS(file.path(DIRS$eli,  "eli_feasibility_pg.rds"))

assert_unique(pg, "gid", "priogrid_static")
assert_unique(eli, "gid", "eli_feasibility_pg")

# months universe from GED aggregates
months <- sort(unique(pm$month))
month_dt <- data.table(month = months)
month_dt[, month_id := as.integer(factor(month, levels = months))]

# grid universe from PRIO
grid_dt <- pg[, .(pg_id = gid)]

# full skeleton
panel <- CJ(pg_id = grid_dt$pg_id, month = months, unique = TRUE)
panel[, month_id := as.integer(factor(month, levels = months))]

# merge GED aggregates (sparse -> full) and fill zeros
setkey(panel, pg_id, month)
setkey(pm, pg_id, month)
panel <- pm[panel]

# fill NAs for event counts and fatalities with 0
zero_cols <- c("events_ged","fat_best_ged","civ_best_ged","events_type3","fat_best_type3","civ_best_type3")
for (cc in zero_cols) panel[is.na(get(cc)), (cc) := 0L]

# attach ELI feasibility (gid == pg_id)
eli2 <- eli[, .(pg_id = gid,
                eli_feas_terrain, min_km_airbase, eli_feas_range, eli_feas)]
setkey(eli2, pg_id)
setkey(panel, pg_id)
panel <- eli2[panel]

# sanity: ELI should exist for all PRIO gids
miss_eli <- panel[is.na(eli_feas), .N]
write_audit_dt(data.table(missing_eli = miss_eli, total_rows = nrow(panel)), "panel_base_eli_missing")
if (miss_eli > 0) stop("Some panel rows missing ELI feasibility. See audit panel_base_eli_missing_*.csv")

# audits
dims <- data.table(
  n_rows = nrow(panel),
  n_pg   = uniqueN(panel$pg_id),
  n_months = uniqueN(panel$month),
  share_zero_events = mean(panel$events_ged == 0),
  share_zero_type3  = mean(panel$events_type3 == 0)
)
write_audit_dt(dims, "panel_base_dims")

month_summ <- panel[, .(
  tot_events = sum(events_ged),
  tot_type3  = sum(events_type3),
  active_cells = sum(events_ged > 0)
), by = .(month)][order(month)]
write_audit_dt(month_summ, "panel_base_months")

out <- file.path(DIRS$base, "panel_pg_month_base.rds")
saveRDS(panel, out)
log_msg("Wrote: ", out)
log_msg("Step 4 complete.")

# --- Audit GED cells not covered by PRIO static (land/covariate subset)
pm_out <- pm[!.(pg_id = pg$gid), on=.(pg_id)]
if (nrow(pm_out) > 0) {
  write_audit_dt(
    pm_out[, .(
      n_rows = .N,
      dropped_events = sum(events_ged),
      dropped_type3  = sum(events_type3)
    )],
    "ged_pg_month_outside_priogrid_static_summary"
  )
  write_audit_dt(pm_out[order(-events_ged)][1:min(200, .N)],
                 "ged_pg_month_outside_priogrid_static_examples")
}

# Restrict pm to PRIO-static-covered cells before building months universe
pm <- pm[.(pg_id = pg$gid), on=.(pg_id)]

