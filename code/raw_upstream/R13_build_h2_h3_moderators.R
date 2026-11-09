# ============================================================
# R13_build_h2_h3_moderators.R
# Build TIME-VARYING H2 (territorial contestation) and
# H3 (rebel capability) moderators from ACLED data.
#
# These REPLACE the static proxies from R05b:
#   H2 old: control_weak = z-scored ttime_mean (static)
#   H2 new: contestation_tv = 12-month rolling battle intensity
#           in the cell (ACLED battles, log1p, z-scored)
#
#   H3 old: rebel_capability = lagged GED type-1 count (static-ish)
#   H3 new: rebel_cap_tv = 12-month rolling count of rebel-initiated
#           events (ACLED), log1p, z-scored
#
# Input:  data/processed/base/acled_pg_month.rds (from R12)
#         data/processed/panels/panel_pg_month_ALLCW_paper1_posteli_v2.rds
# Output: data/processed/panels/panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds
#
# Run AFTER R12, BEFORE R14.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 13: Build time-varying H2 & H3 moderators from ACLED")

# ── LOAD ACLED AGGREGATES ────────────────────────────────────
P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")
if (!file.exists(P_ACLED)) stop("Missing: ", P_ACLED, " — run R12 first.")
agg <- readRDS(P_ACLED); setDT(agg)
agg[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) agg[, month := as_month_idate(month)] else agg[, month := as.IDate(month)]

log_msg("ACLED aggregates: ", nrow(agg), " rows")

# ── LOAD PANEL ───────────────────────────────────────────────
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (!file.exists(P_ALLCW)) stop("Missing: ", P_ALLCW)
dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

log_msg("Panel: ", nrow(dt), " rows, ", uniqueN(dt$pg_id), " cells")


# Verify expected ACLED columns exist before rolling
expected_cols <- c("acled_battles", "acled_rebel_initiated")
missing_cols <- setdiff(expected_cols, names(agg))
if (length(missing_cols)) stop("Missing ACLED columns: ", paste(missing_cols, collapse = ", "),
                               "\nCheck R12 aggregation output.")
# ── BUILD ROLLING WINDOW INDICATORS ─────────────────────────
# We use a 12-month backward-looking rolling sum, LAGGED by 1
# month to avoid simultaneity with the DV.
#
# For each cell-month, sum ACLED events from (t-13) to (t-2),
# i.e., the 12 months ending 1 month before the current month.

setkey(agg, pg_id, month)

# Create month_int for ACLED to enable integer arithmetic
agg[, month_int := year(month) * 12L + month(month)]

# H2: Contestation = rolling battles (12-month, lagged 1)
# H3: Rebel capability = rolling rebel-initiated events (12-month, lagged 1)

# Efficient approach: for each cell, use frollsum
setorder(agg, pg_id, month_int)

agg[, battles_roll12 := frollsum(acled_battles, n = 12L, align = "right",
                                  fill = NA_real_), by = pg_id]
agg[, rebel_roll12   := frollsum(acled_rebel_initiated, n = 12L, align = "right",
                                  fill = NA_real_), by = pg_id]

setorder(agg, pg_id, month_int)
# Lag by 1 month (shift forward: current month gets LAST month's rolling sum)
agg[, battles_roll12_lag := shift(battles_roll12, n = 1L, type = "lag"), by = pg_id]
agg[, rebel_roll12_lag   := shift(rebel_roll12, n = 1L, type = "lag"), by = pg_id]

# Fill NAs with 0 (cells where we don't have 12 months of history yet)
agg[is.na(battles_roll12_lag), battles_roll12_lag := 0]
agg[is.na(rebel_roll12_lag),   rebel_roll12_lag := 0]

# Log transform
agg[, log_battles_roll := log1p(battles_roll12_lag)]
agg[, log_rebel_roll   := log1p(rebel_roll12_lag)]

# Z-score (across full sample)
z_safe <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

agg[, contestation_tv := z_safe(log_battles_roll)]
agg[, rebel_cap_tv    := z_safe(log_rebel_roll)]

# ── ALSO BUILD SPATIAL NEIGHBOURHOOD CONTESTATION ────────────
# For H2, we also want a broader spatial measure: are battles
# happening in ADJACENT cells? This captures the "contested
# territory" concept better than just the focal cell.
# We compute a 3×3 neighbourhood sum (queen contiguity).

# Get coordinates for spatial neighbours
if (all(c("pg_row", "pg_col") %in% names(agg))) {
  # Already have grid coords
} else {
  # Reconstruct from pg_id
  agg[, pg_row := ((pg_id - 1L) %/% 720L) + 1L]
  agg[, pg_col := ((pg_id - 1L) %% 720L) + 1L]
}

# Queen neighbours: for each cell, sum battles in 3×3 window
# This is expensive, so we do it in a vectorised way
log_msg("Computing spatial neighbourhood contestation (3x3 queen)...")

nbr_offsets <- CJ(dr = -1:1, dc = -1:1)
nbr_offsets <- nbr_offsets[!(dr == 0 & dc == 0)]  # exclude self

# Build neighbour lookup
cells_info <- unique(agg[, .(pg_id, pg_row, pg_col)])
setkey(cells_info, pg_row, pg_col)

# For each cell, find its neighbours
nbr_list <- list()
for (i in seq_len(nrow(nbr_offsets))) {
  dr <- nbr_offsets$dr[i]
  dc <- nbr_offsets$dc[i]
  tmp <- copy(cells_info)
  tmp[, `:=`(nbr_row = pg_row + dr, nbr_col = pg_col + dc)]
  tmp[, nbr_id := (nbr_row - 1L) * 720L + nbr_col]
  nbr_list[[i]] <- tmp[, .(pg_id, nbr_id)]
}
nbr_dt <- rbindlist(nbr_list)

# Only keep neighbours that actually exist in our data
nbr_dt <- nbr_dt[nbr_id %in% cells_info$pg_id]

# Merge neighbour battles
nbr_battles <- merge(
  nbr_dt,
  agg[, .(nbr_id = pg_id, month, battles_roll12_lag)],
  by = "nbr_id",
  allow.cartesian = TRUE
)

# Sum across neighbours for each cell-month
nbr_sum <- nbr_battles[, .(
  nbr_battles_sum = sum(battles_roll12_lag, na.rm = TRUE)
), by = .(pg_id, month)]

# Merge back
agg <- merge(agg, nbr_sum, by = c("pg_id", "month"), all.x = TRUE)
agg[is.na(nbr_battles_sum), nbr_battles_sum := 0]
agg[, contest_spatial_tv := z_safe(log1p(nbr_battles_sum))]

log_msg("Spatial contestation computed.")

# ── MERGE INTO PANEL ─────────────────────────────────────────
merge_cols <- c("pg_id", "month", "contestation_tv", "rebel_cap_tv",
                "contest_spatial_tv", "acled_battles", "acled_rebel_events")
merge_cols <- intersect(merge_cols, names(agg))

n_before <- nrow(dt)
dt <- merge(dt, agg[, ..merge_cols], by = c("pg_id", "month"), all.x = TRUE)
stopifnot(nrow(dt) == n_before)

# Fill NAs (cells with no ACLED match = zero conflict)
for (cc in c("contestation_tv", "rebel_cap_tv", "contest_spatial_tv",
             "acled_battles", "acled_rebel_events")) {
  if (cc %in% names(dt)) {
    n_na <- sum(is.na(dt[[cc]]))
    if (n_na > 0) {
      log_msg("  Filling ", n_na, " NAs in ", cc, " with 0 / global mean")
      # For z-scored vars, fill with 0 (= mean); for counts, fill with 0
      fill_val <- if (grepl("_tv$", cc)) 0 else 0L
      dt[is.na(get(cc)), (cc) := fill_val]
    }
  }
}

# ── AUDIT ────────────────────────────────────────────────────
h2h3_audit <- data.table(
  variable = c("contestation_tv", "rebel_cap_tv", "contest_spatial_tv"),
  mean = c(mean(dt$contestation_tv), mean(dt$rebel_cap_tv),
           mean(dt$contest_spatial_tv, na.rm = TRUE)),
  sd   = c(sd(dt$contestation_tv), sd(dt$rebel_cap_tv),
           sd(dt$contest_spatial_tv, na.rm = TRUE)),
  min  = c(min(dt$contestation_tv), min(dt$rebel_cap_tv),
           min(dt$contest_spatial_tv, na.rm = TRUE)),
  max  = c(max(dt$contestation_tv), max(dt$rebel_cap_tv),
           max(dt$contest_spatial_tv, na.rm = TRUE)),
  share_zero = c(mean(dt$contestation_tv == 0), mean(dt$rebel_cap_tv == 0),
                 mean(dt$contest_spatial_tv == 0, na.rm = TRUE))
)
print(h2h3_audit)
write_audit_dt(h2h3_audit, "R13_h2_h3_moderators_summary")

# Correlation with old proxies if they exist
if ("control_weak" %in% names(dt)) {
  r_h2 <- cor(dt$contestation_tv, dt$control_weak, use = "complete.obs")
  log_msg("Cor(contestation_tv, control_weak_old): ", round(r_h2, 3))
}
if ("rebel_capability" %in% names(dt)) {
  r_h3 <- cor(dt$rebel_cap_tv, dt$rebel_capability, use = "complete.obs")
  log_msg("Cor(rebel_cap_tv, rebel_capability_old): ", round(r_h3, 3))
}

# ── SAVE ─────────────────────────────────────────────────────
out_path <- file.path(DIRS$panels,
                      "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds")
saveRDS(dt, out_path)
log_msg("Wrote: ", out_path)

log_msg("Step 13 complete.")
