# ============================================================
# R05b_build_moderators.R — Build H2 and H3 proxy variables
#
# H2 (control_weak): proxy for weak/contested territorial control
#   = z-scored travel time to nearest city (ttime_mean)
#   Rationale: high travel time = remote = weak state reach = contested
#   Source: already in controls_pg_month.rds from PRIO-GRID
#
# H3 (rebel_capability): proxy for rebel organisational capacity
#   = 12-month lagged rolling count of type 1 (state-based battle)
#     events in cell, log1p-transformed and z-scored
#   Rationale: sustained state-based fighting in a cell indicates
#   rebels are actively contesting it and have sufficient organisation
#   to mount repeated engagements. Lagged 12 months to avoid
#   simultaneity with the type 3 outcome.
#   Source: ged_clean.rds
#
# Output:
#   Appends control_weak + rebel_capability to controls_pg_month.rds
#   (overwrites in place — original backed up first)
#
# Run AFTER R05, BEFORE R07.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 5b: build H2 (control_weak) and H3 (rebel_capability) moderators")

# -------------------------
# 1) PATHS
# -------------------------
P_CTRL <- file.path(DIRS$base, "controls_pg_month.rds")
P_GED  <- file.path(DIRS$base, "ged_clean.rds")

if (!file.exists(P_CTRL)) stop("Missing controls file: ", P_CTRL, " — run R05 first.")
if (!file.exists(P_GED))  stop("Missing GED file: ", P_GED,  " — run R03 first.")

# -------------------------
# 2) BACKUP ORIGINAL CONTROLS
# -------------------------
bak <- file.path(DIRS$base, paste0("controls_pg_month_BACKUP_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
file.copy(P_CTRL, bak)
log_msg("Backed up original controls to: ", bak)

# -------------------------
# 3) LOAD CONTROLS
# -------------------------
ctrl <- readRDS(P_CTRL); setDT(ctrl)
assert_cols(ctrl, c("pg_id", "month"), "controls_pg_month")
ctrl[, pg_id := as.integer(pg_id)]
ctrl[, month := as_month_idate(month)]
assert_unique(ctrl, c("pg_id", "month"), "controls_pg_month")

n_ctrl <- nrow(ctrl)
log_msg("Controls loaded: ", n_ctrl, " rows, ", uniqueN(ctrl$pg_id), " cells, ",
        uniqueN(ctrl$month), " months")

# -------------------------
# 4) H2: control_weak = z-scored ttime_mean
# -------------------------
# ttime_mean = mean travel time to nearest city (minutes)
# Higher = more remote = weaker state reach = more likely contested
# Z-scored so coefficient is interpretable in SD units

if (!("ttime_mean" %in% names(ctrl))) {
  stop("ttime_mean not found in controls — cannot build control_weak. ",
       "Check R05 output.")
}

n_na_ttime <- sum(is.na(ctrl$ttime_mean))
if (n_na_ttime > 0) {
  log_msg("ttime_mean has ", n_na_ttime, " NAs (", round(100*n_na_ttime/nrow(ctrl), 2),
          "%) — imputing with median before z-scoring.")
  med_ttime <- median(ctrl$ttime_mean, na.rm = TRUE)
  ctrl[is.na(ttime_mean), ttime_mean := med_ttime]
}

mu_t  <- mean(ctrl$ttime_mean, na.rm = TRUE)
sd_t  <- sd(ctrl$ttime_mean, na.rm = TRUE)
ctrl[, control_weak := (ttime_mean - mu_t) / sd_t]

log_msg("H2 (control_weak): mean=", round(mean(ctrl$control_weak), 4),
        " sd=", round(sd(ctrl$control_weak), 4),
        " min=", round(min(ctrl$control_weak), 2),
        " max=", round(max(ctrl$control_weak), 2))

# -------------------------
# 5) H3: rebel_capability = lagged rolling type-1 event count
# -------------------------
ged <- readRDS(P_GED); setDT(ged)
assert_cols(ged, c("priogrid_gid", "type_of_violence", "date_start"), "ged_clean")

# Keep only type 1 (state-based violence = active rebel contestation)
ged1 <- ged[type_of_violence == 1L]
ged1[, pg_id := as.integer(priogrid_gid)]
ged1[, month := as_month_idate(date_start)]

# Count type-1 events per cell-month
type1_counts <- ged1[, .(type1_count = .N), by = .(pg_id, month)]

# Full skeleton of all cell-months in controls (so zeros are explicit)
skeleton <- ctrl[, .(pg_id, month)]
setkey(skeleton, pg_id, month)
setkey(type1_counts, pg_id, month)

skeleton <- type1_counts[skeleton]  # left join
skeleton[is.na(type1_count), type1_count := 0L]
assert_unique(skeleton, c("pg_id", "month"), "skeleton after type1 join")

# Compute 12-month lagged rolling sum within each cell
# i.e. rebel_capability in month t = sum of type1 events in months t-12 to t-1
setorder(skeleton, pg_id, month)

skeleton[, month_int := year(month) * 12L + month(month)]

# For each cell, for each month, sum type1 in prior 12 months (not including current)
# Using a fast data.table rolling join approach
skeleton[, rebel_capability_raw := {
  mi <- month_int
  tc <- type1_count
  # For each row i, sum tc where mi[j] >= mi[i]-12 and mi[j] < mi[i]
  # Use cumsum trick since data is ordered by month within cell
  cs <- cumsum(tc)
  # lag 1 month
  cs_lag1  <- shift(cs, 1L,  fill = 0L)
  cs_lag13 <- shift(cs, 13L, fill = 0L)
  cs_lag1 - cs_lag13
}, by = pg_id]

# log1p transform (right-skewed count) then z-score
mu_r <- mean(skeleton$rebel_capability_raw, na.rm = TRUE)
sd_r <- sd(skeleton$rebel_capability_raw, na.rm = TRUE)
skeleton[, rebel_capability := (log1p(rebel_capability_raw) - mean(log1p(rebel_capability_raw), na.rm=TRUE)) /
           sd(log1p(rebel_capability_raw), na.rm=TRUE)]

log_msg("H3 (rebel_capability): mean=", round(mean(skeleton$rebel_capability, na.rm=TRUE), 4),
        " sd=", round(sd(skeleton$rebel_capability, na.rm=TRUE), 4),
        " share_zero_raw=", round(mean(skeleton$rebel_capability_raw == 0), 4))

# -------------------------
# 6) MERGE MODERATORS INTO CONTROLS
# -------------------------
# Drop rebel_capability if already exists (re-run safe)
# NOTE: do NOT drop control_weak here — it was computed directly on ctrl above
if ("rebel_capability" %in% names(ctrl)) ctrl[, rebel_capability := NULL]

# Merge rebel_capability from skeleton into ctrl
setkey(ctrl, pg_id, month)
setkey(skeleton, pg_id, month)

ctrl <- skeleton[, .(pg_id, month, rebel_capability)][ctrl]

# control_weak should already be on ctrl from step 4 — verify
if (!("control_weak" %in% names(ctrl))) stop("control_weak missing after merge — unexpected.")
assert_unique(ctrl, c("pg_id", "month"), "controls post-moderator merge")

if (nrow(ctrl) != n_ctrl) stop("Row count changed after merge: before=", n_ctrl, " after=", nrow(ctrl))

# Check no NAs in new vars
na_cw <- sum(is.na(ctrl$control_weak))
na_rc <- sum(is.na(ctrl$rebel_capability))
if (na_cw > 0) stop("control_weak has ", na_cw, " NAs after build — unexpected.")
if (na_rc > 0) stop("rebel_capability has ", na_rc, " NAs after build — unexpected.")

# -------------------------
# 7) AUDIT + SAVE
# -------------------------
audit <- data.table(
  var        = c("control_weak", "rebel_capability"),
  n_rows     = nrow(ctrl),
  n_na       = c(na_cw, na_rc),
  mean       = c(mean(ctrl$control_weak), mean(ctrl$rebel_capability)),
  sd         = c(sd(ctrl$control_weak),   sd(ctrl$rebel_capability)),
  min        = c(min(ctrl$control_weak),  min(ctrl$rebel_capability)),
  p25        = c(quantile(ctrl$control_weak, 0.25), quantile(ctrl$rebel_capability, 0.25)),
  median     = c(median(ctrl$control_weak),  median(ctrl$rebel_capability)),
  p75        = c(quantile(ctrl$control_weak, 0.75), quantile(ctrl$rebel_capability, 0.75)),
  max        = c(max(ctrl$control_weak),  max(ctrl$rebel_capability))
)
write_audit_dt(audit, "R05b_moderators_summary")
log_msg("Audit written.")

saveRDS(ctrl, P_CTRL)
log_msg("Saved updated controls (with control_weak + rebel_capability) to: ", P_CTRL)
log_msg("Step 5b complete.")