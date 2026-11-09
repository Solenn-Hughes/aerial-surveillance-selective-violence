# ============================================================
# R07_build_panel_pg_month.R — Step 7 (Paper 1 master panels)
# Merges controls into the two outcome panels that already include ELI + outcomes.
#
# Inputs:
#   data/processed/panels/panel_pg_month_ALLCW.rds
#   data/processed/panels/panel_pg_month_ALLCW_ACD_STRICT.rds
#   data/processed/base/controls_pg_month.rds
#
# Outputs:
#   data/processed/panels/panel_pg_month_ALLCW_paper1.rds
#   data/processed/panels/panel_pg_month_ALLCW_ACD_STRICT_paper1.rds
#
# Audits (CSV):
#   outputs/paper1/audit/R07_build_panel_pg_month_*.csv
#
# Fail-fast only on TRUE integrity failures:
#   - missing input files
#   - duplicate keys
#   - row-count changes after merge
#   - missing keys after merge
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

# -------------------------
# 0) CONFIG
# -------------------------
if (!file.exists("R00_config.R")) stop("R00_config.R not found in project root.")
source("R00_config.R")

# -------------------------
# 1) PATHS
# -------------------------
P_CTRL   <- file.path(DIRS$base,   "controls_pg_month.rds")
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT.rds")

OUT_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1.rds")
OUT_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1.rds")

for (p in c(P_CTRL, P_ALLCW, P_STRICT)) {
  if (!file.exists(p)) stop(sprintf("Missing required input: %s", p))
}

# -------------------------
# 2) HELPERS
# -------------------------
# Standardize pg_id key
std_pg_id <- function(dt, name) {
  if (!("pg_id" %in% names(dt))) {
    if ("gid" %in% names(dt)) dt[, pg_id := gid]
  }
  assert_cols(dt, "pg_id", name = name)
  dt[, pg_id := as.integer(pg_id)]
  if (anyNA(dt$pg_id)) stop(sprintf("%s: pg_id has NA after coercion.", name))
  dt
}

# Standardize month to IDate (YYYY-MM-01)
std_month <- function(dt, name) {
  if (!("month" %in% names(dt))) {
    # Try common alternatives (only if present)
    for (alt in c("ym", "date", "month_id")) {
      if (alt %in% names(dt)) {
        setnames(dt, alt, "month")
        break
      }
    }
  }
  assert_cols(dt, "month", name = name)
  
  # If already Date/IDate/POSIXct, as_month_idate will handle; if character/numeric we parse safely.
  if (inherits(dt$month, "IDate") || inherits(dt$month, "Date") || inherits(dt$month, "POSIXt")) {
    dt[, month := as_month_idate(month)]
  } else if (is.character(dt$month)) {
    x <- dt$month
    # Accept "YYYY-MM" or "YYYY-MM-DD"
    needs_day <- !is.na(x) & nchar(x) == 7 & grepl("^\\d{4}-\\d{2}$", x)
    x[needs_day] <- paste0(x[needs_day], "-01")
    d <- as.IDate(x)
    if (anyNA(d)) stop(sprintf("%s: month parse produced NA (character).", name))
    dt[, month := as_month_idate(d)]
  } else if (is.numeric(dt$month) || is.integer(dt$month)) {
    xi <- as.integer(round(dt$month))
    med <- suppressWarnings(median(xi, na.rm = TRUE))
    if (!is.finite(med)) stop(sprintf("%s: month column all NA.", name))
    # If YYYYMMDD
    if (med >= 19000101L) {
      y <- xi %/% 10000L
      m <- (xi %/% 100L) %% 100L
      dt[, month := to_month_idate(y, m)]
    } else if (med >= 190001L) { # YYYYMM
      y <- xi %/% 100L
      m <- xi %% 100L
      dt[, month := to_month_idate(y, m)]
    } else {
      stop(sprintf("%s: numeric month format not recognized.", name))
    }
  } else {
    stop(sprintf("%s: month column type not recognized.", name))
  }
  
  if (anyNA(dt$month)) stop(sprintf("%s: month has NA after standardization.", name))
  dt
}

# Merge controls into a panel (all.x) while preventing column collisions
merge_controls <- function(panel_dt, controls_dt, panel_name) {
  # Ensure keys exist and are unique
  assert_unique(panel_dt,   c("pg_id", "month"), name = paste0(panel_name, " (pre-merge)"))
  assert_unique(controls_dt, c("pg_id", "month"), name = "controls (pre-merge)")
  
  n_before <- nrow(panel_dt)
  
  # Avoid collisions: keep panel’s existing columns if same name exists
  key_cols <- c("pg_id", "month")
  overlap <- intersect(setdiff(names(controls_dt), key_cols), setdiff(names(panel_dt), key_cols))
  controls_use <- controls_dt
  if (length(overlap) > 0) controls_use <- controls_use[, setdiff(names(controls_use), overlap), with = FALSE]
  
  merged <- merge(panel_dt, controls_use, by = key_cols, all.x = TRUE, sort = FALSE)
  
  # Integrity checks
  if (nrow(merged) != n_before) {
    stop(sprintf("%s: row count changed after merge (before=%d, after=%d).", panel_name, n_before, nrow(merged)))
  }
  assert_unique(merged, c("pg_id", "month"), name = paste0(panel_name, " (post-merge)"))
  assert_no_na(merged, c("pg_id", "month"), name = paste0(panel_name, " (post-merge keys)"))
  
  list(
    dt = merged,
    overlap_dropped = overlap,
    new_cols = setdiff(names(merged), names(panel_dt))
  )
}

# -------------------------
# 3) LOAD INPUTS
# -------------------------
log_msg("Loading controls: ", P_CTRL)
ctrl <- readRDS(P_CTRL); setDT(ctrl)
ctrl <- std_pg_id(ctrl, "controls")
ctrl <- std_month(ctrl, "controls")
assert_unique(ctrl, c("pg_id", "month"), name = "controls")

log_msg("Loading ALLCW panel: ", P_ALLCW)
allcw <- readRDS(P_ALLCW); setDT(allcw)
allcw <- std_pg_id(allcw, "ALLCW")
allcw <- std_month(allcw, "ALLCW")

log_msg("Loading STRICT panel: ", P_STRICT)
strict <- readRDS(P_STRICT); setDT(strict)
strict <- std_pg_id(strict, "ALLCW_ACD_STRICT")
strict <- std_month(strict, "ALLCW_ACD_STRICT")

# -------------------------
# 4) MERGE + SAVE
# -------------------------
res_allcw  <- merge_controls(allcw,  ctrl, "ALLCW")
res_strict <- merge_controls(strict, ctrl, "ALLCW_ACD_STRICT")

# Save outputs
saveRDS(res_allcw$dt,  OUT_ALLCW)
saveRDS(res_strict$dt, OUT_STRICT)

log_msg("Saved: ", OUT_ALLCW)
log_msg("Saved: ", OUT_STRICT)

# -------------------------
# 5) AUDITS (single CSV)
# -------------------------
mk_meta <- function(panel_name, pre_dt, post_dt, overlap) {
  data.table(
    panel = panel_name,
    pre_rows = nrow(pre_dt),
    post_rows = nrow(post_dt),
    pre_pg = uniqueN(pre_dt$pg_id),
    post_pg = uniqueN(post_dt$pg_id),
    pre_months = uniqueN(pre_dt$month),
    post_months = uniqueN(post_dt$month),
    month_min = as.character(min(post_dt$month)),
    month_max = as.character(max(post_dt$month)),
    overlap_dropped_n = length(overlap)
  )
}

# Missingness of newly added controls only
mk_new_missing <- function(panel_name, post_dt, new_cols) {
  if (length(new_cols) == 0) return(data.table(panel = panel_name, var = character(), n_na = integer(), share_na = numeric()))
  data.table(
    panel = panel_name,
    var = new_cols,
    n_na = vapply(new_cols, function(v) sum(is.na(post_dt[[v]])), integer(1)),
    share_na = vapply(new_cols, function(v) mean(is.na(post_dt[[v]])), numeric(1))
  )[order(-share_na)]
}

audit_meta <- rbindlist(list(
  mk_meta("ALLCW", allcw, res_allcw$dt, res_allcw$overlap_dropped),
  mk_meta("ALLCW_ACD_STRICT", strict, res_strict$dt, res_strict$overlap_dropped)
), fill = TRUE)

audit_overlap <- rbindlist(list(
  data.table(panel = "ALLCW", dropped_col = res_allcw$overlap_dropped),
  data.table(panel = "ALLCW_ACD_STRICT", dropped_col = res_strict$overlap_dropped)
), fill = TRUE)

audit_new_missing <- rbindlist(list(
  mk_new_missing("ALLCW", res_allcw$dt, res_allcw$new_cols),
  mk_new_missing("ALLCW_ACD_STRICT", res_strict$dt, res_strict$new_cols)
), fill = TRUE)

audit_all <- rbindlist(list(
  cbind(section = "panel_meta", audit_meta),
  cbind(section = "overlap_dropped_from_controls", audit_overlap),
  cbind(section = "missingness_new_controls_after_merge", audit_new_missing),
  cbind(section = "missingness_full_post_merge_ALLCW", audit_missing(res_allcw$dt, "ALLCW_post_merge")[1:min(.N, 60)]),
  cbind(section = "missingness_full_post_merge_STRICT", audit_missing(res_strict$dt, "STRICT_post_merge")[1:min(.N, 60)])
), fill = TRUE)

write_audit_dt(audit_all, "R07_build_panel_pg_month", subdir = DIRS$audit)

log_msg("Step 7 complete.")
