# ============================================================
# R10_audits.R — Step 10 (Automated audits, fail-fast on integrity)
#
# Reads final analysis panels (post-ELI merged):
#   data/processed/panels/panel_pg_month_ALLCW_paper1_posteli.rds
#   data/processed/panels/panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds
#
# Writes:
#   outputs/paper1/audit/R10_audits_<timestamp>.csv
#   outputs/paper1/audit/R10_audits_<timestamp>.txt
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

# -------------------------
# 0) CONFIG
# -------------------------
cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found (R00_config.R or R/00_config.R).")
source(cfg_path)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_AUDIT_CSV <- file.path(DIRS$audit, paste0("R10_audits_", TS, ".csv"))
OUT_AUDIT_TXT <- file.path(DIRS$audit, paste0("R10_audits_", TS, ".txt"))

# -------------------------
# 1) INPUTS (final panels)
# -------------------------
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds")

inputs <- data.table(
  panel = c("ALLCW","STRICT"),
  path  = c(P_ALLCW, P_STRICT),
  exists = c(file.exists(P_ALLCW), file.exists(P_STRICT))
)

if (any(!inputs$exists)) {
  stop("Missing required Step 10 input(s):\n",
       paste(inputs[exists == FALSE, paste(panel, path, sep=": ")], collapse = "\n"))
}

# -------------------------
# 2) HELPERS
# -------------------------
std_keys <- function(dt, label) {
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  assert_cols(dt, c("pg_id","month"), name = label)
  dt[, pg_id := as.integer(pg_id)]
  dt[, month := as_month_idate(month)]
  assert_unique(dt, c("pg_id","month"), name = paste0(label, " keys"))
  dt
}

miss_tbl <- function(dt, label) {
  m <- audit_missing(dt, name = label)
  # keep top 40 NA shares (informational)
  m[order(-share_na)][1:min(.N, 40)]
}

range_check <- function(dt, col, lo, hi, label) {
  assert_cols(dt, col, name = label)
  v <- dt[[col]]
  if (anyNA(v)) stop(sprintf("%s: %s has NA (integrity).", label, col))
  if (any(v < lo | v > hi)) {
    badn <- sum(v < lo | v > hi)
    stop(sprintf("%s: %s out of bounds [%s,%s] (n=%d).", label, col, lo, hi, badn))
  }
  TRUE
}

binary_check <- function(dt, col, label) {
  assert_cols(dt, col, name = label)
  v <- dt[[col]]
  if (anyNA(v)) stop(sprintf("%s: %s has NA (integrity).", label, col))
  if (any(!v %in% c(0L,1L))) stop(sprintf("%s: %s not binary 0/1.", label, col))
  TRUE
}

nonneg_check <- function(dt, col, label) {
  assert_cols(dt, col, name = label)
  v <- dt[[col]]
  if (!is.numeric(v)) stop(sprintf("%s: %s not numeric.", label, col))
  if (any(v < 0, na.rm = TRUE)) stop(sprintf("%s: %s has negatives.", label, col))
  TRUE
}

qsum <- function(v) {
  if (!is.numeric(v)) return(data.table())
  qs <- as.numeric(quantile(v, probs = c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE))
  data.table(p = c("0","1","5","25","50","75","95","99","100"), value = qs)
}

# -------------------------
# 3) LOAD + AUDIT EACH PANEL
# -------------------------
audit_parts <- list()
txt_lines <- c(
  paste0("R10 audits run: ", Sys.time()),
  paste0("PROJECT_ROOT: ", PROJECT_ROOT),
  ""
)

for (i in 1:nrow(inputs)) {
  panel_name <- inputs$panel[i]
  path <- inputs$path[i]
  log_msg("Loading: ", panel_name, " -> ", path)
  
  dt <- readRDS(path); setDT(dt)
  dt <- std_keys(dt, panel_name)
  
  # Required vars
  required <- c("eli_feas","post_eli","event_time","treated")
  assert_cols(dt, required, name = panel_name)
  
  # Bounds/integrity
  range_check(dt, "eli_feas", 0, 1, panel_name)
  binary_check(dt, "post_eli", panel_name)
  binary_check(dt, "treated", panel_name)
  
  # event_time should exist for treated rows
  if (dt[treated == 1L, sum(is.na(event_time))] > 0) {
    stop(panel_name, ": treated==1 rows have NA event_time (integrity).")
  }
  
  # Outcomes: audit if present
  outcomes <- intersect(c(
    "events_type3","events_ged",
    "fat_best_type3","fat_best_ged",
    "civ_best_type3","civ_best_ged"
  ), names(dt))
  
  for (y in outcomes) nonneg_check(dt, y, panel_name)
  
  # Controls presence: informational only (don’t fail if extra columns)
  ctrl_like <- setdiff(names(dt), c("pg_id","month"))
  # Missingness top table
  audit_parts[[paste0(panel_name, "_missingness_top")]] <-
    cbind(section = "missingness_top", panel = panel_name, miss_tbl(dt, panel_name))
  
  # Key summaries
  summ <- data.table(
    section = "panel_summary",
    panel = panel_name,
    rows = nrow(dt),
    pg = uniqueN(dt$pg_id),
    months = uniqueN(dt$month),
    countries = if ("country_id" %in% names(dt)) uniqueN(dt$country_id) else NA_integer_,
    share_post = mean(dt$post_eli == 1L),
    eli_min = min(dt$eli_feas),
    eli_max = max(dt$eli_feas),
    event_time_min = min(dt$event_time, na.rm = TRUE),
    event_time_max = max(dt$event_time, na.rm = TRUE)
  )
  audit_parts[[paste0(panel_name, "_summary")]] <- summ
  
  # Quantile summaries for ELI + outcomes
  q_eli <- qsum(dt$eli_feas)
  if (nrow(q_eli)) audit_parts[[paste0(panel_name, "_q_eli")]] <-
    cbind(section = "quantiles", panel = panel_name, var = "eli_feas", q_eli)
  
  for (y in outcomes) {
    qy <- qsum(dt[[y]])
    if (nrow(qy)) audit_parts[[paste0(panel_name, "_q_", y)]] <-
        cbind(section = "quantiles", panel = panel_name, var = y, qy)
  }
  
  # Informational: are there any events pre vs post?
  if ("events_type3" %in% names(dt)) {
    pre_sum  <- dt[post_eli == 0L, sum(events_type3, na.rm = TRUE)]
    post_sum <- dt[post_eli == 1L, sum(events_type3, na.rm = TRUE)]
    audit_parts[[paste0(panel_name, "_events_prepost")]] <- data.table(
      section="prepost_sums",
      panel=panel_name,
      var="events_type3",
      sum_pre=pre_sum,
      sum_post=post_sum
    )
  }
  if ("events_ged" %in% names(dt)) {
    pre_sum  <- dt[post_eli == 0L, sum(events_ged, na.rm = TRUE)]
    post_sum <- dt[post_eli == 1L, sum(events_ged, na.rm = TRUE)]
    audit_parts[[paste0(panel_name, "_ged_prepost")]] <- data.table(
      section="prepost_sums",
      panel=panel_name,
      var="events_ged",
      sum_pre=pre_sum,
      sum_post=post_sum
    )
  }
  
  txt_lines <- c(txt_lines,
                 paste0("Panel: ", panel_name),
                 paste0("  rows=", nrow(dt), " pg=", uniqueN(dt$pg_id), " months=", uniqueN(dt$month)),
                 paste0("  eli_feas in [", min(dt$eli_feas), ", ", max(dt$eli_feas), "]"),
                 paste0("  share_post=", round(mean(dt$post_eli==1L), 4)),
                 paste0("  event_time range [", min(dt$event_time, na.rm=TRUE), ", ", max(dt$event_time, na.rm=TRUE), "]"),
                 ""
  )
}

# -------------------------
# 4) WRITE OUTPUTS
# -------------------------
audit_dt <- rbindlist(audit_parts, fill = TRUE)
fwrite(audit_dt, OUT_AUDIT_CSV)
writeLines(txt_lines, con = OUT_AUDIT_TXT, useBytes = TRUE)

log_msg("Wrote audit CSV: ", OUT_AUDIT_CSV)
log_msg("Wrote audit TXT: ", OUT_AUDIT_TXT)
log_msg("Step 10 complete.")
