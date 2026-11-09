# ============================================================
# R08_define_post_eli.R — Step 8 (GLOBAL ELI onset, one-go)
#
# LOCKED RULE:
#   theatre = country_id (from ALLCW paper1 panel)
#   treated = 1 for all theatres in sample
#   eli_onset_month = GLOBAL_ONSET (single date for all theatres)
#   post_eli = month >= eli_onset_month
#   event_time = months_since_onset
#
# Inputs:
#   data/processed/panels/panel_pg_month_ALLCW_paper1.rds
#   data/processed/panels/panel_pg_month_ALLCW_ACD_STRICT_paper1.rds
#
# Outputs:
#   data/processed/eli/eli_onsets.csv
#   data/processed/eli/post_eli_theatre_month.rds
#   data/processed/panels/panel_pg_month_ALLCW_paper1_posteli.rds
#   data/processed/panels/panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds
#
# Audit:
#   outputs/paper1/audit/R08_define_post_eli_*.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

# -------------------------
# 0) CONFIG
# -------------------------
cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found (R00_config.R or R/00_config.R).")
source(cfg_path)

# -------------------------
# 1) SETTINGS (locked baseline)
# -------------------------
GLOBAL_ONSET_STR <- "2001-10-01"
# (or "2001-07-01")
GLOBAL_ONSET <- as.IDate(GLOBAL_ONSET_STR)

# If an older eli_onsets.csv exists (e.g. your conflict-start version), archive it automatically
ARCHIVE_OLD_ONSETS <- TRUE

# -------------------------
# 2) PATHS
# -------------------------
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1.rds")

OUT_ONSETS_CSV    <- file.path(DIRS$eli,   "eli_onsets.csv")
OUT_THEATRE_MONTH <- file.path(DIRS$eli,   "post_eli_theatre_month.rds")
OUT_ALLCW_POST    <- file.path(DIRS$panels,"panel_pg_month_ALLCW_paper1_posteli.rds")
OUT_STRICT_POST   <- file.path(DIRS$panels,"panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds")

for (p in c(P_ALLCW, P_STRICT)) if (!file.exists(p)) stop(sprintf("Missing required input: %s", p))

# -------------------------
# 3) HELPERS
# -------------------------
std_pg_id <- function(dt, name) {
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  assert_cols(dt, "pg_id", name = name)
  dt[, pg_id := as.integer(pg_id)]
  if (anyNA(dt$pg_id)) stop(sprintf("%s: pg_id has NA after coercion.", name))
  dt
}

std_month <- function(dt, name) {
  if (!("month" %in% names(dt))) {
    for (alt in c("ym","date","month_id")) {
      if (alt %in% names(dt)) { setnames(dt, alt, "month"); break }
    }
  }
  assert_cols(dt, "month", name = name)
  dt[, month := as_month_idate(month)]
  if (anyNA(dt$month)) stop(sprintf("%s: month has NA after as_month_idate().", name))
  dt
}

month_diff <- function(a, b) (year(a) * 12L + month(a)) - (year(b) * 12L + month(b))

# -------------------------
# 4) LOAD PANELS
# -------------------------
log_msg("Loading ALLCW paper1: ", P_ALLCW)
p1 <- readRDS(P_ALLCW); setDT(p1)
p1 <- std_pg_id(p1, "ALLCW_paper1")
p1 <- std_month(p1, "ALLCW_paper1")
assert_unique(p1, c("pg_id","month"), name = "ALLCW_paper1")
assert_cols(p1, c("country_id","country"), name = "ALLCW_paper1")

p1[, country_id := as.integer(country_id)]
if (anyNA(p1$country_id)) stop("ALLCW_paper1: country_id has NA after coercion.")

log_msg("Loading STRICT paper1: ", P_STRICT)
p2 <- readRDS(P_STRICT); setDT(p2)
p2 <- std_pg_id(p2, "STRICT_paper1")
p2 <- std_month(p2, "STRICT_paper1")
assert_unique(p2, c("pg_id","month"), name = "STRICT_paper1")
assert_cols(p2, c("country_id","country"), name = "STRICT_paper1")

p2[, country_id := as.integer(country_id)]
if (anyNA(p2$country_id)) stop("STRICT_paper1: country_id has NA after coercion.")

months_univ <- sort(unique(p1$month))
if (length(months_univ) == 0) stop("ALLCW_paper1: no months found.")
panel_month_min <- min(months_univ)
panel_month_max <- max(months_univ)

# Integrity: onset must be within panel time window
if (GLOBAL_ONSET < panel_month_min || GLOBAL_ONSET > panel_month_max) {
  stop(sprintf("GLOBAL_ONSET=%s lies outside panel month range [%s, %s].",
               GLOBAL_ONSET_STR, as.character(panel_month_min), as.character(panel_month_max)))
}

theatres <- unique(p1[, .(country_id, country)])
assert_unique(theatres, "country_id", name = "theatres(country_id)")

# -------------------------
# 5) ARCHIVE OLD eli_onsets.csv (if present)
# -------------------------
if (ARCHIVE_OLD_ONSETS && file.exists(OUT_ONSETS_CSV)) {
  bak <- file.path(DIRS$eli, paste0("eli_onsets_ARCHIVE_", .now_stamp(), ".csv"))
  ok <- file.copy(OUT_ONSETS_CSV, bak, overwrite = FALSE)
  if (ok) log_msg("Archived old eli_onsets.csv to: ", bak)
}

# -------------------------
# 6) WRITE NEW eli_onsets.csv (GLOBAL)
# -------------------------
onsets <- theatres[, .(
  country_id,
  country,
  treated = 1L,
  eli_onset_month = as.character(GLOBAL_ONSET),
  onset_source = "GLOBAL",
  notes = paste0("LOCKED: global ELI onset = ", GLOBAL_ONSET_STR)
)]
fwrite(onsets, OUT_ONSETS_CSV)
log_msg("Wrote: ", OUT_ONSETS_CSV)

# -------------------------
# 7) BUILD THEATRE×MONTH (post_eli + event_time)
# -------------------------
theatre_month <- CJ(country_id = theatres$country_id, month = months_univ, unique = TRUE)
theatre_month <- merge(theatre_month, onsets[, .(country_id, country, treated, eli_onset_month = as.IDate(eli_onset_month), onset_source)],
                       by = "country_id", all.x = TRUE, sort = FALSE)
assert_unique(theatre_month, c("country_id","month"), name = "theatre_month")

theatre_month[, post_eli := as.integer(month >= eli_onset_month)]
theatre_month[, event_time := month_diff(month, eli_onset_month)]

# hard integrity
if (anyNA(theatre_month$post_eli)) stop("post_eli has NA — integrity failure.")
if (any(!theatre_month$post_eli %in% c(0L,1L))) stop("post_eli not binary — integrity failure.")
if (anyNA(theatre_month$event_time)) stop("event_time has NA — integrity failure.")

saveRDS(theatre_month, OUT_THEATRE_MONTH)
log_msg("Wrote: ", OUT_THEATRE_MONTH)

# -------------------------
# 8) MERGE INTO PANELS + SAVE
# -------------------------
merge_posteli <- function(panel_dt, label) {
  assert_unique(panel_dt, c("pg_id","month"), name = paste0(label, " pre-merge"))
  n_before <- nrow(panel_dt)
  
  merged <- merge(panel_dt,
                  theatre_month[, .(country_id, month, treated, eli_onset_month, onset_source, post_eli, event_time)],
                  by = c("country_id","month"),
                  all.x = TRUE, sort = FALSE)
  
  if (nrow(merged) != n_before) stop(sprintf("%s: row count changed after merge.", label))
  assert_unique(merged, c("pg_id","month"), name = paste0(label, " post-merge"))
  assert_no_na(merged, c("pg_id","month"), name = paste0(label, " keys post-merge"))
  
  # should never be NA
  if (anyNA(merged$post_eli)) stop(sprintf("%s: post_eli NA after merge.", label))
  if (anyNA(merged$event_time)) stop(sprintf("%s: event_time NA after merge.", label))
  merged
}

log_msg("Merging into ALLCW...")
p1_post <- merge_posteli(p1, "ALLCW_paper1")
saveRDS(p1_post, OUT_ALLCW_POST)
log_msg("Wrote: ", OUT_ALLCW_POST)

log_msg("Merging into STRICT...")
p2_post <- merge_posteli(p2, "STRICT_paper1")
saveRDS(p2_post, OUT_STRICT_POST)
log_msg("Wrote: ", OUT_STRICT_POST)

# -------------------------
# 9) AUDIT (coverage + overlap)
# -------------------------
# How many theatres span the onset within ALLCW panel? (min<onset<=max)
span <- p1[, .(min_m = min(month), max_m = max(month)), by = .(country_id, country)]
span[, spans_onset := as.integer(min_m < GLOBAL_ONSET & max_m >= GLOBAL_ONSET)]
span_summary <- span[, .(
  n_theatres = .N,
  n_spanning = sum(spans_onset),
  share_spanning = mean(spans_onset)
)]

# Pre/post row shares (ALLCW)
pre_post_rows <- p1_post[, .(
  rows_pre = sum(post_eli == 0L),
  rows_post = sum(post_eli == 1L),
  share_post = mean(post_eli == 1L)
)]

meta <- data.table(
  section = "meta",
  global_onset = GLOBAL_ONSET_STR,
  panel_month_min = as.character(panel_month_min),
  panel_month_max = as.character(panel_month_max),
  out_onsets_csv = OUT_ONSETS_CSV,
  out_theatre_month = OUT_THEATRE_MONTH,
  out_allcw_post = OUT_ALLCW_POST,
  out_strict_post = OUT_STRICT_POST
)

audit <- rbindlist(list(
  meta,
  cbind(section = "theatre_span_summary", span_summary),
  cbind(section = "pre_post_row_summary_ALLCW", pre_post_rows),
  cbind(section = "theatre_span_details_head", span[order(-spans_onset, country)][1:min(.N, 120)])
), fill = TRUE)

write_audit_dt(audit, "R08_define_post_eli", subdir = DIRS$audit)
log_msg("Step 8 (GLOBAL) complete. Audit written to: ", DIRS$audit)
