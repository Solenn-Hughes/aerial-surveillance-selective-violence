# ============================================================
# R01_ingest_priogrid.R — Ingest PRIO-GRID static + yearly (ROBUST)
# Outputs:
#   data/processed/base/priogrid_static.rds
#   data/processed/base/priogrid_yearly.rds
# ============================================================

source("R00_config.R")
log_msg("Step 1: ingest PRIO-GRID (static + yearly)")

# -------------------------
# 1) Locate input files (auto-pick latest by mtime)
# -------------------------
prio_dir <- file.path(PROJECT_ROOT, "data", "raw", "priogrid")
if (!dir.exists(prio_dir)) stop("PRIO-GRID raw dir not found: ", prio_dir)

static_candidates <- list.files(prio_dir, pattern = "PRIO[- ]?GRID.*Static.*\\.csv$", ignore.case = TRUE, full.names = TRUE)
yearly_candidates <- list.files(prio_dir, pattern = "PRIO[- ]?GRID.*Yearly.*\\.csv$", ignore.case = TRUE, full.names = TRUE)

if (!length(static_candidates)) stop("No PRIO-GRID Static CSV found in: ", prio_dir)
if (!length(yearly_candidates)) stop("No PRIO-GRID Yearly CSV found in: ", prio_dir)

pick_latest <- function(paths) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  paths[order(file.info(paths)$mtime, decreasing = TRUE)][1]
}

STATIC_CSV <- pick_latest(static_candidates)
YEARLY_CSV <- pick_latest(yearly_candidates)

log_msg("Using STATIC: ", STATIC_CSV)
log_msg("Using YEARLY: ", YEARLY_CSV)

write_audit_dt(
  data.table(
    kind = c("static", "yearly"),
    path = c(STATIC_CSV, YEARLY_CSV),
    mtime = c(file.info(STATIC_CSV)$mtime, file.info(YEARLY_CSV)$mtime)
  ),
  "priogrid_inputs"
)

# -------------------------
# 2) Helpers
# -------------------------
rename_if_present <- function(dt, old, new) {
  if (old %in% names(dt) && !(new %in% names(dt))) setnames(dt, old, new)
  invisible(dt)
}

# robust numeric parser:
# - trims
# - removes non-breaking spaces
# - converts decimal comma to dot when appropriate
# - strips harmless trailing stuff
to_num_locale <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\u00A0", "", x, fixed = TRUE)   # NBSP
  x <- gsub(" ", "", x, fixed = TRUE)
  
  # If it has comma but no dot, treat comma as decimal separator
  x <- ifelse(grepl(",", x, fixed = TRUE) & !grepl("\\.", x), gsub(",", ".", x, fixed = TRUE), x)
  
  # remove anything that's not part of a standard number (keeps - + e E .)
  x <- gsub("[^0-9eE\\+\\-\\.]", "", x)
  
  suppressWarnings(as.numeric(x))
}

to_int_locale <- function(x) suppressWarnings(as.integer(round(to_num_locale(x))))

# Rescale suffix vars that look like 0..100 (%)
rescale_suffix_pct_if_needed <- function(dt, suffix_regex) {
  cols <- grep(suffix_regex, names(dt), value = TRUE, ignore.case = TRUE)
  
  empty_audit <- data.table(
    var = character(),
    min_before = numeric(),
    max_before = numeric(),
    scaled = logical(),
    min_after = numeric(),
    max_after = numeric()
  )
  if (!length(cols)) return(list(dt = dt, audit = empty_audit))
  
  aud_list <- lapply(cols, function(cc) {
    x <- dt[[cc]]
    if (!is.numeric(x)) return(NULL)
    x_ok <- x[!is.na(x)]
    if (!length(x_ok)) return(NULL)
    
    mn <- suppressWarnings(min(x_ok))
    mx <- suppressWarnings(max(x_ok))
    
    do_scale <- is.finite(mn) && is.finite(mx) && mn >= 0 && mx > 1.5 && mx <= 100.5
    if (do_scale) dt[, (cc) := get(cc) / 100]
    
    data.table(
      var = cc,
      min_before = mn,
      max_before = mx,
      scaled = do_scale,
      min_after = if (do_scale) mn/100 else mn,
      max_after = if (do_scale) mx/100 else mx
    )
  })
  
  aud <- rbindlist(aud_list, fill = TRUE)
  if (is.null(aud) || nrow(aud) == 0) aud <- empty_audit
  list(dt = dt, audit = aud[order(-scaled, var)])
}

# -------------------------
# 3) Read STATIC
# -------------------------
dt_static <- fread(STATIC_CSV, showProgress = FALSE, fill = TRUE, quote = "")

# Normalize key names
dt_static <- rename_if_present(dt_static, "GID", "gid")
dt_static <- rename_if_present(dt_static, "ROW", "row")
dt_static <- rename_if_present(dt_static, "COL", "col")
dt_static <- rename_if_present(dt_static, "XCOORD", "xcoord")
dt_static <- rename_if_present(dt_static, "YCOORD", "ycoord")

assert_cols(dt_static, c("gid","row","col","xcoord","ycoord"), "dt_static")

# Coerce keys robustly
dt_static[, gid := to_int_locale(gid)]
dt_static[, row := to_int_locale(row)]
dt_static[, col := to_int_locale(col)]
dt_static[, xcoord := to_num_locale(xcoord)]
dt_static[, ycoord := to_num_locale(ycoord)]

bad_static <- dt_static[
  is.na(gid) | is.na(row) | is.na(col) | is.na(xcoord) | is.na(ycoord),
  .(gid,row,col,xcoord,ycoord)
]

if (nrow(bad_static)) {
  # save examples
  write_audit_dt(bad_static[1:min(200, .N)], "priogrid_static_bad_keys")
  
  # If only a handful of rows are bad, drop them (often blank/footer rows)
  share_bad <- nrow(bad_static) / nrow(dt_static)
  if (share_bad <= 0.001 || nrow(bad_static) <= 10) {
    log_msg("STATIC: dropping ", nrow(bad_static), " bad-key rows (share=", round(share_bad, 6), ").")
    dt_static <- dt_static[!(is.na(gid) | is.na(row) | is.na(col) | is.na(xcoord) | is.na(ycoord))]
  } else {
    stop("STATIC: too many bad-key rows (share=", round(share_bad, 4),
         "). Open the audit priogrid_static_bad_keys_*.csv; this indicates a systematic parse issue.")
  }
}

assert_unique(dt_static, c("gid"), "dt_static")

# Rescale *_gc vars if they look like percentages
res_static <- rescale_suffix_pct_if_needed(dt_static, "_gc$")
dt_static <- res_static$dt
audit_static_scale <- res_static$audit

write_audit_dt(
  data.table(n_rows = nrow(dt_static), n_cols = ncol(dt_static), n_scaled_gc = sum(audit_static_scale$scaled, na.rm = TRUE)),
  "priogrid_static_dims"
)
if (nrow(audit_static_scale)) write_audit_dt(audit_static_scale, "priogrid_static_gc_rescale")

miss_static <- audit_missing(dt_static, "priogrid_static")[1:min(30, .N)]
write_audit_dt(miss_static, "priogrid_static_missing_top30")

out_static <- file.path(DIRS$base, "priogrid_static.rds")
saveRDS(dt_static, out_static)
log_msg("Wrote: ", out_static)

# -------------------------
# 4) Read YEARLY
# -------------------------
dt_yearly <- fread(YEARLY_CSV, showProgress = FALSE, fill = TRUE, quote = "")

dt_yearly <- rename_if_present(dt_yearly, "GID", "gid")
dt_yearly <- rename_if_present(dt_yearly, "YEAR", "year")

assert_cols(dt_yearly, c("gid","year"), "dt_yearly")

dt_yearly[, gid := to_int_locale(gid)]
dt_yearly[, year := to_int_locale(year)]

bad_yearly <- dt_yearly[is.na(gid) | is.na(year), .(gid, year)]
if (nrow(bad_yearly)) {
  write_audit_dt(bad_yearly[1:min(200, .N)], "priogrid_yearly_bad_keys")
  share_bad <- nrow(bad_yearly) / nrow(dt_yearly)
  if (share_bad <= 0.001 || nrow(bad_yearly) <= 10) {
    log_msg("YEARLY: dropping ", nrow(bad_yearly), " bad-key rows (share=", round(share_bad, 6), ").")
    dt_yearly <- dt_yearly[!(is.na(gid) | is.na(year))]
  } else {
    stop("YEARLY: too many bad-key rows (share=", round(share_bad, 4),
         "). See priogrid_yearly_bad_keys_*.csv")
  }
}

assert_unique(dt_yearly, c("gid","year"), "dt_yearly")

# Rescale *_ih vars if they look like percentages
res_yearly_ih <- rescale_suffix_pct_if_needed(dt_yearly, "_ih$")
dt_yearly <- res_yearly_ih$dt
audit_yearly_ih <- res_yearly_ih$audit

write_audit_dt(
  data.table(
    n_rows = nrow(dt_yearly),
    n_cols = ncol(dt_yearly),
    years_min = min(dt_yearly$year, na.rm = TRUE),
    years_max = max(dt_yearly$year, na.rm = TRUE),
    n_scaled_ih = sum(audit_yearly_ih$scaled, na.rm = TRUE)
  ),
  "priogrid_yearly_dims"
)
if (nrow(audit_yearly_ih)) write_audit_dt(audit_yearly_ih, "priogrid_yearly_ih_rescale")

miss_yearly <- audit_missing(dt_yearly, "priogrid_yearly")[1:min(30, .N)]
write_audit_dt(miss_yearly, "priogrid_yearly_missing_top30")

out_yearly <- file.path(DIRS$base, "priogrid_yearly.rds")
saveRDS(dt_yearly, out_yearly)
log_msg("Wrote: ", out_yearly)

log_msg("STATIC: rows=", nrow(dt_static), " cols=", ncol(dt_static), " scaled_gc=", sum(audit_static_scale$scaled, na.rm = TRUE))
log_msg("YEARLY: rows=", nrow(dt_yearly), " cols=", ncol(dt_yearly),
        " years=", min(dt_yearly$year, na.rm = TRUE), "–", max(dt_yearly$year, na.rm = TRUE),
        " scaled_ih=", sum(audit_yearly_ih$scaled, na.rm = TRUE))

log_msg("Step 1 complete.")
