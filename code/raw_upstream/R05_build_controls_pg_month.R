# ============================================================
# Step 5 — Build controls panel at (pg_id, month)
# Output: data/processed/base/controls_pg_month.rds
# Audit: outputs/paper1/audit/R05_build_controls_pg_month_<timestamp>.csv
# Fail-fast ONLY on integrity failures (duplicates, key loss, impossible ranges)
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

# -------------------------
# 0) CONFIG
# -------------------------
if (!file.exists("R00_config.R")) stop("R00_config.R not found in project root.")
source("R00_config.R")

# Lightweight fallbacks (won't override your config helpers if present)
if (!exists("assert_cols")) {
  assert_cols <- function(dt, cols, label = "") {
    miss <- setdiff(cols, names(dt))
    if (length(miss) > 0) stop(sprintf("Missing columns in %s: %s", label, paste(miss, collapse = ", ")))
  }
}
if (!exists("assert_unique")) {
  assert_unique <- function(dt, cols, label = "") {
    if (dt[duplicated(dt, by = cols), .N] > 0) stop(sprintf("Duplicate keys in %s on (%s).", label, paste(cols, collapse = ",")))
  }
}
if (!exists("log_msg")) {
  log_msg <- function(...) message(sprintf(...))
}

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Always create audit dir (no assumptions)
PROJECT_ROOT <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
AUDIT_DIR <- file.path(PROJECT_ROOT, "outputs", "paper1", "audit")
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

# Paths
P_STATIC  <- file.path(PROJECT_ROOT, "data", "processed", "base", "priogrid_static.rds")
P_YEARLY  <- file.path(PROJECT_ROOT, "data", "processed", "base", "priogrid_yearly.rds")

P_SKEL_ALLCW <- file.path(PROJECT_ROOT, "data", "processed", "panels", "panel_pg_month_ALLCW.rds")
P_SKEL_BASE  <- file.path(PROJECT_ROOT, "data", "processed", "base", "panel_pg_month_base.rds")

OUT_PATH <- file.path(PROJECT_ROOT, "data", "processed", "base", "controls_pg_month.rds")

# Fail fast on missing required inputs
for (p in c(P_STATIC, P_YEARLY)) {
  if (!file.exists(p)) stop(sprintf("Required input missing: %s", p))
}

# -------------------------
# 1) AUDIT WRITER (robust to your helper signatures)
# -------------------------
emit_audit <- function(name_no_ts, x_dt) {
  # Preferred: your config helper write_audit(name, dt_or_list) which adds timestamp & dir
  if (exists("write_audit")) {
    ok <- tryCatch({
      get("write_audit")(name_no_ts, x_dt)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
  }
  # Secondary: your config helper write_audit_dt(name, dt) (also adds dir+timestamp typically)
  if (exists("write_audit_dt")) {
    ok <- tryCatch({
      get("write_audit_dt")(name_no_ts, x_dt)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
  }
  
  # Fallback: manual file write with exactly one timestamp and correct directory
  out <- file.path(AUDIT_DIR, sprintf("%s_%s.csv", name_no_ts, ts))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  fwrite(x_dt, out)
  invisible(TRUE)
}

# -------------------------
# 2) KEY + MONTH STANDARDISATION
# -------------------------
std_pg_key <- function(dt, label = "") {
  if (!("pg_id" %in% names(dt))) {
    if ("gid" %in% names(dt)) dt[, pg_id := gid]
  }
  assert_cols(dt, "pg_id", label)
  dt[, pg_id := as.integer(pg_id)]
  if (anyNA(dt$pg_id)) stop(sprintf("NA pg_id after coercion in %s", label))
  dt
}

detect_month_col <- function(dt) {
  cand <- c("month", "ym", "date", "month_id")
  hit <- cand[cand %in% names(dt)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

to_month_idate <- function(x) {
  if (inherits(x, "IDate")) return(x)
  if (inherits(x, "Date"))  return(as.IDate(x))
  
  if (is.character(x)) {
    xx <- x
    needs_day <- !is.na(xx) & nchar(xx) == 7 & grepl("^\\d{4}-\\d{2}$", xx)
    xx[needs_day] <- paste0(xx[needs_day], "-01")
    return(as.IDate(xx))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    xi <- as.integer(round(x))
    med <- suppressWarnings(median(xi, na.rm = TRUE))
    if (is.finite(med) && med >= 19000101L) {
      y <- xi %/% 10000L
      m <- (xi %/% 100L) %% 100L
      return(as.IDate(sprintf("%04d-%02d-01", y, m)))
    }
    if (is.finite(med) && med >= 190001L) {
      y <- xi %/% 100L
      m <- xi %% 100L
      return(as.IDate(sprintf("%04d-%02d-01", y, m)))
    }
  }
  
  stop("Unrecognised month format (expected Date/IDate, 'YYYY-MM[-DD]', or YYYYMM[DD]).")
}

# -------------------------
# 3) QUOTED NUMERIC + SHARE RESCALE
# -------------------------
strip_quotes <- function(v) {
  if (!is.character(v)) return(v)
  out <- gsub('"', "", v, fixed = TRUE)
  out[out %in% c("", "NA", "NaN")] <- NA_character_
  out
}

maybe_numericize <- function(v) {
  if (!is.character(v)) return(v)
  vv <- strip_quotes(v)
  num <- suppressWarnings(as.numeric(vv))
  denom <- sum(!is.na(vv))
  if (denom == 0) return(v)
  share_num <- sum(!is.na(num)) / denom
  if (share_num >= 0.80) return(num)
  v
}

pct_name <- function(nm) grepl("(_gc$|pct|perc|percent|share|prop)", nm, ignore.case = TRUE)

rescale_pct_if_needed <- function(dt, cols) {
  for (cc in cols) {
    if (!(cc %in% names(dt))) next
    if (!is.numeric(dt[[cc]])) next
    if (!pct_name(cc)) next
    v <- dt[[cc]]
    mx <- suppressWarnings(max(v, na.rm = TRUE))
    mn <- suppressWarnings(min(v, na.rm = TRUE))
    if (is.finite(mx) && is.finite(mn) && mx > 1.01 && mx <= 100.01 && mn >= -1e-9) {
      dt[, (cc) := get(cc) / 100]
    }
  }
  invisible(dt)
}

range_check_pct <- function(dt, cols, label = "") {
  bad <- character()
  for (cc in cols) {
    if (!(cc %in% names(dt))) next
    if (!is.numeric(dt[[cc]])) next
    if (!pct_name(cc)) next
    v <- dt[[cc]]
    if (any(v < -1e-6, na.rm = TRUE) || any(v > 1 + 1e-6, na.rm = TRUE)) bad <- c(bad, cc)
  }
  if (length(bad) > 0) stop(sprintf("Pct/share columns out of [0,1] in %s: %s", label, paste(unique(bad), collapse = ", ")))
}

median_impute_numeric <- function(dt, cols) {
  imp <- data.table(var = character(), n_imputed = integer(), median_value = numeric())
  for (cc in cols) {
    if (!(cc %in% names(dt))) next
    if (!is.numeric(dt[[cc]])) next
    n_na <- sum(is.na(dt[[cc]]))
    if (n_na == 0) next
    med <- suppressWarnings(median(dt[[cc]], na.rm = TRUE))
    if (!is.finite(med)) {
      imp <- rbind(imp, data.table(var = cc, n_imputed = 0L, median_value = NA_real_), fill = TRUE)
      next
    }
    dt[is.na(get(cc)), (cc) := med]
    imp <- rbind(imp, data.table(var = cc, n_imputed = n_na, median_value = med), fill = TRUE)
  }
  list(dt = dt, impute = imp)
}

# -------------------------
# 4) SKELETON (prefer ALLCW)
# -------------------------
skel_path <- if (file.exists(P_SKEL_ALLCW)) P_SKEL_ALLCW else P_SKEL_BASE
if (!file.exists(skel_path)) stop("No skeleton available (ALLCW missing and base panel missing).")

log_msg("Using skeleton: %s", skel_path)
skel <- readRDS(skel_path); setDT(skel)
skel <- std_pg_key(skel, "skeleton")

mcol <- detect_month_col(skel)
if (is.na(mcol)) stop("Skeleton missing a recognised month column (month/ym/date/month_id).")
if (mcol != "month") setnames(skel, mcol, "month")
skel[, month := to_month_idate(month)]
if (anyNA(skel$month)) stop("NA month after standardisation in skeleton.")

skel <- skel[, .(pg_id, month)]
assert_unique(skel, c("pg_id", "month"), "skeleton (pg_id,month)")
skel[, year := as.integer(format(month, "%Y"))]
if (anyNA(skel$year)) stop("NA year derived from month in skeleton.")

pg_year <- unique(skel[, .(pg_id, year)])
assert_unique(pg_year, c("pg_id", "year"), "pg_year skeleton")

# -------------------------
# 5) YEARLY CONTROLS
# -------------------------
prio_y <- readRDS(P_YEARLY); setDT(prio_y)
prio_y <- std_pg_key(prio_y, "priogrid_yearly")
assert_cols(prio_y, "year", "priogrid_yearly")
prio_y[, year := as.integer(year)]
if (anyNA(prio_y$year)) stop("NA year in priogrid_yearly.")

if (prio_y[duplicated(prio_y, by = c("pg_id", "year")), .N] > 0) {
  stop("priogrid_yearly has duplicate (pg_id,year). Integrity failure.")
}

for (nm in names(prio_y)) if (is.character(prio_y[[nm]])) prio_y[, (nm) := maybe_numericize(get(nm))]
rescale_pct_if_needed(prio_y, setdiff(names(prio_y), c("pg_id", "year")))
range_check_pct(prio_y, names(prio_y), "priogrid_yearly")

y_ctrl <- merge(pg_year, prio_y, by = c("pg_id", "year"), all.x = TRUE, sort = FALSE)
year_covars <- setdiff(names(y_ctrl), c("pg_id", "year"))

miss_before_year <- data.table(
  var = year_covars,
  missing_before = vapply(year_covars, function(v) sum(is.na(y_ctrl[[v]])), integer(1))
)

setorder(y_ctrl, pg_id, year)
for (cc in year_covars) {
  if (!is.numeric(y_ctrl[[cc]])) next
  y_ctrl[, (cc) := nafill(get(cc), type = "locf"), by = pg_id]
  y_ctrl[, (cc) := nafill(get(cc), type = "nocb"), by = pg_id]
}

miss_after_locf <- data.table(
  var = year_covars,
  missing_after_locf_nocb = vapply(year_covars, function(v) sum(is.na(y_ctrl[[v]])), integer(1))
)

med_imp_year <- median_impute_numeric(y_ctrl, year_covars)
y_ctrl <- med_imp_year$dt
imp_year <- med_imp_year$impute

miss_after_year <- data.table(
  var = year_covars,
  missing_after_median = vapply(year_covars, function(v) sum(is.na(y_ctrl[[v]])), integer(1))
)

range_check_pct(y_ctrl, names(y_ctrl), "yearly controls after fill+median")

# -------------------------
# 6) STATIC CONTROLS
# -------------------------
prio_s <- readRDS(P_STATIC); setDT(prio_s)
prio_s <- std_pg_key(prio_s, "priogrid_static")

if (prio_s[duplicated(prio_s, by = "pg_id"), .N] > 0) {
  stop("priogrid_static has duplicate pg_id. Integrity failure.")
}

for (nm in names(prio_s)) if (is.character(prio_s[[nm]])) prio_s[, (nm) := maybe_numericize(get(nm))]

static_covars <- setdiff(names(prio_s), "pg_id")
rescale_pct_if_needed(prio_s, static_covars)
range_check_pct(prio_s, names(prio_s), "priogrid_static")

miss_before_static <- data.table(
  var = static_covars,
  missing_before = vapply(static_covars, function(v) sum(is.na(prio_s[[v]])), integer(1))
)

med_imp_static <- median_impute_numeric(prio_s, static_covars)
prio_s <- med_imp_static$dt
imp_static <- med_imp_static$impute

miss_after_static <- data.table(
  var = static_covars,
  missing_after_median = vapply(static_covars, function(v) sum(is.na(prio_s[[v]])), integer(1))
)

range_check_pct(prio_s, names(prio_s), "static controls after median")

# -------------------------
# 7) BUILD MONTHLY CONTROLS: (pg_id, month) + yearly + static
# -------------------------
ctrl <- merge(skel, y_ctrl, by = c("pg_id", "year"), all.x = TRUE, sort = FALSE)
ctrl <- merge(ctrl, prio_s, by = "pg_id", all.x = TRUE, sort = FALSE)

# Integrity checks
if (nrow(ctrl) != nrow(skel)) stop("Row count changed after merging controls. Integrity failure.")
assert_unique(ctrl, c("pg_id", "month"), "controls (pg_id,month)")
if (anyNA(ctrl$pg_id) || anyNA(ctrl$month)) stop("Missing key values in final controls. Integrity failure.")

# No NaN/Inf in numeric columns (integrity)
num_cols <- names(ctrl)[vapply(ctrl, is.numeric, logical(1))]
bad_inf <- vapply(num_cols, function(v) any(!is.finite(ctrl[[v]]), na.rm = TRUE), logical(1))
if (any(bad_inf)) stop(sprintf("Non-finite values found in numeric controls: %s", paste(names(bad_inf)[bad_inf], collapse = ", ")))

setorder(ctrl, pg_id, month)

# -------------------------
# 8) AUDIT + SAVE
# -------------------------
meta <- data.table(
  section = "meta",
  timestamp = ts,
  skeleton_used = skel_path,
  n_rows = nrow(ctrl),
  n_pg = uniqueN(ctrl$pg_id),
  n_months = uniqueN(ctrl$month),
  month_min = as.character(min(ctrl$month)),
  month_max = as.character(max(ctrl$month)),
  out_path = OUT_PATH
)

ctrl_cols <- setdiff(names(ctrl), c("pg_id", "month"))
final_missing <- data.table(
  var = ctrl_cols,
  missing_final = vapply(ctrl_cols, function(v) sum(is.na(ctrl[[v]])), integer(1))
)

audit <- rbindlist(list(
  meta,
  cbind(section = "yearly_missing_before", miss_before_year),
  cbind(section = "yearly_missing_after_locf_nocb", miss_after_locf),
  cbind(section = "yearly_impute_median", imp_year),
  cbind(section = "yearly_missing_after_median", miss_after_year),
  cbind(section = "static_missing_before", miss_before_static),
  cbind(section = "static_impute_median", imp_static),
  cbind(section = "static_missing_after_median", miss_after_static),
  cbind(section = "final_missing_controls", final_missing)
), fill = TRUE)

emit_audit("R05_build_controls_pg_month", audit)

dir.create(dirname(OUT_PATH), recursive = TRUE, showWarnings = FALSE)
saveRDS(ctrl, OUT_PATH)

log_msg("Saved: %s", OUT_PATH)
log_msg("Audit written under: %s", AUDIT_DIR)

# ============================================================
# Checkpoint — Step 5 controls integrity
# Verifies:
#   - unique (pg_id, month)
#   - controls cover ALLCW and STRICT panel keys
#   - no NaN/Inf in numeric columns
#   - pct/share variables in [0,1]
# Writes audit: outputs/paper1/audit/R05_checkpoint_controls_<timestamp>.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

if (!file.exists("R00_config.R")) stop("R00_config.R not found in project root.")
source("R00_config.R")

if (!exists("assert_cols")) {
  assert_cols <- function(dt, cols, label = "") {
    miss <- setdiff(cols, names(dt))
    if (length(miss) > 0) stop(sprintf("Missing columns in %s: %s", label, paste(miss, collapse = ", ")))
  }
}
if (!exists("assert_unique")) {
  assert_unique <- function(dt, cols, label = "") {
    if (dt[duplicated(dt, by = cols), .N] > 0) stop(sprintf("Duplicate keys in %s on (%s).", label, paste(cols, collapse = ",")))
  }
}

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
PROJECT_ROOT <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
AUDIT_DIR <- file.path(PROJECT_ROOT, "outputs", "paper1", "audit")
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

emit_audit <- function(name_no_ts, x_dt) {
  if (exists("write_audit")) {
    ok <- tryCatch({ get("write_audit")(name_no_ts, x_dt); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
  }
  if (exists("write_audit_dt")) {
    ok <- tryCatch({ get("write_audit_dt")(name_no_ts, x_dt); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(TRUE))
  }
  out <- file.path(AUDIT_DIR, sprintf("%s_%s.csv", name_no_ts, ts))
  fwrite(x_dt, out)
  invisible(TRUE)
}

pct_name <- function(nm) grepl("(_gc$|pct|perc|percent|share|prop)", nm, ignore.case = TRUE)

# Paths
P_CTRL   <- file.path(PROJECT_ROOT, "data", "processed", "base", "controls_pg_month.rds")
P_ALLCW  <- file.path(PROJECT_ROOT, "data", "processed", "panels", "panel_pg_month_ALLCW.rds")
P_STRICT <- file.path(PROJECT_ROOT, "data", "processed", "panels", "panel_pg_month_ALLCW_ACD_STRICT.rds")

for (p in c(P_CTRL, P_ALLCW, P_STRICT)) if (!file.exists(p)) stop(sprintf("Missing required file: %s", p))

# Load
ctrl <- readRDS(P_CTRL); setDT(ctrl)
allcw <- readRDS(P_ALLCW); setDT(allcw)
strict <- readRDS(P_STRICT); setDT(strict)

# Standardise keys
if (!("pg_id" %in% names(ctrl)) && "gid" %in% names(ctrl)) ctrl[, pg_id := gid]
if (!("pg_id" %in% names(allcw)) && "gid" %in% names(allcw)) allcw[, pg_id := gid]
if (!("pg_id" %in% names(strict)) && "gid" %in% names(strict)) strict[, pg_id := gid]

assert_cols(ctrl, c("pg_id","month"), "controls")
assert_cols(allcw, c("pg_id","month"), "ALLCW panel")
assert_cols(strict, c("pg_id","month"), "STRICT panel")

# Month to IDate if needed
to_idate <- function(x) {
  if (inherits(x, "IDate")) return(x)
  if (inherits(x, "Date"))  return(as.IDate(x))
  if (is.character(x)) {
    xx <- x
    needs_day <- !is.na(xx) & nchar(xx) == 7 & grepl("^\\d{4}-\\d{2}$", xx)
    xx[needs_day] <- paste0(xx[needs_day], "-01")
    return(as.IDate(xx))
  }
  if (is.numeric(x) || is.integer(x)) {
    xi <- as.integer(round(x))
    med <- suppressWarnings(median(xi, na.rm = TRUE))
    if (is.finite(med) && med >= 190001L) {
      y <- xi %/% 100L
      m <- xi %% 100L
      return(as.IDate(sprintf("%04d-%02d-01", y, m)))
    }
  }
  stop("Unrecognised month format in checkpoint.")
}

ctrl[, month := to_idate(month)]
allcw[, month := to_idate(month)]
strict[, month := to_idate(month)]

# Uniqueness
assert_unique(ctrl, c("pg_id","month"), "controls (pg_id,month)")
assert_unique(allcw, c("pg_id","month"), "ALLCW (pg_id,month)")
assert_unique(strict, c("pg_id","month"), "STRICT (pg_id,month)")

# Coverage: controls must cover ALLCW and STRICT keys
setkey(ctrl, pg_id, month)
miss_allcw  <- allcw[!ctrl, on = .(pg_id, month), .N]
miss_strict <- strict[!ctrl, on = .(pg_id, month), .N]
if (miss_allcw > 0) stop(sprintf("Integrity failure: controls missing %d ALLCW keys.", miss_allcw))
if (miss_strict > 0) stop(sprintf("Integrity failure: controls missing %d STRICT keys.", miss_strict))

# Numeric finiteness
num_cols <- names(ctrl)[vapply(ctrl, is.numeric, logical(1))]
bad_inf <- vapply(num_cols, function(v) any(!is.finite(ctrl[[v]]), na.rm = TRUE), logical(1))
if (any(bad_inf)) stop(sprintf("Integrity failure: non-finite in numeric cols: %s", paste(names(bad_inf)[bad_inf], collapse = ", ")))

# Pct/share bounds
pct_cols <- num_cols[pct_name(num_cols)]
bad_pct <- character()
for (cc in pct_cols) {
  v <- ctrl[[cc]]
  if (any(v < -1e-6, na.rm = TRUE) || any(v > 1 + 1e-6, na.rm = TRUE)) bad_pct <- c(bad_pct, cc)
}
if (length(bad_pct) > 0) stop(sprintf("Integrity failure: pct/share out of [0,1]: %s", paste(unique(bad_pct), collapse = ", ")))

# Missingness summary (not a stop; written to audit)
missing_tbl <- data.table(
  var = setdiff(names(ctrl), c("pg_id","month")),
  missing = vapply(setdiff(names(ctrl), c("pg_id","month")), function(v) sum(is.na(ctrl[[v]])), integer(1)),
  missing_rate = vapply(setdiff(names(ctrl), c("pg_id","month")), function(v) mean(is.na(ctrl[[v]])), numeric(1))
)[order(-missing)]

meta <- data.table(
  section = "meta",
  timestamp = ts,
  ctrl_rows = nrow(ctrl),
  ctrl_pg = uniqueN(ctrl$pg_id),
  ctrl_months = uniqueN(ctrl$month),
  month_min = as.character(min(ctrl$month)),
  month_max = as.character(max(ctrl$month)),
  allcw_rows = nrow(allcw),
  strict_rows = nrow(strict),
  missing_allcw_keys = miss_allcw,
  missing_strict_keys = miss_strict
)

audit <- rbindlist(list(
  meta,
  cbind(section = "missingness_controls", missing_tbl[1:min(.N, 100)])
), fill = TRUE)

emit_audit("R05_checkpoint_controls", audit)

message("✅ Step 5 checkpoint PASSED. Audit written to outputs/paper1/audit/.")

