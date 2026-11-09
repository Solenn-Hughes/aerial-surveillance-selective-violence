# ============================================================
# 00_config.R  — Paper 1 pipeline config + invariants
# Usage:
#   source("R/00_config.R")
# ============================================================

# -------------------------
# 0) Clean workspace hazards
# -------------------------
# data.table uses .N internally; if you've created .N in GlobalEnv, it can break things.
if (exists(".N", envir = .GlobalEnv, inherits = FALSE)) {
  rm(list = ".N", envir = .GlobalEnv)
  message("Removed .N from .GlobalEnv (it masks data.table's .N).")
}

# -------------------------
# 1) Packages + options
# -------------------------
suppressMessages({
  library(data.table)
  # library(fixest)   # uncomment when you start modeling
})

options(stringsAsFactors = FALSE)
options(scipen = 999)   # avoid scientific notation surprises
# options(warn = 1)     # uncomment if you want warnings printed as they occur

# -------------------------
# 2) Project root (EDIT THIS)
# -------------------------
PROJECT_ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
setwd(PROJECT_ROOT)

# -------------------------
# 3) Directory layout (created if missing)
# -------------------------
DIRS <- list(
  data_raw       = file.path(PROJECT_ROOT, "data", "raw"),
  data_processed = file.path(PROJECT_ROOT, "data", "processed"),
  base           = file.path(PROJECT_ROOT, "data", "processed", "base"),
  events         = file.path(PROJECT_ROOT, "data", "processed", "events"),
  outcomes       = file.path(PROJECT_ROOT, "data", "processed", "outcomes"),
  controls       = file.path(PROJECT_ROOT, "data", "processed", "controls"),
  eli            = file.path(PROJECT_ROOT, "data", "processed", "eli"),
  panels         = file.path(PROJECT_ROOT, "data", "processed", "panels"),
  
  out_paper1     = file.path(PROJECT_ROOT, "outputs", "paper1"),
  audit          = file.path(PROJECT_ROOT, "outputs", "paper1", "audit"),
  figures        = file.path(PROJECT_ROOT, "outputs", "paper1", "figures"),
  tables         = file.path(PROJECT_ROOT, "outputs", "paper1", "tables")
)

.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
invisible(lapply(DIRS, .ensure_dir))

# -------------------------
# 4) Logging helpers
# -------------------------
.now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(msg)
  invisible(msg)
}

# -------------------------
# 5) Canonical month handling (use IDate = fast, stable)
# -------------------------
# Always represent month as first day of month: YYYY-MM-01 (IDate)
to_month_idate <- function(year, month) {
  if (length(year) != length(month)) stop("to_month_idate: year and month must have same length.")
  y <- as.integer(year)
  m <- as.integer(month)
  if (anyNA(y) || anyNA(m)) stop("to_month_idate: year/month contain NA after coercion.")
  if (any(m < 1 | m > 12)) stop("to_month_idate: month must be 1..12.")
  data.table::as.IDate(sprintf("%04d-%02d-01", y, m))
}

# Convert arbitrary Date/IDate/datetime to month IDate (YYYY-MM-01)
as_month_idate <- function(x) {
  # Accept Date, IDate, POSIXct
  d <- as.IDate(x)
  data.table::as.IDate(sprintf("%04d-%02d-01", year(d), month(d)))
}

# -------------------------
# 6) Assertions (fail fast)
# -------------------------
assert_cols <- function(dt, cols, name = deparse(substitute(dt))) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(sprintf("Missing cols in %s: %s", name, paste(miss, collapse = ", ")))
  invisible(TRUE)
}

assert_unique <- function(dt, key_cols, name = deparse(substitute(dt))) {
  assert_cols(dt, key_cols, name = name)
  dup_n <- dt[, .N, by = key_cols][N > 1L, .N]
  if (dup_n > 0L) {
    example <- dt[, .N, by = key_cols][N > 1L][1:min(10, .N)]
    stop(sprintf(
      "%s is NOT unique on [%s]. Duplicate keys found: %d. Example:\n%s",
      name, paste(key_cols, collapse = ", "), dup_n,
      paste(capture.output(print(example)), collapse = "\n")
    ))
  }
  invisible(TRUE)
}

assert_no_na <- function(dt, cols, name = deparse(substitute(dt))) {
  assert_cols(dt, cols, name = name)
  bad <- sapply(cols, function(cc) sum(is.na(dt[[cc]])))
  bad <- bad[bad > 0]
  if (length(bad)) stop(sprintf("%s has NA in: %s", name,
                                paste(names(bad), bad, sep = "=", collapse = ", ")))
  invisible(TRUE)
}

# -------------------------
# 7) Audit writers
# -------------------------
audit_missing <- function(dt, name = deparse(substitute(dt))) {
  data.table(
    dataset = name,
    var = names(dt),
    n = nrow(dt),
    n_na = vapply(dt, function(x) sum(is.na(x)), integer(1)),
    share_na = vapply(dt, function(x) mean(is.na(x)), numeric(1))
  )[order(-share_na)]
}

write_audit_dt <- function(dt, file_stub, subdir = DIRS$audit) {
  .ensure_dir(subdir)
  out <- file.path(subdir, paste0(file_stub, "_", .now_stamp(), ".csv"))
  fwrite(dt, out)
  log_msg("Wrote audit: ", out)
  invisible(out)
}

write_audit_text <- function(lines, file_stub, subdir = DIRS$audit) {
  .ensure_dir(subdir)
  out <- file.path(subdir, paste0(file_stub, "_", .now_stamp(), ".txt"))
  writeLines(lines, con = out, useBytes = TRUE)
  log_msg("Wrote audit: ", out)
  invisible(out)
}

# -------------------------
# 8) Session stamp (helps debugging reproducibility)
# -------------------------
write_session_info <- function() {
  si <- c(
    paste0("PROJECT_ROOT: ", PROJECT_ROOT),
    paste0("R.version: ", R.version.string),
    paste0("data.table: ", as.character(packageVersion("data.table"))),
    paste0("Threads (data.table): ", getDTthreads())
  )
  write_audit_text(si, "session_info")
}

# Run once per pipeline start (comment out if you don’t want this every time)
write_session_info()

log_msg("Config loaded. Working directory set to: ", PROJECT_ROOT)
