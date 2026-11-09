# ============================================================
# R09h_m0boot_pretrend_controlled.R
# Feasibility-differential trend robustness under PRIMARY (m0_boot) spec
# - FE: pg_id + theatre_id^month_int
# - VCOV: ~theatre_id
# - Controls: 6 surviving controls
# Outputs:
#   outputs/paper1/audit/R09h_m0boot_pretrend_controlled_*.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

# -------------------------
# 0) Config + helpers
# -------------------------
cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!exists("write_audit_dt")) {
  write_audit_dt <- function(dt, stem) {
    out_dir <- file.path("outputs", "paper1", "audit")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    f <- file.path(out_dir, paste0(stem, "_", ts, ".csv"))
    fwrite(as.data.table(dt), f)
    message("Wrote audit: ", f)
    invisible(f)
  }
}

assert_file <- function(p) if (!file.exists(p)) stop("INTEGRITY FAIL: Missing file: ", p)
assert_cols2 <- function(dt, cols, nm = "dt") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("INTEGRITY FAIL: ", nm, " missing cols: ", paste(miss, collapse = ", "))
}

# -------------------------
# 1) Load panel
# -------------------------
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
assert_file(P_PANEL)

dt <- readRDS(P_PANEL)
setDT(dt)

if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]
dt[, theatre_id := as.character(theatre_id)]

req <- c("pg_id", "theatre_id", "month_int", "post_eli", "eli_feas", "events_type3")
assert_cols2(dt, req, "ALLCW_posteli_v2")

# Controls (m0_boot)
CTRL_M0 <- c("nlights_calib_mean", "nlights_mean", "nlights_min", "prec_gpcp", "temp", "water_ih")
missing_ctrl <- setdiff(CTRL_M0, names(dt))
if (length(missing_ctrl)) stop("INTEGRITY FAIL: missing m0_boot controls in panel: ", paste(missing_ctrl, collapse = ", "))

# trend term
dt[, month_int2 := month_int^2]

# -------------------------
# 2) Fit models
# -------------------------
fit_one <- function(extra_rhs, label) {
  rhs <- paste(c("post_eli * eli_feas", CTRL_M0, extra_rhs), collapse = " + ")
  fml <- as.formula(paste0("log1p(events_type3) ~ ", rhs, " | pg_id + theatre_id^month_int"))
  m <- feols(fml, data = dt, vcov = ~ theatre_id)
  
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  did <- ct[term %in% c("post_eli:eli_feas", "eli_feas:post_eli")][1]
  if (is.null(did$Estimate)) stop("INTEGRITY FAIL: did term not found.")
  
  out <- data.table(
    spec = label,
    coef = as.numeric(did$Estimate),
    se   = as.numeric(did$`Std. Error`),
    t    = as.numeric(did$`t value`),
    p    = as.numeric(did$`Pr(>|t|)`),
    nobs = as.integer(nobs(m))
  )
  
  if ("eli_feas:month_int" %in% ct$term) {
    tr <- ct[term == "eli_feas:month_int"][1]
    out[, trend_lin := as.numeric(tr$Estimate)]
    out[, trend_lin_p := as.numeric(tr$`Pr(>|t|)`)]
  }
  
  if ("eli_feas:month_int2" %in% ct$term) {
    tr2 <- ct[term == "eli_feas:month_int2"][1]
    out[, trend_quad := as.numeric(tr2$Estimate)]
    out[, trend_quad_p := as.numeric(tr2$`Pr(>|t|)`)]
  }
  
  out
}

res <- rbindlist(list(
  fit_one(extra_rhs = character(), label = "M1 baseline"),
  fit_one(extra_rhs = c("eli_feas:month_int"), label = "M2 + eli_feas×month_int"),
  fit_one(extra_rhs = c("eli_feas:month_int", "eli_feas:month_int2"), label = "M3 + lin+quad")
), fill = TRUE)

write_audit_dt(res, "R09h_m0boot_pretrend_controlled")
print(res)

cat("\n✅ Done. Use the printed table to replace Table 5 in your results document.\n")