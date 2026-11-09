# ============================================================
# R22_cross_validation_acled.R
# Cross-dataset validation: replicate core β using ACLED
#
# If β replicates across UCDP-GED and ACLED, the result is
# not an artefact of one dataset's coding decisions.
#
# Runs the core DiD specification with:
#   (a) UCDP-GED type-3 (main — already done, repeated for comparison)
#   (b) ACLED violence against civilians (acled_vac)
#   (c) ACLED state-initiated VAC (acled_state_vac)
#   (d) ACLED air/drone strikes (acled_airdrone_strike)
#
# Also runs the three core quasi-experiments with ACLED outcome.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
#         acled_pg_month.rds
# Output: cross-validation comparison table
#
# Run AFTER R16 + R12.
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 22: Cross-dataset validation (ACLED)")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")

if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL)
if (!file.exists(P_ACLED)) stop("Missing: ", P_ACLED, " — run R12 first.")

dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

acled <- readRDS(P_ACLED); setDT(acled)
acled[, pg_id := as.integer(pg_id)]

# Merge ACLED outcomes
acled_cols <- c("acled_vac", "acled_state_vac", "acled_airdrone_strike")
acled_merge <- acled[, c("pg_id", "month", acled_cols), with = FALSE]

for (cc in acled_cols) {
  if (cc %in% names(dt)) dt[, (cc) := NULL]
}
dt <- merge(dt, acled_merge, by = c("pg_id", "month"), all.x = TRUE)
for (cc in acled_cols) dt[is.na(get(cc)), (cc) := 0L]

# Create log outcomes
dt[, log1p_ged_type3     := log1p(events_type3)]
dt[, log1p_acled_vac     := log1p(acled_vac)]
dt[, log1p_acled_state   := log1p(acled_state_vac)]
dt[, log1p_acled_drone   := log1p(acled_airdrone_strike)]

CTRL <- c("nlights_calib_mean", "nlights_mean", "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── HELPER ───────────────────────────────────────────────────
run_spec <- function(outcome_col, label, data) {
  f <- as.formula(paste0(outcome_col, " ~ post_eli:eli_feas + ", ctrl_str, " | ", FE))
  m <- tryCatch(feols(f, data = data, cluster = ~theatre_id, notes = FALSE),
                error = function(e) { log_msg("ERROR: ", label, " — ", e$message); NULL })
  if (is.null(m)) return(data.table(outcome = label, coef = NA, se = NA, p = NA, nobs = NA))
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  hit <- ct[grepl("eli_feas", term)]
  if (nrow(hit) == 0) return(data.table(outcome = label, coef = NA, se = NA, p = NA, nobs = nobs(m)))
  data.table(outcome = label,
             coef = round(hit$Estimate[1], 5),
             se   = round(hit$`Std. Error`[1], 5),
             p    = round(hit$`Pr(>|t|)`[1], 4),
             nobs = nobs(m))
}

# ══════════════════════════════════════════════════════════════
# PART 1: CORE DiD WITH DIFFERENT OUTCOMES
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 1: Core DiD across datasets\n")
cat("============================================================\n")

core_results <- rbindlist(list(
  run_spec("log1p_ged_type3",   "UCDP-GED type-3 (main)", dt),
  run_spec("log1p_acled_vac",   "ACLED: all VAC", dt),
  run_spec("log1p_acled_state", "ACLED: state-initiated VAC", dt),
  run_spec("log1p_acled_drone", "ACLED: air/drone strikes", dt)
))

cat("\nCore DiD comparison:\n")
print(core_results)

# Check replication
if (all(!is.na(core_results$coef))) {
  signs_agree <- all(sign(core_results$coef) == sign(core_results$coef[1]))
  if (signs_agree) {
    cat("\n\u2713  All outcomes have the same sign. Cross-validation passed.\n")
  } else {
    cat("\n\u26A0\uFE0F  Signs differ across datasets. Investigate.\n")
  }
}

# ══════════════════════════════════════════════════════════════
# PART 2: QUASI-EXPERIMENTS WITH ACLED OUTCOME
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 2: Quasi-experiments with ACLED outcome\n")
cat("============================================================\n")

run_shock_acled <- function(shock_col, shock_name, theatres) {
  sub <- dt[theatre_id %in% theatres]
  f <- as.formula(paste0(
    "log1p_acled_vac ~ post_eli:eli_feas + post_eli:eli_feas:", shock_col,
    " + ", ctrl_str, " | ", FE
  ))
  m <- tryCatch(feols(f, data = sub, cluster = ~theatre_id, notes = FALSE),
                error = function(e) { log_msg("ERROR: ", e$message); NULL })
  if (is.null(m)) return(data.table(shock = shock_name, coef = NA, se = NA, p = NA, nobs = NA))
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  hit <- ct[grepl(shock_col, term)]
  if (nrow(hit) == 0) return(data.table(shock = shock_name, coef = NA, se = NA, p = NA, nobs = nobs(m)))
  data.table(shock = shock_name,
             coef = round(hit$Estimate[1], 5),
             se   = round(hit$`Std. Error`[1], 5),
             p    = round(hit$`Pr(>|t|)`[1], 4),
             nobs = nobs(m))
}

shock_acled <- rbindlist(list(
  run_shock_acled("shock_bagram", "Bagram (ACLED)", "Afghanistan"),
  run_shock_acled("shock_somalia_wd", "Somalia withdrawal (ACLED)", "Somalia"),
  run_shock_acled("shock_biden_pause", "Biden pause (ACLED)", c("Somalia", "Yemen (North Yemen)"))
))

cat("\nQuasi-experiments with ACLED outcome:\n")
print(shock_acled)

# ══════════════════════════════════════════════════════════════
# COMPILE
# ══════════════════════════════════════════════════════════════

all_out <- rbindlist(list(
  cbind(section = "Core DiD", core_results),
  cbind(section = "Quasi-experiments (ACLED)", shock_acled)
), fill = TRUE)

cat("\n============================================================\n")
cat("ALL CROSS-VALIDATION RESULTS\n")
cat("============================================================\n")
print(all_out)

write_audit_dt(all_out, "R22_cross_validation_acled")

# Save models
out_dir <- file.path(DIRS$out_paper1, "models")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (oc in c("log1p_ged_type3", "log1p_acled_vac", "log1p_acled_state", "log1p_acled_drone")) {
  f <- as.formula(paste0(oc, " ~ post_eli:eli_feas + ", ctrl_str, " | ", FE))
  m <- tryCatch(feols(f, data = dt, cluster = ~theatre_id, notes = FALSE), error = function(e) NULL)
  if (!is.null(m)) saveRDS(m, file.path(out_dir, paste0("R22_", gsub("log1p_", "", oc), ".rds")))
}

log_msg("Step 22 complete.")
