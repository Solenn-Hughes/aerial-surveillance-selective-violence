# ============================================================
# R18_placebos.R — Placebo tests
#
# Two placebos (outcomes/environments drones should NOT affect):
#   (a) Riots + protests as outcome (drones cannot cause riots)
#   (b) Dense canopy sample split (drones cannot see through forest)
#
# Each placebo is run for:
#   - Core DiD spec (post_eli × eli_feas)
#   - Each of the three core quasi-experiments
#
# If any placebo β is significant, the design has a problem.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds (from R16)
#         acled_pg_month.rds (from R12)
# Output: placebo results table
#
# Run AFTER R16 + R12.
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 18: Placebo tests")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")

dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

# Merge ACLED riots
if (file.exists(P_ACLED)) {
  acled <- readRDS(P_ACLED); setDT(acled)
  acled[, pg_id := as.integer(pg_id)]
  acled[, acled_riots_protests := acled_riots + acled_protests]
  acled_m <- acled[, .(pg_id, month, acled_riots_protests)]
  if ("acled_riots_protests" %in% names(dt)) dt[, acled_riots_protests := NULL]
  dt <- merge(dt, acled_m, by = c("pg_id", "month"), all.x = TRUE)
  dt[is.na(acled_riots_protests), acled_riots_protests := 0L]
} else {
  stop("ACLED not found. Download and run R12 first.")
}

dt[, log1p_y := log1p(events_type3)]
dt[, log1p_riots := log1p(acled_riots_protests)]

CTRL <- c("nlights_calib_mean", "nlights_mean", "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── HELPER ───────────────────────────────────────────────────
run_placebo <- function(label, outcome, fml_rhs, data) {
  f <- as.formula(paste0(outcome, " ~ ", fml_rhs, " + ", ctrl_str, " | ", FE))
  m <- tryCatch(feols(f, data = data, cluster = ~theatre_id, notes = FALSE),
                error = function(e) { log_msg("ERROR: ", label, " — ", e$message); NULL })
  if (is.null(m)) return(data.table(test = label, term = NA, coef = NA, se = NA, p = NA, nobs = NA))
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  hit <- ct[grepl("eli_feas", term)]
  if (nrow(hit) == 0) return(data.table(test = label, term = "none found", coef = NA, se = NA, p = NA, nobs = nobs(m)))
  data.table(test = label, term = hit$term[1],
             coef = round(hit$Estimate[1], 5), se = round(hit$`Std. Error`[1], 5),
             p = round(hit$`Pr(>|t|)`[1], 4), nobs = nobs(m))
}

results <- list()

# ══════════════════════════════════════════════════════════════
# PLACEBO A: RIOTS AS OUTCOME
# ══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("PLACEBO A: Riots/protests as outcome\n")
cat("============================================================\n")

# Core DiD with riots
results$riots_core <- run_placebo("Riots: Core DiD", "log1p_riots",
                                   "post_eli:eli_feas", dt)

# Bagram shock with riots (Afghanistan only)
results$riots_bagram <- run_placebo("Riots: Bagram shock", "log1p_riots",
                                     "post_eli:eli_feas + post_eli:eli_feas:shock_bagram",
                                     dt[theatre_id == "Afghanistan"])

# Somalia withdrawal with riots
results$riots_somalia <- run_placebo("Riots: Somalia withdrawal", "log1p_riots",
                                      "post_eli:eli_feas + post_eli:eli_feas:shock_somalia_wd",
                                      dt[theatre_id == "Somalia"])

# Biden pause with riots
results$riots_biden <- run_placebo("Riots: Biden pause", "log1p_riots",
                                    "post_eli:eli_feas + post_eli:eli_feas:shock_biden_pause",
                                    dt[theatre_id %in% c("Somalia", "Yemen (North Yemen)")])

# ══════════════════════════════════════════════════════════════
# PLACEBO B: CANOPY SPLIT (>40% forest = should be null)
# ══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("PLACEBO B: Canopy split (dense forest >40% vs open <20%)\n")
cat("============================================================\n")

# Identify canopy variable
canopy_col <- NULL
if ("forest_gc" %in% names(dt)) canopy_col <- "forest_gc"
if ("forest_ih" %in% names(dt) && is.null(canopy_col)) canopy_col <- "forest_ih"

if (!is.null(canopy_col)) {
  # Dense forest (>40% canopy) — β should be null
  dense <- dt[get(canopy_col) > 0.40]
  # Open terrain (<20% canopy) — β should be positive
  open  <- dt[get(canopy_col) < 0.20]
  
  log_msg("Dense forest cells: ", uniqueN(dense$pg_id), " (", nrow(dense), " cell-months)")
  log_msg("Open terrain cells: ", uniqueN(open$pg_id), " (", nrow(open), " cell-months)")
  
  results$canopy_dense <- run_placebo("Canopy: Dense (>40%) — should be NULL", "log1p_y",
                                       "post_eli:eli_feas", dense)
  results$canopy_open  <- run_placebo("Canopy: Open (<20%) — should be POSITIVE", "log1p_y",
                                       "post_eli:eli_feas", open)
  
  # Also run canopy split for Bagram specifically
  afg_dense <- dt[theatre_id == "Afghanistan" & get(canopy_col) > 0.40]
  afg_open  <- dt[theatre_id == "Afghanistan" & get(canopy_col) < 0.20]
  
  if (nrow(afg_dense) > 100) {
    results$canopy_bagram_dense <- run_placebo("Canopy: Bagram + Dense", "log1p_y",
                                                "post_eli:eli_feas + post_eli:eli_feas:shock_bagram", afg_dense)
  }
  if (nrow(afg_open) > 100) {
    results$canopy_bagram_open <- run_placebo("Canopy: Bagram + Open", "log1p_y",
                                               "post_eli:eli_feas + post_eli:eli_feas:shock_bagram", afg_open)
  }
} else {
  log_msg("WARNING: No canopy variable found. Skipping canopy placebo.")
}

# ══════════════════════════════════════════════════════════════
# COMPILE
# ══════════════════════════════════════════════════════════════
out <- rbindlist(results, fill = TRUE)

cat("\n============================================================\n")
cat("ALL PLACEBO RESULTS\n")
cat("============================================================\n")
print(out)

# Flag any significant placebos (p < 0.05)
sig <- out[!is.na(p) & p < 0.05]
if (nrow(sig) > 0) {
  cat("\n⚠️  WARNING: Significant placebos detected!\n")
  print(sig)
  cat("Investigate these before proceeding.\n")
} else {
  cat("\n✓  No significant placebos. Design looks clean.\n")
}

write_audit_dt(out, "R18_placebo_results")

log_msg("Step 18 complete.")
