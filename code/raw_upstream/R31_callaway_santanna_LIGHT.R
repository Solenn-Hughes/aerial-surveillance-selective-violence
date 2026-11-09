# ============================================================================
# R31_callaway_santanna_LIGHT.R
# Lightweight version to avoid memory crash
# ============================================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

# Skip if 'did' package not installed
if (!requireNamespace("did", quietly = TRUE)) {
  cat("CS2021 package not installed. Skipping R31.\n")
  cat("Install with: install.packages('did')\n")
  quit(save = "no")
}

library(did)
source("R00_config.R")
log_msg("R31: Callaway-Sant'Anna (lightweight)")

OUT <- file.path("outputs", "paper1", "callaway_santanna")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── LOAD & SUBSET ──
dt <- readRDS(file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds"))
setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_type3 := log1p(events_type3)]

# LIGHTER: Just Afghanistan (single shock, cleanest)
dt_afg <- dt[theatre_id == "Afghanistan"]
log_msg("Afghanistan only: ", nrow(dt_afg), " obs, ", uniqueN(dt_afg$pg_id), " cells")

# Binary treatment: above-median eli_feas
dt_afg[, median_feas := median(eli_feas, na.rm = TRUE)]
dt_afg[, treated := eli_feas > median_feas]

# Group = shock month for treated, 0 for never-treated
SHOCK_MONTH <- 2021L * 12L + 7L  # July 2021
dt_afg[, G := ifelse(treated, SHOCK_MONTH, 0L)]

# LIGHTER: ±12 month window only (not 24)
dt_afg[, months_from_shock := month_int - SHOCK_MONTH]
dt_afg <- dt_afg[months_from_shock >= -12L & months_from_shock <= 12L]
log_msg("12-month window: ", nrow(dt_afg), " obs")

# Ensure unique
dt_afg <- unique(dt_afg, by = c("pg_id", "month_int"))

cat("\n── Running att_gt() (light spec) ──\n")
att_res <- tryCatch({
  att_gt(
    yname         = "log1p_type3",
    tname         = "month_int",
    idname        = "pg_id",
    gname         = "G",
    data          = dt_afg,
    control_group = "nevertreated",
    bstrap        = FALSE,
    biters        = 50,  # Reduced from 200
    clustervars   = "pg_id",
    allow_unbalanced_panel = TRUE
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

if (is.null(att_res)) {
  cat("\nCS2021 failed. This is optional - your main results stand.\n")
  log_msg("R31 skipped (computational limits)")
  quit(save = "no")
}

# ── AGGREGATIONS ──
cat("\n── Overall ATT ──\n")
agg_overall <- aggte(att_res, type = "simple", na.rm = TRUE)
print(agg_overall)

cat(sprintf("\nCS2021 ATT = %.4f (SE = %.4f)\n",
            agg_overall$overall.att, agg_overall$overall.se))

# Event study
cat("\n── Event-study ──\n")
agg_es <- aggte(att_res, type = "dynamic", na.rm = TRUE)
print(summary(agg_es))

es_dt <- data.table(
  event_time = agg_es$egt,
  ATT = agg_es$att.egt,
  SE = agg_es$se.egt,
  p = 2 * pnorm(-abs(agg_es$att.egt / agg_es$se.egt))
)

fwrite(es_dt, file.path(OUT, "CS2021_event_study.csv"))

cat("\nComparison to your estimates:\n")
cat("  TWFE (8 theatres):      β = 0.045\n")
cat("  Stacked shocks:         β = -0.037\n")
cat("  CS2021 (Afghanistan):   ATT =", round(agg_overall$overall.att, 3), "\n")

log_msg("R31 complete: ", OUT)
