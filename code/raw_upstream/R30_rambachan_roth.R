# ============================================================================
# R30_rambachan_roth_WORKING.R
# Honest DiD sensitivity - FIXED reference period issue
# ============================================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

# Check for HonestDiD
if (!requireNamespace("HonestDiD", quietly = TRUE)) {
  cat("HonestDiD package not found. Skipping R30.\n")
  cat("Optional test - your main results stand without it.\n")
  quit(save = "no")
}
library(HonestDiD)

source("R00_config.R")
log_msg("R30: Honest DiD sensitivity")

OUT <- file.path("outputs", "paper1", "honest_did")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── LOAD ──
dt <- readRDS(file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds"))
setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_type3 := log1p(events_type3)]

# ── AFGHANISTAN ONLY (Bagram shock) ──
cat("\nBagram closure (Afghanistan) - July 2021\n")

BAGRAM_DATE <- as.IDate("2021-07-01")
dt_afg <- dt[theatre_id == "Afghanistan"]
dt_afg[, event_time := as.integer(difftime(month, BAGRAM_DATE, units = "days") / 30.44)]

# Window: ±12 months
dt_afg <- dt_afg[event_time >= -12 & event_time <= 12]

# BINNED PERIODS (fewer bins = more observations per bin)
dt_afg[, et_bin := fcase(
  event_time %in% c(-12:-6), "pre_long",    # Reference period
  event_time %in% c(-5:-3),  "pre_mid",
  event_time %in% c(-2:-1),  "pre_short",
  event_time %in% c(0:2),    "post_short",
  event_time %in% c(3:5),    "post_mid",
  event_time >= 6,           "post_long"
)]

# Check we have observations in each bin
cat("\nObservations per bin:\n")
print(dt_afg[, .N, by = et_bin][order(et_bin)])

dt_afg[, et_bin := factor(et_bin, levels = c(
  "pre_long", "pre_mid", "pre_short",
  "post_short", "post_mid", "post_long"
))]

cat("\nTotal obs:", nrow(dt_afg), "\n")

# ── EVENT STUDY ──
CTRL <- c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(CTRL, names(dt_afg)), collapse = " + ")

cat("\n── Event study ──\n")

f_es <- as.formula(paste0(
  "log1p_type3 ~ i(et_bin, eli_feas, ref = 'pre_long') + ",
  ctrl_str, " | pg_id + month_int"
))

m_es <- feols(f_es, data = dt_afg, vcov = ~pg_id, notes = FALSE)

cat("\nEvent study coefficients:\n")
ct <- coeftable(m_es)
print(ct[grepl("et_bin", rownames(ct)), ])

# ── HONEST DID ──
coefs <- coef(m_es)
vcov_mat <- vcov(m_es)

et_terms <- grep("et_bin.*eli_feas", names(coefs), value = TRUE)
beta <- coefs[et_terms]
sigma <- vcov_mat[et_terms, et_terms]

# Map bins to approximate event time
et_map <- data.table(
  term = et_terms,
  period = c(-4, -1.5, 1, 4, 7)  # Midpoints
)

pre_idx <- which(et_map$period < 0)
post_idx <- which(et_map$period >= 0)

cat("\n── Honest DiD Sensitivity ──\n")
cat("Testing robustness to parallel trends violations\n\n")

M_grid <- c(0, 0.01, 0.02, 0.05, 0.1)

results <- rbindlist(lapply(M_grid, function(M_val) {
  
  robust <- tryCatch({
    createSensitivityResults(
      betahat = beta,
      sigma = sigma,
      numPrePeriods = length(pre_idx),
      numPostPeriods = length(post_idx),
      method = "C-LF",
      Mvec = M_val
    )
  }, error = function(e) NULL)
  
  if (is.null(robust)) {
    return(data.table(M = M_val, lb = NA, ub = NA, note = "failed"))
  }
  
  data.table(
    M = M_val,
    lb = robust$lb,
    ub = robust$ub,
    significant = !(robust$lb <= 0 & robust$ub >= 0)
  )
}))

cat("\nSensitivity Results:\n")
cat("M = max allowed violation of parallel trends\n\n")
print(results)

fwrite(results, file.path(OUT, "honest_did_sensitivity.csv"))

# ── INTERPRETATION ──
cat("\n── INTERPRETATION ──\n")
if (all(results$significant, na.rm = TRUE)) {
  cat("✓ Effect robust across all M values\n")
  cat("  → Result withstands PT violations\n")
} else {
  first_null <- min(results[significant == FALSE, M], na.rm = TRUE)
  cat("✗ Effect becomes null at M =", first_null, "\n")
  cat("  → Sensitive to PT violations\n")
  cat("\nNote: This suggests the Bagram effect may be driven by\n")
  cat("pre-existing trends. Consider this when interpreting.\n")
}

log_msg("R30 complete: ", OUT)