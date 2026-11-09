# ============================================================
# R21_robustness_shocks.R — Additional quasi-experiments (robustness)
#
# Runs the four additional shocks from Section 3.4:
#   1. Incirlik restrictions (Syria/Iraq, Oct 2019)
#   2. Starlink geofencing (Ukraine, 2022)
#   3. Trump active hostilities expansion (Somalia/Yemen, Mar 2017)
#   4. Solar storm (Ukraine, May 2024) — monthly proxy
#
# Same specification as R17: triple interaction + riots placebo.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
#         acled_pg_month.rds
# Output: robustness shock results table
#
# Run AFTER R16 + R12. Can run after R17 (core results first).
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 21: Robustness shocks")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")

dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

# Merge ACLED riots
if (file.exists(P_ACLED)) {
  acled <- readRDS(P_ACLED); setDT(acled)
  acled[, pg_id := as.integer(pg_id)]
  acled[, acled_riots_protests := acled_riots + acled_protests]
  acled_m <- acled[, .(pg_id, month, acled_riots_protests)]
  if ("acled_riots_protests" %in% names(dt)) dt[, acled_riots_protests := NULL]
  dt <- merge(dt, acled_m, by = c("pg_id", "month"), all.x = TRUE)
  dt[is.na(acled_riots_protests), acled_riots_protests := 0L]
  dt[, log1p_riots := log1p(acled_riots_protests)]
} else {
  dt[, log1p_riots := NA_real_]
}

CTRL <- c("nlights_calib_mean", "nlights_mean", "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── HELPER ───────────────────────────────────────────────────
run_robustness_shock <- function(name, shock_col, theatres, extra_spatial = NULL) {
  cat("\n", strrep("-", 50), "\n")
  cat("Robustness shock: ", name, "\n")
  cat(strrep("-", 50), "\n")
  
  sub <- dt[theatre_id %in% theatres]
  if (nrow(sub) == 0) {
    log_msg("No data for theatres: ", paste(theatres, collapse = ", "))
    return(NULL)
  }
  
  results <- list()
  
  # Triple interaction: post_eli × eli_feas × shock
  f_main <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:", shock_col,
    " + ", ctrl_str, " | ", FE
  ))
  m_main <- tryCatch(feols(f_main, data = sub, cluster = ~theatre_id, notes = FALSE),
                     error = function(e) { log_msg("ERROR: ", e$message); NULL })
  
  # Riots placebo
  m_riots <- NULL
  if (!all(is.na(sub$log1p_riots))) {
    f_riots <- as.formula(paste0(
      "log1p_riots ~ post_eli:eli_feas + post_eli:eli_feas:", shock_col,
      " + ", ctrl_str, " | ", FE
    ))
    m_riots <- tryCatch(feols(f_riots, data = sub, cluster = ~theatre_id, notes = FALSE),
                        error = function(e) NULL)
  }
  
  # Spatial interaction (if provided)
  m_spatial <- NULL
  if (!is.null(extra_spatial) && extra_spatial %in% names(sub)) {
    f_sp <- as.formula(paste0(
      "log1p_y ~ ", shock_col, ":", extra_spatial, " + ", shock_col, " + ", extra_spatial,
      " + ", ctrl_str, " | pg_id + month_int"
    ))
    m_spatial <- tryCatch(feols(f_sp, data = sub, cluster = ~pg_id, notes = FALSE),
                          error = function(e) NULL)
  }
  
  extract <- function(m, label) {
    if (is.null(m)) return(data.table(shock = name, spec = label, coef = NA, se = NA, p = NA, nobs = NA))
    ct <- as.data.table(coeftable(m), keep.rownames = "term")
    hit <- ct[grepl(shock_col, term) | grepl("eli_feas", term)]
    if (nrow(hit) == 0) return(data.table(shock = name, spec = label, coef = NA, se = NA, p = NA, nobs = nobs(m)))
    # Return the most informative row (triple interaction if present)
    triple <- hit[grepl(shock_col, term)]
    if (nrow(triple) > 0) hit <- triple
    data.table(shock = name, spec = label,
               coef = round(hit$Estimate[1], 5), se = round(hit$`Std. Error`[1], 5),
               p = round(hit$`Pr(>|t|)`[1], 4), nobs = nobs(m))
  }
  
  out <- rbindlist(list(
    extract(m_main, "Triple interaction"),
    extract(m_riots, "Riots placebo"),
    if (!is.null(m_spatial)) extract(m_spatial, paste0("Spatial: ", extra_spatial)) else NULL
  ), fill = TRUE)
  
  print(out)
  
  # Save models
  out_dir <- file.path(DIRS$out_paper1, "models")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  slug <- gsub("[^a-zA-Z0-9]", "_", name)
  if (!is.null(m_main)) saveRDS(m_main, file.path(out_dir, paste0("R21_", slug, "_main.rds")))
  if (!is.null(m_riots)) saveRDS(m_riots, file.path(out_dir, paste0("R21_", slug, "_riots.rds")))
  
  out
}

# ══════════════════════════════════════════════════════════════
# RUN ROBUSTNESS SHOCKS
# ══════════════════════════════════════════════════════════════

all_results <- list()

# 1. INCIRLIK
all_results$incirlik <- run_robustness_shock(
  "Incirlik_restrictions", "shock_incirlik",
  theatres = c("Syria", "Iraq"),
  extra_spatial = "incirlik_was_closer"
)

# 2. STARLINK
all_results$starlink <- run_robustness_shock(
  "Starlink_geofencing", "shock_starlink",
  theatres = "Ukraine"
)

# 3. TRUMP 2017
# Note: This tests the "ON" direction — expanding strike authority
# Prediction: β2 > 0 (violence INCREASES when authority expands)
all_results$trump2017 <- run_robustness_shock(
  "Trump_active_hostilities", "shock_trump2017",
  theatres = c("Somalia", "Yemen (North Yemen)")
)

# 4. SOLAR STORM (monthly proxy — for a proper test, need ACLED daily)
all_results$solar <- run_robustness_shock(
  "Solar_storm_May2024", "shock_solar",
  theatres = "Ukraine"
)

# ══════════════════════════════════════════════════════════════
# COMPILE
# ══════════════════════════════════════════════════════════════
final <- rbindlist(all_results, fill = TRUE)

cat("\n============================================================\n")
cat("ALL ROBUSTNESS SHOCK RESULTS\n")
cat("============================================================\n")
print(final)

# Interpretation notes
cat("\n── NOTES ──\n")
cat("Incirlik: β2 < 0 expected (ISR removed from N Syria).\n")
cat("Starlink: β2 < 0 expected (connectivity severed).\n")
cat("Trump 2017: β2 > 0 expected (authority EXPANDED — 'on' direction).\n")
cat("Solar: β2 < 0 expected (connectivity disrupted). Monthly proxy — noisy.\n")
cat("All riots placebos should be null.\n")

write_audit_dt(final, "R21_robustness_shocks")
log_msg("Step 21 complete.")
