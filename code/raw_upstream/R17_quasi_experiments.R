# ============================================================
# R17_quasi_experiments.R — Three core quasi-experiments
#
# 1. Bagram closure (Afghanistan, July 2021)
#    - Triple interaction: post_eli × eli_feas × shock_bagram
#    - Diff-in-RD at MQ-9 range boundary
# 2. Somalia troop withdrawal (Dec 2020 – Jan 2021)
#    - Triple interaction: post_eli × eli_feas × shock_somalia_wd
#    - Interior-to-interior comparison (near-FOB vs far-from-FOB)
# 3. Biden drone pause (Somalia + Yemen, Jan 2021)
#    - Triple interaction: post_eli × eli_feas × shock_biden_pause
#    - Tests information-only vs information+strike mechanism
#
# For each: reports β with/without controls, riots placebo.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds (from R16)
#         acled_pg_month.rds (from R12, for riots placebo)
# Output: audit tables + model objects
#
# Run AFTER R16.
# ============================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 17: Core quasi-experiments")

# ── LOAD DATA ────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")

if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL, " — run R16 first.")
dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

# Merge ACLED riots for placebo
if (file.exists(P_ACLED)) {
  acled <- readRDS(P_ACLED); setDT(acled)
  acled[, pg_id := as.integer(pg_id)]
  # Compute placebo outcome: riots + protests
  acled[, acled_riots_protests := acled_riots + acled_protests]
  acled_merge <- acled[, .(pg_id, month, acled_riots_protests, acled_state_vac, acled_vac)]
  
  # Remove if already present
  for (cc in c("acled_riots_protests", "acled_state_vac", "acled_vac")) {
    if (cc %in% names(dt)) dt[, (cc) := NULL]
  }
  dt <- merge(dt, acled_merge, by = c("pg_id", "month"), all.x = TRUE)
  # Fill NAs with 0
  for (cc in c("acled_riots_protests", "acled_state_vac", "acled_vac")) {
    dt[is.na(get(cc)), (cc) := 0L]
  }
  dt[, log1p_riots := log1p(acled_riots_protests)]
  log_msg("ACLED riots merged. Non-zero months: ", dt[acled_riots_protests > 0, .N])
} else {
  log_msg("WARNING: ACLED not found. Riots placebo will be skipped.")
  dt[, log1p_riots := NA_real_]
}

# Outcome
dt[, log1p_y := log1p(events_type3)]

# Controls
CTRL <- c("nlights_calib_mean", "nlights_mean", "prec_gpcp", "temp", "water_ih",
          "pop_gpw_sum", "forest_ih", "urban_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

# FE
FE <- "pg_id + theatre_id^month_int"

log_msg("Controls: ", ctrl_str)
log_msg("FE: ", FE)

# ── HELPER ───────────────────────────────────────────────────
run_shock <- function(shock_name, shock_col, theatre_filter = NULL, extra_note = "") {
  cat("\n", strrep("=", 60), "\n")
  cat("SHOCK: ", shock_name, "\n")
  cat(strrep("=", 60), "\n")
  
  # Subset if needed
  if (!is.null(theatre_filter)) {
    sub <- dt[theatre_id %in% theatre_filter]
  } else {
    sub <- dt
  }
  
  results <- list()
  
  # ── Model 1: Triple interaction (no extra controls) ──
  f1 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:", shock_col,
    " | ", FE
  ))
  m1 <- tryCatch(feols(f1, data = sub, cluster = ~theatre_id, notes = FALSE),
                 error = function(e) { log_msg("ERROR m1: ", e$message); NULL })
  
  # ── Model 2: Triple interaction + controls ──
  f2 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:", shock_col,
    " + ", ctrl_str, " | ", FE
  ))
  m2 <- tryCatch(feols(f2, data = sub, cluster = ~theatre_id, notes = FALSE),
                 error = function(e) { log_msg("ERROR m2: ", e$message); NULL })
  
  # ── Model 3: Riots placebo ──
  m3 <- NULL
  if (!all(is.na(sub$log1p_riots))) {
    f3 <- as.formula(paste0(
      "log1p_riots ~ post_eli:eli_feas + post_eli:eli_feas:", shock_col,
      " + ", ctrl_str, " | ", FE
    ))
    m3 <- tryCatch(feols(f3, data = sub, cluster = ~theatre_id, notes = FALSE),
                   error = function(e) { log_msg("ERROR m3 riots: ", e$message); NULL })
  }
  
  # ── Extract results ──
  extract <- function(m, label) {
    if (is.null(m)) return(data.table(shock = shock_name, spec = label,
                                       term = NA, coef = NA, se = NA, p = NA, nobs = NA))
    ct <- as.data.table(coeftable(m), keep.rownames = "term")
    # Find the triple interaction term
    triple <- ct[grepl(shock_col, term)]
    if (nrow(triple) == 0) {
      # Try the base interaction
      triple <- ct[grepl("eli_feas", term)]
    }
    rbindlist(lapply(seq_len(nrow(triple)), function(i) {
      data.table(shock = shock_name, spec = label,
                 term = triple$term[i],
                 coef = round(triple$Estimate[i], 5),
                 se   = round(triple$`Std. Error`[i], 5),
                 p    = round(triple$`Pr(>|t|)`[i], 4),
                 nobs = nobs(m))
    }))
  }
  
  results <- rbindlist(list(
    extract(m1, "No controls"),
    extract(m2, "With controls"),
    extract(m3, "Riots placebo")
  ), fill = TRUE)
  
  cat("\nResults:\n")
  print(results)
  
  # Save models
  out_dir <- file.path(DIRS$out_paper1, "models")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(m1)) saveRDS(m1, file.path(out_dir, paste0("R17_", shock_name, "_m1.rds")))
  if (!is.null(m2)) saveRDS(m2, file.path(out_dir, paste0("R17_", shock_name, "_m2.rds")))
  if (!is.null(m3)) saveRDS(m3, file.path(out_dir, paste0("R17_", shock_name, "_riots.rds")))
  
  results
}

# ══════════════════════════════════════════════════════════════
# RUN CORE QUASI-EXPERIMENTS
# ══════════════════════════════════════════════════════════════

all_results <- list()

# 1. BAGRAM
all_results$bagram <- run_shock(
  "Bagram", "shock_bagram", theatre_filter = "Afghanistan"
)

# 1b. BAGRAM DIFF-IN-RD
cat("\n", strrep("=", 60), "\n")
cat("BAGRAM DIFF-IN-RD (at MQ-9 range boundary)\n")
cat(strrep("=", 60), "\n")

afg <- dt[theatre_id == "Afghanistan"]
rd_results <- list()

for (cutoff_km in c(800, 1000, 1200)) {
  afg[, in_range := as.integer(dist_bagram_km <= cutoff_km)]
  
  f_rd <- as.formula(paste0(
    "log1p_y ~ shock_bagram:in_range + shock_bagram + in_range + ",
    ctrl_str, " | pg_id + month_int"
  ))
  m_rd <- tryCatch(
    feols(f_rd, data = afg, cluster = ~pg_id, notes = FALSE),
    error = function(e) { log_msg("ERROR RD ", cutoff_km, ": ", e$message); NULL }
  )
  
  if (!is.null(m_rd)) {
    ct <- as.data.table(coeftable(m_rd), keep.rownames = "term")
    hit <- ct[grepl("shock_bagram.*in_range|in_range.*shock_bagram", term)]
    if (nrow(hit) > 0) {
      rd_results[[as.character(cutoff_km)]] <- data.table(
        cutoff_km = cutoff_km,
        term = hit$term[1],
        coef = round(hit$Estimate[1], 5),
        se   = round(hit$`Std. Error`[1], 5),
        p    = round(hit$`Pr(>|t|)`[1], 4),
        nobs = nobs(m_rd)
      )
    }
  }
}

if (length(rd_results)) {
  rd_out <- rbindlist(rd_results)
  cat("\nDiff-in-RD results:\n")
  print(rd_out)
  write_audit_dt(rd_out, "R17_bagram_diff_in_RD")
}

# 2. SOMALIA WITHDRAWAL
all_results$somalia <- run_shock(
  "Somalia_withdrawal", "shock_somalia_wd", theatre_filter = "Somalia"
)

# 2b. SOMALIA INTERIOR-TO-INTERIOR
cat("\n", strrep("=", 60), "\n")
cat("SOMALIA INTERIOR-TO-INTERIOR (near-FOB vs far-from-FOB)\n")
cat(strrep("=", 60), "\n")

som <- dt[theatre_id == "Somalia"]
f_som_int <- as.formula(paste0(
  "log1p_y ~ shock_somalia_wd:somalia_near_fob + shock_somalia_wd + somalia_near_fob + ",
  ctrl_str, " | pg_id + month_int"
))
m_som_int <- tryCatch(
  feols(f_som_int, data = som, cluster = ~pg_id, notes = FALSE),
  error = function(e) { log_msg("ERROR Somalia int: ", e$message); NULL }
)
if (!is.null(m_som_int)) {
  cat("\nSomalia interior-to-interior:\n")
  print(coeftable(m_som_int))
  saveRDS(m_som_int, file.path(DIRS$out_paper1, "models", "R17_Somalia_interior.rds"))
}

# 3. BIDEN PAUSE
all_results$biden <- run_shock(
  "Biden_pause", "shock_biden_pause",
  theatre_filter = c("Somalia", "Yemen (North Yemen)")
)

# ══════════════════════════════════════════════════════════════
# COMPILE AND SAVE
# ══════════════════════════════════════════════════════════════

final <- rbindlist(all_results, fill = TRUE)
cat("\n", strrep("=", 60), "\n")
cat("ALL QUASI-EXPERIMENT RESULTS\n")
cat(strrep("=", 60), "\n")
print(final)

write_audit_dt(final, "R17_quasi_experiments_all")

log_msg("Step 17 complete.")
