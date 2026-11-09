# ============================================================
# R14b_info_substitution.R
# Test the INFORMATION SUBSTITUTION mechanism (Group E)
#
# Theory: drones substitute for absent informant networks.
#   → effect should be LARGEST where the state was most
#     information-blind pre-acquisition.
#
# Information deficit index constructed from:
#   (a) Travel time to nearest city (ttime_mean) — PRIO-GRID
#   (b) Nightlights (nlights_mean) — inverse = dark = blind
#   (c) Pre-treatment ACLED event density — low = poorly covered
#
# Test: post_eli × eli_feas × info_deficit
#   Positive coefficient = effect strongest in info deserts
#   = confirms substitution mechanism
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
#         acled_pg_month.rds (for pre-treatment density)
# Output: model results + audit
#
# Run AFTER R16 + R12.
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 14b: Information substitution test (Group E)")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")

if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL, " — run R16 first.")
dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

# ══════════════════════════════════════════════════════════════
# PART 1: BUILD INFORMATION DEFICIT INDEX
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 1: Building information deficit index\n")
cat("============================================================\n")

# ── Component A: Travel time (remoteness) ────────────────────
# High ttime = remote = state is blind = high deficit
# Use cell-level time-invariant ttime_mean from PRIO-GRID
ttime_col <- NULL
for (cc in c("ttime_mean", "ttime_max")) {
  if (cc %in% names(dt) && !all(is.na(dt[[cc]]))) { ttime_col <- cc; break }
}

if (!is.null(ttime_col)) {
  # z-score within dataset (higher = more remote = more deficit)
  cell_ttime <- dt[, .(ttime = mean(get(ttime_col), na.rm = TRUE)), by = pg_id]
  cell_ttime[, z_ttime := (ttime - mean(ttime, na.rm = TRUE)) / sd(ttime, na.rm = TRUE)]
  log_msg("Travel time component: ", ttime_col, " (", sum(!is.na(cell_ttime$z_ttime)), " cells)")
} else {
  log_msg("WARNING: No travel time variable found. Using 0.")
  cell_ttime <- unique(dt[, .(pg_id)])[, `:=`(ttime = NA, z_ttime = 0)]
}

# ── Component B: Nightlights (inverse — dark = blind) ────────
nl_col <- NULL
for (cc in c("nlights_mean", "nlights_calib_mean")) {
  if (cc %in% names(dt) && !all(is.na(dt[[cc]]))) { nl_col <- cc; break }
}

if (!is.null(nl_col)) {
  # Use pre-treatment average (before any theatre's onset)
  # Pre-treatment = before the EARLIEST onset (2001-10-01)
  pre_nl <- dt[month < as.IDate("2001-10-01"),
               .(nl_pre = mean(get(nl_col), na.rm = TRUE)), by = pg_id]
  # If cell has no pre-2001 data, use all data before its own onset
  cells_missing <- setdiff(unique(dt$pg_id), pre_nl$pg_id)
  if (length(cells_missing) > 0) {
    extra <- dt[pg_id %in% cells_missing & post_eli == 0L,
                .(nl_pre = mean(get(nl_col), na.rm = TRUE)), by = pg_id]
    pre_nl <- rbindlist(list(pre_nl, extra), fill = TRUE)
  }
  # Invert: low light = high deficit
  pre_nl[, z_nl_deficit := -(nl_pre - mean(nl_pre, na.rm = TRUE)) / sd(nl_pre, na.rm = TRUE)]
  log_msg("Nightlights deficit component: ", nl_col, " (", sum(!is.na(pre_nl$z_nl_deficit)), " cells)")
} else {
  log_msg("WARNING: No nightlights variable found. Using 0.")
  pre_nl <- unique(dt[, .(pg_id)])[, `:=`(nl_pre = NA, z_nl_deficit = 0)]
}

# ── Component C: Pre-treatment ACLED density (low = poorly monitored) ──
if (file.exists(P_ACLED)) {
  acled <- readRDS(P_ACLED); setDT(acled)
  acled[, pg_id := as.integer(pg_id)]
  
  # Merge post_eli to identify pre-treatment
  onset_map <- unique(dt[, .(pg_id, eli_onset_month)])
  acled <- merge(acled, onset_map, by = "pg_id", all.x = TRUE)
  
  pre_acled <- acled[is.na(eli_onset_month) | month < eli_onset_month,
                     .(acled_pre_density = sum(acled_all_events, na.rm = TRUE) / .N),
                     by = pg_id]
  # Low density = high deficit (invert)
  pre_acled[, z_acled_deficit := -(acled_pre_density - mean(acled_pre_density, na.rm = TRUE)) / 
              sd(acled_pre_density, na.rm = TRUE)]
  log_msg("ACLED pre-treatment density: ", sum(!is.na(pre_acled$z_acled_deficit)), " cells")
} else {
  log_msg("WARNING: ACLED not found. Using 2-component index.")
  pre_acled <- unique(dt[, .(pg_id)])[, `:=`(acled_pre_density = NA, z_acled_deficit = 0)]
}

# ── Combine into index ───────────────────────────────────────
idx <- merge(cell_ttime[, .(pg_id, z_ttime)],
             pre_nl[, .(pg_id, z_nl_deficit)], by = "pg_id", all = TRUE)
idx <- merge(idx, pre_acled[, .(pg_id, z_acled_deficit)], by = "pg_id", all = TRUE)

# Replace NAs with 0 (conservative — doesn't inflate the index)
for (cc in c("z_ttime", "z_nl_deficit", "z_acled_deficit")) {
  idx[is.na(get(cc)), (cc) := 0]
}

# Equal-weighted average
idx[, info_deficit := (z_ttime + z_nl_deficit + z_acled_deficit) / 3]

# Standardise to mean-0, sd-1
idx[, info_deficit := (info_deficit - mean(info_deficit, na.rm = TRUE)) / sd(info_deficit, na.rm = TRUE)]

log_msg("Information deficit index built. Range: ",
        round(min(idx$info_deficit), 2), " to ", round(max(idx$info_deficit), 2))

# Audit
idx_audit <- idx[, .(n = .N,
                     deficit_p10 = quantile(info_deficit, 0.10),
                     deficit_med = median(info_deficit),
                     deficit_p90 = quantile(info_deficit, 0.90),
                     cor_ttime_nl = cor(z_ttime, z_nl_deficit, use = "complete.obs"),
                     cor_ttime_acled = cor(z_ttime, z_acled_deficit, use = "complete.obs"))]
cat("\nIndex summary:\n")
print(idx_audit)
write_audit_dt(idx_audit, "R14b_info_deficit_index_summary")

# ── Merge into panel ─────────────────────────────────────────
if ("info_deficit" %in% names(dt)) dt[, info_deficit := NULL]
dt <- merge(dt, idx[, .(pg_id, info_deficit)], by = "pg_id", all.x = TRUE)
dt[is.na(info_deficit), info_deficit := 0]

# ══════════════════════════════════════════════════════════════
# PART 2: TEST INFORMATION SUBSTITUTION
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 2: Information substitution test\n")
cat("============================================================\n")

CTRL <- c("nlights_calib_mean", "nlights_mean", "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

FE_SPECS <- list(
  theatre_month = "pg_id + theatre_id^month_int"
)

results <- list()

for (fe_name in names(FE_SPECS)) {
  fe_part <- FE_SPECS[[fe_name]]
  
  # ── M0: Baseline (no moderation) ──
  f0 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + ", ctrl_str, " | ", fe_part
  ))
  m0 <- feols(f0, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── M1: + info_deficit moderation ──
  f1 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:info_deficit + ",
    ctrl_str, " | ", fe_part
  ))
  m1 <- feols(f1, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  # Extract
  ct0 <- as.data.table(coeftable(m0), keep.rownames = "term")
  ct1 <- as.data.table(coeftable(m1), keep.rownames = "term")
  
  base_hit <- ct0[grepl("eli_feas", term)]
  triple_hit <- ct1[grepl("info_deficit", term)]
  main_hit <- ct1[grepl("eli_feas", term) & !grepl("info_deficit", term)]
  
  if (nrow(base_hit) > 0) {
    results$baseline <- data.table(
      fe = fe_name, model = "M0: Baseline",
      term = base_hit$term[1],
      coef = round(base_hit$Estimate[1], 5),
      se = round(base_hit$`Std. Error`[1], 5),
      p = round(base_hit$`Pr(>|t|)`[1], 4),
      nobs = nobs(m0)
    )
  }
  
  if (nrow(triple_hit) > 0) {
    results$triple <- data.table(
      fe = fe_name, model = "M1: + info_deficit triple",
      term = triple_hit$term[1],
      coef = round(triple_hit$Estimate[1], 5),
      se = round(triple_hit$`Std. Error`[1], 5),
      p = round(triple_hit$`Pr(>|t|)`[1], 4),
      nobs = nobs(m1)
    )
  }
  
  if (nrow(main_hit) > 0) {
    results$main_after <- data.table(
      fe = fe_name, model = "M1: main effect (after adding triple)",
      term = main_hit$term[1],
      coef = round(main_hit$Estimate[1], 5),
      se = round(main_hit$`Std. Error`[1], 5),
      p = round(main_hit$`Pr(>|t|)`[1], 4),
      nobs = nobs(m1)
    )
  }
  
  # Save models
  out_dir <- file.path(DIRS$out_paper1, "models")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(m0, file.path(out_dir, paste0("R14b_baseline_", fe_name, ".rds")))
  saveRDS(m1, file.path(out_dir, paste0("R14b_info_sub_", fe_name, ".rds")))
}

out <- rbindlist(results, fill = TRUE)

cat("\nResults:\n")
print(out)

cat("\n\u2500\u2500 INTERPRETATION \u2500\u2500\n")
cat("If triple interaction (info_deficit) is POSITIVE:\n")
cat("  \u2192 ELI has a STRONGER effect in information deserts.\n")
cat("  \u2192 Confirms drones substitute for absent local intelligence.\n")
cat("If triple interaction is NULL or NEGATIVE:\n")
cat("  \u2192 ELI effect does not depend on pre-existing info environment.\n")
cat("  \u2192 Mechanism may be strike reach rather than info substitution.\n")

write_audit_dt(out, "R14b_info_substitution")

# Also save the index itself for potential use in other scripts
idx_out <- file.path(DIRS$eli, "info_deficit_index.rds")
saveRDS(idx, idx_out)
log_msg("Saved info deficit index: ", idx_out)

log_msg("Step 14b complete.")
