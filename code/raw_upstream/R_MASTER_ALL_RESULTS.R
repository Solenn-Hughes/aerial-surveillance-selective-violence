# ============================================================================
# MASTER RESULTS — CORRECT VERSION
# Based on HANDOFF_COMPLETE_20260415.md
#
# MAIN RESULT: Stacked quasi-experiment β = -0.037 (ISR removal reduces violence)
# MECHANISMS: post_eli:eli_feas tests (ISR presence increases violence)
# OUTCOME: events_type3 (UCDP-GED one-sided violence)
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})

source("R00_config.R")

cat("\n════════════════════════════════════════════════════════════\n")
cat("MASTER RESULTS — COMPREHENSIVE VERIFICATION\n")
cat("════════════════════════════════════════════════════════════\n\n")

OUT <- file.path("outputs", "paper1", "MASTER_RESULTS")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════
# LOAD PANEL
# ══════════════════════════════════════════════════════════════════════

cat("Loading panel...\n")
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(panel_path); setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]

if (!all(c("lat", "lon") %in% names(dt))) {
  if (all(c("ycoord", "xcoord") %in% names(dt))) {
    dt[, `:=`(lat = ycoord, lon = xcoord)]
  }
}

# Build outcome
if (!"log1p_type3" %in% names(dt)) {
  dt[, log1p_type3 := log1p(events_type3)]
}

cat("Panel: ", nrow(dt), " rows\n")
cat("Outcome: events_type3 (", sum(dt$events_type3), " total events)\n\n")

# ══════════════════════════════════════════════════════════════════════
# CONTROLS
# ══════════════════════════════════════════════════════════════════════

CTRL <- intersect(
  c("nlights_calib_mean", "prec_gpcp", "temp", "water_gc", "forest_gc"),
  names(dt)
)
ctrl_str <- paste(CTRL, collapse = " + ")

# ══════════════════════════════════════════════════════════════════════
# PART 1: MAIN CAUSAL EFFECT (STACKED QUASI-EXPERIMENT)
# ══════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════════\n")
cat("MAIN CAUSAL EFFECT: STACKED QUASI-EXPERIMENT\n")
cat("══════════════════════════════════════════════════════════\n\n")

# Build shock indicators
if (!"any_isr_removed" %in% names(dt)) {
  dt[, any_isr_removed := 0L]
  dt[theatre_id == "Afghanistan" & month >= as.IDate("2021-07-01"), any_isr_removed := 1L]
  dt[theatre_id %in% c("Syria", "Iraq") & month >= as.IDate("2019-10-01"), any_isr_removed := 1L]
}

# Stacked QE specification
fml_main <- as.formula(paste0(
  "log1p_type3 ~ any_isr_removed:eli_feas + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))

m_main <- feols(fml_main, data = dt,
                vcov = conley(cutoff = 200, distance = "spherical"))

main_beta <- coef(m_main)["any_isr_removed:eli_feas"]
main_se <- se(m_main)["any_isr_removed:eli_feas"]
main_p <- summary(m_main)$coeftable["any_isr_removed:eli_feas", "Pr(>|t|)"]

main_result <- data.table(
  section = "MAIN CAUSAL EFFECT",
  test = "Stacked quasi-experiment (ISR removal)",
  outcome = "events_type3",
  specification = "any_isr_removed × eli_feas",
  beta = main_beta,
  se = main_se,
  p_value = main_p,
  prediction = "β < 0 (violence decreases)",
  interpretation = ifelse(
    main_beta < 0 & main_p < 0.05,
    "✓ CAUSAL EFFECT CONFIRMED",
    ifelse(main_beta < 0, "○ Correct sign, marginal", "✗ WRONG SIGN")
  )
)

cat("Result:\n")
print(main_result[, .(test, beta, se, p_value, interpretation)])
cat("\n")

if (main_beta < 0 & main_p < 0.05) {
  cat("✓ ISR removal REDUCES violence in high-feasibility cells\n")
  cat("  This is your PRIMARY CAUSAL CLAIM\n\n")
} else {
  cat("⚠️ WARNING: Main effect not replicating expected result!\n")
  cat("   Expected: β ≈ -0.037, p < 0.01\n")
  cat("   Got: β = ", sprintf("%.4f", main_beta), ", p = ", sprintf("%.4f", main_p), "\n\n")
}

# ══════════════════════════════════════════════════════════════════════
# PART 2: LOAD ACLED FOR ACTOR DECOMPOSITION
# ══════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════════\n")
cat("ACTOR DECOMPOSITION (CORRECTED)\n")
cat("══════════════════════════════════════════════════════════\n\n")

cat("Loading ACLED...\n")
ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT, "data", "raw", "acled", "acled_full.csv")
)
acled_path <- Filter(file.exists, ACLED_PATHS)[1]
ac <- fread(acled_path, showProgress = FALSE, encoding = "UTF-8")
setnames(ac, tolower(names(ac)))
ac[, month := as.IDate(format(as.Date(event_date), "%Y-%m-01"))]
ac[, fatalities := as.integer(fatalities)]
ac[is.na(fatalities), fatalities := 0L]
ac[, pg_row := floor((latitude + 90) / 0.5) + 1L]
ac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
ac[, pg_id := (pg_row - 1L) * 720L + pg_col]

# Actor classification
ac[, state_inv := (inter1 == 1L)]
ac[, rebel_inv := (inter1 == 2L)]
ac[, is_vac := grepl("Violence against civilians", event_type, ignore.case = TRUE)]
ac[, is_remote := grepl("Explosions|Remote violence", event_type, ignore.case = TRUE)]
ac[, civ_target := (inter1 == 7L | inter2 == 7L)]

# CORRECTED state VAC (includes remote violence with civilian targets)
state_vac_events <- ac[state_inv == TRUE & (is_vac == TRUE | 
                                              (is_remote == TRUE & civ_target == TRUE))]
rebel_vac_events <- ac[rebel_inv == TRUE & is_vac == TRUE]

agg_state <- state_vac_events[, .(state_vac = .N), by = .(pg_id, month)]
agg_rebel <- rebel_vac_events[, .(rebel_vac = .N), by = .(pg_id, month)]

dt <- merge(dt, agg_state, by = c("pg_id", "month"), all.x = TRUE)
dt <- merge(dt, agg_rebel, by = c("pg_id", "month"), all.x = TRUE)
dt[is.na(state_vac), state_vac := 0L]
dt[is.na(rebel_vac), rebel_vac := 0L]
dt[, log1p_sv := log1p(state_vac)]
dt[, log1p_rv := log1p(rebel_vac)]

cat("State VAC: ", sum(dt$state_vac), " events\n")
cat("Rebel VAC: ", sum(dt$rebel_vac), " events\n")
cat("Rebel/State ratio: ", round(sum(dt$rebel_vac) / sum(dt$state_vac), 2), "\n\n")

# Actor decomposition tests
actor_results <- list()

fml_state <- as.formula(paste0("log1p_sv ~ post_eli:eli_feas + ", ctrl_str,
                               " | pg_id + theatre_id^month_int"))
fml_rebel <- as.formula(paste0("log1p_rv ~ post_eli:eli_feas + ", ctrl_str,
                               " | pg_id + theatre_id^month_int"))

m_state <- feols(fml_state, data = dt, 
                 vcov = conley(cutoff = 200, distance = "spherical"))
m_rebel <- feols(fml_rebel, data = dt,
                 vcov = conley(cutoff = 200, distance = "spherical"))

actor_results[[1]] <- data.table(
  section = "ACTOR DECOMPOSITION",
  test = "State VAC (substitution test)",
  outcome = "State civilian harm",
  specification = "post_eli × eli_feas",
  beta = coef(m_state)["post_eli:eli_feas"],
  se = se(m_state)["post_eli:eli_feas"],
  p_value = summary(m_state)$coeftable["post_eli:eli_feas", "Pr(>|t|)"],
  prediction = "β ≈ 0 (substitution rejected)",
  interpretation = ifelse(
    summary(m_state)$coeftable["post_eli:eli_feas", "Pr(>|t|)"] > 0.10,
    "✓ NULL (as predicted)", "✗ Unexpected significance"
  )
)

actor_results[[2]] <- data.table(
  section = "ACTOR DECOMPOSITION",
  test = "Rebel VAC (displacement test)",
  outcome = "Rebel civilian harm",
  specification = "post_eli × eli_feas",
  beta = coef(m_rebel)["post_eli:eli_feas"],
  se = se(m_rebel)["post_eli:eli_feas"],
  p_value = summary(m_rebel)$coeftable["post_eli:eli_feas", "Pr(>|t|)"],
  prediction = "β > 0 (displacement confirmed)",
  interpretation = ifelse(
    coef(m_rebel)["post_eli:eli_feas"] > 0 & 
      summary(m_rebel)$coeftable["post_eli:eli_feas", "Pr(>|t|)"] < 0.05,
    "✓ DISPLACEMENT CONFIRMED", "○ Weak/null"
  )
)

actor_final <- rbindlist(actor_results)
cat("Results:\n")
print(actor_final[, .(test, beta, se, p_value, interpretation)])
cat("\n")

# ══════════════════════════════════════════════════════════════════════
# PART 3: CONSOLIDATE AND SAVE
# ══════════════════════════════════════════════════════════════════════

all_results <- rbindlist(list(main_result, actor_final), fill = TRUE)
fwrite(all_results, file.path(OUT, "MASTER_RESULTS_CORRECT.csv"))

cat("══════════════════════════════════════════════════════════\n")
cat("SUMMARY\n")
cat("══════════════════════════════════════════════════════════\n\n")

cat("MAIN CAUSAL EFFECT (ISR removal):\n")
cat(sprintf("  β = %.4f (SE = %.4f, p = %.4f) %s\n\n",
            main_result$beta, main_result$se, main_result$p_value,
            main_result$interpretation))

cat("ACTOR DECOMPOSITION:\n")
cat(sprintf("  State VAC:  β = %.4f (p = %.4f) %s\n",
            actor_final[1]$beta, actor_final[1]$p_value, actor_final[1]$interpretation))
cat(sprintf("  Rebel VAC:  β = %.4f (p = %.4f) %s\n\n",
            actor_final[2]$beta, actor_final[2]$p_value, actor_final[2]$interpretation))

cat("INTERPRETATION:\n")
if (main_result$beta < 0 & main_result$p_value < 0.05 &
    actor_final[2]$beta > 0 & actor_final[2]$p_value < 0.05) {
  cat("✓ PAPER STORY CONFIRMED:\n")
  cat("  - ISR removal REDUCES violence (main effect)\n")
  cat("  - Rebels DRIVE the result (actor decomp)\n")
  cat("  - States show NULL effect (substitution rejected)\n")
  cat("  → DISPLACEMENT MECHANISM\n")
} else {
  cat("⚠️ RESULTS NOT FULLY REPLICATING:\n")
  cat("  Check panel file, shock dates, and ACLED version\n")
}

cat("\n══════════════════════════════════════════════════════════\n")
cat("SAVED TO: ", file.path(OUT, "MASTER_RESULTS_CORRECT.csv"), "\n")
cat("══════════════════════════════════════════════════════════\n\n")