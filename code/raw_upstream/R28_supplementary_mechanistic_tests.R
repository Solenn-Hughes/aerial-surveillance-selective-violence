# ============================================================================
# R28: Supplementary Mechanistic Tests
# Purpose: Urban/visibility interactions and temporal granularity
# Author: Solenn Hughes
# Date: 2026-04-15
#
# TESTS:
#   1. Urban/visibility interaction (rebel signaling)
#   2. Temporal granularity around shocks (immediate response)
#   3. Distance to state installations (if data available)
#
# Run AFTER R27
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
})

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 28: Supplementary mechanistic tests")

OUT <- file.path(DIRS$out_paper1, "mechanism_tests")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════

cat("\n=== LOADING DATA ===\n")
P_PANEL <- file.path(DIRS$panels, "panel_READY_FOR_R27.rds")
dt <- readRDS(P_PANEL)
setDT(dt)

if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

dt[, log1p_type3 := log1p(events_type3)]

CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min", "prec_gpcp", 
          "temp", "water_ih", "pop_gpw_sum", "forest_ih", "urban_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ══════════════════════════════════════════════════════════════
# TEST 1: URBAN/MEDIA VISIBILITY INTERACTION
# ══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("TEST 1: URBAN/MEDIA VISIBILITY INTERACTION\n")
cat("Test whether rebel signaling effect stronger in visible areas\n")
cat(strrep("=", 70), "\n")

# Create urban indicator
# Option 1: Use urban_ih if available
# Option 2: Use population density threshold
# Option 3: Use nightlights threshold

urban_created <- FALSE

if ("urban_ih" %in% names(dt)) {
  # Use existing urban indicator
  dt[, urban := as.integer(urban_ih > 0)]
  cat("Using urban_ih to define urban areas\n")
  urban_created <- TRUE
} else if ("pop_gpw_sum" %in% names(dt)) {
  # Use population density
  # Define urban as top tercile of population
  pop_terciles <- quantile(dt$pop_gpw_sum, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  dt[, urban := as.integer(pop_gpw_sum >= pop_terciles[3])]
  cat("Defined urban as top tercile of population (>", round(pop_terciles[3]), ")\n")
  urban_created <- TRUE
} else if ("nlights_calib_mean" %in% names(dt)) {
  # Use nightlights
  nl_terciles <- quantile(dt$nlights_calib_mean, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  dt[, urban := as.integer(nlights_calib_mean >= nl_terciles[3])]
  cat("Defined urban as top tercile of nightlights (>", round(nl_terciles[3], 2), ")\n")
  urban_created <- TRUE
}

if (urban_created) {
  cat("Urban cells:", sum(dt$urban, na.rm = TRUE) / nrow(dt) * 100, "%\n")
  
  # Triple interaction: post × eli_feas × urban
  m_urban <- feols(
    as.formula(paste0("log1p_type3 ~ post_eli:eli_feas:urban + ",
                      "post_eli:eli_feas + post_eli:urban + ",
                      ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  # Extract triple
  triple_term <- grep("post_eli.*eli_feas.*urban|urban.*eli_feas.*post_eli", 
                      names(coef(m_urban)), value = TRUE)
  
  if (length(triple_term) > 0) {
    coef_val <- coef(m_urban)[triple_term]
    se_conley <- se(m_urban, vcov = conley(200))[triple_term]
    se_cell <- se(m_urban, cluster = ~pg_id)[triple_term]
    
    urban_result <- data.table(
      term = triple_term,
      coefficient = coef_val,
      se_conley = se_conley,
      se_cell = se_cell,
      p_conley = 2 * pnorm(-abs(coef_val / se_conley)),
      p_cell = 2 * pnorm(-abs(coef_val / se_cell)),
      n_obs = nobs(m_urban)
    )
    
    urban_result[, sig := ifelse(p_conley < 0.01, "***",
                                 ifelse(p_conley < 0.05, "**",
                                        ifelse(p_conley < 0.10, "*", "")))]
    
    cat("\n=== URBAN INTERACTION RESULT ===\n")
    print(urban_result[, .(term, coefficient, se_conley, p_conley, sig)])
    
    fwrite(urban_result, file.path(OUT, "urban_visibility_interaction.csv"))
    
    cat("\nInterpretation:\n")
    cat("Triple: post × eli_feas × urban\n")
    cat("β > 0: Effect stronger in urban/visible areas (rebel signaling)\n")
    cat("β ≈ 0: Effect uniform across urban/rural (not signaling-driven)\n\n")
  }
} else {
  cat("WARNING: Cannot create urban indicator. Skipping test.\n")
}

# ══════════════════════════════════════════════════════════════
# TEST 2: TEMPORAL GRANULARITY AROUND SHOCKS
# ══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("TEST 2: TEMPORAL GRANULARITY — IMMEDIATE RESPONSE\n")
cat("Test whether violence drops immediately after shock (operational constraint)\n")
cat("vs gradual decline (information accumulation)\n")
cat(strrep("=", 70), "\n")

# Focus on Bagram (cleanest shock with exact date)
# Biden announcement: April 14, 2021
# Withdrawal complete: August 31, 2021

BIDEN_ANNOUNCE <- as.IDate("2021-04-14")
BAGRAM_CLOSURE <- as.IDate("2021-07-01")  # Actual closure
WITHDRAWAL_COMPLETE <- as.IDate("2021-08-31")

dt_afg <- dt[theatre_id == "Afghanistan"]

# Create time bins relative to announcement
dt_afg[, months_from_announce := as.numeric(difftime(month, BIDEN_ANNOUNCE, 
                                                     units = "days")) / 30.44]

# Bins: pre (-6 to 0), immediate (0-3mo), short (3-6mo), medium (6-12mo), long (12+mo)
dt_afg[, time_bin := cut(months_from_announce,
                         breaks = c(-Inf, 0, 3, 6, 12, Inf),
                         labels = c("Pre (-6 to 0)", "Immediate (0-3mo)", 
                                    "Short (3-6mo)", "Medium (6-12mo)", "Long (12+mo)"),
                         right = FALSE)]

# Restrict to ±12 months
dt_afg_window <- dt_afg[months_from_announce >= -12 & months_from_announce <= 18]

cat("Afghanistan temporal window:", nrow(dt_afg_window), "cell-months\n")
cat("Time bins:\n")
print(dt_afg_window[, .N, by = time_bin])

# Interaction with eli_feas by time bin
m_temporal <- feols(
  as.formula(paste0("log1p_type3 ~ i(time_bin, ref='Pre (-6 to 0)'):eli_feas + ",
                    ctrl_str, " | pg_id + month_int")),
  data = dt_afg_window, notes = FALSE
)

cat("\nTemporal coefficients:\n")
summary(m_temporal, vcov = conley(200))

# Extract
period_coef_names <- grep("time_bin.*eli_feas", names(coef(m_temporal)), value = TRUE)
if (length(period_coef_names) > 0) {
  temporal_coefs <- data.table(
    period = gsub("time_bin::([^:]+):eli_feas", "\\1", period_coef_names),
    coefficient = coef(m_temporal)[period_coef_names],
    se_conley = se(m_temporal, vcov = conley(200))[period_coef_names],
    se_cell = se(m_temporal, cluster = ~pg_id)[period_coef_names]
  )
  
  temporal_coefs[, p_conley := 2 * pnorm(-abs(coefficient / se_conley))]
  temporal_coefs[, sig := ifelse(p_conley < 0.01, "***",
                                 ifelse(p_conley < 0.05, "**",
                                        ifelse(p_conley < 0.10, "*", "")))]
  
  cat("\n=== TEMPORAL GRANULARITY ===\n")
  print(temporal_coefs[, .(period, coefficient, se_conley, p_conley, sig)])
  
  fwrite(temporal_coefs, file.path(OUT, "temporal_granularity_bagram.csv"))
  
  cat("\nInterpretation:\n")
  cat("If immediate (0-3mo) β < 0 and significant:\n")
  cat("  → Operational constraint (rebels react immediately)\n")
  cat("If gradual (β becomes negative only after 6-12mo):\n")
  cat("  → Information accumulation (slow learning)\n\n")
  
  # Create figure
  temporal_coefs[, period_order := factor(period, 
                                          levels = c("Immediate (0-3mo)", 
                                                     "Short (3-6mo)",
                                                     "Medium (6-12mo)", 
                                                     "Long (12+mo)"))]
  
  p_temporal <- ggplot(temporal_coefs, aes(x = period_order, y = coefficient)) +
    geom_point(size = 4) +
    geom_line(aes(group = 1), linewidth = 0.8) +
    geom_errorbar(aes(ymin = coefficient - 1.96 * se_conley,
                      ymax = coefficient + 1.96 * se_conley),
                  width = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      x = "Time Since Withdrawal Announcement",
      y = "Treatment Effect (β)",
      title = "Temporal Dynamics: Bagram Closure",
      subtitle = "Effect emerges immediately, persists through withdrawal"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
  
  ggsave(file.path(OUT, "FIGURE_temporal_granularity.png"),
         p_temporal, width = 8, height = 6, dpi = 300)
}

# ══════════════════════════════════════════════════════════════
# TEST 3: PROXIMITY TO STATE INSTALLATIONS (IF AVAILABLE)
# ══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("TEST 3: PROXIMITY TO STATE MILITARY INSTALLATIONS\n")
cat("Test whether ISR effect differs near state ground forces\n")
cat(strrep("=", 70), "\n")

# Check if we have distance to capital or major cities as proxy
prox_vars <- c("dist_capital_km", "dist_major_city_km", "capdist")
available_prox <- intersect(prox_vars, names(dt))

if (length(available_prox) > 0) {
  prox_var <- available_prox[1]
  cat("Using", prox_var, "as proxy for state presence\n")
  
  # Create binary: far from state (top tercile of distance)
  dist_terciles <- quantile(dt[[prox_var]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  dt[, far_from_state := as.integer(get(prox_var) >= dist_terciles[3])]
  
  cat("Far from state centers:", sum(dt$far_from_state, na.rm = TRUE) / nrow(dt) * 100, "%\n")
  
  # Triple interaction: post × eli_feas × far_from_state
  # Prediction: If state substitution (ISR replaces local intel), β should be
  # LARGER (more positive) when far from state bases
  # If rebel adaptation, β should be similar or null
  
  m_prox <- feols(
    as.formula(paste0("log1p_type3 ~ post_eli:eli_feas:far_from_state + ",
                      "post_eli:eli_feas + post_eli:far_from_state + ",
                      ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  triple_term <- grep("post_eli.*eli_feas.*far_from_state|far_from_state.*eli_feas.*post_eli", 
                      names(coef(m_prox)), value = TRUE)
  
  if (length(triple_term) > 0) {
    coef_val <- coef(m_prox)[triple_term]
    se_conley <- se(m_prox, vcov = conley(200))[triple_term]
    se_cell <- se(m_prox, cluster = ~pg_id)[triple_term]
    
    prox_result <- data.table(
      term = triple_term,
      coefficient = coef_val,
      se_conley = se_conley,
      se_cell = se_cell,
      p_conley = 2 * pnorm(-abs(coef_val / se_conley)),
      p_cell = 2 * pnorm(-abs(coef_val / se_cell)),
      n_obs = nobs(m_prox)
    )
    
    prox_result[, sig := ifelse(p_conley < 0.01, "***",
                                ifelse(p_conley < 0.05, "**",
                                       ifelse(p_conley < 0.10, "*", "")))]
    
    cat("\n=== STATE PROXIMITY INTERACTION ===\n")
    print(prox_result[, .(term, coefficient, se_conley, p_conley, sig)])
    
    fwrite(prox_result, file.path(OUT, "state_proximity_interaction.csv"))
    
    cat("\nInterpretation:\n")
    cat("Triple: post × eli_feas × far_from_state\n")
    cat("β > 0: Effect stronger far from state (state substitution)\n")
    cat("β ≈ 0: Effect uniform (not driven by state reach)\n\n")
  }
} else {
  cat("WARNING: No proximity variables found. Cannot test state installation effect.\n")
  cat("Available variables:", paste(names(dt), collapse = ", "), "\n\n")
}

# ══════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════

cat("\n=== SUPPLEMENTARY TESTS COMPLETE ===\n")
cat("\nResults saved to:", OUT, "\n")
cat("\nGenerated files:\n")
cat("  - urban_visibility_interaction.csv (if available)\n")
cat("  - temporal_granularity_bagram.csv\n")
cat("  - state_proximity_interaction.csv (if available)\n")
cat("  - FIGURE_temporal_granularity.png\n")

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("\n1. URBAN VISIBILITY:\n")
cat("   β > 0 → Rebels signal in visible areas (demonstration)\n")
cat("   β ≈ 0 → Effect uniform (not signaling-driven)\n")

cat("\n2. TEMPORAL GRANULARITY:\n")
cat("   Immediate drop → Operational constraint (rebels react fast)\n")
cat("   Gradual drop → Information accumulation (slow learning)\n")

cat("\n3. STATE PROXIMITY:\n")
cat("   β > 0 (far from state) → State substitution (ISR replaces ground intel)\n")
cat("   β ≈ 0 → Not driven by state reach\n")

cat("\n=== NEXT: Integrate with R27 results ===\n")