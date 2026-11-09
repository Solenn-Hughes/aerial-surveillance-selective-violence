# ============================================================================
# R39_asymmetry_cannot_reciprocate.R
#
# DESIGN 4 — ASYMMETRY SIGNATURE (PROPER TEST)
#
# Theory prediction: Effects should be LARGER where rebels CANNOT reciprocate
# (i.e., rebels WITHOUT drones hurt MORE by ISR removal)
#
# SPECIFICATION:
#   y_igt = α_{ig} + γ_{gt} + λ_{it} + β₁(shock × eli_feas × no_drones)_igt
#                                     + β₂(shock × eli_feas × has_drones)_igt
#                                     + ε_igt
#
# Where no_drones = 1 - group_has_drones (inverted)
#
# HYPOTHESIS: β₁ > β₂ (groups without drones hurt more)
#
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})

source("R00_config.R")
log_msg("R39: Asymmetry signature — rebels cannot reciprocate")

OUT <- file.path("outputs", "paper1", "design4_asymmetry_disadvantage")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1) LOAD DYADIC PANEL ───────────────────────────────────────────────────
dyad_path <- file.path(DIRS$panels, "panel_dyadic_cell_group_month.rds")
if (!file.exists(dyad_path)) {
  stop("Dyadic panel not found. Run R34_build_dyadic_panel.R first.")
}

dyad <- readRDS(dyad_path); setDT(dyad)
log_msg("Loaded dyadic panel: ", nrow(dyad), " dyad-months")

# Add year for FE
dyad[, year := month_int %/% 12L]

# ── 2) LOAD GROUP DRONE ACQUISITION DATA ──────────────────────────────────
# From R23/R29 drone acquisition coding
drone_path <- file.path(DIRS$base, "group_drone_acquisitions.rds")

if (!file.exists(drone_path)) {
  cat("\n⚠ WARNING: group_drone_acquisitions.rds not found.\n")
  cat("Creating placeholder from actor names.\n")
  cat("Replace with actual SIPRI/CTC/UN Panel of Experts coding.\n\n")
  
  # Placeholder: flag known drone-operating groups
  drone_groups <- c(
    "Houthi", "Ansar Allah",  # Yemen
    "Hezbollah",              # Lebanon/Syria  
    "Hamas",                  # Gaza
    "Islamic State", "IS",    # Iraq/Syria
    "PKK", "YPG"              # Turkey/Syria
  )
  
  dyad[, group_has_drones := as.integer(
    grepl(paste(drone_groups, collapse = "|"), actor1, ignore.case = TRUE)
  )]
  
} else {
  drone_data <- readRDS(drone_path); setDT(drone_data)
  log_msg("Loaded drone acquisition data: ", nrow(drone_data), " group-months")
  
  # Merge into dyadic panel
  dyad <- merge(dyad, 
                drone_data[, .(actor1, month_int, group_has_drones)],
                by = c("actor1", "month_int"), 
                all.x = TRUE)
  
  dyad[is.na(group_has_drones), group_has_drones := 0L]
}

# Create INVERSE: groups WITHOUT drones
dyad[, group_no_drones := 1L - group_has_drones]

# ── 3) PRE-FLIGHT DIAGNOSTICS ──────────────────────────────────────────────
cat("\n== PRE-FLIGHT DIAGNOSTICS ==\n")
cat(sprintf("Total dyad-months: %d\n", nrow(dyad)))
cat(sprintf("Active-presence dyad-months: %d\n", sum(dyad$group_present == 1)))
cat(sprintf("Groups WITH drones: %d\n", uniqueN(dyad[group_has_drones == 1, actor1])))
cat(sprintf("Groups WITHOUT drones: %d\n", uniqueN(dyad[group_no_drones == 1, actor1])))

# Baseline violence by drone status
baseline_vac <- dyad[any_isr_removed == 0 & group_present == 1, .(
  mean_vac = mean(log1p(vac_events)),
  n = .N
), by = group_has_drones]

cat("\nBaseline VAC by drone status:\n")
print(baseline_vac)

# ── 4) SPEC 1 — STRATIFIED BY DRONE STATUS ────────────────────────────────
cat("\n== SPEC 1 — Stratified Analysis ==\n")

CTRL <- intersect(
  c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih"),
  names(dyad)
)
ctrl_str <- paste(CTRL, collapse = " + ")

# Sample 1: Groups WITHOUT drones (should show LARGE effect)
cat("\n--- Groups WITHOUT Drones ---\n")
dyad_no_drones <- dyad[group_no_drones == 1L]

fml_no <- as.formula(paste0(
  "log1p_vac ~ eli_feas:any_isr_removed:group_present + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))

m_no <- feols(fml_no, data = dyad_no_drones, notes = FALSE, mem.clean = TRUE)
cv_no <- vcov_conley(m_no, lat = "lat", lon = "lon", cutoff = 200)
ct_no <- as.data.table(coeftable(m_no, vcov = cv_no), keep.rownames = "term")

print(ct_no[grepl("eli_feas.*any_isr_removed.*group_present", term)])

# Sample 2: Groups WITH drones (should show SMALL or null effect)
cat("\n--- Groups WITH Drones ---\n")
dyad_has_drones <- dyad[group_has_drones == 1L]

if (nrow(dyad_has_drones[group_present == 1]) > 1000) {
  fml_has <- as.formula(paste0(
    "log1p_vac ~ eli_feas:any_isr_removed:group_present + ",
    ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
  ))
  
  m_has <- feols(fml_has, data = dyad_has_drones, notes = FALSE, mem.clean = TRUE)
  cv_has <- vcov_conley(m_has, lat = "lat", lon = "lon", cutoff = 200)
  ct_has <- as.data.table(coeftable(m_has, vcov = cv_has), keep.rownames = "term")
  
  print(ct_has[grepl("eli_feas.*any_isr_removed.*group_present", term)])
} else {
  cat("Insufficient observations for groups WITH drones (n =", 
      nrow(dyad_has_drones[group_present == 1]), ")\n")
  m_has <- NULL
  ct_has <- data.table()
}

# ── 5) SPEC 2 — TRIPLE INTERACTION (POOLED) ───────────────────────────────
cat("\n== SPEC 2 — Triple Interaction (Pooled Sample) ==\n")

# Test both groups in one model with triple interaction
fml_pool <- as.formula(paste0(
  "log1p_vac ~ eli_feas:any_isr_removed:group_present:group_no_drones + ",
  "eli_feas:any_isr_removed:group_present:group_has_drones + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))

m_pool <- feols(fml_pool, data = dyad, notes = FALSE, mem.clean = TRUE)
cv_pool <- vcov_conley(m_pool, lat = "lat", lon = "lon", cutoff = 200)
ct_pool <- as.data.table(coeftable(m_pool, vcov = cv_pool), keep.rownames = "term")

print(ct_pool[grepl("group_no_drones|group_has_drones", term)])

# ── 6) TEST DIFFERENCE ─────────────────────────────────────────────────────
cat("\n== HYPOTHESIS TEST: β_no_drones > β_has_drones ==\n")

b_no <- ct_no[grepl("eli_feas.*any_isr_removed.*group_present", term), Estimate]
se_no <- ct_no[grepl("eli_feas.*any_isr_removed.*group_present", term), `Std. Error`]
p_no <- ct_no[grepl("eli_feas.*any_isr_removed.*group_present", term), `Pr(>|t|)`]

if (nrow(ct_has) > 0) {
  b_has <- ct_has[grepl("eli_feas.*any_isr_removed.*group_present", term), Estimate]
  se_has <- ct_has[grepl("eli_feas.*any_isr_removed.*group_present", term), `Std. Error`]
  p_has <- ct_has[grepl("eli_feas.*any_isr_removed.*group_present", term), `Pr(>|t|)`]
  
  diff <- b_no - b_has
  se_diff <- sqrt(se_no^2 + se_has^2)
  t_diff <- diff / se_diff
  p_diff <- 2 * pt(-abs(t_diff), df = m_no$nobs - m_no$nparams)
  
  cat(sprintf("\nβ_no_drones   = %.5f (SE = %.5f, p = %.4f)\n", b_no, se_no, p_no))
  cat(sprintf("β_has_drones  = %.5f (SE = %.5f, p = %.4f)\n", b_has, se_has, p_has))
  cat(sprintf("Difference    = %.5f (SE = %.5f)\n", diff, se_diff))
  cat(sprintf("t-statistic   = %.3f\n", t_diff))
  cat(sprintf("p-value       = %.4f\n", p_diff))
  
  if (diff > 0 && p_diff < 0.10) {
    cat("\n✓ ASYMMETRY SIGNATURE CONFIRMED\n")
    cat("  Groups WITHOUT drones show significantly larger effects.\n")
    cat("  Supports displacement theory: asymmetric disadvantage matters.\n")
  } else {
    cat("\n○ No clear asymmetry pattern detected.\n")
  }
} else {
  cat(sprintf("\nβ_no_drones   = %.5f (SE = %.5f, p = %.4f)\n", b_no, se_no, p_no))
  cat("\nCannot test difference (insufficient drone-equipped group observations)\n")
  
  if (!is.na(b_no) && p_no < 0.10) {
    cat("\n○ Groups WITHOUT drones show significant effects.\n")
    cat("  But cannot compare to groups WITH drones.\n")
  }
  
  diff <- NA_real_
  p_diff <- NA_real_
}

# ── 7) RESULTS TABLE ───────────────────────────────────────────────────────
results <- data.table(
  spec = c("No drones (disadvantaged)", "Has drones (capable)", "Difference"),
  beta = c(b_no, 
           ifelse(nrow(ct_has) > 0, b_has, NA_real_),
           diff),
  se = c(se_no,
         ifelse(nrow(ct_has) > 0, se_has, NA_real_),
         ifelse(!is.na(diff), se_diff, NA_real_)),
  p_value = c(p_no,
              ifelse(nrow(ct_has) > 0, p_has, NA_real_),
              p_diff),
  n_groups = c(uniqueN(dyad_no_drones[group_present == 1, actor1]),
               ifelse(nrow(ct_has) > 0, uniqueN(dyad_has_drones[group_present == 1, actor1]), NA_integer_),
               NA_integer_)
)

print(results)
fwrite(results, file.path(OUT, "R39_asymmetry_disadvantage_results.csv"))

# ── 8) INTERPRETATION ──────────────────────────────────────────────────────
cat("\n== INTERPRETATION ==\n")

if (!is.na(diff) && diff > 0 && p_diff < 0.10) {
  cat("✓ ASYMMETRY SIGNATURE CONFIRMED\n")
  cat(sprintf("  Groups without drones: β = %.5f, p = %.4f\n", b_no, p_no))
  cat(sprintf("  Groups with drones:    β = %.5f, p = %.4f\n", b_has, p_has))
  cat(sprintf("  Difference:            β = %.5f, p = %.4f\n", diff, p_diff))
  cat("\n  Interpretation: Rebels lacking drone capability suffer MORE from ISR removal.\n")
  cat("  This supports displacement theory's asymmetry prediction:\n")
  cat("  when one side has ISR advantage, the disadvantaged side is hurt more.\n")
  
} else if (!is.na(b_no) && p_no < 0.10) {
  cat("○ PARTIAL SUPPORT\n")
  cat(sprintf("  Groups without drones: β = %.5f, p = %.4f (significant)\n", b_no, p_no))
  if (nrow(ct_has) > 0) {
    cat(sprintf("  Groups with drones:    β = %.5f, p = %.4f\n", b_has, p_has))
    cat("  But difference not statistically significant.\n")
  } else {
    cat("  Cannot compare to drone-equipped groups (insufficient data).\n")
  }
  
} else {
  cat("✗ ASYMMETRY SIGNATURE NOT CONFIRMED\n")
  cat("  No clear pattern of asymmetric disadvantage.\n")
  cat("  Effects are similar regardless of drone capability.\n")
}

log_msg("R39 complete. Outputs in: ", OUT)
