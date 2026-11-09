# ============================================================
# R09g_leave_one_out.R
# Leave-one-out theatre sensitivity analysis
#
# Re-runs the primary DiD model 8 times, dropping one treated
# theatre each time. Tests whether the result is driven by a
# single theatre or is genuinely cross-theatre consistent.
#
# Requires: morning startup code already run (dt_clean, fml, m0_boot)
# ============================================================

source("R00_config.R")
library(fixest)
library(data.table)

TREATED <- c("Afghanistan", "Iraq", "Libya", "Pakistan",
             "Somalia", "Syria", "Ukraine", "Yemen (North Yemen)")

# ── BASELINE ──────────────────────────────────────────────
baseline_coef <- coef(m0_boot)["post_eli:eli_feas"]
baseline_se   <- se(m0_boot)["post_eli:eli_feas"]
baseline_p    <- pvalue(m0_boot)["post_eli:eli_feas"]

cat(sprintf("Baseline (all 8 theatres): coef=%.5f  SE=%.5f  p=%.4f\n\n",
            baseline_coef, baseline_se, baseline_p))

# ── LEAVE-ONE-OUT LOOP ────────────────────────────────────
results <- list()

for (drop_t in TREATED) {
  cat(sprintf("Dropping: %-25s ... ", drop_t))

  dt_loo <- dt_clean[theatre_id != drop_t]

  # Recompute theatre_month_fe factor without the dropped theatre
  dt_loo[, theatre_month_fe := as.factor(paste0(theatre_id, "__", month_int))]

  m_loo <- tryCatch(
    feols(fml, data = dt_loo, cluster = ~theatre_id, warn = FALSE, notes = FALSE),
    error = function(e) { cat("ERROR:", e$message, "\n"); NULL }
  )

  if (!is.null(m_loo)) {
    c_val <- coef(m_loo)["post_eli:eli_feas"]
    s_val <- se(m_loo)["post_eli:eli_feas"]
    p_val <- pvalue(m_loo)["post_eli:eli_feas"]
    cat(sprintf("coef=%.5f  SE=%.5f  p=%.4f\n", c_val, s_val, p_val))
    results[[drop_t]] <- data.table(
      dropped   = drop_t,
      coef      = c_val,
      se        = s_val,
      p_value   = p_val,
      n_theatre = uniqueN(dt_loo$theatre_id)
    )
  }
}

# ── COMPILE RESULTS ───────────────────────────────────────
out <- rbindlist(results)

# Add baseline row
out <- rbindlist(list(
  data.table(dropped="[Baseline: all 8]", coef=baseline_coef,
             se=baseline_se, p_value=baseline_p, n_theatre=8L),
  out
))

cat("\n====== LEAVE-ONE-OUT RESULTS ======\n")
print(out[, .(dropped, coef=round(coef,5), se=round(se,5),
              p_value=round(p_value,4), n_theatre)])
cat("===================================\n\n")

# Flag any sign changes or large shifts
sign_changes <- out[dropped != "[Baseline: all 8]" & sign(coef) != sign(baseline_coef)]
large_shifts  <- out[dropped != "[Baseline: all 8]" & abs(coef - baseline_coef) > 0.02]

if (nrow(sign_changes) > 0) {
  cat("⚠ SIGN CHANGES when dropping:", paste(sign_changes$dropped, collapse=", "), "\n")
} else {
  cat("✓ No sign changes — coefficient direction stable across all leave-one-out runs.\n")
}
if (nrow(large_shifts) > 0) {
  cat("⚠ LARGE SHIFTS (>0.02) when dropping:", paste(large_shifts$dropped, collapse=", "), "\n")
} else {
  cat("✓ No large shifts — no single theatre drives the result.\n")
}

# ── SAVE ──────────────────────────────────────────────────
out_path <- file.path("outputs/paper1",
  paste0("leave_one_out_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
fwrite(out, out_path)
log_msg("Saved: ", out_path)
write_audit_dt(out, "R09g_leave_one_out", subdir = DIRS$audit)
log_msg("Leave-one-out complete.")
