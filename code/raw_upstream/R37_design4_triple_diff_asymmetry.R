# ============================================================================
# R37_design4_triple_diff_asymmetry.R
#
# DESIGN 4 — ASYMMETRY SIGNATURE, TRIPLE-DIFFERENCE
#
# Purpose: Disentangle two channels that R32 conflated:
#   (A) CAPACITY — rebel drone acquisition enables more sophisticated
#       targeting → more VAC (positive β on group_has_drones regardless
#       of ISR exposure).
#   (B) ASYMMETRY — rebel drone acquisition reduces the state's
#       informational advantage → less coercion pressure in surveilled
#       cells (negative β on eli_feas × group_has_drones, NET of the
#       general capacity effect).
#
#   The asymmetry prediction is theoretically cleaner: it says that
#   the differential effect of surveillance (eli_feas) on rebel VAC
#   should attenuate when the rebel has their own surveillance.
#
# SPECIFICATION (DDD):
#   y_igt = α_{ig} + γ_{gt} + λ_{it} +
#            β₁(eli_feas × group_has_drones)_igt +     ← the asymmetry test
#            β₂(group_has_drones)_gt          [absorbed by γ_{gt}] +
#            β₃(eli_feas)_i                   [absorbed by α_{ig}] +
#            ε_igt
#
#   β₁ is the triple interaction. Under:
#     — Pure capacity:   β₁ = 0   (drones don't affect the surveillance
#                                    differential; only shift levels)
#     — Pure asymmetry:  β₁ < 0  (drones attenuate the surveillance-
#                                    driven coercion pressure)
#     — Both:            β₁ < 0  but smaller in magnitude than expected
#                                  under pure asymmetry
#
# IDENTIFYING ASSUMPTION:
#   Absent drone acquisition, the difference in VAC between surveilled
#   and non-surveilled cells would have evolved identically for groups
#   that did vs didn't acquire drones. Cell × group FE absorb baseline
#   group presence. Group × month FE absorb all group-specific shocks
#   (organizational change, battlefield developments, leadership).
#
# INPUT:  panel_dyadic_cell_group_month.rds (from R34)
# OUTPUT: Design 4 regression tables
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})
source("R00_config.R")
log_msg("R37: Design 4 — DDD for asymmetry signature")

OUT <- file.path("outputs", "paper1", "design4_asymmetry_ddd")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1) LOAD DYADIC PANEL ──────────────────────────────────────────────
dyad_path <- file.path(DIRS$panels, "panel_dyadic_cell_group_month.rds")
if (!file.exists(dyad_path)) {
  stop("Dyadic panel not found. Run R34_build_dyadic_panel.R first.")
}
dyad <- readRDS(dyad_path); setDT(dyad)
log_msg("Loaded: ", nrow(dyad), " dyad-months")

# Add year
dyad[, year := month_int %/% 12L]

# Restrict to groups that had activity in the panel
# (group_present == 1 at some point in the cell)
# This reduces noise from dyads where the group never operated in that cell
dyad_group_cell_active <- dyad[, .(was_ever_present = max(group_present)),
                                by = .(pg_id, actor1)]
active_combos <- dyad_group_cell_active[was_ever_present == 1L,
                                         .(pg_id, actor1)]
dyad_active <- merge(dyad, active_combos, by = c("pg_id", "actor1"))
log_msg("Dyad-months with active cell × group: ", nrow(dyad_active))

# ── 2) DIAGNOSTICS ────────────────────────────────────────────────────
cat("\n== DIAGNOSTICS ==\n")
cat(sprintf("Groups that acquired drones: %d / %d\n",
            uniqueN(dyad_active[group_has_drones == 1L, actor1]),
            uniqueN(dyad_active$actor1)))
cat("Acquisition summary by group:\n")
acq_summary <- dyad_active[, .(
  has_drones_share = mean(group_has_drones),
  n_dyad_months    = .N,
  n_vac_events     = sum(vac_events)
), by = actor1]
print(acq_summary[order(-has_drones_share)][1:15])

CTRL <- intersect(
  c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih"),
  names(dyad_active)
)
ctrl_str <- paste(CTRL, collapse = " + ")

# ── 3) SPEC 1 — NAIVE CAPACITY-ONLY TEST ──────────────────────────────
# Replicates R32 conceptually: just eli_feas × group_has_drones
# without cell × group × time FE. This is what R32 did.
cat("\n== SPEC 1 — Naive capacity test (R32 equivalent) ==\n")
fml_1 <- as.formula(paste0(
  "log1p_vac ~ eli_feas:group_has_drones + ",
  ctrl_str, " | pg_id + actor1^month_int"
))
m1 <- feols(fml_1, data = dyad_active, notes = FALSE, mem.clean = TRUE)
cv1 <- vcov_conley(m1, lat = "lat", lon = "lon", cutoff = 200)
ct1 <- as.data.table(coeftable(m1, vcov = cv1), keep.rownames = "term")
print(ct1[grepl("eli_feas|has_drones", term)])

# ── 4) SPEC 2 — FULL DDD ──────────────────────────────────────────────
# Cell × group FE + group × month FE + cell × year FE
# The β₁ coefficient on the triple interaction is the asymmetry test
cat("\n== SPEC 2 — Full DDD (asymmetry test) ==\n")
fml_2 <- as.formula(paste0(
  "log1p_vac ~ eli_feas:group_has_drones + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))
m2 <- feols(fml_2, data = dyad_active, notes = FALSE, mem.clean = TRUE)
cv2 <- vcov_conley(m2, lat = "lat", lon = "lon", cutoff = 200)
ct2 <- as.data.table(coeftable(m2, vcov = cv2), keep.rownames = "term")
print(ct2[grepl("eli_feas|has_drones", term)])

# ── 5) SPEC 3 — DDD WITH PANEL TIMING AS WELL ─────────────────────────
# Isolate asymmetry net of capacity + panel ISR onset
cat("\n== SPEC 3 — DDD with panel ISR interaction ==\n")
fml_3 <- as.formula(paste0(
  "log1p_vac ~ eli_feas:post_eli + ",
  "eli_feas:group_has_drones + ",
  "eli_feas:post_eli:group_has_drones + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))
m3 <- feols(fml_3, data = dyad_active, notes = FALSE, mem.clean = TRUE)
cv3 <- vcov_conley(m3, lat = "lat", lon = "lon", cutoff = 200)
ct3 <- as.data.table(coeftable(m3, vcov = cv3), keep.rownames = "term")
print(ct3[grepl("eli_feas|has_drones|post_eli", term)])

# ── 6) SPEC 4 — SELECTIVE VAC OUTCOME ─────────────────────────────────
cat("\n== SPEC 4 — Selective VAC (asymmetry + targeting) ==\n")
fml_4 <- as.formula(paste0(
  "log1p_vac_sel ~ eli_feas:group_has_drones + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))
m4 <- feols(fml_4, data = dyad_active, notes = FALSE, mem.clean = TRUE)
cv4 <- vcov_conley(m4, lat = "lat", lon = "lon", cutoff = 200)
ct4 <- as.data.table(coeftable(m4, vcov = cv4), keep.rownames = "term")
print(ct4[grepl("eli_feas|has_drones", term)])

# ── 7) SPEC 5 — INDISCRIMINATE VAC OUTCOME ────────────────────────────
cat("\n== SPEC 5 — Indiscriminate VAC ==\n")
fml_5 <- as.formula(paste0(
  "log1p_vac_ind ~ eli_feas:group_has_drones + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))
m5 <- feols(fml_5, data = dyad_active, notes = FALSE, mem.clean = TRUE)
cv5 <- vcov_conley(m5, lat = "lat", lon = "lon", cutoff = 200)
ct5 <- as.data.table(coeftable(m5, vcov = cv5), keep.rownames = "term")
print(ct5[grepl("eli_feas|has_drones", term)])

# ── 8) EVENT-STUDY AROUND GROUP'S DRONE-ACQUISITION DATE ──────────────
# For each group, define event time as months since acquisition.
# Run event study: does eli_feas × post_acq change smoothly or jump?
cat("\n== SPEC 6 — Event study around acquisition ==\n")

# Compute group-specific acquisition month
group_acq <- dyad_active[group_has_drones == 1L,
                          .(acq_month_int = min(month_int)),
                          by = actor1]
dyad_es <- merge(dyad_active, group_acq, by = "actor1", all.x = TRUE)
dyad_es[, event_time := month_int - acq_month_int]
dyad_es <- dyad_es[!is.na(event_time) & event_time >= -24L & event_time <= 24L]
dyad_es[, et_cat := factor(event_time,
  levels = setdiff(-24L:24L, -1L))]

fml_es <- as.formula(paste0(
  "log1p_vac ~ i(et_cat, eli_feas, ref = NULL) + ",
  ctrl_str, " | pg_id^actor1 + actor1^month_int + pg_id^year"
))
m_es <- tryCatch({
  feols(fml_es, data = dyad_es, notes = FALSE, mem.clean = TRUE)
}, error = function(e) {
  message("Event-study failed: ", e$message); NULL
})

if (!is.null(m_es)) {
  es_ct <- as.data.table(coeftable(m_es),
                          keep.rownames = "term")
  es_ct <- es_ct[grepl("et_cat::", term)]
  es_ct[, et := as.integer(sub(".*et_cat::(-?\\d+).*", "\\1", term))]
  setorder(es_ct, et)
  cat("Event-study coefficients around acquisition:\n")
  print(head(es_ct[, .(et, Estimate, `Std. Error`, `Pr(>|t|)`)], 10))
  fwrite(es_ct, file.path(OUT, "R37_event_study_acquisition.csv"))
}

# ── 9) RESULTS SUMMARY ────────────────────────────────────────────────
results <- rbindlist(list(
  data.table(
    spec = "Spec 1 — naive (R32 eq.)",
    outcome = "log1p_vac",
    term    = "eli_feas:group_has_drones",
    beta    = round(ct1[grepl("eli_feas:group_has_drones", term), Estimate][1], 5),
    p_value = round(ct1[grepl("eli_feas:group_has_drones", term), `Pr(>|t|)`][1], 4)
  ),
  data.table(
    spec = "Spec 2 — full DDD",
    outcome = "log1p_vac",
    term    = "eli_feas:group_has_drones",
    beta    = round(ct2[grepl("eli_feas:group_has_drones", term), Estimate][1], 5),
    p_value = round(ct2[grepl("eli_feas:group_has_drones", term), `Pr(>|t|)`][1], 4)
  ),
  data.table(
    spec = "Spec 3 — triple w/ post_eli",
    outcome = "log1p_vac",
    term    = "eli_feas:post_eli:group_has_drones",
    beta    = round(ct3[grepl("post_eli:group_has_drones", term), Estimate][1], 5),
    p_value = round(ct3[grepl("post_eli:group_has_drones", term), `Pr(>|t|)`][1], 4)
  ),
  data.table(
    spec = "Spec 4 — selective",
    outcome = "log1p_vac_sel",
    term    = "eli_feas:group_has_drones",
    beta    = round(ct4[grepl("eli_feas:group_has_drones", term), Estimate][1], 5),
    p_value = round(ct4[grepl("eli_feas:group_has_drones", term), `Pr(>|t|)`][1], 4)
  ),
  data.table(
    spec = "Spec 5 — indiscriminate",
    outcome = "log1p_vac_ind",
    term    = "eli_feas:group_has_drones",
    beta    = round(ct5[grepl("eli_feas:group_has_drones", term), Estimate][1], 5),
    p_value = round(ct5[grepl("eli_feas:group_has_drones", term), `Pr(>|t|)`][1], 4)
  )
), fill = TRUE)
print(results)
fwrite(results, file.path(OUT, "R37_design4_main_results.csv"))

# ── 10) INTERPRETATION ────────────────────────────────────────────────
cat("\n== INTERPRETATION ==\n")
b_naive <- results[spec == "Spec 1 — naive (R32 eq.)", beta]
b_ddd   <- results[spec == "Spec 2 — full DDD",         beta]
p_ddd   <- results[spec == "Spec 2 — full DDD",         p_value]

cat(sprintf("Naive β (R32-equivalent):  %.5f\n", b_naive))
cat(sprintf("Full-DDD β (asymmetry net of capacity): %.5f  p = %.4f\n",
            b_ddd, p_ddd))

if (!is.na(b_naive) && !is.na(b_ddd)) {
  cat("\nThe DIFFERENCE between the naive and DDD coefficients is\n")
  cat("the amount absorbed by the FE structure — i.e., the part of\n")
  cat("the 'capacity' effect explained by group-specific narrative\n")
  cat("trends and baseline cell presence.\n\n")

  if (b_naive > 0 && b_ddd < 0) {
    cat("✓ EXPECTED RESULT: naive positive sign was capacity-driven.\n")
    cat("  Full DDD reveals the ASYMMETRY channel is negative, as theory\n")
    cat("  predicts: acquiring drones reduces the surveillance-driven\n")
    cat("  coercion pressure on this group.\n")
  } else if (b_naive > 0 && b_ddd > 0) {
    cat("→ Capacity dominates. Drones increase VAC even after netting\n")
    cat("  out group-specific trends. The asymmetry mechanism is\n")
    cat("  either absent or overwhelmed by capacity.\n")
  } else if (b_naive < 0 && b_ddd < 0) {
    cat("✓ Both coefficients negative — asymmetry is the sole channel.\n")
    cat("  Drones REDUCE VAC even naively, with stronger effect in DDD.\n")
  } else {
    cat("→ Mixed results — interpret carefully and report in appendix.\n")
  }
}

log_msg("R37 complete. Outputs in: ", OUT)
