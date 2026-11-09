# ============================================================================
# R27_civilian_distinctions.R
#
# Purpose: Close the theory–method gap flagged in pre-presentation audit.
#
#   The existing analysis conflates (i) selective, low-fatality, targeted
#   rebel violence with (ii) indiscriminate, high-fatality collective
#   violence. Both are counted equally as UCDP type-3 events. The Kalyvas
#   framework the theory inherits depends on this distinction.
#
#   This script adds three tests — all using data already in the project.
#
#   TEST A. Selective vs indiscriminate (fatalities-per-event split)
#     Under collaboration-punishment → selective (low-fatality) dominates.
#     Under escalatory pressure / operational displacement → indiscriminate
#     (high-fatality) dominates. Both together → mixed positive signs.
#
#   TEST B. Territorial control interaction (Kalyvas moderator)
#     Kalyvas predicts selective violence concentrates in zones of PARTIAL
#     control. We proxy control using the Herfindahl of battle actors in a
#     cell-year: high-HHI = contested, low-HHI = single-actor dominated.
#
#   TEST C. Mechanism signature table
#     Writes a publication-ready table mapping every sub-mechanism to its
#     distinctive empirical signature, with the β and p-value from the
#     corresponding regression. Makes the theory-method fit explicit.
#
# Author: Solenn Hughes
# Date:   (fill in on run)
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})

source("R00_config.R")
log_msg("R27: civilian distinctions — selective vs indiscriminate + control")

OUT <- file.path("outputs", "paper1", "civilian_distinctions")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── LOAD PANEL + RAW ACLED ─────────────────────────────────────────────
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(panel_path); setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]

ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT,  "data", "raw", "acled", "acled_full.csv")
)
acled_path <- Filter(file.exists, ACLED_PATHS)[1]
if (is.na(acled_path)) stop("ACLED raw CSV not found.")
ac <- fread(acled_path, showProgress = FALSE, encoding = "UTF-8")
setnames(ac, tolower(names(ac)))
ac[, month := as.IDate(format(as.Date(event_date), "%Y-%m-01"))]
ac[, fatalities := as.integer(fatalities)]
ac[is.na(fatalities), fatalities := 0L]

# Map each event to a PRIO-GRID cell (0.5-degree grid)
ac[, pg_row := floor((latitude  + 90)  / 0.5) + 1L]
ac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
ac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# Restrict to rebel-initiated VAC events (the mechanism object)
rebel_vac <- ac[
  inter1 == 2L &
    grepl("Violence against civilians", event_type, ignore.case = TRUE)
]
log_msg("Rebel-initiated VAC events: ", nrow(rebel_vac))

# ═══════════════════════════════════════════════════════════════════════
# TEST A — SELECTIVE VS INDISCRIMINATE (fatalities-per-event)
# ═══════════════════════════════════════════════════════════════════════
# Definition:
#   SELECTIVE      = 1 or 2 fatalities per event, attack or abduction
#                    (targeted, named-individual killing in Kalyvas' sense)
#   INDISCRIMINATE = 5+ fatalities per event
#                    (mass-casualty collective violence)
#   Events with 3–4 fatalities are ambiguous; we assign them to "mid"
#   and report separately so the reader can see the boundary case.
#
# This mirrors Kalyvas (2006: 142–143) on selective targeting and
# Eck & Hultman (2007) on one-sided mass killing.

rebel_vac[, fat_class := fcase(
  fatalities <= 2, "selective",
  fatalities <= 4, "mid",
  fatalities >= 5, "indiscriminate",
  default         = NA_character_
)]

cat("\nRebel VAC event distribution by fatality class:\n")
print(rebel_vac[, .(events = .N, fatalities = sum(fatalities)),
                by = fat_class][order(-events)])

agg_fat <- rebel_vac[, .(
  rebel_vac_selective      = sum(fat_class == "selective",       na.rm = TRUE),
  rebel_vac_indiscriminate = sum(fat_class == "indiscriminate",  na.rm = TRUE),
  rebel_vac_mid            = sum(fat_class == "mid",             na.rm = TRUE),
  rebel_vac_allfat         = sum(fatalities,                     na.rm = TRUE)
), by = .(pg_id, month)]

# Merge into main panel (zero-fill cells with no rebel VAC)
for (cc in names(agg_fat)[-(1:2)]) if (cc %in% names(dt)) dt[, (cc) := NULL]
dt <- merge(dt, agg_fat, by = c("pg_id", "month"), all.x = TRUE)
for (cc in c("rebel_vac_selective", "rebel_vac_indiscriminate",
             "rebel_vac_mid", "rebel_vac_allfat")) {
  dt[is.na(get(cc)), (cc) := 0L]
}

# Add coordinates for Conley (if not already present)
if (!all(c("lat", "lon") %in% names(dt))) {
  if (all(c("ycoord", "xcoord") %in% names(dt))) {
    dt[, `:=`(lat = ycoord, lon = xcoord)]
  } else {
    stop("Need coordinates for Conley SEs")
  }
}

# Log-transform outcomes
dt[, log1p_rv_sel     := log1p(rebel_vac_selective)]
dt[, log1p_rv_indis   := log1p(rebel_vac_indiscriminate)]
dt[, log1p_rv_mid     := log1p(rebel_vac_mid)]

CTRL <- c("nlights_calib_mean", "prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(CTRL, names(dt)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

run_model <- function(yvar, label) {
  f <- as.formula(paste0(yvar, " ~ post_eli:eli_feas + ", ctrl_str, " | ", FE))
  m <- feols(f, data = dt, notes = FALSE)
  cv <- vcov_conley(m, lat = "lat", lon = "lon", cutoff = 200)
  ct <- as.data.table(coeftable(m, vcov = cv), keep.rownames = "term")
  r <- ct[grepl("eli_feas", term)]
  data.table(
    outcome  = label,
    beta     = round(r$Estimate[1], 4),
    se       = round(r$`Std. Error`[1], 4),
    p_value  = round(r$`Pr(>|t|)`[1], 4),
    n_nonzero = sum(dt[[gsub("log1p_rv_", "rebel_vac_", yvar)]] > 0, na.rm = TRUE)
  )
}

cat("\n── TEST A: Selective vs Indiscriminate rebel VAC (Conley 200km) ──\n")
testA <- rbindlist(list(
  run_model("log1p_rv_sel",   "Rebel VAC — SELECTIVE (≤2 fatal.)"),
  run_model("log1p_rv_mid",   "Rebel VAC — mid (3–4 fatal.)"),
  run_model("log1p_rv_indis", "Rebel VAC — INDISCRIMINATE (≥5 fatal.)")
))
print(testA)
fwrite(testA, file.path(OUT, "testA_selective_vs_indiscriminate.csv"))

# Ratio diagnostic
b_sel   <- testA[grepl("SELECTIVE",      outcome), beta]
b_indis <- testA[grepl("INDISCRIMINATE", outcome), beta]
cat(sprintf("\n  β(indiscriminate) / β(selective) = %.2f\n", b_indis / b_sel))
cat("  Interpretation:\n")
if (b_indis > b_sel & b_indis > 0) {
  cat("    ✓ Indiscriminate component dominates — consistent with\n")
  cat("      escalatory pressure / operational displacement,\n")
  cat("      NOT narrow collaboration-punishment.\n")
} else if (b_sel > b_indis & b_sel > 0) {
  cat("    ✓ Selective component dominates — consistent with\n")
  cat("      collaboration-punishment (Kalyvas reverse).\n")
} else {
  cat("    → Both components positive; mechanism is mixed.\n")
}


# ═══════════════════════════════════════════════════════════════════════
# TEST B — TERRITORIAL-CONTROL INTERACTION
# ═══════════════════════════════════════════════════════════════════════
# Proxy for Kalyvas' control dimension: for each cell-year, compute the
# Herfindahl index (HHI) of battle-fatality shares across armed actors
# that engaged in battles in that cell. Low HHI ≈ contested (multiple
# actors fighting) → Kalyvas' zone of PARTIAL control where selective
# violence should concentrate. High HHI ≈ one actor dominant.

# Build from raw ACLED battle events
ac_batt <- ac[grepl("^Battles", event_type, ignore.case = TRUE) &
                !is.na(actor1) & !is.na(latitude) & !is.na(longitude)]
ac_batt[, year := year(month)]

# Fatalities per actor per cell-year
actor_share <- ac_batt[, .(fat = sum(fatalities, na.rm = TRUE)),
                       by = .(pg_id, year, actor1)]
actor_share[, total := sum(fat), by = .(pg_id, year)]
actor_share[, share := fat / total]
control_hhi <- actor_share[total > 0, .(
  hhi       = sum(share^2),
  n_actors  = .N,
  total_fat = total[1]
), by = .(pg_id, year)]
# Inverted HHI: 1 = maximally contested (partial control), 0 = single dominant actor
control_hhi[, contested := 1 - hhi]

# If a cell-year has no battles → mark as "uncontested/uninformative"
# We treat missing as contested = 0 (no observed contestation)
dt[, year := year(month)]
dt <- merge(dt, control_hhi[, .(pg_id, year, contested, n_actors)],
            by = c("pg_id", "year"), all.x = TRUE)
dt[is.na(contested), contested := 0]
dt[is.na(n_actors),  n_actors  := 0L]

# Z-standardise contested for interpretable coefficients
dt[, contested_z := (contested - mean(contested, na.rm = TRUE)) /
     sd(contested, na.rm = TRUE)]

cat("\n── TEST B: Kalyvas territorial-control interaction ──\n")

# Model 1: simple interaction with treatment
fB1 <- as.formula(paste0(
  "log1p_rv_sel ~ post_eli:eli_feas + post_eli:eli_feas:contested_z + ",
  ctrl_str, " | ", FE
))
mB1 <- feols(fB1, data = dt, notes = FALSE)
cvB1 <- vcov_conley(mB1, lat = "lat", lon = "lon", cutoff = 200)
ctB1 <- as.data.table(coeftable(mB1, vcov = cvB1), keep.rownames = "term")
cat("Selective rebel VAC × contested control:\n")
print(ctB1[grepl("eli_feas|contested", term)])

# Model 2: same for indiscriminate outcome
fB2 <- as.formula(paste0(
  "log1p_rv_indis ~ post_eli:eli_feas + post_eli:eli_feas:contested_z + ",
  ctrl_str, " | ", FE
))
mB2 <- feols(fB2, data = dt, notes = FALSE)
cvB2 <- vcov_conley(mB2, lat = "lat", lon = "lon", cutoff = 200)
ctB2 <- as.data.table(coeftable(mB2, vcov = cvB2), keep.rownames = "term")
cat("\nIndiscriminate rebel VAC × contested control:\n")
print(ctB2[grepl("eli_feas|contested", term)])

# Clean summary
extract_iaB <- function(ct, outcome_label) {
  r1 <- ct[grepl("post_eli:eli_feas$", term)]
  r2 <- ct[grepl("contested_z", term)]
  data.table(
    outcome = outcome_label,
    beta_base   = round(r1$Estimate[1], 4),
    p_base      = round(r1$`Pr(>|t|)`[1], 4),
    beta_inter  = round(r2$Estimate[1], 4),
    p_inter     = round(r2$`Pr(>|t|)`[1], 4)
  )
}
testB <- rbindlist(list(
  extract_iaB(ctB1, "Selective rebel VAC"),
  extract_iaB(ctB2, "Indiscriminate rebel VAC")
))
fwrite(testB, file.path(OUT, "testB_territorial_control.csv"))

cat("\nKalyvas prediction:\n")
cat("  If β_inter > 0 for SELECTIVE → effect concentrates in contested cells.\n")
cat("  If β_inter ~ 0 (or negative) → surveillance bypasses the control logic,\n")
cat("    consistent with the claim that sensor-based intelligence does NOT\n")
cat("    need local informants and thus decouples from Kalyvas' control moderator.\n")


# ═══════════════════════════════════════════════════════════════════════
# TEST C — MECHANISM SIGNATURE TABLE
# ═══════════════════════════════════════════════════════════════════════
# Maps each theoretical sub-mechanism to its predicted signature and the
# empirical finding. This is the publication-ready theory-method fit
# table for the dissertation.

signature <- data.table(
  sub_mechanism = c(
    "Substitution (state targeting)",
    "Displacement: collaboration-punishment",
    "Displacement: escalatory pressure",
    "Displacement: operational displacement",
    "Observer-effect (infrastructure, not optics)",
    "Asymmetric surveillance"
  ),
  empirical_signature = c(
    "β(state VAC) > 0",
    "β(rebel abductions) > β(rebel attacks); β(selective) > β(indiscriminate)",
    "β(rebel attacks) > β(rebel abductions); β(indiscriminate) > β(selective)",
    "β(rebel VAC) rises in civilian-dense cells; high-fatality dominates",
    "Cloud-cover does not moderate β",
    "β(rebels without drones) > β(rebels with drones)"
  ),
  finding = c(
    "NULL (β = 0.022, p = 0.387) — rejected",
    "Partial: abductions smaller than attacks; selective vs indiscriminate → see TEST A",
    "Supported: attacks > abductions; indiscriminate ratio → see TEST A",
    "See TEST A indiscriminate component",
    "Supported: cloud-cover β = –0.002, p = 0.874",
    "Supported: no-drone β / drone β = 2.23"
  ),
  status = c("REJECTED", "PARTIAL", "SUPPORTED", "SEE TEST A", "SUPPORTED", "SUPPORTED")
)

fwrite(signature, file.path(OUT, "testC_mechanism_signature.csv"))

cat("\n── TEST C: Mechanism signature table written. ──\n")
print(signature)

# ═══════════════════════════════════════════════════════════════════════
# SUMMARY — what each test contributes to the dissertation
# ═══════════════════════════════════════════════════════════════════════
sink(file.path(OUT, "R27_summary.txt"))
cat("R27 — CIVILIAN DISTINCTIONS CLOSURE\n")
cat(strrep("=", 60), "\n\n")

cat("TEST A — Selective vs Indiscriminate (NEW KEY RESULT)\n")
cat(strrep("-", 60), "\n")
print(testA); cat("\n")

cat("TEST B — Territorial control moderator (Kalyvas proxy)\n")
cat(strrep("-", 60), "\n")
print(testB); cat("\n")

cat("TEST C — Mechanism signature table\n")
cat(strrep("-", 60), "\n")
print(signature); cat("\n")

cat("\nINTERPRETATION FOR DISSERTATION:\n")
cat(strrep("-", 60), "\n")
cat("If TEST A shows β(indiscriminate) > β(selective), that is the\n")
cat("empirically sharp version of the attack-over-abduction finding.\n")
cat("It means the displacement mechanism operates through BROAD coercion\n")
cat("of the observed population, not through narrow Kalyvas-style\n")
cat("informant-targeting. This UPDATES the Kalyvas framework rather\n")
cat("than simply reverse-applying it.\n\n")
cat("If TEST B shows that contested-control does NOT moderate the\n")
cat("selective-violence effect, that is FURTHER evidence that sensor-\n")
cat("based intelligence bypasses the denunciation-frontier logic.\n")
cat("The geography of surveillance replaces the geography of control\n")
cat("as the organising principle of civilian victimisation under\n")
cat("persistent aerial ISR.\n")
sink()

log_msg("R27 complete. Outputs in: ", OUT)