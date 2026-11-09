# ============================================================
# R23_rebel_drone_decomp.R
# Rebel ISR capability decomposition
#
# THEORETICAL PREDICTION:
#   The displacement mechanism predicts that rebel groups facing
#   a PURELY ASYMMETRIC surveillance environment (state has ISR,
#   rebels do not) will show the strongest reactive violence.
#   Rebel groups that ALSO acquired drone capability can partially
#   reciprocate observation — reducing the asymmetry and therefore
#   the displacement incentive.
#
#   Prediction: β(rebel_WITHOUT_drones) >> β(rebel_WITH_drones)
#   If confirmed: this is the third layer of mechanism discrimination.
#     Layer 1: Rebel > State (displacement over substitution)
#     Layer 2: During ISR removal → rebel VAC falls (causal)
#     Layer 3: Strongest for rebels who CANNOT reciprocate (asymmetry)
#
# REBEL DRONE ACQUISITION DATES (sourced from SIPRI, Airwars,
# conflict documentation):
#
#   Afghanistan:
#     Taliban: basic commercial surveillance drones ~Jul 2019+
#              (pre-2019 = NO capability)
#
#   Iraq:
#     ISIS/ISIL: weaponised commercial drones ~Jul 2016+
#     Kata'ib Hezbollah / PMF: Iranian-supplied UAVs ~Jan 2018+
#
#   Libya:
#     LNA (Haftar): UAE-supplied Wing Loong II ~Jun 2019+
#     GNA militias: Turkey-supplied Bayraktar TB2 ~Oct 2019+
#
#   Pakistan:
#     TTP: NO drone capability throughout panel period
#
#   Somalia:
#     Al-Shabaab: NO drone capability throughout panel period
#                 (clearest asymmetric case)
#
#   Syria:
#     ISIS/ISIL: same as Iraq ~Jul 2016+
#     HTS (Hayat Tahrir al-Sham): limited drones ~Jun 2018+
#
#   Ukraine:
#     DPR/LPR separatists: Russian-supplied drones ~Apr 2014+
#                          (have drone capability throughout)
#
#   Yemen:
#     Houthis (Ansar Allah): Iranian-supplied Qasef-1/Shahed
#                            ~Jan 2017+ (documented strikes)
#                            Pre-2017 = marginal capability
#
# APPROACH:
#   Go back to raw ACLED, tag rebel VAC events by whether the
#   perpetrating actor (actor1) had drone capability at event time.
#   Aggregate to pg_id × month.
#   Run core DiD separately for each group.
#
# INPUT:  data/raw/acled/acled_full.csv
#         panel_pg_month_ALLCW_paper1_shocks.rds
# OUTPUT: R23_rebel_drone_decomp_*.csv
#         R23_rebel_drone_decomp_audit_*.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!requireNamespace("fixest", quietly = TRUE))
  stop("Package 'fixest' required.")
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_msg("R23: Rebel drone capability decomposition — ", TS)

# ── PATHS ─────────────────────────────────────────────────────
ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT, "data", "raw", "acled", "acled_full.csv")
)
acled_path <- Filter(file.exists, ACLED_PATHS)[1]
if (is.na(acled_path)) stop("ACLED raw CSV not found.")

P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL)

# ── 1) REBEL DRONE ACQUISITION TABLE ─────────────────────────
# Each row: actor name pattern (grepl), theatre, acquisition date
# after which the actor has drone capability.
# NA acquisition_date = never had capability in our panel period.

drone_acq <- rbindlist(list(
  # ── GROUPS WITH DRONE CAPABILITY ──────────────────────────
  # Afghanistan
  data.table(
    actor_pattern = "Taliban",
    theatre       = "Afghanistan",
    acq_date      = as.IDate("2019-07-01"),
    has_drones    = TRUE,
    source_note   = "Commercial surveillance drones documented 2019+ (Long War Journal)"
  ),
  # Iraq — ISIS
  data.table(
    actor_pattern = "Islamic State|ISIS|ISIL|Daesh",
    theatre       = "Iraq",
    acq_date      = as.IDate("2016-07-01"),
    has_drones    = TRUE,
    source_note   = "Weaponised DJI Phantom at Mosul Jul 2016 (West Point CTC)"
  ),
  # Iraq — Iranian-backed PMF
  data.table(
    actor_pattern = "Kata.ib Hezbollah|Popular Mobilization|PMF|Hashd",
    theatre       = "Iraq",
    acq_date      = as.IDate("2018-01-01"),
    has_drones    = TRUE,
    source_note   = "Iranian UAV transfers to PMF ~2018 (IISS, Airwars)"
  ),
  # Libya — LNA (Haftar)
  data.table(
    actor_pattern = "Libyan National Army|LNA|Haftar",
    theatre       = "Libya",
    acq_date      = as.IDate("2019-06-01"),
    has_drones    = TRUE,
    source_note   = "UAE Wing Loong II from June 2019 (SIPRI, UN Panel of Experts)"
  ),
  # Libya — GNA/Turkey-backed
  data.table(
    actor_pattern = "Government of National Accord|GNA|Tripoli|National Salvation",
    theatre       = "Libya",
    acq_date      = as.IDate("2019-10-01"),
    has_drones    = TRUE,
    source_note   = "Turkey Bayraktar TB2 from Oct 2019 (SIPRI)"
  ),
  # Syria — ISIS (same acquisition as Iraq)
  data.table(
    actor_pattern = "Islamic State|ISIS|ISIL|Daesh",
    theatre       = "Syria",
    acq_date      = as.IDate("2016-07-01"),
    has_drones    = TRUE,
    source_note   = "Same as Iraq — ISIS operated across both theatres"
  ),
  # Syria — HTS
  data.table(
    actor_pattern = "Hayat Tahrir|HTS|Jabhat|Nusra",
    theatre       = "Syria",
    acq_date      = as.IDate("2018-06-01"),
    has_drones    = TRUE,
    source_note   = "HTS limited commercial drone use from ~2018 (Airwars, PAX)"
  ),
  # Ukraine — separatists (have capability throughout our Ukraine period)
  data.table(
    actor_pattern = "DPR|LPR|Donetsk|Luhansk|Separatists|Russian-backed",
    theatre       = "Ukraine",
    acq_date      = as.IDate("2014-04-01"),
    has_drones    = TRUE,
    source_note   = "Russian-supplied Orlan-10 from start of conflict 2014"
  ),
  # Yemen — Houthis (pre-2017 only marginal)
  data.table(
    actor_pattern = "Houthi|Ansar Allah|Ansarallah",
    theatre       = "Yemen (North Yemen)",
    acq_date      = as.IDate("2017-01-01"),
    has_drones    = TRUE,
    source_note   = "Qasef-1 first documented strike Jan 2017 (Airwars, ACLED)"
  ),

  # ── GROUPS WITHOUT DRONE CAPABILITY (never in panel) ──────
  # Somalia — al-Shabaab (clearest case)
  data.table(
    actor_pattern = "Al.?Shabaab|Al Shabaab|Harakat",
    theatre       = "Somalia",
    acq_date      = NA_real_,
    has_drones    = FALSE,
    source_note   = "Al-Shabaab: no drone capability documented in any source"
  ),
  # Pakistan — TTP
  data.table(
    actor_pattern = "Tehrik.i.Taliban|TTP|Pakistani Taliban",
    theatre       = "Pakistan",
    acq_date      = NA_real_,
    has_drones    = FALSE,
    source_note   = "TTP: no drone capability documented (SIPRI, ICG)"
  ),
  # Pakistan — al-Qaeda
  data.table(
    actor_pattern = "Al.?Qaeda|AQIS|Al Qaida",
    theatre       = "Pakistan",
    acq_date      = NA_real_,
    has_drones    = FALSE,
    source_note   = "AQIS: no drone capability documented"
  )
))

drone_acq[, acq_date := as.IDate(acq_date)]

cat("\nRebel drone acquisition table:\n")
print(drone_acq[, .(actor_pattern, theatre, acq_date, has_drones)])

# ── 2) LOAD AND PROCESS RAW ACLED ────────────────────────────
log_msg("Reading raw ACLED...")
ac <- fread(acled_path, showProgress = FALSE, encoding = "UTF-8")
setnames(ac, tolower(names(ac)))
log_msg("  Raw ACLED rows: ", nrow(ac))

# Standardise date
if ("event_date" %in% names(ac)) {
  ac[, month := as.IDate(format(as.Date(event_date), "%Y-%m-01"))]
} else stop("event_date column not found in ACLED.")

# Filter to VAC events with rebel actor1
req <- c("event_type", "sub_event_type", "inter1", "inter2",
         "actor1", "latitude", "longitude", "fatalities")
miss <- setdiff(req, names(ac))
if (length(miss)) stop("ACLED missing columns: ", paste(miss, collapse = ", "))

# Rebel-initiated VAC
ac[, is_vac    := grepl("Violence against civilians", event_type, ignore.case = TRUE)]
ac[, is_rebel1 := inter1 == 2L]

rebel_vac <- ac[is_vac == TRUE & is_rebel1 == TRUE]
log_msg("  Rebel-initiated VAC events: ", nrow(rebel_vac))

# Map to PRIO-GRID
rebel_vac[, pg_row := floor((latitude  + 90)  / 0.5) + 1L]
rebel_vac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
rebel_vac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# ── 3) TAG EACH EVENT: DRONE-GROUP vs NO-DRONE-GROUP ─────────
log_msg("Tagging events by rebel drone capability...")

# Restrict to our 8 theatres by country name
theatre_country_map <- c(
  "Afghanistan"         = "Afghanistan",
  "Iraq"                = "Iraq",
  "Libya"               = "Libya",
  "Pakistan"            = "Pakistan",
  "Somalia"             = "Somalia",
  "Syria"               = "Syria",
  "Ukraine"             = "Ukraine",
  "Yemen"               = "Yemen (North Yemen)"
)

# ACLED country column name
country_col <- intersect(c("country", "country_name"), names(rebel_vac))[1]
if (is.na(country_col)) stop("No country column in ACLED.")

# Keep only events in our theatres
rebel_vac <- rebel_vac[get(country_col) %in% names(theatre_country_map)]
log_msg("  After theatre filter: ", nrow(rebel_vac), " rebel VAC events")

# Assign theatre_id
country_to_theatre <- names(theatre_country_map)
names(country_to_theatre) <- names(theatre_country_map)
# Reverse map: country name → theatre_id (same here)
rebel_vac[, theatre_id := get(country_col)]
# Fix Yemen name
rebel_vac[theatre_id == "Yemen", theatre_id := "Yemen (North Yemen)"]

# Tag each event
rebel_vac[, rebel_has_drone := FALSE]
rebel_vac[, drone_group_name := NA_character_]

for (i in seq_len(nrow(drone_acq))) {
  row      <- drone_acq[i]
  pattern  <- row$actor_pattern
  th       <- row$theatre
  acq      <- row$acq_date
  has_it   <- row$has_drones

  if (has_it) {
    # Actor matches pattern AND in correct theatre AND after acquisition date
    if (!is.na(acq)) {
      idx <- rebel_vac[
        grepl(pattern, actor1, ignore.case = TRUE, perl = TRUE) &
        theatre_id == th &
        month >= acq,
        which = TRUE
      ]
    } else {
      # Never acquired — skip (will remain FALSE)
      idx <- integer(0)
    }
    if (length(idx) > 0) {
      rebel_vac[idx, rebel_has_drone := TRUE]
      rebel_vac[idx, drone_group_name := sub("\\|.*", "", pattern)]
    }
  }
}

n_drone    <- sum(rebel_vac$rebel_has_drone)
n_nodrone  <- sum(!rebel_vac$rebel_has_drone)
log_msg("  Rebel VAC events — drone-group: ", n_drone,
        "  no-drone-group: ", n_nodrone)

# Audit: breakdown by theatre and group
audit_groups <- rebel_vac[, .(
  n_events      = .N,
  n_drone       = sum(rebel_has_drone),
  n_nodrone     = sum(!rebel_has_drone),
  pct_drone     = round(100 * mean(rebel_has_drone), 1)
), by = .(theatre_id)]
cat("\nEvents by theatre:\n")
print(audit_groups)

audit_actors <- rebel_vac[rebel_has_drone == TRUE,
                           .N, by = .(theatre_id, drone_group_name)][order(theatre_id)]
cat("\nDrone-group events by named actor:\n")
print(audit_actors)

# ── 4) AGGREGATE TO PG × MONTH ────────────────────────────────
log_msg("Aggregating to pg_id × month...")

panel_ref <- readRDS(P_PANEL); setDT(panel_ref)
if (!("pg_id" %in% names(panel_ref)) && "gid" %in% names(panel_ref))
  panel_ref[, pg_id := gid]
panel_cells  <- unique(as.integer(panel_ref$pg_id))
panel_months <- sort(unique(as.IDate(panel_ref$month)))

# Filter to panel cells
rebel_vac <- rebel_vac[pg_id %in% panel_cells]
rebel_vac[, month := as.IDate(format(month, "%Y-%m-01"))]

# Aggregate
agg <- rebel_vac[, .(
  rebel_vac_drone    = sum(rebel_has_drone  == TRUE),
  rebel_vac_nodrone  = sum(rebel_has_drone  == FALSE),
  rebel_vac_total    = .N
), by = .(pg_id, month)]

# Complete grid with zeros
grid <- CJ(pg_id = panel_cells, month = panel_months)
agg  <- merge(grid, agg, by = c("pg_id", "month"), all.x = TRUE)
for (cc in c("rebel_vac_drone","rebel_vac_nodrone","rebel_vac_total")) {
  agg[is.na(get(cc)), (cc) := 0L]
}

log_msg("Aggregated: ", nrow(agg), " cell-months")
log_msg("  Non-zero rebel_vac_drone:   ", agg[rebel_vac_drone   > 0, .N])
log_msg("  Non-zero rebel_vac_nodrone: ", agg[rebel_vac_nodrone > 0, .N])

# Save the aggregated series
out_agg <- file.path(DIRS$base,
                     paste0("rebel_vac_drone_decomp_pg_month.rds"))
saveRDS(agg, out_agg)
log_msg("Saved: ", out_agg)

# ── 5) MERGE INTO PANEL ───────────────────────────────────────
log_msg("Merging into panel...")
dt <- copy(panel_ref)
dt[, pg_id     := as.integer(pg_id)]
dt[, month     := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

# Drop if already present
for (cc in c("rebel_vac_drone","rebel_vac_nodrone","rebel_vac_total")) {
  if (cc %in% names(dt)) dt[, (cc) := NULL]
}
dt <- merge(dt, agg, by = c("pg_id","month"), all.x = TRUE)
for (cc in c("rebel_vac_drone","rebel_vac_nodrone","rebel_vac_total")) {
  dt[is.na(get(cc)), (cc) := 0L]
}

# Log transforms
dt[, log1p_rebel_drone   := log1p(rebel_vac_drone)]
dt[, log1p_rebel_nodrone := log1p(rebel_vac_nodrone)]
dt[, log1p_rebel_total   := log1p(rebel_vac_total)]

# Coordinates for Conley
if ("xcoord" %in% names(dt)) {
  dt[, longitude := xcoord]; dt[, latitude := ycoord]
} else stop("xcoord/ycoord not found.")

# Controls
CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str     <- paste(ctrl_present, collapse = " + ")
FE           <- "pg_id + theatre_id^month_int"

# Conley helper
safe_conley <- function(m, cut = 200) {
  r <- tryCatch(summary(m, vcov = conley(cutoff = cut)), error = function(e) NULL)
  if (!is.null(r)) return(list(m = r, vcov = paste0("Conley_", cut, "km")))
  r <- tryCatch(summary(m, vcov = vcov_conley(cutoff = cut)), error = function(e) NULL)
  if (!is.null(r)) return(list(m = r, vcov = paste0("Conley_", cut, "km")))
  list(m = summary(m, cluster = ~pg_id + month_int), vcov = "TwoWay_fallback")
}

extract_row <- function(m_summary, pattern, label, vcov_type) {
  ct  <- as.data.table(coeftable(m_summary), keep.rownames = "term")
  hit <- ct[grepl(pattern, term, perl = TRUE)][1]
  if (nrow(hit) == 0L || is.na(hit$Estimate)) {
    return(data.table(outcome = label, vcov_type = vcov_type,
                      coef = NA_real_, se = NA_real_, p = NA_real_,
                      n_nonzero = NA_integer_, nobs = NA_integer_))
  }
  # Count non-zero cells for the outcome variable
  yvar <- gsub("log1p_", "", label)
  nz   <- if (yvar %in% names(dt)) dt[get(yvar) > 0, .N] else NA_integer_
  data.table(
    outcome   = label,
    vcov_type = vcov_type,
    coef      = round(hit$Estimate[1],     5),
    se        = round(hit$`Std. Error`[1], 5),
    p         = round(hit$`Pr(>|t|)`[1],   4),
    n_nonzero = nz,
    nobs      = nobs(m_summary)
  )
}

# ── 6) RUN DECOMPOSED DiD ────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("REBEL DRONE CAPABILITY DECOMPOSITION — DiD\n")
cat("Prediction: β(nodrone) >> β(drone)\n")
cat("Displacement mechanism: strongest where asymmetry is greatest\n")
cat(strrep("=", 60), "\n\n")

OUTCOMES <- list(
  list(ylog = "log1p_rebel_nodrone", label = "rebel_vac_nodrone",
       desc  = "Rebel VAC — NO drone capability (pure asymmetry)"),
  list(ylog = "log1p_rebel_drone",   label = "rebel_vac_drone",
       desc  = "Rebel VAC — HAS drone capability (reduced asymmetry)"),
  list(ylog = "log1p_rebel_total",   label = "rebel_vac_total",
       desc  = "Rebel VAC — total (benchmark)")
)

all_rows <- list()

for (oc in OUTCOMES) {
  cat("── ", oc$desc, "──\n")

  fml <- as.formula(paste0(
    oc$ylog, " ~ post_eli:eli_feas + ", ctrl_str, " | ", FE
  ))

  m_raw <- tryCatch(
    feols(fml, data = dt, notes = FALSE),
    error = function(e) { log_msg("ERROR: ", e$message); NULL }
  )
  if (is.null(m_raw)) next

  # Primary: Conley 200km
  rc <- safe_conley(m_raw, 200)
  r  <- extract_row(rc$m, "post_eli.*eli_feas|eli_feas.*post_eli",
                    oc$label, rc$vcov)
  cat(sprintf("  β = %.5f  SE = %.5f  p = %.4f  (%s)\n",
              r$coef, r$se, r$p, rc$vcov))
  all_rows[[length(all_rows) + 1]] <- r

  # Secondary: cell clustering
  m_cl <- tryCatch(summary(m_raw, cluster = ~pg_id), error = function(e) NULL)
  if (!is.null(m_cl)) {
    r_cl <- extract_row(m_cl, "post_eli.*eli_feas|eli_feas.*post_eli",
                        oc$label, "cell_cluster")
    all_rows[[length(all_rows) + 1]] <- r_cl
  }
}

results_dt <- rbindlist(all_rows, fill = TRUE)

cat("\n", strrep("=", 60), "\n")
cat("FULL RESULTS TABLE (Conley 200km)\n")
cat(strrep("=", 60), "\n")
conley_only <- results_dt[grepl("Conley", vcov_type)]
print(conley_only[, .(outcome, coef, se, p, n_nonzero)])

# ── 7) COMPUTE THE RATIO ─────────────────────────────────────
b_nodrone <- conley_only[outcome == "rebel_vac_nodrone", coef]
b_drone   <- conley_only[outcome == "rebel_vac_drone",   coef]
b_total   <- conley_only[outcome == "rebel_vac_total",   coef]

cat("\n── Mechanism diagnostic ──\n")
cat(sprintf("β(rebel NO drones):   %.5f  p = %.4f\n",
            b_nodrone, conley_only[outcome == "rebel_vac_nodrone", p]))
cat(sprintf("β(rebel HAS drones):  %.5f  p = %.4f\n",
            b_drone,   conley_only[outcome == "rebel_vac_drone",   p]))
cat(sprintf("β(rebel total):       %.5f  p = %.4f\n",
            b_total,   conley_only[outcome == "rebel_vac_total",   p]))

if (!is.na(b_nodrone) && !is.na(b_drone) && b_drone != 0) {
  ratio <- round(b_nodrone / b_drone, 2)
  cat(sprintf("\nRatio β(nodrone) / β(drone): %.2f\n", ratio))
  if (ratio > 1) {
    cat("✓ SUPPORTS displacement via asymmetric surveillance.\n")
    cat("  Effect is concentrated in rebel groups that cannot reciprocate ISR.\n")
  } else {
    cat("✗ Does NOT support asymmetric displacement prediction.\n")
    cat("  Effect is similar (or larger) for drone-equipped rebels.\n")
  }
} else if (!is.na(b_drone) && b_drone == 0) {
  cat("\nβ(drone) = 0: ratio undefined. No effect for drone-equipped rebels.\n")
  cat("✓ This STRONGLY supports asymmetric displacement: only groups\n")
  cat("  that cannot reciprocate show reactive violence.\n")
}

# ── 8) THEATRE-LEVEL BREAKDOWN ────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("THEATRE BREAKDOWN — rebel_vac_nodrone (pure asymmetry cases)\n")
cat(strrep("=", 60), "\n")
cat("Theatres where rebels have NO drone capability:\n")
cat("  Somalia (al-Shabaab), Pakistan (TTP)\n")
cat("  + Afghanistan pre-2019, Yemen pre-2017, Iraq/Syria pre-2016\n\n")

# Split: theatre-periods where NO rebel group had drones
dt[, theatre_nodrone_period := FALSE]

# Somalia and Pakistan: always no drones
dt[theatre_id %in% c("Somalia", "Pakistan"), theatre_nodrone_period := TRUE]

# Afghanistan pre-2019
dt[theatre_id == "Afghanistan" & month < as.IDate("2019-07-01"),
   theatre_nodrone_period := TRUE]

# Yemen pre-2017 (Houthis acquired Jan 2017)
dt[theatre_id == "Yemen (North Yemen)" & month < as.IDate("2017-01-01"),
   theatre_nodrone_period := TRUE]

n_nodrone_period <- dt[theatre_nodrone_period == TRUE, .N]
cat("Cell-months in pure-asymmetry theatre-periods: ", n_nodrone_period, "\n\n")

if (n_nodrone_period > 500) {
  dt_na <- dt[theatre_nodrone_period == TRUE]
  fe_na <- if (uniqueN(dt_na$theatre_id) > 1) FE else "pg_id + month_int"

  dt_na[, log1p_y := log1p(events_type3)]
  fml_na <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + ", ctrl_str, " | ", fe_na
  ))
  m_na_raw <- tryCatch(feols(fml_na, data = dt_na, notes = FALSE),
                       error = function(e) { log_msg("ERROR: ", e$message); NULL })

  if (!is.null(m_na_raw)) {
    rc_na <- safe_conley(m_na_raw, 200)
    hit_na <- extract_row(rc_na$m,
                          "post_eli.*eli_feas|eli_feas.*post_eli",
                          "pure_asymmetry_theatres", rc_na$vcov)
    cat("DiD in pure-asymmetry theatre-periods (Conley 200km):\n")
    cat(sprintf("  β = %.5f  SE = %.5f  p = %.4f\n",
                hit_na$coef, hit_na$se, hit_na$p))
    cat(sprintf("  N = %d cell-months, %d theatres, FE = %s\n",
                n_nodrone_period, uniqueN(dt_na$theatre_id), fe_na))
    results_dt <- rbind(results_dt, hit_na, fill = TRUE)
  }
}

# ── 9) SAVE ───────────────────────────────────────────────────
out_res   <- file.path(DIRS$out_paper1,
                       paste0("R23_rebel_drone_decomp_", TS, ".csv"))
out_audit <- file.path(DIRS$audit,
                       paste0("R23_rebel_drone_decomp_audit_", TS, ".csv"))

fwrite(results_dt, out_res)
log_msg("Wrote results: ", out_res)

audit_out <- rbindlist(list(
  cbind(section = "theatre_breakdown", audit_groups),
  cbind(section = "actor_breakdown",
        if (nrow(audit_actors) > 0) audit_actors
        else data.table(theatre_id = NA, drone_group_name = NA, N = NA))
), fill = TRUE)
fwrite(audit_out, out_audit)
log_msg("Wrote audit: ", out_audit)

# ── 10) FINAL INTERPRETATION GUIDE ────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("INTERPRETATION GUIDE\n")
cat(strrep("=", 60), "\n")
cat("
Three possible outcomes and what they mean:

(A) β(nodrone) >> β(drone), both positive:
    STRONGLY supports asymmetric displacement.
    Groups that cannot reciprocate observation respond most violently.
    This is the third mechanism layer in the paper.

(B) β(nodrone) >> β(drone) ≈ 0:
    MOST STRONGLY supports asymmetric displacement.
    No effect at all for groups with own ISR — the surveillance
    asymmetry is the specific driver, not drone warfare generally.

(C) β(nodrone) ≈ β(drone), both positive:
    Displacement is not asymmetric — both types respond equally.
    Could indicate the mechanism is deterrence-related rather than
    surveillance-asymmetry-related, or that drone-group classification
    is too noisy at this spatial resolution.

(D) Both null:
    The ACLED rebel VAC disaggregation by specific actor names is too
    sparse to identify effects. Fall back to the aggregate rebel β.

In case (C) or (D): the original aggregate actor decomposition
(rebel β = 0.131, state β = 0.022) remains valid and the drone-group
split does not add to the paper. Report in appendix.
In case (A) or (B): this is a major addition to the paper.
")

log_msg("R23 complete.")
