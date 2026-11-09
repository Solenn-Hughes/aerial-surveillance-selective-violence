# ============================================================
# R24_subevent_decomp.R
# Rebel VAC sub-event type decomposition
#
# THEORETICAL PURPOSE:
#   The displacement mechanism has three named sub-mechanisms:
#     (1) Collaboration punishment — rebels punish civilians for
#         being passively observed (suspected informants).
#         Predicted signature: ABDUCTIONS / FORCED DISAPPEARANCES
#         rise in high-feasibility cells after ISR onset.
#         Rationale: targeted removal of suspected informants is
#         the most direct response to surveillance — it mirrors
#         Kalyvas's denunciation logic but with the information
#         channel reversed (punished for being watched, not for
#         actively informing).
#
#     (2) Operational displacement — persistent ISR pushes rebel
#         operations into civilian-dense areas where targeting is
#         harder (engagement rule constraints).
#         Predicted signature: ATTACKS rise broadly, especially
#         in dense-populated cells, as rebels embed in civilian
#         environments to complicate drone targeting.
#
#     (3) Escalatory pressure — loss of operational freedom under
#         ISR prompts rebels to intensify violence against soft
#         targets to demonstrate continued capability.
#         Predicted signature: ATTACKS with high fatality rates
#         rise, or SEXUAL VIOLENCE as collective coercion.
#
#   If ABDUCTIONS specifically drive the result, that is the
#   most theoretically precise finding: collaboration punishment
#   is the cleanest expression of the displacement mechanism and
#   maps directly onto Kalyvas (2006) run in reverse.
#
# ACLED sub-event types within Violence against civilians:
#   "Attack"                        — direct violence
#   "Abduction/forced disappearance" — targeted removal
#   "Sexual violence"               — collective coercion
#
# INPUT:  data/raw/acled/acled_full.csv
#         panel_pg_month_ALLCW_paper1_shocks.rds
# OUTPUT: R24_subevent_decomp_*.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!requireNamespace("fixest", quietly = TRUE)) stop("fixest required.")
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_msg("R24: Rebel VAC sub-event type decomposition — ", TS)

# ── PATHS ─────────────────────────────────────────────────────
ACLED_PATHS <- c(
  file.path(DIRS$data_raw, "acled", "acled_full.csv"),
  file.path(PROJECT_ROOT, "data", "raw", "acled", "acled_full.csv")
)
acled_path <- Filter(file.exists, ACLED_PATHS)[1]
if (is.na(acled_path)) stop("ACLED raw CSV not found.")

P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL)

# ── 1) LOAD ACLED ─────────────────────────────────────────────
log_msg("Reading raw ACLED...")
ac <- fread(acled_path, showProgress = FALSE, encoding = "UTF-8")
setnames(ac, tolower(names(ac)))
log_msg("  Raw rows: ", nrow(ac))

# Date
ac[, month := as.IDate(format(as.Date(event_date), "%Y-%m-01"))]

# ── 2) FILTER TO REBEL VAC ────────────────────────────────────
req <- c("event_type", "sub_event_type", "inter1",
         "actor1", "latitude", "longitude", "fatalities")
miss <- setdiff(req, names(ac))
if (length(miss)) stop("ACLED missing: ", paste(miss, collapse = ", "))

ac[, is_vac    := grepl("Violence against civilians", event_type, ignore.case = TRUE)]
ac[, is_rebel1 := inter1 == 2L]

rebel_vac <- ac[is_vac == TRUE & is_rebel1 == TRUE]
log_msg("  Rebel VAC events total: ", nrow(rebel_vac))

# Sub-event type inventory
cat("\nSub-event type distribution in rebel VAC:\n")
sub_tab <- rebel_vac[, .N, by = sub_event_type][order(-N)]
print(sub_tab)

# ── 3) CODE SUB-EVENT CATEGORIES ─────────────────────────────
rebel_vac[, sub_abduction := grepl(
  "Abduction|forced disappearance|kidnap",
  sub_event_type, ignore.case = TRUE
)]

rebel_vac[, sub_attack := grepl(
  "^Attack$|^Attack ",
  sub_event_type, ignore.case = TRUE
) & !sub_abduction]

rebel_vac[, sub_sexual := grepl(
  "Sexual violence|Sexual",
  sub_event_type, ignore.case = TRUE
)]

rebel_vac[, sub_other := !sub_abduction & !sub_attack & !sub_sexual]

cat("\nSub-event category counts:\n")
cat("  Abduction/forced disappearance: ", sum(rebel_vac$sub_abduction), "\n")
cat("  Direct attack:                  ", sum(rebel_vac$sub_attack),    "\n")
cat("  Sexual violence:                ", sum(rebel_vac$sub_sexual),    "\n")
cat("  Other:                          ", sum(rebel_vac$sub_other),     "\n")

# ── 4) MAP TO PRIO-GRID ───────────────────────────────────────
rebel_vac[, pg_row := floor((latitude  + 90)  / 0.5) + 1L]
rebel_vac[, pg_col := floor((longitude + 180) / 0.5) + 1L]
rebel_vac[, pg_id  := (pg_row - 1L) * 720L + pg_col]

# Restrict to panel theatres
country_col <- intersect(c("country", "country_name"), names(rebel_vac))[1]
theatre_countries <- c("Afghanistan", "Iraq", "Libya", "Pakistan",
                       "Somalia", "Syria", "Ukraine", "Yemen")
rebel_vac <- rebel_vac[get(country_col) %in% theatre_countries]
log_msg("  After theatre filter: ", nrow(rebel_vac), " rebel VAC events")

# ── 5) AGGREGATE TO PG × MONTH ────────────────────────────────
log_msg("Aggregating to pg_id × month...")

panel_ref <- readRDS(P_PANEL); setDT(panel_ref)
if (!("pg_id" %in% names(panel_ref)) && "gid" %in% names(panel_ref))
  panel_ref[, pg_id := gid]
panel_cells  <- unique(as.integer(panel_ref$pg_id))
panel_months <- sort(unique(as.IDate(panel_ref$month)))

rebel_vac <- rebel_vac[pg_id %in% panel_cells]

agg <- rebel_vac[, .(
  rebel_vac_abduction = sum(sub_abduction),
  rebel_vac_attack    = sum(sub_attack),
  rebel_vac_sexual    = sum(sub_sexual),
  rebel_vac_other     = sum(sub_other),
  rebel_vac_total     = .N,
  rebel_fatalities    = sum(fatalities, na.rm = TRUE)
), by = .(pg_id, month)]

# Complete grid
grid <- CJ(pg_id = panel_cells, month = panel_months)
agg  <- merge(grid, agg, by = c("pg_id", "month"), all.x = TRUE)
fill_cols <- c("rebel_vac_abduction","rebel_vac_attack",
               "rebel_vac_sexual","rebel_vac_other",
               "rebel_vac_total","rebel_fatalities")
for (cc in fill_cols) agg[is.na(get(cc)), (cc) := 0L]

log_msg("Aggregated: ", nrow(agg), " cell-months")
for (cc in fill_cols) {
  log_msg("  Non-zero ", cc, ": ", agg[get(cc) > 0, .N])
}

# ── 6) MERGE INTO PANEL ───────────────────────────────────────
dt <- copy(panel_ref)
dt[, pg_id     := as.integer(pg_id)]
dt[, month     := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

for (cc in fill_cols) if (cc %in% names(dt)) dt[, (cc) := NULL]
dt <- merge(dt, agg, by = c("pg_id","month"), all.x = TRUE)
for (cc in fill_cols) dt[is.na(get(cc)), (cc) := 0L]

# Log transforms
for (cc in fill_cols) {
  dt[, paste0("log1p_", cc) := log1p(get(cc))]
}

# Coordinates
if ("xcoord" %in% names(dt)) {
  dt[, longitude := xcoord]; dt[, latitude := ycoord]
} else stop("xcoord/ycoord not found.")

# Controls
CTRL     <- c("nlights_calib_mean","nlights_mean","nlights_min",
              "prec_gpcp","temp","water_ih")
ctrl_str <- paste(intersect(CTRL, names(dt)), collapse = " + ")
FE       <- "pg_id + theatre_id^month_int"

# Conley helper
safe_conley <- function(m, cut = 200) {
  r <- tryCatch(summary(m, vcov = conley(cutoff = cut)), error = function(e) NULL)
  if (!is.null(r)) return(list(m = r, vcov = paste0("Conley_", cut, "km")))
  r <- tryCatch(summary(m, vcov = vcov_conley(cutoff = cut)), error = function(e) NULL)
  if (!is.null(r)) return(list(m = r, vcov = paste0("Conley_", cut, "km")))
  list(m = summary(m, cluster = ~pg_id + month_int), vcov = "TwoWay_fallback")
}

extract_row <- function(m_sum, label, vcov_type, yvar) {
  ct  <- as.data.table(coeftable(m_sum), keep.rownames = "term")
  hit <- ct[grepl("post_eli.*eli_feas|eli_feas.*post_eli", term, perl = TRUE)][1]
  nz  <- if (yvar %in% names(dt)) dt[get(yvar) > 0, .N] else NA_integer_
  if (nrow(hit) == 0L || is.na(hit$Estimate)) {
    return(data.table(outcome = label, vcov_type = vcov_type,
                      coef = NA_real_, se = NA_real_, p = NA_real_,
                      n_nonzero = nz, nobs = NA_integer_))
  }
  data.table(
    outcome   = label,
    vcov_type = vcov_type,
    coef      = round(hit$Estimate[1],     5),
    se        = round(hit$`Std. Error`[1], 5),
    p         = round(hit$`Pr(>|t|)`[1],   4),
    n_nonzero = nz,
    nobs      = nobs(m_sum)
  )
}

# ── 7) RUN DiD FOR EACH SUB-EVENT TYPE ───────────────────────
cat("\n", strrep("=", 60), "\n")
cat("SUB-EVENT TYPE DECOMPOSITION — DiD (Conley 200km)\n")
cat("Theoretical predictions:\n")
cat("  Abduction >> Attack >> Sexual/Other: collaboration punishment\n")
cat("  Attack >> Abduction: escalatory pressure\n")
cat("  All equal: cannot discriminate sub-mechanisms\n")
cat(strrep("=", 60), "\n\n")

OUTCOMES <- list(
  list(yvar = "rebel_vac_abduction", label = "Abduction/forced disappearance",
       mechanism = "Collaboration punishment — targeted removal of suspected informants"),
  list(yvar = "rebel_vac_attack",    label = "Direct attack",
       mechanism = "Escalatory pressure / operational displacement"),
  list(yvar = "rebel_vac_sexual",    label = "Sexual violence",
       mechanism = "Collective coercion"),
  list(yvar = "rebel_vac_total",     label = "Total rebel VAC",
       mechanism = "Benchmark (should match R19/R23)")
)

all_rows <- list()

for (oc in OUTCOMES) {
  ylog <- paste0("log1p_", oc$yvar)
  cat(sprintf("── %s ──\n   (%s)\n", oc$label, oc$mechanism))

  fml <- as.formula(paste0(
    ylog, " ~ post_eli:eli_feas + ", ctrl_str, " | ", FE
  ))

  m_raw <- tryCatch(
    feols(fml, data = dt, notes = FALSE),
    error = function(e) { log_msg("  ERROR: ", e$message); NULL }
  )
  if (is.null(m_raw)) { cat("  SKIPPED (model error)\n\n"); next }

  # Primary: Conley 200km
  rc <- safe_conley(m_raw, 200)
  r  <- extract_row(rc$m, oc$label, rc$vcov, oc$yvar)
  cat(sprintf("  β = %7.5f  SE = %7.5f  p = %.4f  n_nonzero = %d\n\n",
              r$coef, r$se, r$p, r$n_nonzero))
  all_rows[[length(all_rows) + 1]] <- r

  # Secondary: cell clustering
  m_cl <- tryCatch(summary(m_raw, cluster = ~pg_id), error = function(e) NULL)
  if (!is.null(m_cl)) {
    r_cl <- extract_row(m_cl, oc$label, "cell_cluster", oc$yvar)
    all_rows[[length(all_rows) + 1]] <- r_cl
  }
}

results_dt <- rbindlist(all_rows, fill = TRUE)
conley_res <- results_dt[grepl("Conley", vcov_type)]

cat(strrep("=", 60), "\n")
cat("RESULTS TABLE (Conley 200km)\n")
cat(strrep("=", 60), "\n")
print(conley_res[, .(outcome, coef, se, p, n_nonzero)])

# ── 8) SUB-MECHANISM DIAGNOSTICS ─────────────────────────────
cat("\n── Sub-mechanism diagnostic ──\n")

b_abduct <- conley_res[outcome == "Abduction/forced disappearance", coef]
b_attack <- conley_res[outcome == "Direct attack",                  coef]
b_sexual <- conley_res[outcome == "Sexual violence",                coef]
p_abduct <- conley_res[outcome == "Abduction/forced disappearance", p]
p_attack <- conley_res[outcome == "Direct attack",                  p]

cat(sprintf("β(abduction): %.5f  p = %.4f\n", b_abduct, p_abduct))
cat(sprintf("β(attack):    %.5f  p = %.4f\n", b_attack, p_attack))
cat(sprintf("β(sexual):    %.5f  p = %.4f\n", b_sexual,
            conley_res[outcome == "Sexual violence", p]))

cat("\n")
if (!is.na(b_abduct) && !is.na(b_attack)) {
  if (b_abduct > b_attack && p_abduct < 0.10) {
    cat("✓ ABDUCTIONS dominate:\n")
    cat("  Collaboration punishment is the primary sub-mechanism.\n")
    cat("  Rebels are responding to surveillance by targeting civilians\n")
    cat("  suspected of providing (passive or active) intelligence.\n")
    cat("  This is Kalyvas (2006) denunciation logic run in reverse:\n")
    cat("  the informational trigger is being OBSERVED, not informing.\n")
  } else if (b_attack > b_abduct && p_attack < 0.10) {
    cat("✓ DIRECT ATTACKS dominate:\n")
    cat("  Escalatory pressure or operational displacement is primary.\n")
    cat("  Rebels intensify broad violence rather than targeted removal.\n")
    cat("  Consistent with: embedding in civilian areas to complicate\n")
    cat("  drone targeting, or demonstrating capability under pressure.\n")
  } else if (b_abduct > 0 && b_attack > 0) {
    cat("→ Both abduction and attack are positive.\n")
    cat("  Multiple sub-mechanisms may operate simultaneously.\n")
    cat("  Cannot discriminate cleanly at this spatial resolution.\n")
  } else {
    cat("→ Neither abduction nor attack is clearly dominant.\n")
    cat("  Sub-mechanism decomposition inconclusive — report in appendix.\n")
  }
}

# ── 9) ABDUCTION QUASI-EXPERIMENT (Bagram reversal) ──────────
# Does abduction specifically reverse when ISR is removed?
# This is the sharpest test of collaboration punishment.
cat("\n", strrep("=", 60), "\n")
cat("ABDUCTION QUASI-EXPERIMENT: Bagram reversal\n")
cat("If collaboration punishment drives the effect:\n")
cat("β(abduction × shock_bagram) should be NEGATIVE (reversal)\n")
cat(strrep("=", 60), "\n\n")

if ("shock_bagram" %in% names(dt)) {
  afg <- dt[theatre_id == "Afghanistan"]
  cat("Afghanistan cell-months: ", nrow(afg), "\n")
  cat("Non-zero abductions: ", afg[rebel_vac_abduction > 0, .N], "\n\n")

  if (afg[rebel_vac_abduction > 0, .N] > 10) {
    fe_afg <- "pg_id + month_int"   # single theatre — no theatre_id^month FE

    f_abd_bag <- as.formula(paste0(
      "log1p_rebel_vac_abduction ~ post_eli:eli_feas +",
      " post_eli:eli_feas:shock_bagram + ",
      ctrl_str, " | ", fe_afg
    ))

    m_abd_bag <- tryCatch(
      feols(f_abd_bag, data = afg, notes = FALSE),
      error = function(e) { log_msg("ERROR: ", e$message); NULL }
    )

    if (!is.null(m_abd_bag)) {
      rc_ab <- safe_conley(m_abd_bag, 200)
      ct_ab <- as.data.table(coeftable(rc_ab$m), keep.rownames = "term")
      trip  <- ct_ab[grepl("shock_bagram", term)]
      if (nrow(trip) > 0) {
        cat(sprintf("Triple interaction (Bagram × eli_feas on abductions):\n"))
        cat(sprintf("  β = %.5f  SE = %.5f  p = %.4f  (%s)\n",
                    trip$Estimate[1], trip$`Std. Error`[1],
                    trip$`Pr(>|t|)`[1], rc_ab$vcov))
        if (!is.na(trip$Estimate[1]) && trip$Estimate[1] < 0) {
          cat("  → NEGATIVE: abductions FALL when Bagram closes.\n")
          cat("    Collaboration punishment reversal confirmed.\n")
        } else if (!is.na(trip$Estimate[1]) && trip$Estimate[1] > 0) {
          cat("  → POSITIVE: abductions rise (unexpected for this mechanism).\n")
        }
      }
    }
  } else {
    cat("Too few abduction events in Afghanistan for quasi-experiment.\n")
  }
}

# ── 10) SAVE ──────────────────────────────────────────────────
out_res <- file.path(DIRS$out_paper1,
                     paste0("R24_subevent_decomp_", TS, ".csv"))
fwrite(results_dt, out_res)
log_msg("Wrote: ", out_res)

out_audit <- file.path(DIRS$audit,
                       paste0("R24_subevent_decomp_audit_", TS, ".csv"))
fwrite(sub_tab, out_audit)
log_msg("Wrote audit: ", out_audit)

# ── 11) INTERPRETATION ────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("WHAT THIS MEANS FOR THE PAPER\n")
cat(strrep("=", 60), "\n")
cat("
Layer 1 (rebel vs state):    displacement over substitution
Layer 2 (drone asymmetry):   effect concentrated where rebels cannot reciprocate
Layer 3 (sub-event type):    [see diagnostic above]

If abductions dominate → collaboration punishment sub-mechanism.
  Write: 'The sub-event decomposition reveals that the displacement
  effect operates primarily through abductions and forced disappearances
  (β = X, p = Y) rather than direct attacks (β = Z, p = W). This is
  consistent with collaboration punishment: rebels respond to aerial
  observation by removing civilians suspected of providing intelligence,
  whether actively or passively. This is Kalyvas's (2006) denunciation
  dynamic run in reverse — the informational trigger is being watched,
  not informing.'

If attacks dominate → escalatory pressure / operational displacement.
  Write: 'The sub-event decomposition reveals that direct attacks drive
  the effect (β = X, p = Y), while targeted abductions are smaller
  (β = Z, p = W). This is consistent with rebels embedding in
  civilian-dense environments to complicate drone targeting, or
  intensifying violence to demonstrate operational capacity under
  persistent surveillance pressure.'

If both → all three sub-mechanisms may operate; report descriptively.
")

log_msg("R24 complete.")
