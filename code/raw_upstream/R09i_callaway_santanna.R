# ============================================================
# R09i_callaway_santanna.R  — FIXED (self-contained)
#
# FIX: No longer depends on dt_clean / log1p_events_type3
# existing in GlobalEnv. Loads the panel itself.
#
# Callaway-Sant'Anna (2021) heterogeneity-robust DiD.
# Addresses the concern that staggered treatment onsets across
# 8 theatres can bias TWFE estimates when treatment effects
# are heterogeneous (Callaway & Sant'Anna 2021; Roth et al. 2023).
#
# Setup:
#   Unit:    pg_id (grid cell)
#   Time:    month_int
#   Group:   theatre ELI onset month_int (0 = never treated)
#   Outcome: log1p(events_type3)
#
# Two control groups:
#   "nevertreated"  — cells whose theatre never receives ELI
#   "notyettreated" — cells not yet treated at each calendar time
#
# Runtime: ~10–30 min depending on machine.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!requireNamespace("did", quietly = TRUE))
  stop("Package 'did' required. Install via install.packages('did').")
if (!requireNamespace("fixest", quietly = TRUE))
  stop("Package 'fixest' required. Install via install.packages('fixest').")

suppressMessages(library(did))
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_msg("Step 9i: Callaway-Sant'Anna DiD (memory-safe version)")

# ── 1) LOAD PANEL ─────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
if (!file.exists(P_PANEL)) {
  P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
  if (!file.exists(P_PANEL)) stop("No panel found. Run R16 (or R08v2) first.")
  log_msg("Using posteli_v2 panel (shocks panel not found).")
}

dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) {
  dt[, month := as_month_idate(month)]
} else {
  dt[, month := as.IDate(month)]
}
dt[, month_int := year(month) * 12L + month(month)]

log_msg("Panel loaded: ", nrow(dt), " rows, ",
        uniqueN(dt$pg_id), " cells, ",
        uniqueN(dt$theatre_id), " theatres")

# ── 2) OUTCOME ─────────────────────────────────────────────────
if (!("events_type3" %in% names(dt)))
  stop("events_type3 not found. Check R03/R04 pipeline.")

dt[, log1p_events_type3 := log1p(events_type3)]

# ── 3) BUILD ONSET MAP ─────────────────────────────────────────
# ELI onset month_int = year * 12 + month (1-indexed)
onset_map <- c(
  "Afghanistan"         = 2001L * 12L + 10L,  # Oct 2001
  "Iraq"                = 2014L * 12L + 8L,   # Aug 2014
  "Libya"               = 2011L * 12L + 3L,   # Mar 2011
  "Pakistan"            = 2004L * 12L + 6L,   # Jun 2004
  "Somalia"             = 2007L * 12L + 1L,   # Jan 2007
  "Syria"               = 2014L * 12L + 9L,   # Sep 2014
  "Ukraine"             = 2022L * 12L + 2L,   # Feb 2022
  "Yemen (North Yemen)" = 2009L * 12L + 12L   # Dec 2009
)

cat("ELI onset month_int values:\n")
print(onset_map)
cat("\n")

# ── 4) ASSIGN GNAME ────────────────────────────────────────────
# gname = ELI onset month_int for treated theatre; 0 = never treated
dt[, gname := 0L]
for (th in names(onset_map)) {
  n_matched <- dt[theatre_id == th, .N]
  if (n_matched > 0L) {
    dt[theatre_id == th, gname := onset_map[th]]
    log_msg("  Theatre '", th, "': gname = ", onset_map[th],
            "  (", n_matched, " cell-months)")
  } else {
    log_msg("  WARNING: theatre '", th, "' not found in panel.")
  }
}

cat("\nTreatment group distribution:\n")
print(dt[, .N, by = .(gname)][order(gname)])
cat("\n")

# ── 5) SAMPLE FOR SPEED ────────────────────────────────────────
# CS on 300k+ rows is slow. Keep ALL treated-theatre cells + a
# 3× sample of never-treated control cells (balanced across time).
treated_cells  <- unique(dt[gname > 0L, pg_id])
control_cells  <- unique(dt[gname == 0L, pg_id])

set.seed(42)
n_ctrl_sample <- min(length(control_cells), length(treated_cells) * 3L)
ctrl_sample   <- sample(control_cells, n_ctrl_sample, replace = FALSE)

dt_cs <- dt[pg_id %in% c(treated_cells, ctrl_sample)]

cat(sprintf("CS sample: %d treated cells + %d control cells = %d total cells\n",
            length(treated_cells), length(ctrl_sample),
            uniqueN(dt_cs$pg_id)))
cat(sprintf("Rows: %d\n\n", nrow(dt_cs)))

# ── 6) CONTROLS ────────────────────────────────────────────────
cs_ctrl_cand <- c("nlights_calib_mean", "nlights_mean",
                  "prec_gpcp", "temp", "water_ih")
cs_controls  <- intersect(cs_ctrl_cand, names(dt_cs))
xformla      <- as.formula(paste("~", paste(cs_controls, collapse = " + ")))
cat("CS controls:", paste(cs_controls, collapse = ", "), "\n\n")

# Check for any remaining NAs in outcome or controls
bad_cols <- c("log1p_events_type3", cs_controls)
for (cc in bad_cols) {
  n_na <- sum(is.na(dt_cs[[cc]]))
  if (n_na > 0L) {
    log_msg("  WARNING: ", n_na, " NAs in '", cc, "' — filling with column median.")
    med_val <- median(dt_cs[[cc]], na.rm = TRUE)
    dt_cs[is.na(get(cc)), (cc) := med_val]
  }
}

# CS requires a strictly balanced or explicitly unbalanced flag; use
# allow_unbalanced_panel = TRUE throughout.

# ── 7) RUN CS: NEVER-TREATED ───────────────────────────────────
log_msg("Running CS with never-treated control group (bstrap=TRUE, biters=999)...")
log_msg("  Expected runtime: 5-20 min depending on machine.")

cs_never <- tryCatch(
  att_gt(
    yname        = "log1p_events_type3",
    tname        = "month_int",
    idname       = "pg_id",
    gname        = "gname",
    xformla      = xformla,
    data         = dt_cs,
    control_group          = "nevertreated",
    anticipation           = 0L,
    panel                  = TRUE,
    allow_unbalanced_panel = TRUE,
    bstrap       = TRUE,
    biters       = 999L,
    clustervars  = "pg_id",   # cluster on cell, not theatre (8-cluster problem)
    print_details = FALSE
  ),
  error = function(e) { log_msg("ERROR (never-treated): ", e$message); NULL }
)

# ── 8) RUN CS: NOT-YET-TREATED ─────────────────────────────────
log_msg("Running CS with not-yet-treated control group...")

cs_notyet <- tryCatch(
  att_gt(
    yname        = "log1p_events_type3",
    tname        = "month_int",
    idname       = "pg_id",
    gname        = "gname",
    xformla      = xformla,
    data         = dt_cs,
    control_group          = "notyettreated",
    anticipation           = 0L,
    panel                  = TRUE,
    allow_unbalanced_panel = TRUE,
    bstrap       = TRUE,
    biters       = 999L,
    clustervars  = "pg_id",
    print_details = FALSE
  ),
  error = function(e) { log_msg("ERROR (not-yet-treated): ", e$message); NULL }
)

# ── 9) AGGREGATE ───────────────────────────────────────────────
cat("\n====== CALLAWAY-SANT'ANNA RESULTS ======\n\n")

agg_never  <- NULL
agg_notyet <- NULL

if (!is.null(cs_never)) {
  agg_never <- aggte(cs_never, type = "simple")
  cat("--- Never-treated control group ---\n")
  print(summary(agg_never))
  cat("\n")
  
  cat("--- ATT by onset cohort (never-treated control) ---\n")
  agg_group <- aggte(cs_never, type = "group")
  print(summary(agg_group))
  cat("\n")
}

if (!is.null(cs_notyet)) {
  agg_notyet <- aggte(cs_notyet, type = "simple")
  cat("--- Not-yet-treated control group ---\n")
  print(summary(agg_notyet))
  cat("\n")
}

# ── 10) COMPARE WITH TWFE ──────────────────────────────────────
cat("=========================================\n")
cat("Comparison with TWFE baseline:\n")
cat(sprintf("  TWFE (Conley 200km):   coef = 0.045  p = 0.026\n"))
if (!is.null(agg_never)) {
  cat(sprintf("  CS (never-treated):    ATT  = %.5f  SE = %.5f  p ≈ %.4f\n",
              agg_never$overall.att,
              agg_never$overall.se,
              2 * pnorm(-abs(agg_never$overall.att / agg_never$overall.se))))
}
if (!is.null(agg_notyet)) {
  cat(sprintf("  CS (not-yet-treated):  ATT  = %.5f  SE = %.5f  p ≈ %.4f\n",
              agg_notyet$overall.att,
              agg_notyet$overall.se,
              2 * pnorm(-abs(agg_notyet$overall.att / agg_notyet$overall.se))))
}
cat("=========================================\n\n")
cat("Interpretation:\n")
cat("  If CS ATT ≈ TWFE β (same sign, overlapping CI) → staggered-treatment\n")
cat("  heterogeneity bias is not driving the main result.\n")
cat("  If CS ATT diverges → heterogeneous treatment effects across cohorts;\n")
cat("  report CS as the preferred estimator in robustness table.\n\n")

# ── 11) SAVE ───────────────────────────────────────────────────
out_rds <- file.path(DIRS$out_paper1,
                     paste0("R09i_cs_models_", TS, ".rds"))
saveRDS(
  list(
    cs_never   = cs_never,
    cs_notyet  = cs_notyet,
    agg_never  = agg_never,
    agg_notyet = agg_notyet
  ),
  out_rds
)
log_msg("Saved models: ", out_rds)

# Write summary CSV for dissertation table
summary_rows <- list()
if (!is.null(agg_never)) {
  summary_rows[[1]] <- data.table(
    estimator    = "CS_never_treated",
    ATT          = round(agg_never$overall.att,  5),
    SE           = round(agg_never$overall.se,   5),
    p_approx     = round(2 * pnorm(-abs(agg_never$overall.att / agg_never$overall.se)), 4),
    control_group = "nevertreated",
    clustervars  = "pg_id",
    n_cells      = uniqueN(dt_cs$pg_id),
    n_obs        = nrow(dt_cs)
  )
}
if (!is.null(agg_notyet)) {
  summary_rows[[2]] <- data.table(
    estimator    = "CS_not_yet_treated",
    ATT          = round(agg_notyet$overall.att, 5),
    SE           = round(agg_notyet$overall.se,  5),
    p_approx     = round(2 * pnorm(-abs(agg_notyet$overall.att / agg_notyet$overall.se)), 4),
    control_group = "notyettreated",
    clustervars  = "pg_id",
    n_cells      = uniqueN(dt_cs$pg_id),
    n_obs        = nrow(dt_cs)
  )
}

if (length(summary_rows)) {
  summary_dt <- rbindlist(summary_rows)
  out_csv <- file.path(DIRS$audit,
                       paste0("R09i_cs_summary_", TS, ".csv"))
  fwrite(summary_dt, out_csv)
  log_msg("Wrote summary: ", out_csv)
  cat("\nSummary table:\n")
  print(summary_dt)
}

log_msg("Step 9i complete.")