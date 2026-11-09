# ============================================================================
# R31_callaway_santanna.R — CS2021 on 3-theatre shock subset
#
# Purpose: Callaway & Sant'Anna (2021, Journal of Econometrics) is the
#   standard staggered-adoption-robust estimator for DiD designs. The
#   earlier attempt (R09i) failed on memory using the full 8-theatre
#   panel. Here we restrict to the 3 shock theatres, which both:
#     (a) fits in memory (~3× smaller sample)
#     (b) matches our causal design (shocks identify the effect)
#
# Design:
#   — Subset: Afghanistan, Syria, Iraq
#   — Treatment: high-eli_feas cells (above theatre-median)
#   — Never-treated control: low-eli_feas cells in the same theatre
#   — Group variable: first month treated (i.e. shock date, by theatre)
#   — Estimator: ATT(g,t) aggregated to:
#       (a) event-study aggregation
#       (b) overall ATT
#
# Install:
#   # install.packages("did")
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})

if (!requireNamespace("did", quietly = TRUE)) install.packages("did")
library(did)

source("R00_config.R")
log_msg("R31: Callaway-Sant'Anna on 3-theatre subset")

OUT <- file.path("outputs", "paper1", "callaway_santanna")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ── 1) LOAD PANEL ──────────────────────────────────────────────────────
panel_path <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(panel_path); setDT(dt)
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_type3 := log1p(events_type3)]

# ── 2) RESTRICT TO SHOCK THEATRES ──────────────────────────────────────
SHOCK_THEATRES <- c("Afghanistan", "Syria", "Iraq")
dt_s <- dt[theatre_id %in% SHOCK_THEATRES]
log_msg("Shock subset: ", nrow(dt_s), " obs, ",
        uniqueN(dt_s$pg_id), " cells, ",
        length(SHOCK_THEATRES), " theatres")

# ── 3) DEFINE BINARY TREATMENT & GROUP ─────────────────────────────────
# CS2021 requires binary treatment assignment.
# Treated = cells with above-theatre-median eli_feas.
# Group variable G = first treatment month (shock month) for treated cells,
#                    0 for never-treated.

dt_s[, theatre_median_feas := median(eli_feas, na.rm = TRUE),
     by = theatre_id]
dt_s[, treated_cell := eli_feas > theatre_median_feas]
log_msg("Treated cells: ",
        length(unique(dt_s[treated_cell == TRUE, pg_id])),
        " / ", length(unique(dt_s$pg_id)))

# Group = first treated month_int for treated cells (theatre-specific shock)
shock_month <- c(
  "Afghanistan" = as.integer(as.IDate("2021-07-01")) %/% 30 + 1L,
  "Syria"       = as.integer(as.IDate("2019-10-01")) %/% 30 + 1L,
  "Iraq"        = as.integer(as.IDate("2019-10-01")) %/% 30 + 1L
)
# Use month_int directly for CS2021 (year*12 + month)
shock_month_int <- c(
  "Afghanistan" = 2021L * 12L + 7L,
  "Syria"       = 2019L * 12L + 10L,
  "Iraq"        = 2019L * 12L + 10L
)

dt_s[, G := 0L]
for (th in names(shock_month_int)) {
  dt_s[theatre_id == th & treated_cell == TRUE,
       G := shock_month_int[[th]]]
}

# Filter: keep only observations near the shock windows to manage size
# Window: 24 months before, 24 months after the shock
dt_s[, shock_m := shock_month_int[theatre_id]]
dt_s[, months_from_shock := month_int - shock_m]
dt_s <- dt_s[months_from_shock >= -24L & months_from_shock <= 24L]
log_msg("Windowed subset: ", nrow(dt_s), " obs")

# ── 4) RUN CS2021 ──────────────────────────────────────────────────────
# att_gt expects:
#   yname       : outcome
#   tname       : time variable
#   idname      : unit id
#   gname       : group (first treated period); 0 = never-treated
#   control_group: "nevertreated" (low-eli_feas cells in same theatres)
#   xformla     : covariates (pre-treatment)
#   panel       : TRUE (balanced panel)
#   allow_unbalanced_panel: TRUE if cells enter/exit
#   bstrap      : TRUE
#   biters      : 200 (enough for top-journal, not too slow)
#   clustervars : cluster SEs at cell level

# Ensure unique cell-month observations
dt_s <- unique(dt_s, by = c("pg_id", "month_int"))

cat("\n── Running att_gt() — may take 3-8 minutes ──\n")
att_res <- tryCatch({
  att_gt(
    yname         = "log1p_type3",
    tname         = "month_int",
    idname        = "pg_id",
    gname         = "G",
    data          = dt_s,
    xformla       = ~ eli_feas,
    control_group = "nevertreated",
    bstrap        = TRUE,
    biters        = 200,
    cband         = TRUE,
    clustervars   = "pg_id",
    allow_unbalanced_panel = TRUE
  )
}, error = function(e) {
  message("att_gt failed: ", e$message)
  NULL
})

if (is.null(att_res)) {
  cat("\nATT(g,t) estimation failed — likely memory or panel imbalance.\n")
  cat("Fallback: run without covariates.\n")
  att_res <- tryCatch({
    att_gt(
      yname         = "log1p_type3",
      tname         = "month_int",
      idname        = "pg_id",
      gname         = "G",
      data          = dt_s,
      control_group = "nevertreated",
      bstrap        = TRUE,
      biters        = 100,
      clustervars   = "pg_id",
      allow_unbalanced_panel = TRUE
    )
  }, error = function(e) {
    message("Fallback also failed: ", e$message); NULL
  })
}

if (is.null(att_res)) stop("CS2021 could not be estimated even in reduced spec.")

# ── 5) AGGREGATIONS ────────────────────────────────────────────────────

cat("\n── Event-study aggregation (dynamic effects) ──\n")
agg_es <- tryCatch({
  aggte(att_res, type = "dynamic", na.rm = TRUE)
}, error = function(e) { message("Event-study agg failed: ", e$message); NULL })
if (!is.null(agg_es)) {
  print(agg_es)
  es_dt <- data.table(
    event_time = agg_es$egt,
    ATT        = agg_es$att.egt,
    SE         = agg_es$se.egt
  )
  es_dt[, p_value := 2 * pnorm(-abs(ATT / SE))]
  fwrite(es_dt, file.path(OUT, "R31_CS2021_event_study.csv"))
}

cat("\n── Overall ATT aggregation ──\n")
agg_overall <- tryCatch({
  aggte(att_res, type = "simple", na.rm = TRUE)
}, error = function(e) { message("Overall agg failed: ", e$message); NULL })
if (!is.null(agg_overall)) {
  print(agg_overall)
  cat(sprintf("\nCS2021 OVERALL ATT = %.4f  SE = %.4f\n",
              agg_overall$overall.att, agg_overall$overall.se))
}

cat("\n── Group-specific ATT (by shock theatre) ──\n")
agg_group <- tryCatch({
  aggte(att_res, type = "group", na.rm = TRUE)
}, error = function(e) { message("Group agg failed: ", e$message); NULL })
if (!is.null(agg_group)) {
  print(agg_group)
  gr_dt <- data.table(
    group       = agg_group$egt,
    ATT         = agg_group$att.egt,
    SE          = agg_group$se.egt
  )
  fwrite(gr_dt, file.path(OUT, "R31_CS2021_by_group.csv"))
}

# ── 6) COMPARISON WITH TWFE ESTIMATE ───────────────────────────────────
cat("\n── COMPARISON ──\n")
cat("Your TWFE estimate on the 8-theatre panel: β(post × eli_feas) ≈ 0.045\n")
cat("Your stacked shock-reversal estimate       : β ≈ -0.037\n")
cat("CS2021 overall ATT (binary treatment)      : see above\n")
cat("\nCS2021 differs conceptually because:\n")
cat("  — It uses binary treatment (high- vs low-feas), not continuous\n")
cat("  — It avoids negative-weighting problems in staggered TWFE\n")
cat("  — It's robust to heterogeneous treatment effects across cohorts\n")

log_msg("R31 complete. Outputs in: ", OUT)
