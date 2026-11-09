# ============================================================
# R09k_ratio_outcome.R  (FIXED)
# Tests proposal DV (iii): share of selective vs indiscriminate violence
# FIX: removed summary(m, "compact") — not valid in this fixest version
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

try(setFixest_nthreads(getDTthreads()), silent = TRUE)
log_msg("Step 9k: Ratio outcome — selective/indiscriminate violence share")

# ── LOAD PANEL ───────────────────────────────────────────────
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (!file.exists(P_ALLCW)) stop("Missing: ", P_ALLCW)

dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

req <- c("pg_id", "month", "month_int", "post_eli", "eli_feas",
         "theatre_id", "events_type3", "events_ged")
miss <- setdiff(req, names(dt))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

# ── APPROACH A: CELL-MONTH RATIO ─────────────────────────────
log_msg("Approach A: Cell-month ratio (conditional on events_ged > 0)")

dt[, share_osv := fifelse(events_ged > 0, events_type3 / events_ged, NA_real_)]

n_with_events <- dt[events_ged > 0, .N]
n_total       <- nrow(dt)
log_msg("  Cell-months with events_ged > 0: ", n_with_events,
        " (", round(100 * n_with_events / n_total, 2), "%)")
log_msg("  Mean share_osv (conditional): ", round(mean(dt$share_osv, na.rm = TRUE), 4))

dt_pos <- dt[events_ged > 0]

CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt_pos))
ctrl_str <- paste(ctrl_present, collapse = " + ")

fml_a1 <- as.formula(paste0(
  "share_osv ~ post_eli:eli_feas + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))

m_a1 <- feols(fml_a1, data = dt_pos, cluster = ~theatre_id, notes = FALSE)

log_msg("Model A1 (cell-month ratio, conditional on events):")
# FIX: use print(summary()) instead of summary(m, "compact")
print(summary(m_a1))

ct_a1 <- as.data.table(coeftable(m_a1), keep.rownames = "term")
hit_a1 <- ct_a1[grepl("post_eli.*eli_feas|eli_feas.*post_eli", term)]

# ── APPROACH B: THEATRE-MONTH RATIO ──────────────────────────
log_msg("\nApproach B: Theatre-month ratio")

tm <- dt[, .(
  sum_type3 = sum(events_type3, na.rm = TRUE),
  sum_ged   = sum(events_ged, na.rm = TRUE),
  post_eli  = max(as.integer(post_eli), na.rm = TRUE)
), by = .(theatre_id, month, month_int)]

tm[, share_osv_tm := fifelse(sum_ged > 0, sum_type3 / sum_ged, NA_real_)]

log_msg("  Theatre-months with sum_ged > 0: ", tm[sum_ged > 0, .N], " / ", nrow(tm))
log_msg("  Mean share_osv_tm: ", round(mean(tm$share_osv_tm, na.rm = TRUE), 4))

tm_pos <- tm[sum_ged > 0]

fml_b1 <- as.formula("share_osv_tm ~ post_eli | theatre_id + month_int")
m_b1 <- feols(fml_b1, data = tm_pos, cluster = ~theatre_id,
              ssc = ssc(adj = TRUE, cluster.adj = TRUE, cluster.df = "min"),
              notes = FALSE)

log_msg("\nModel B1 (theatre-month ratio):")
print(summary(m_b1))

ct_b1 <- as.data.table(coeftable(m_b1), keep.rownames = "term")
hit_b1 <- ct_b1[grepl("post_eli", term)]

# ── RESULTS TABLE ─────────────────────────────────────────────
results <- rbindlist(list(
  data.table(
    approach  = "A: Cell-month (conditional)",
    outcome   = "share_osv = type3/ged",
    coef      = if (nrow(hit_a1)) round(hit_a1$Estimate[1], 5) else NA,
    se        = if (nrow(hit_a1)) round(hit_a1$`Std. Error`[1], 5) else NA,
    t         = if (nrow(hit_a1)) round(hit_a1$`t value`[1], 3) else NA,
    p         = if (nrow(hit_a1)) round(hit_a1$`Pr(>|t|)`[1], 4) else NA,
    n         = nobs(m_a1),
    fe        = "pg_id + theatre*month",
    vcov      = "theatre-clustered"
  ),
  data.table(
    approach  = "B: Theatre-month",
    outcome   = "share_osv_tm = Sum_type3/Sum_ged",
    coef      = if (nrow(hit_b1)) round(hit_b1$Estimate[1], 5) else NA,
    se        = if (nrow(hit_b1)) round(hit_b1$`Std. Error`[1], 5) else NA,
    t         = if (nrow(hit_b1)) round(hit_b1$`t value`[1], 3) else NA,
    p         = if (nrow(hit_b1)) round(hit_b1$`Pr(>|t|)`[1], 4) else NA,
    n         = nobs(m_b1),
    fe        = "theatre + month",
    vcov      = "theatre-clustered + ssc"
  )
))

cat("\n============================================\n")
cat("RATIO OUTCOME: Selective / Total Violence Share\n")
cat("============================================\n")
print(results)

cat("\n── Interpretation ──\n")
if (!is.na(results$coef[1])) {
  if (results$coef[1] > 0) {
    cat("Cell-level: ELI onset associated with HIGHER share of selective violence\n")
    cat("  in ISR-accessible cells. Consistent with H1 substitution story.\n")
  } else if (results$coef[1] < 0) {
    cat("Cell-level: ELI onset associated with LOWER share of selective violence\n")
    cat("  in ISR-accessible cells. Could indicate more battle events or\n")
    cat("  reclassification of OSV as battle-related in ISR-intense areas.\n")
  }
}

write_audit_dt(results, "R09k_ratio_outcome")
log_msg("Step 9k complete.")