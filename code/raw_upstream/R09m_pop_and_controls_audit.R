# ============================================================
# R09m_pop_and_controls_audit.R
# 1) Audits what variables are already in the panel
# 2) Adds population (pop_gpw_sum) as control
# 3) Shows coefficient stability
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 9m: Population robustness + controls audit")

# ── LOAD PANEL ───────────────────────────────────────────────
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (!file.exists(P_ALLCW)) stop("Missing: ", P_ALLCW)

dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

# ══════════════════════════════════════════════════════════════
# PART 1: AUDIT ALL VARIABLES IN PANEL
# ══════════════════════════════════════════════════════════════
cat("\n============================================\n")
cat("PANEL VARIABLE AUDIT\n")
cat("============================================\n")
cat("Total columns:", ncol(dt), "\n")
cat("Total rows:", nrow(dt), "\n\n")

# Identify all numeric columns
num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]

audit <- rbindlist(lapply(num_cols, function(cc) {
  x <- dt[[cc]]
  data.table(
    variable = cc,
    class    = class(x)[1],
    n_na     = sum(is.na(x)),
    pct_na   = round(100 * mean(is.na(x)), 2),
    min      = round(min(x, na.rm = TRUE), 4),
    median   = round(median(x, na.rm = TRUE), 4),
    max      = round(max(x, na.rm = TRUE), 4),
    sd       = round(sd(x, na.rm = TRUE), 6),
    n_unique = uniqueN(x, na.rm = TRUE)
  )
}))

# Flag time-varying vs absorbed-by-FE
# A variable is "time-varying within cell" if it has >1 unique value per pg_id
tv_check <- dt[, lapply(.SD, uniqueN), by = pg_id, .SDcols = intersect(num_cols, names(dt))]
tv_max <- tv_check[, lapply(.SD, max), .SDcols = setdiff(names(tv_check), "pg_id")]
tv_vec <- unlist(tv_max)

audit[, time_varying := fifelse(variable %in% names(tv_vec) & tv_vec[variable] > 1,
                                 "YES", "no (absorbed by cell FE)")]

cat("── Numeric variables in panel ──\n")
print(audit[order(-time_varying, variable)], nrows = 100)

# ── POPULATION CHECK ─────────────────────────────────────────
pop_candidates <- grep("pop|gpw|worldpop|popd", names(dt), value = TRUE, ignore.case = TRUE)
cat("\n── Population-related columns found: ──\n")
if (length(pop_candidates)) {
  cat(paste(pop_candidates, collapse = ", "), "\n")
} else {
  cat("NONE found in panel.\n")
  cat("Checking controls_pg_month.rds...\n")
  
  P_CTRL <- file.path(DIRS$controls %||% file.path("data", "processed", "base"),
                       "controls_pg_month.rds")
  if (!file.exists(P_CTRL)) P_CTRL <- file.path("data", "processed", "base", "controls_pg_month.rds")
  
  if (file.exists(P_CTRL)) {
    ctrl <- readRDS(P_CTRL); setDT(ctrl)
    pop_in_ctrl <- grep("pop|gpw|worldpop|popd", names(ctrl), value = TRUE, ignore.case = TRUE)
    cat("In controls_pg_month.rds:", paste(pop_in_ctrl, collapse = ", "), "\n")
    
    if (length(pop_in_ctrl)) {
      cat("Merging population into panel...\n")
      if (!("pg_id" %in% names(ctrl)) && "gid" %in% names(ctrl)) ctrl[, pg_id := gid]
      ctrl[, pg_id := as.integer(pg_id)]
      
      mcol <- intersect(c("month", "ym", "date"), names(ctrl))[1]
      if (!is.na(mcol)) {
        if (mcol != "month") setnames(ctrl, mcol, "month")
        if (exists("as_month_idate")) ctrl[, month := as_month_idate(month)] else ctrl[, month := as.IDate(month)]
        
        pop_cols <- c("pg_id", "month", pop_in_ctrl)
        pop_dt <- ctrl[, ..pop_cols]
        dt <- merge(dt, pop_dt, by = c("pg_id", "month"), all.x = TRUE)
        pop_candidates <- pop_in_ctrl
        cat("Merged successfully.\n")
      }
    }
  }
  
  # Also try priogrid_yearly directly
  P_YEARLY <- file.path(DIRS$base %||% file.path("data", "processed", "base"),
                         "priogrid_yearly.rds")
  if (!file.exists(P_YEARLY)) P_YEARLY <- file.path("data", "processed", "base", "priogrid_yearly.rds")
  
  if (file.exists(P_YEARLY) && !length(pop_candidates)) {
    yearly <- readRDS(P_YEARLY); setDT(yearly)
    pop_in_yearly <- grep("pop|gpw|worldpop", names(yearly), value = TRUE, ignore.case = TRUE)
    cat("In priogrid_yearly.rds:", paste(pop_in_yearly, collapse = ", "), "\n")
    
    if (length(pop_in_yearly)) {
      if (!("pg_id" %in% names(yearly)) && "gid" %in% names(yearly)) yearly[, pg_id := as.integer(gid)]
      yearly[, year := as.integer(year)]
      dt[, year := as.integer(year(month))]
      
      pop_cols <- c("pg_id", "year", pop_in_yearly)
      pop_dt <- unique(yearly[, ..pop_cols])
      dt <- merge(dt, pop_dt, by = c("pg_id", "year"), all.x = TRUE)
      pop_candidates <- pop_in_yearly
      cat("Merged from priogrid_yearly.\n")
    }
  }
}

# ══════════════════════════════════════════════════════════════
# PART 2: POPULATION ROBUSTNESS
# ══════════════════════════════════════════════════════════════

# m0_boot baseline controls
CTRL_BASE <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
               "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL_BASE, names(dt))

# M0: Baseline (no population)
fml_base <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + ", paste(ctrl_present, collapse = " + "),
  " | pg_id + theatre_id^month_int"
))
m0 <- feols(fml_base, data = dt, cluster = ~theatre_id, notes = FALSE)

cat("\n============================================\n")
cat("M0: BASELINE (m0_boot, no population)\n")
cat("============================================\n")
ct0 <- as.data.table(coeftable(m0), keep.rownames = "term")
hit0 <- ct0[grepl("post_eli.*eli_feas|eli_feas.*post_eli", term)]
if (nrow(hit0)) cat(sprintf("  Coef: %.5f  SE: %.5f  p: %.4f\n",
  hit0$Estimate[1], hit0$`Std. Error`[1], hit0$`Pr(>|t|)`[1]))

# Now try adding population
pop_col <- NULL
for (pc in c("pop_gpw_sum", "pop_gpw_mean", "pop_total", "pop_gpw",
             pop_candidates)) {
  if (pc %in% names(dt) && is.numeric(dt[[pc]])) {
    # Check it actually varies within cells
    test_var <- dt[, .(uq = uniqueN(get(pc))), by = pg_id][uq > 1, .N]
    if (test_var > 0) {
      pop_col <- pc
      break
    }
  }
}

results <- data.table(
  model = "M0: Baseline (6 controls)",
  controls = paste(ctrl_present, collapse = ", "),
  coef = if (nrow(hit0)) round(hit0$Estimate[1], 5) else NA,
  se   = if (nrow(hit0)) round(hit0$`Std. Error`[1], 5) else NA,
  p    = if (nrow(hit0)) round(hit0$`Pr(>|t|)`[1], 4) else NA,
  n    = nobs(m0)
)

if (!is.null(pop_col)) {
  cat("\n── Population variable found:", pop_col, "──\n")
  
  # Log-transform population (right-skewed, many zeros possible)
  dt[, log1p_pop := log1p(get(pop_col))]
  
  cat("  Median:", round(median(dt[[pop_col]], na.rm = TRUE), 2), "\n")
  cat("  Mean:", round(mean(dt[[pop_col]], na.rm = TRUE), 2), "\n")
  cat("  NA:", sum(is.na(dt[[pop_col]])), "\n")
  
  # M1: Add log population
  ctrl_with_pop <- c(ctrl_present, "log1p_pop")
  fml_pop <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + ", paste(ctrl_with_pop, collapse = " + "),
    " | pg_id + theatre_id^month_int"
  ))
  m1 <- feols(fml_pop, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  cat("\n============================================\n")
  cat("M1: + log(1 + population)\n")
  cat("============================================\n")
  ct1 <- as.data.table(coeftable(m1), keep.rownames = "term")
  hit1 <- ct1[grepl("post_eli.*eli_feas|eli_feas.*post_eli", term)]
  if (nrow(hit1)) cat(sprintf("  Coef: %.5f  SE: %.5f  p: %.4f\n",
    hit1$Estimate[1], hit1$`Std. Error`[1], hit1$`Pr(>|t|)`[1]))
  
  results <- rbindlist(list(results, data.table(
    model = "M1: + log(1+pop)",
    controls = paste(ctrl_with_pop, collapse = ", "),
    coef = if (nrow(hit1)) round(hit1$Estimate[1], 5) else NA,
    se   = if (nrow(hit1)) round(hit1$`Std. Error`[1], 5) else NA,
    p    = if (nrow(hit1)) round(hit1$`Pr(>|t|)`[1], 4) else NA,
    n    = nobs(m1)
  )), fill = TRUE)
  
  # M2: Add population + all other available time-varying controls
  # Find all time-varying numeric columns not already in controls
  potential_extra <- setdiff(num_cols, c(ctrl_with_pop, "pg_id", "month", "month_int",
    "year", "post_eli", "eli_feas", "treated", "event_time", "eli_onset_month",
    "log1p_y", "log1p_pop", "events_type3", "events_ged", "fat_best_type3",
    "civ_best_type3", "eli_feas_terrain", "eli_feas_range", "min_km_airbase",
    "theatre_id", "country_id"))
  
  # Keep only those that vary within cells
  tv_extras <- character()
  for (v in potential_extra) {
    if (is.numeric(dt[[v]]) && sum(is.na(dt[[v]])) < nrow(dt) * 0.5) {
      uq <- dt[, uniqueN(get(v)), by = pg_id][V1 > 1, .N]
      if (uq > 100) tv_extras <- c(tv_extras, v)
    }
  }
  
  if (length(tv_extras)) {
    cat("\n── Additional time-varying controls available: ──\n")
    cat(paste(tv_extras, collapse = ", "), "\n")
    
    ctrl_full <- c(ctrl_with_pop, tv_extras)
    fml_full <- as.formula(paste0(
      "log1p_y ~ post_eli:eli_feas + ", paste(ctrl_full, collapse = " + "),
      " | pg_id + theatre_id^month_int"
    ))
    
    m2 <- tryCatch(
      feols(fml_full, data = dt, cluster = ~theatre_id, notes = FALSE),
      error = function(e) { cat("M2 failed:", e$message, "\n"); NULL }
    )
    
    if (!is.null(m2)) {
      ct2 <- as.data.table(coeftable(m2), keep.rownames = "term")
      hit2 <- ct2[grepl("post_eli.*eli_feas|eli_feas.*post_eli", term)]
      
      cat("\n============================================\n")
      cat("M2: + all available time-varying controls\n")
      cat("============================================\n")
      if (nrow(hit2)) cat(sprintf("  Coef: %.5f  SE: %.5f  p: %.4f\n",
        hit2$Estimate[1], hit2$`Std. Error`[1], hit2$`Pr(>|t|)`[1]))
      
      results <- rbindlist(list(results, data.table(
        model = paste0("M2: + all TV controls (", length(ctrl_full), " vars)"),
        controls = paste(ctrl_full, collapse = ", "),
        coef = if (nrow(hit2)) round(hit2$Estimate[1], 5) else NA,
        se   = if (nrow(hit2)) round(hit2$`Std. Error`[1], 5) else NA,
        p    = if (nrow(hit2)) round(hit2$`Pr(>|t|)`[1], 4) else NA,
        n    = nobs(m2)
      )), fill = TRUE)
    }
  }
  
} else {
  cat("\n⚠ No time-varying population variable found in panel or controls.\n")
  cat("PRIO-GRID yearly should have pop_gpw_sum. Check your priogrid_yearly.rds.\n")
  cat("If it's missing, download WorldPop or use GPWv4 (see download instructions).\n")
}

# ══════════════════════════════════════════════════════════════
# PART 3: RESULTS COMPARISON
# ══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("CONTROL ROBUSTNESS: Coefficient stability\n")
cat("============================================================\n")
print(results[, .(model, coef, se, p, n)])

if (nrow(results) > 1) {
  base_coef <- results$coef[1]
  for (i in 2:nrow(results)) {
    pct_change <- round(100 * (results$coef[i] - base_coef) / abs(base_coef), 1)
    cat(sprintf("  %s → %s: coefficient changes by %+.1f%%\n",
                results$model[1], results$model[i], pct_change))
  }
  cat("\nIf changes are < 10%, coefficient is robust to control set.\n")
  cat("If changes are > 20%, investigate which control is absorbing variation.\n")
}

write_audit_dt(results, "R09m_pop_robustness")
write_audit_dt(audit, "R09m_panel_variable_audit")

log_msg("Step 9m complete.")
