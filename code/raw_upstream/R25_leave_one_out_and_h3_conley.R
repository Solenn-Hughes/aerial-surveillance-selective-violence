# ============================================================
# R25_leave_one_out_and_h3_conley.R
#
# PART 1: Leave-one-out — drop each theatre, re-estimate β.
#   Shows the result doesn't depend on any single theatre.
#
# PART 2: H3 rebel capability with Conley SEs.
#   The p=0.004 result used 8-cluster VCOV with warnings.
#   Need to confirm or kill it with proper inference.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
#         panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds (from R13)
# Output: LOO table + H3 Conley result
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 25: Leave-one-out + H3 Conley")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

# Coordinates for Conley
eli <- readRDS(file.path(DIRS$eli, "eli_feasibility_pg.rds")); setDT(eli)
coords <- eli[, .(pg_id = as.integer(gid), lon, lat)]
if ("lon" %in% names(dt)) dt[, c("lon","lat") := NULL]
dt <- merge(dt, coords, by = "pg_id", all.x = TRUE)

ctrl <- c("prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(ctrl, names(dt)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

theatres <- sort(unique(dt$theatre_id))

# ══════════════════════════════════════════════════════════════
# PART 1: LEAVE-ONE-OUT
# ══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("PART 1: LEAVE-ONE-OUT (drop each theatre, re-estimate)\n")
cat("============================================================\n")

loo_results <- list()

# Full-sample reference
f <- as.formula(paste0("log1p_y ~ post_eli:eli_feas + ", ctrl_str, " | ", FE))
m_full <- feols(f, data = dt, notes = FALSE)
cv_full <- vcov_conley(m_full, lat = "lat", lon = "lon", cutoff = 200)
ct_full <- as.data.table(coeftable(m_full, vcov = cv_full), keep.rownames = "term")
hit_full <- ct_full[grepl("eli_feas", term)]

loo_results[["Full sample (all 8)"]] <- data.table(
  dropped = "None (full)",
  coef = round(hit_full$Estimate[1], 5),
  se = round(hit_full$`Std. Error`[1], 5),
  p = round(hit_full$`Pr(>|t|)`[1], 4),
  nobs = nobs(m_full),
  n_theatres = 8L
)

for (th in theatres) {
  sub <- dt[theatre_id != th]
  n_th <- uniqueN(sub$theatre_id)
  
  m <- tryCatch(feols(f, data = sub, notes = FALSE), error = function(e) NULL)
  if (!is.null(m)) {
    cv <- tryCatch(vcov_conley(m, lat = "lat", lon = "lon", cutoff = 200),
                   error = function(e) NULL)
    if (!is.null(cv)) {
      ct <- as.data.table(coeftable(m, vcov = cv), keep.rownames = "term")
    } else {
      ct <- as.data.table(coeftable(m, vcov = ~pg_id), keep.rownames = "term")
    }
    hit <- ct[grepl("eli_feas", term)]
    if (nrow(hit) > 0) {
      loo_results[[paste0("Drop ", th)]] <- data.table(
        dropped = th,
        coef = round(hit$Estimate[1], 5),
        se = round(hit$`Std. Error`[1], 5),
        p = round(hit$`Pr(>|t|)`[1], 4),
        nobs = nobs(m),
        n_theatres = n_th
      )
    }
  }
}

loo_out <- rbindlist(loo_results)
cat("\nLeave-one-out results (Conley 200km):\n")
print(loo_out)

# Check if any theatre drives the result
if (nrow(loo_out) > 1) {
  signs_ok <- all(loo_out$coef > 0)
  cat(sprintf("\nAll LOO coefficients positive: %s\n", ifelse(signs_ok, "YES", "NO — investigate!")))
  cat(sprintf("Coefficient range: [%.5f, %.5f]\n", min(loo_out$coef), max(loo_out$coef)))
  
  # Which theatre matters most?
  full_beta <- loo_out[dropped == "None (full)", coef]
  loo_out[dropped != "None (full)", change_pct := round(100 * (coef - full_beta) / abs(full_beta), 1)]
  cat("\nChange from full sample when each theatre is dropped:\n")
  print(loo_out[dropped != "None (full)", .(dropped, coef, p, change_pct)])
}

write_audit_dt(loo_out, "R25_leave_one_out")


# ══════════════════════════════════════════════════════════════
# PART 2: H3 REBEL CAPABILITY WITH CONLEY SEs
# ══════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("PART 2: H3 rebel capability — Conley SEs\n")
cat("============================================================\n")

P_H2H3 <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds")

if (file.exists(P_H2H3)) {
  dt_h <- readRDS(P_H2H3); setDT(dt_h)
  if (!("pg_id" %in% names(dt_h)) && "gid" %in% names(dt_h)) dt_h[, pg_id := gid]
  dt_h[, pg_id := as.integer(pg_id)]
  if (!("month_int" %in% names(dt_h))) dt_h[, month_int := year(month) * 12L + month(month)]
  dt_h[, log1p_y := log1p(events_type3)]
  
  # Restrict to 8 theatres
  THEATRES_8 <- c("Afghanistan", "Pakistan", "Somalia", "Yemen (North Yemen)",
                  "Libya", "Iraq", "Syria", "Ukraine")
  dt_h <- dt_h[theatre_id %in% THEATRES_8]
  
  # Add coordinates
  if ("lon" %in% names(dt_h)) dt_h[, c("lon","lat") := NULL]
  dt_h <- merge(dt_h, coords, by = "pg_id", all.x = TRUE)
  
  if ("rebel_cap_tv" %in% names(dt_h)) {
    cat("rebel_cap_tv found. Running H3 with Conley...\n")
    
    f_h3 <- as.formula(paste0(
      "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:rebel_cap_tv + ",
      ctrl_str, " | ", FE
    ))
    m_h3 <- tryCatch(feols(f_h3, data = dt_h, notes = FALSE),
                     error = function(e) { log_msg("ERROR: ", e$message); NULL })
    
    if (!is.null(m_h3)) {
      # Conley
      cv_h3 <- tryCatch(vcov_conley(m_h3, lat = "lat", lon = "lon", cutoff = 200),
                        error = function(e) NULL)
      
      if (!is.null(cv_h3)) {
        cat("\nH3 rebel capability (Conley 200km):\n")
        print(coeftable(m_h3, vcov = cv_h3))
      }
      
      cat("\nH3 rebel capability (cell clustering):\n")
      m_h3_cell <- feols(f_h3, data = dt_h, cluster = ~pg_id, notes = FALSE)
      print(coeftable(m_h3_cell))
      
      # Also H2 for completeness
      if ("contestation_tv" %in% names(dt_h)) {
        f_h2 <- as.formula(paste0(
          "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:contestation_tv + ",
          ctrl_str, " | ", FE
        ))
        m_h2 <- tryCatch(feols(f_h2, data = dt_h, notes = FALSE), error = function(e) NULL)
        if (!is.null(m_h2)) {
          cv_h2 <- tryCatch(vcov_conley(m_h2, lat = "lat", lon = "lon", cutoff = 200),
                           error = function(e) NULL)
          if (!is.null(cv_h2)) {
            cat("\nH2 contestation (Conley 200km):\n")
            print(coeftable(m_h2, vcov = cv_h2))
          }
        }
      }
    }
  } else {
    cat("rebel_cap_tv not found in H2H3 panel. Run R13 first.\n")
  }
} else {
  cat("H2H3 panel not found: ", P_H2H3, "\n")
  cat("Run R13_build_h2_h3_moderators.R first.\n")
}

log_msg("Step 25 complete.")
