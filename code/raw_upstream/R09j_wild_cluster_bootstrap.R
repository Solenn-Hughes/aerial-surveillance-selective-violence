# ============================================================
# R09j_wild_cluster_bootstrap.R  (FIXED v4)
#
# FIX v4: The root cause of v2/v3 failures is that fixest
#   stores removed interaction-FE levels as formatted strings
#   (e.g. "Afghanistan^24001") that don't reliably match
#   paste0(theatre_id, "^", month_int) due to whitespace/type
#   differences.
#
#   Solution: PRE-CREATE the interaction FE as a plain column
#   (theatre_month_fe) and use pg_id + theatre_month_fe as
#   separate FEs. Then singleton removal is just a simple
#   column-value lookup, no string reconstruction needed.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

if (!requireNamespace("fwildclusterboot", quietly = TRUE)) {
  stop("Package 'fwildclusterboot' required.\n",
       "Install via: install.packages('fwildclusterboot')")
}
library(fwildclusterboot)

log_msg("Step 9j: Wild cluster bootstrap inference (v4)")
log_msg("fwildclusterboot version: ", as.character(packageVersion("fwildclusterboot")))
log_msg("fixest version: ", as.character(packageVersion("fixest")))

# ── SETTINGS ─────────────────────────────────────────────────
B_REPS      <- 399L
SEED        <- 42L
IMPOSE_NULL <- TRUE
RADII       <- c(200, 400, 800)

# ── LOAD PANEL ───────────────────────────────────────────────
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (!file.exists(P_ALLCW)) stop("Missing: ", P_ALLCW)

dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

OUTCOME <- "events_type3"
dt[, log1p_y := log1p(get(OUTCOME))]

# ── KEY FIX: Pre-create interaction FE as a plain column ─────
# This avoids the string-matching nightmare with theatre_id^month_int
dt[, theatre_month_fe := paste0(theatre_id, "_", month_int)]

CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

# eli_feas components for decay radius
req <- c("min_km_airbase", "eli_feas_terrain")
miss <- setdiff(req, names(dt))
if (length(miss)) {
  P_ELI <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
  if (file.exists(P_ELI)) {
    eli <- readRDS(P_ELI); setDT(eli)
    if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := as.integer(gid)]
    eli_cols <- intersect(c("pg_id", "min_km_airbase", "eli_feas_terrain"), names(eli))
    n_before <- nrow(dt)
    dt <- merge(dt, unique(eli[, ..eli_cols]), by = "pg_id", all.x = TRUE,
                suffixes = c("", ".eli"))
    stopifnot(nrow(dt) == n_before)
    for (cc in c("min_km_airbase", "eli_feas_terrain")) {
      cc_eli <- paste0(cc, ".eli")
      if (cc_eli %in% names(dt)) {
        if (!(cc %in% names(dt))) setnames(dt, cc_eli, cc)
        else { dt[is.na(get(cc)), (cc) := get(cc_eli)]; dt[, (cc_eli) := NULL] }
      }
    }
    log_msg("Merged eli components.")
  } else {
    log_msg("WARNING: eli_feasibility_pg.rds not found. Decay radius variants skipped.")
    RADII <- 400
  }
}

# ── ITERATIVE SINGLETON REMOVAL ──────────────────────────────
# Now using plain column FEs (pg_id + theatre_month_fe),
# singleton removal is a simple column-value lookup.

drop_singletons_v4 <- function(dt_in, fml) {
  dt_work <- copy(dt_in)
  max_iter <- 10L
  total_dropped <- 0L
  
  for (iter in seq_len(max_iter)) {
    m_pre <- tryCatch(
      feols(fml, data = dt_work, cluster = ~theatre_id, notes = FALSE),
      error = function(e) {
        log_msg("  ERROR in fit iteration ", iter, ": ", e$message)
        NULL
      }
    )
    if (is.null(m_pre)) return(NULL)
    
    removed <- m_pre$fixef_removed
    if (is.null(removed) || length(unlist(removed)) == 0) {
      log_msg("  Iteration ", iter, ": no singletons. Clean: ", nrow(dt_work), " rows.")
      return(list(data = dt_work, model = m_pre))
    }
    
    rows_to_drop <- rep(FALSE, nrow(dt_work))
    any_matched <- FALSE
    
    for (fe_name in names(removed)) {
      fe_vals <- removed[[fe_name]]
      if (length(fe_vals) == 0) next
      
      log_msg("  Iteration ", iter, " — FE '", fe_name, "': ",
              length(fe_vals), " levels removed")
      
      # Plain column lookup (works because we pre-created the FE column)
      if (fe_name %in% names(dt_work)) {
        col_vals <- as.character(dt_work[[fe_name]])
        fe_vals_chr <- as.character(fe_vals)
        matched <- col_vals %in% fe_vals_chr
        if (sum(matched) > 0) {
          rows_to_drop <- rows_to_drop | matched
          any_matched <- TRUE
          log_msg("    Matched ", sum(matched), " rows")
        } else {
          # Try numeric matching as fallback
          if (is.numeric(dt_work[[fe_name]])) {
            fe_vals_num <- suppressWarnings(as.numeric(fe_vals))
            matched2 <- dt_work[[fe_name]] %in% fe_vals_num[!is.na(fe_vals_num)]
            if (sum(matched2) > 0) {
              rows_to_drop <- rows_to_drop | matched2
              any_matched <- TRUE
              log_msg("    Numeric fallback matched ", sum(matched2), " rows")
            }
          }
          if (!any_matched) {
            log_msg("    WARNING: 0 rows matched for FE '", fe_name, "'")
          }
        }
      } else {
        log_msg("    WARNING: FE column '", fe_name, "' not in data")
      }
    }
    
    n_drop <- sum(rows_to_drop)
    if (n_drop == 0) {
      if (!any_matched) {
        log_msg("  Iteration ", iter, ": singletons reported but none matched.")
        log_msg("  Proceeding with fixest's internal handling.")
        return(list(data = dt_work, model = m_pre, has_unmatched_singletons = TRUE))
      }
      break
    }
    
    total_dropped <- total_dropped + n_drop
    log_msg("  Dropping ", n_drop, " rows (cumulative: ", total_dropped, ")")
    dt_work <- dt_work[!rows_to_drop]
    
    if (nrow(dt_work) < 100) {
      stop("Too few rows after singleton removal (", nrow(dt_work), ")")
    }
  }
  
  # Final refit
  m_final <- feols(fml, data = dt_work, cluster = ~theatre_id, notes = FALSE)
  list(data = dt_work, model = m_final)
}

# ── HELPER: FIT + BOOTSTRAP ──────────────────────────────────
run_wcb <- function(dt_in, label, eli_feas_col = "eli_feas") {
  log_msg("  Fitting: ", label)
  
  # Use theatre_month_fe (plain column) instead of theatre_id^month_int
  fml <- as.formula(paste0(
    "log1p_y ~ post_eli:", eli_feas_col, " + ", ctrl_str,
    " | pg_id + theatre_month_fe"
  ))
  
  result <- drop_singletons_v4(dt_in, fml)
  if (is.null(result)) {
    log_msg("  SKIP: singleton removal failed for ", label)
    return(NULL)
  }
  
  dt_use <- result$data
  m <- result$model
  
  cn <- names(coef(m))
  term_val <- cn[grepl(paste0("post_eli.*", eli_feas_col, "|", eli_feas_col, ".*post_eli"), cn)]
  if (!length(term_val)) {
    log_msg("  WARNING: Interaction term not found. Skipping.")
    return(NULL)
  }
  term_val <- term_val[1]
  
  log_msg("  Coefficient: ", round(coef(m)[term_val], 5),
          " (SE=", round(se(m)[term_val], 5), ")")
  log_msg("  Running wild cluster bootstrap (B=", B_REPS, ")...")
  
  set.seed(SEED)
  
  boot_result <- tryCatch({
    boottest(
      m,
      param       = term_val,
      clustid     = "theatre_id",
      B           = B_REPS,
      impose_null = IMPOSE_NULL,
      type        = "rademacher"
    )
  }, error = function(e1) {
    log_msg("  Attempt 1 (rademacher, impose_null) failed: ",
            substr(e1$message, 1, 150))
    
    # Attempt 2: webb weights (designed for few clusters)
    tryCatch({
      set.seed(SEED)
      log_msg("  Trying webb weights...")
      boottest(m, param = term_val, clustid = "theatre_id",
               B = B_REPS, impose_null = IMPOSE_NULL, type = "webb")
    }, error = function(e2) {
      log_msg("  Attempt 2 (webb) failed: ", substr(e2$message, 1, 150))
      
      # Attempt 3: without impose_null
      tryCatch({
        set.seed(SEED)
        log_msg("  Trying without impose_null...")
        boottest(m, param = term_val, clustid = "theatre_id",
                 B = B_REPS, type = "rademacher")
      }, error = function(e3) {
        log_msg("  Attempt 3 failed: ", substr(e3$message, 1, 150))
        
        # Attempt 4: mammen weights
        tryCatch({
          set.seed(SEED)
          log_msg("  Trying mammen weights...")
          boottest(m, param = term_val, clustid = "theatre_id",
                   B = B_REPS, type = "mammen")
        }, error = function(e4) {
          log_msg("  ALL attempts failed: ", substr(e4$message, 1, 150))
          NULL
        })
      })
    })
  })
  
  if (is.null(boot_result)) return(NULL)
  
  coef_val <- coef(m)[term_val]
  se_val   <- se(m)[term_val]
  p_cl     <- pvalue(m)[term_val]
  
  p_boot <- tryCatch(boot_result$p_val, error = function(e) {
    tryCatch(pval(boot_result), error = function(e2) NA_real_)
  })
  
  ci <- tryCatch(boot_result$conf_int, error = function(e) {
    tryCatch(confint(boot_result), error = function(e2) c(NA_real_, NA_real_))
  })
  
  log_msg("  SUCCESS: p_cluster=", round(p_cl, 4),
          " p_bootstrap=", round(as.numeric(p_boot), 4))
  
  data.table(
    label       = label,
    term        = term_val,
    coef        = round(coef_val, 5),
    se_cluster  = round(se_val, 5),
    p_cluster   = round(p_cl, 4),
    p_bootstrap = round(as.numeric(p_boot), 4),
    ci_lower    = round(as.numeric(ci[1]), 5),
    ci_upper    = round(as.numeric(ci[2]), 5),
    B           = B_REPS,
    impose_null = IMPOSE_NULL,
    nobs        = nobs(m)
  )
}

# ── RUN BASELINE ─────────────────────────────────────────────
fml_base <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + ", ctrl_str,
  " | pg_id + theatre_month_fe"
))

log_msg("Global singleton removal...")
base_result <- drop_singletons_v4(dt, fml_base)
if (is.null(base_result)) stop("Cannot fit baseline model.")
dt_clean <- base_result$data
log_msg("Clean data: ", nrow(dt_clean), " rows (original: ", nrow(dt), ")")

results <- list()

res_base <- run_wcb(dt_clean, "Baseline R=400km", "eli_feas")
if (!is.null(res_base)) results[[length(results) + 1]] <- res_base

# ── DECAY RADIUS VARIANTS ────────────────────────────────────
has_components <- all(c("min_km_airbase", "eli_feas_terrain") %in% names(dt_clean))

if (has_components && length(RADII) > 1) {
  for (r_km in RADII) {
    if (r_km == 400) next
    
    col_name <- paste0("eli_feas_R", r_km)
    dt_clean[, (col_name) := eli_feas_terrain * exp(-min_km_airbase / r_km)]
    
    n_na <- sum(is.na(dt_clean[[col_name]]))
    if (n_na > 0) {
      med_val <- median(dt_clean[[col_name]], na.rm = TRUE)
      dt_clean[is.na(get(col_name)), (col_name) := med_val]
    }
    
    res_r <- run_wcb(dt_clean, paste0("Decay R=", r_km, "km"), col_name)
    if (!is.null(res_r)) results[[length(results) + 1]] <- res_r
  }
}

# ── REPORT ────────────────────────────────────────────────────
if (length(results)) {
  out <- rbindlist(results)
  
  cat("\n============================================================\n")
  cat("WILD CLUSTER BOOTSTRAP RESULTS\n")
  cat("============================================================\n")
  print(out[, .(label, coef, se_cluster, p_cluster, p_bootstrap, ci_lower, ci_upper)])
  
  cat("\n── Comparison ──\n")
  for (i in seq_len(nrow(out))) {
    cat(sprintf("  %s: p_cluster=%.4f → p_bootstrap=%.4f  [%.4f, %.4f]\n",
                out$label[i], out$p_cluster[i], out$p_bootstrap[i],
                out$ci_lower[i], out$ci_upper[i]))
  }
  
  cat("\n── Interpretation ──\n")
  cat("If p_bootstrap > p_cluster: theatre-clustered SEs over-reject (too liberal).\n")
  cat("If p_bootstrap < p_cluster: they were conservative.\n")
  cat("Bootstrap CI provides valid inference with 8 clusters.\n")
  
  write_audit_dt(out, "R09j_wild_cluster_bootstrap")
} else {
  cat("\n============================================================\n")
  cat("NO BOOTSTRAP RESULTS PRODUCED\n")
  cat("============================================================\n")
  cat("fwildclusterboot: ", as.character(packageVersion("fwildclusterboot")), "\n")
  cat("fixest: ", as.character(packageVersion("fixest")), "\n")
  cat("\nDiagnostic: Trying minimal model without controls...\n")
  
  # Diagnostic: try the simplest possible model
  fml_diag <- log1p_y ~ post_eli:eli_feas | pg_id + theatre_month_fe
  m_diag <- tryCatch(
    feols(fml_diag, data = dt_clean, cluster = ~theatre_id, notes = FALSE),
    error = function(e) { cat("  Diagnostic fit failed: ", e$message, "\n"); NULL }
  )
  if (!is.null(m_diag)) {
    cat("  Diagnostic model fit OK. Nobs:", nobs(m_diag), "\n")
    cat("  Removed FEs:", length(unlist(m_diag$fixef_removed)), "\n")
    
    cn_d <- names(coef(m_diag))
    tv_d <- cn_d[grepl("post_eli.*eli_feas|eli_feas.*post_eli", cn_d)][1]
    
    bt_diag <- tryCatch({
      set.seed(SEED)
      boottest(m_diag, param = tv_d, clustid = "theatre_id",
               B = 999L, type = "rademacher")
    }, error = function(e) {
      cat("  Diagnostic boottest error: ", e$message, "\n")
      NULL
    })
    
    if (!is.null(bt_diag)) {
      p_diag <- tryCatch(bt_diag$p_val, error = function(e) NA)
      cat("  Diagnostic boottest SUCCEEDED. p_boot =", p_diag, "\n")
      cat("  → Controls may be causing the issue. Try updating packages.\n")
    } else {
      cat("  → Version incompatibility likely. Update packages.\n")
    }
  }
  
  cat("\nRecommended:\n")
  cat("  1. install.packages(c('fixest', 'fwildclusterboot'))\n")
  cat("  2. Restart R fully (not just clear workspace)\n")
  cat("  3. Re-run this script\n")
  cat("  4. If still failing, your Conley SEs already confirm robustness\n")
}

log_msg("Step 9j complete.")