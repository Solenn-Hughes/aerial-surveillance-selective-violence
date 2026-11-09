# ============================================================
# R09p_negative_weights.R
# de Chaisemartin & D'Haultfoeuille (2020) TWFE weight diagnostic
#
# With staggered treatment across 8 theatres (2001–2022),
# standard TWFE can assign negative weights to some group-time
# ATTs, biasing the pooled estimate. This script:
#   (a) Computes the share of negative weights
#   (b) Reports the distribution of weights
#   (c) If negative share is large (>10%), flags concern
#
# Uses the DIDmultiplegt package (or twowayfeweights).
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
# Output: audit table with weight diagnostics
#
# Run AFTER R16.
# ============================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
})

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 9p: de Chaisemartin & D'Haultfoeuille negative weights diagnostic")

# ── CHECK PACKAGE ────────────────────────────────────────────
# Try TwoWayFEWeights first, then DIDmultiplegt
pkg_available <- FALSE

if (requireNamespace("TwoWayFEWeights", quietly = TRUE)) {
  library(TwoWayFEWeights)
  pkg_available <- TRUE
  pkg_used <- "TwoWayFEWeights"
} else if (requireNamespace("DIDmultiplegt", quietly = TRUE)) {
  library(DIDmultiplegt)
  pkg_available <- TRUE
  pkg_used <- "DIDmultiplegt"
} else {
  cat("\n============================================================\n")
  cat("INSTALL ONE OF THESE PACKAGES:\n")
  cat("  install.packages('TwoWayFEWeights')\n")
  cat("  OR\n")
  cat("  install.packages('DIDmultiplegt')\n")
  cat("============================================================\n")
  cat("\nSkipping negative weight diagnostic.\n")
  cat("This is a robustness check, not blocking.\n")
}

if (pkg_available) {
  
  # ── LOAD DATA ──────────────────────────────────────────────
  P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
  if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL)
  
  dt <- readRDS(P_PANEL); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  dt[, pg_id := as.integer(pg_id)]
  if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
  dt[, log1p_y := log1p(events_type3)]
  
  log_msg("Panel loaded: ", nrow(dt), " rows")
  
  if (pkg_used == "TwoWayFEWeights") {
    # ── TwoWayFEWeights approach ───────────────────────────────
    # This directly computes the weights from a TWFE regression
    
    cat("\nUsing TwoWayFEWeights package\n")
    
    # Need a binary treatment for the weight decomposition
    # Use post_eli as the treatment indicator
    # The function needs: outcome, group, time, treatment
    
    tw <- tryCatch({
      twowayfeweights(
        df = as.data.frame(dt),
        Y = "log1p_y",
        G = "pg_id",
        T = "month_int",
        D = "post_eli",
        cmd_type = "feTR"
      )
    }, error = function(e) {
      log_msg("TwoWayFEWeights error: ", e$message)
      NULL
    })
    
    if (!is.null(tw)) {
      cat("\n============================================================\n")
      cat("TWFE WEIGHT DECOMPOSITION (de Chaisemartin & D'Haultfoeuille)\n")
      cat("============================================================\n")
      print(summary(tw))
      
      # Extract key statistics
      if (is.data.frame(tw) || is.data.table(tw)) {
        setDT(tw)
        n_neg <- sum(tw$weight < 0, na.rm = TRUE)
        n_pos <- sum(tw$weight > 0, na.rm = TRUE)
        n_total <- n_neg + n_pos
        share_neg <- n_neg / n_total
        sum_neg <- sum(tw$weight[tw$weight < 0], na.rm = TRUE)
        
        diag_out <- data.table(
          metric = c("N positive weights", "N negative weights", "Total",
                     "Share negative", "Sum negative weights"),
          value = c(n_pos, n_neg, n_total,
                    round(share_neg, 4), round(sum_neg, 4))
        )
        
        cat("\nDiagnostic summary:\n")
        print(diag_out)
        write_audit_dt(diag_out, "R09p_negative_weights")
        
        if (share_neg > 0.10) {
          cat("\n\u26A0\uFE0F  WARNING: >10% negative weights. TWFE may be biased.\n")
          cat("Rely on Callaway & Sant'Anna (R09i) as primary estimator.\n")
        } else {
          cat("\n\u2713  Negative weight share is small (<10%). TWFE likely unbiased.\n")
        }
      }
    }
    
  } else {
    # ── DIDmultiplegt approach ─────────────────────────────────
    cat("\nUsing DIDmultiplegt package\n")
    cat("Running did_multiplegt... (this may take a while)\n")
    
    # Sample for speed if panel is very large
    sample_cells <- unique(dt$pg_id)
    if (length(sample_cells) > 5000) {
      set.seed(42)
      sample_cells <- sample(sample_cells, 5000)
      dt_sub <- dt[pg_id %in% sample_cells]
      log_msg("Sampled ", length(sample_cells), " cells for speed")
    } else {
      dt_sub <- dt
    }
    
    dm <- tryCatch({
      did_multiplegt(
        df = as.data.frame(dt_sub),
        Y = "log1p_y",
        G = "pg_id",
        T = "month_int",
        D = "post_eli",
        brep = 0  # no bootstrap for speed; just weights
      )
    }, error = function(e) {
      log_msg("DIDmultiplegt error: ", e$message)
      NULL
    })
    
    if (!is.null(dm)) {
      cat("\n============================================================\n")
      cat("DIDmultiplegt RESULTS\n")
      cat("============================================================\n")
      print(dm)
      
      # Save whatever output we get
      if (is.list(dm)) {
        diag_out <- data.table(
          metric = names(dm),
          value = as.character(unlist(dm))
        )
        write_audit_dt(diag_out, "R09p_negative_weights_didm")
      }
    }
  }
}

log_msg("Step 9p complete.")
