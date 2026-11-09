# ============================================================================
# R27_FIXED: Mechanistic Tests with Proper Actor Variable Loading
# Purpose: Fixes actor variable loading issue
# ============================================================================

rm(list = ls())
suppressMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
})

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 27 FIXED: Comprehensive mechanistic validation tests")

OUT <- file.path(DIRS$out_paper1, "mechanism_tests")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════
# LOAD DATA — WITH ACTOR VARIABLES
# ══════════════════════════════════════════════════════════════

cat("\n=== LOADING DATA WITH ACTOR VARIABLES ===\n")
P_PANEL <- file.path(DIRS$panels, "panel_READY_FOR_R27.rds")
if (!file.exists(P_PANEL)) {
  stop("Run FIX_add_actors.R first to create panel_READY_FOR_R27.rds")
}
dt <- readRDS(P_PANEL)
setDT(dt)

# Standardize variable names
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

# Verify actor variables are present
actor_vars <- c("state_civ_harm_n", "rebel_civ_harm_n", "nonstate_civ_harm_n", "state_drone_n")
if (!all(actor_vars %in% names(dt))) {
  stop("Actor variables missing. Re-run FIX_add_actors.R")
}

# Create outcomes
dt[, log1p_type3 := log1p(events_type3)]
dt[, binary_violence := as.integer(events_type3 > 0)]

# Verify coordinates
if (!all(c("latitude", "longitude") %in% names(dt))) {
  stop("Coordinates missing. Re-run FIX_add_actors.R")
}

# Controls
CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min", "prec_gpcp", 
          "temp", "water_ih", "pop_gpw_sum", "forest_ih", "urban_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

FE <- "pg_id + theatre_id^month_int"

cat("Panel:", nrow(dt), "obs,", uniqueN(dt$pg_id), "cells,", 
    uniqueN(dt$theatre_id), "theatres\n")
cat("Controls:", ctrl_str, "\n")
cat("Actor variables available:", all(actor_vars %in% names(dt)), "\n")

# ══════════════════════════════════════════════════════════════
# TEST 1: INCIRLIK ACTOR DECOMPOSITION ⭐⭐⭐
# ══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("TEST 1: INCIRLIK SHOCK — ACTOR DECOMPOSITION\n")
cat(strrep("=", 70), "\n")

if (!all(actor_vars %in% names(dt))) {
  cat("WARNING: Missing actor variables. Skipping.\n")
  incirlik_actor_results <- NULL
} else {
  
  # Create log outcomes
  dt[, log1p_state_harm := log1p(state_civ_harm_n)]
  dt[, log1p_rebel_harm := log1p(rebel_civ_harm_n)]
  dt[, log1p_nonstate_harm := log1p(nonstate_civ_harm_n)]
  dt[, log1p_state_drone := log1p(state_drone_n)]
  
  # Check for Incirlik variables
  if (!all(c("shock_incirlik", "dist_incirlik_km") %in% names(dt))) {
    cat("WARNING: Missing Incirlik shock variables. Skipping.\n")
    incirlik_actor_results <- NULL
  } else {
    
    dt_incirlik <- dt[theatre_id %in% c("Syria", "Iraq")]
    cat("Incirlik subset: Syria + Iraq,", nrow(dt_incirlik), "obs\n")
    
    dt_incirlik[, dist_incirlik_norm := (dist_incirlik_km - mean(dist_incirlik_km, na.rm=TRUE)) / 
                  sd(dist_incirlik_km, na.rm=TRUE)]
    
    incirlik_models <- list()
    
    # Models for each outcome
    outcomes <- list(
      type3 = "log1p_type3",
      rebel = "log1p_rebel_harm",
      state = "log1p_state_harm",
      nonstate = "log1p_nonstate_harm",
      drone = "log1p_state_drone"
    )
    
    for (nm in names(outcomes)) {
      outcome_var <- outcomes[[nm]]
      
      incirlik_models[[nm]] <- feols(
        as.formula(paste0(outcome_var, " ~ shock_incirlik:dist_incirlik_norm:eli_feas + ",
                          "shock_incirlik:dist_incirlik_norm + shock_incirlik:eli_feas + ",
                          ctrl_str, " | pg_id + month_int")),
        data = dt_incirlik, notes = FALSE
      )
    }
    
    # Extract results
    incirlik_actor_results <- rbindlist(lapply(names(incirlik_models), function(nm) {
      m <- incirlik_models[[nm]]
      
      triple_term <- grep("shock_incirlik.*dist_incirlik_norm.*eli_feas", 
                          names(coef(m)), value = TRUE)
      
      if (length(triple_term) == 0) {
        return(data.table(outcome = nm, coefficient = NA, se_conley = NA, 
                          se_cell = NA, p_conley = NA, p_cell = NA))
      }
      
      coef_val <- coef(m)[triple_term]
      se_conley <- tryCatch(
        se(m, vcov = conley(200))[triple_term],
        error = function(e) NA_real_
      )
      se_cell <- se(m, cluster = ~pg_id)[triple_term]
      
      data.table(
        outcome = nm,
        coefficient = coef_val,
        se_conley = se_conley,
        se_cell = se_cell,
        p_conley = if (!is.na(se_conley)) 2 * pnorm(-abs(coef_val / se_conley)) else NA_real_,
        p_cell = 2 * pnorm(-abs(coef_val / se_cell)),
        n_obs = nobs(m)
      )
    }))
    
    incirlik_actor_results[, sig_conley := ifelse(is.na(p_conley), "",
                                                  ifelse(p_conley < 0.01, "***",
                                                         ifelse(p_conley < 0.05, "**",
                                                                ifelse(p_conley < 0.10, "*", ""))))]
    
    cat("\n=== INCIRLIK ACTOR DECOMPOSITION ===\n")
    print(incirlik_actor_results[, .(outcome, coefficient, se_conley, p_conley, sig_conley)])
    
    fwrite(incirlik_actor_results, file.path(OUT, "incirlik_actor_decomposition.csv"))
  }
}

# ══════════════════════════════════════════════════════════════
# TEST 2: DISTANCE GRADIENT BY ACTOR TYPE ⭐⭐⭐
# ══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("TEST 2: DISTANCE GRADIENT BY ACTOR TYPE\n")
cat(strrep("=", 70), "\n")

if (!all(actor_vars %in% names(dt))) {
  cat("WARNING: Missing actor variables. Skipping.\n")
  distance_gradient_results <- NULL
} else {
  
  distance_models <- list()
  
  distance_models$type3 <- feols(
    as.formula(paste0("log1p_type3 ~ post_eli:eli_feas + ", ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  distance_models$rebel <- feols(
    as.formula(paste0("log1p_rebel_harm ~ post_eli:eli_feas + ", ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  distance_models$state <- feols(
    as.formula(paste0("log1p_state_harm ~ post_eli:eli_feas + ", ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  distance_models$nonstate <- feols(
    as.formula(paste0("log1p_nonstate_harm ~ post_eli:eli_feas + ", ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  distance_models$drone <- feols(
    as.formula(paste0("log1p_state_drone ~ post_eli:eli_feas + ", ctrl_str, " | ", FE)),
    data = dt, notes = FALSE
  )
  
  distance_gradient_results <- rbindlist(lapply(names(distance_models), function(nm) {
    m <- distance_models[[nm]]
    
    coef_val <- coef(m)["post_eli:eli_feas"]
    se_conley <- tryCatch(
      se(m, vcov = conley(200))["post_eli:eli_feas"],
      error = function(e) NA_real_
    )
    se_cell <- se(m, cluster = ~pg_id)["post_eli:eli_feas"]
    
    data.table(
      outcome = nm,
      coefficient = coef_val,
      se_conley = se_conley,
      se_cell = se_cell,
      p_conley = if (!is.na(se_conley)) 2 * pnorm(-abs(coef_val / se_conley)) else NA_real_,
      p_cell = 2 * pnorm(-abs(coef_val / se_cell)),
      n_obs = nobs(m)
    )
  }))
  
  distance_gradient_results[, sig_conley := ifelse(is.na(p_conley), "",
                                                   ifelse(p_conley < 0.01, "***",
                                                          ifelse(p_conley < 0.05, "**",
                                                                 ifelse(p_conley < 0.10, "*", ""))))]
  
  cat("\n=== DISTANCE GRADIENT BY ACTOR ===\n")
  print(distance_gradient_results[, .(outcome, coefficient, se_conley, p_conley, sig_conley)])
  
  fwrite(distance_gradient_results, file.path(OUT, "distance_gradient_by_actor.csv"))
  
  # Create figure
  p_gradient <- ggplot(distance_gradient_results, 
                       aes(x = factor(outcome, levels = c("type3", "rebel", "state", 
                                                          "nonstate", "drone")),
                           y = coefficient)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = coefficient - 1.96 * se_conley,
                      ymax = coefficient + 1.96 * se_conley),
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      x = "",
      y = "Treatment Effect (β)",
      title = "ISR Effect by Actor Type",
      subtitle = "Rebel violence drives aggregate effect; state violence null"
    ) +
    scale_x_discrete(labels = c("type3" = "All (UCDP)",
                                "rebel" = "Rebel harm",
                                "state" = "State harm",
                                "nonstate" = "Non-state",
                                "drone" = "State drone")) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(OUT, "FIGURE_actor_gradient.png"),
         p_gradient, width = 8, height = 6, dpi = 300)
}

cat("\n=== R27 FIXED COMPLETE ===\n")
cat("Results saved to:", OUT, "\n")
cat("\nIf actor variables loaded successfully, you should see:\n")
cat("  - Incirlik actor decomposition results\n")
cat("  - Distance gradient by actor results\n")
cat("  - Actor gradient figure\n")