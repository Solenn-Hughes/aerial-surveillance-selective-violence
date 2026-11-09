# ============================================================================
# R27_PREP: Data Preparation for Mechanistic Tests
# Purpose: Ensures all required variables are in the main panel
# Run this ONCE before R27 and R28
# ============================================================================

rm(list = ls())
suppressMessages({ library(data.table) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
source(cfg_path)

cat("\n", strrep("=", 70), "\n")
cat("PREPARING DATA FOR MECHANISTIC TESTS\n")
cat(strrep("=", 70), "\n\n")

# ══════════════════════════════════════════════════════════════
# STEP 1: LOAD MAIN PANEL
# ══════════════════════════════════════════════════════════════

P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
cat("Loading main panel:", P_PANEL, "\n")

if (!file.exists(P_PANEL)) {
  stop("Main panel not found. Run pipeline through R16 first.")
}

dt <- readRDS(P_PANEL)
setDT(dt)

if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]

cat("Panel loaded:", nrow(dt), "obs,", uniqueN(dt$pg_id), "cells\n\n")

# ══════════════════════════════════════════════════════════════
# STEP 2: CHECK AND ADD ACTOR VARIABLES
# ══════════════════════════════════════════════════════════════

cat("Checking actor variables...\n")

actor_vars <- c("state_civ_harm_n", "rebel_civ_harm_n", "nonstate_civ_harm_n", "state_drone_n")
missing_actors <- setdiff(actor_vars, names(dt))

if (length(missing_actors) > 0) {
  cat("Missing from panel:", paste(missing_actors, collapse = ", "), "\n")
  cat("Attempting to load from acled_pg_month.rds...\n")
  
  P_ACLED <- file.path(DIRS$base, "acled_pg_month.rds")
  
  if (file.exists(P_ACLED)) {
    acled <- readRDS(P_ACLED)
    setDT(acled)
    acled[, pg_id := as.integer(pg_id)]
    
    still_missing <- setdiff(actor_vars, names(acled))
    
    if (length(still_missing) > 0) {
      cat("\nERROR: Actor variables also missing from ACLED file.\n")
      cat("Missing:", paste(still_missing, collapse = ", "), "\n")
      cat("\nSOLUTION: Run R23_audit_and_fix.R first to create corrected actor variables.\n")
      cat("Then re-run this prep script.\n\n")
      stop("Cannot proceed without actor variables.")
    }
    
    cat("Found in ACLED. Merging...\n")
    
    acled_merge <- acled[, c("pg_id", "month", actor_vars), with = FALSE]
    
    # Remove if partially present
    for (av in actor_vars) {
      if (av %in% names(dt)) dt[, (av) := NULL]
    }
    
    n_before <- nrow(dt)
    dt <- merge(dt, acled_merge, by = c("pg_id", "month"), all.x = TRUE)
    
    if (nrow(dt) != n_before) {
      stop("ERROR: Merge changed row count. Check for duplicate keys.")
    }
    
    for (av in actor_vars) {
      dt[is.na(get(av)), (av) := 0L]
    }
    
    cat("✓ Actor variables merged successfully\n\n")
  } else {
    cat("\nERROR: acled_pg_month.rds not found.\n")
    cat("SOLUTION: Run R23_audit_and_fix.R to create it.\n\n")
    stop("Cannot proceed without ACLED file.")
  }
} else {
  cat("✓ All actor variables already in panel\n\n")
}

# ══════════════════════════════════════════════════════════════
# STEP 3: CHECK AND ADD COORDINATES
# ══════════════════════════════════════════════════════════════

cat("Checking coordinates...\n")

coord_pairs <- list(
  c("latitude", "longitude"),
  c("lat", "lon"),
  c("ycoord", "xcoord")
)

coords_found <- FALSE
for (pair in coord_pairs) {
  if (all(pair %in% names(dt))) {
    coords_found <- TRUE
    cat("✓ Found coordinates:", paste(pair, collapse = ", "), "\n")
    
    # Ensure we have latitude/longitude (Conley requirement)
    if (!("latitude" %in% names(dt))) {
      if ("lat" %in% names(dt)) {
        dt[, latitude := lat]
        dt[, longitude := lon]
      } else if ("ycoord" %in% names(dt)) {
        dt[, latitude := ycoord]
        dt[, longitude := xcoord]
      }
      cat("✓ Created latitude/longitude from", paste(pair, collapse = ", "), "\n")
    }
    break
  }
}

if (!coords_found) {
  cat("\nWARNING: No coordinate variables found.\n")
  cat("Attempting to load from eli_feasibility_pg.rds...\n")
  
  P_ELI <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
  
  if (file.exists(P_ELI)) {
    eli <- readRDS(P_ELI)
    setDT(eli)
    
    if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := as.integer(gid)]
    
    if (all(c("lon", "lat") %in% names(eli))) {
      coords <- unique(eli[, .(pg_id, latitude = lat, longitude = lon)])
      
      n_before <- nrow(dt)
      dt <- merge(dt, coords, by = "pg_id", all.x = TRUE)
      
      if (nrow(dt) != n_before) {
        stop("ERROR: Coordinate merge changed row count.")
      }
      
      cat("✓ Coordinates merged from ELI file\n")
    } else {
      cat("ERROR: No coordinates in ELI file either.\n")
      stop("Cannot proceed without coordinates.")
    }
  } else {
    cat("ERROR: eli_feasibility_pg.rds not found.\n")
    stop("Cannot proceed without coordinates.")
  }
}

cat("\n")

# ══════════════════════════════════════════════════════════════
# STEP 4: FINAL VALIDATION
# ══════════════════════════════════════════════════════════════

cat("Final validation...\n\n")

required_vars <- c(
  "pg_id", "month", "month_int", "theatre_id",
  "events_type3", "post_eli", "eli_feas",
  "state_civ_harm_n", "rebel_civ_harm_n", "nonstate_civ_harm_n", "state_drone_n",
  "latitude", "longitude"
)

missing_final <- setdiff(required_vars, names(dt))

if (length(missing_final) > 0) {
  cat("ERROR: Still missing required variables:\n")
  print(missing_final)
  stop("Cannot proceed. Check pipeline.")
}

cat("✓ All required variables present\n\n")

cat("Variable summary:\n")
cat("  Events (type-3):", sum(dt$events_type3, na.rm=TRUE), "\n")
cat("  State civilian harm:", sum(dt$state_civ_harm_n, na.rm=TRUE), "\n")
cat("  Rebel civilian harm:", sum(dt$rebel_civ_harm_n, na.rm=TRUE), "\n")
cat("  Non-state harm:", sum(dt$nonstate_civ_harm_n, na.rm=TRUE), "\n")
cat("  State drone strikes:", sum(dt$state_drone_n, na.rm=TRUE), "\n\n")

# ══════════════════════════════════════════════════════════════
# STEP 5: SAVE PREPARED PANEL
# ══════════════════════════════════════════════════════════════

P_OUT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks_PREPARED.rds")

saveRDS(dt, P_OUT)

cat("✓ Prepared panel saved to:\n")
cat(" ", P_OUT, "\n\n")

cat(strrep("=", 70), "\n")
cat("DATA PREPARATION COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Next steps:\n")
cat("1. Source R27_mechanistic_tests_FIXED.R (updated to use prepared panel)\n")
cat("2. Source R28_supplementary_mechanistic_tests.R\n")
cat("3. Review results in outputs/paper1/mechanism_tests/\n\n")

cat("✓ Ready for mechanistic tests!\n")