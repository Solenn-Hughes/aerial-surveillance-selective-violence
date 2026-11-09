# Just run R28-R33 in order
# Skipping R31 (too memory intensive)

cat("\n=== Running R28_FIXED ===\n")
source("R28_FIXED.R")

cat("\n=== Running R29 ===\n")
if (file.exists("R29_shock_actor_decomp.R")) {
  source("R29_shock_actor_decomp.R")
} else {
  cat("R29 not found, skipping\n")
}

cat("\n=== Running R30 (Rambachan-Roth) ===\n")
if (file.exists("R30_rambachan_roth.R")) {
  source("R30_rambachan_roth.R")
} else {
  cat("R30 not found, skipping\n")
}

cat("\n=== Skipping R31 (memory issues) ===\n")

cat("\n=== Running R32 ===\n")
if (file.exists("R32_rebel_drone_acquisition_shock.R")) {
  source("R32_rebel_drone_acquisition_shock.R")
} else {
  cat("R32 not found, skipping\n")
}

cat("\n=== Running R33 ===\n")
if (file.exists("R33_spillover_check.R")) {
  source("R33_spillover_check.R")
} else {
  cat("R33 not found, skipping\n")
}

cat("\n=== DONE ===\n")