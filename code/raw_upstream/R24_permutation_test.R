# ============================================================
# R24_permutation_test.R — Randomisation inference
#
# Permute theatre-level onset dates 500 times, re-estimate β
# each time, and compute the share of permuted β that exceed
# the real β. This is a non-parametric p-value that doesn't
# rely on distributional assumptions or cluster counts.
#
# If the real β is larger than 95% of permuted βs, the result
# is significant at 5% by randomisation inference.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
# Output: permutation distribution + rank-based p-value
#
# Run time: ~15-30 min (500 iterations)
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 24: Permutation / randomisation inference")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

ctrl <- c("prec_gpcp", "temp", "water_ih", "forest_ih")
ctrl_str <- paste(intersect(ctrl, names(dt)), collapse = " + ")
FE <- "pg_id + theatre_id^month_int"

# ── REAL β ───────────────────────────────────────────────────
f_real <- as.formula(paste0("log1p_y ~ post_eli:eli_feas + ", ctrl_str, " | ", FE))
m_real <- feols(f_real, data = dt, cluster = ~pg_id, notes = FALSE)
ct_real <- as.data.table(coeftable(m_real), keep.rownames = "term")
beta_real <- ct_real[grepl("eli_feas", term), Estimate][1]
cat(sprintf("\nReal \u03B2: %.5f\n", beta_real))

# ── ONSET DATES TO PERMUTE ───────────────────────────────────
# Get unique theatre onset dates
onset_map <- unique(dt[treated == 1L, .(theatre_id, eli_onset_month)])
setorder(onset_map, theatre_id)
n_theatres <- nrow(onset_map)
cat("Theatres with onset dates:\n")
print(onset_map)

# All unique months in the panel (potential onset candidates)
all_months <- sort(unique(dt$month))
# Restrict to reasonable range (not too early, not too late)
month_range <- all_months[all_months >= as.IDate("2000-01-01") & 
                          all_months <= as.IDate("2022-06-01")]

# ── PERMUTATION LOOP ─────────────────────────────────────────
N_PERM <- 500
set.seed(42)

cat(sprintf("\nRunning %d permutations...\n", N_PERM))

perm_betas <- numeric(N_PERM)

for (i in seq_len(N_PERM)) {
  if (i %% 50 == 0) cat(sprintf("  Permutation %d / %d\n", i, N_PERM))
  
  # Randomly assign onset dates to theatres
  perm_onsets <- sample(month_range, n_theatres, replace = FALSE)
  perm_map <- data.table(theatre_id = onset_map$theatre_id,
                         perm_onset = perm_onsets)
  
  # Create permuted post_eli
  dt_perm <- copy(dt)
  dt_perm[, perm_onset := NULL]
  dt_perm <- merge(dt_perm, perm_map, by = "theatre_id", all.x = TRUE)
  dt_perm[, post_eli_perm := as.integer(!is.na(perm_onset) & month >= perm_onset)]
  
  # Estimate
  f_perm <- as.formula(paste0("log1p_y ~ post_eli_perm:eli_feas + ", ctrl_str, " | ", FE))
  m_perm <- tryCatch(
    feols(f_perm, data = dt_perm, cluster = ~pg_id, notes = FALSE),
    error = function(e) NULL
  )
  
  if (!is.null(m_perm)) {
    ct_perm <- as.data.table(coeftable(m_perm), keep.rownames = "term")
    hit <- ct_perm[grepl("eli_feas", term)]
    perm_betas[i] <- if (nrow(hit) > 0) hit$Estimate[1] else NA_real_
  } else {
    perm_betas[i] <- NA_real_
  }
}

# ── RESULTS ──────────────────────────────────────────────────
perm_betas <- perm_betas[!is.na(perm_betas)]
n_valid <- length(perm_betas)

# One-sided p-value: share of permuted β ≥ real β
p_perm <- mean(perm_betas >= beta_real)

cat(sprintf("\n============================================================\n"))
cat(sprintf("PERMUTATION INFERENCE RESULTS\n"))
cat(sprintf("============================================================\n"))
cat(sprintf("Real \u03B2:                    %.5f\n", beta_real))
cat(sprintf("Permutation mean:           %.5f\n", mean(perm_betas)))
cat(sprintf("Permutation SD:             %.5f\n", sd(perm_betas)))
cat(sprintf("Permutation 95th pctile:    %.5f\n", quantile(perm_betas, 0.95)))
cat(sprintf("Valid permutations:         %d / %d\n", n_valid, N_PERM))
cat(sprintf("Randomisation p-value:      %.4f\n", p_perm))

if (p_perm < 0.05) {
  cat("\n\u2713  Real \u03B2 exceeds 95%% of permuted \u03B2s. Significant at 5%%.\n")
} else {
  cat(sprintf("\n\u26A0\uFE0F  Real \u03B2 does not exceed 95%% of permuted \u03B2s (p = %.3f).\n", p_perm))
}

# Save
perm_out <- data.table(
  beta_real = beta_real,
  perm_mean = mean(perm_betas),
  perm_sd = sd(perm_betas),
  perm_p95 = quantile(perm_betas, 0.95),
  p_value = p_perm,
  n_perms = n_valid
)
write_audit_dt(perm_out, "R24_permutation_inference")

# Save full distribution for potential histogram
perm_dist <- data.table(iteration = seq_along(perm_betas), beta_perm = perm_betas)
write_audit_dt(perm_dist, "R24_permutation_distribution")

log_msg("Step 24 complete.")
