# ============================================================
# R09k_joint_post_test.R
# 
# Computes a formal joint Wald test that all post-treatment
# event-study coefficients = 0. More powerful than eyeballing
# individual period t-tests.
#
# Also computes average post-treatment effect with SE.
#
# Requires: saved models RDS from R09_models_event_study_v2.R
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 9k: Joint post-treatment tests on event study models")

# ── FIND SAVED MODELS ──────────────────────────────────────
models_rds <- list.files(DIRS$out_paper1, pattern = "R09v2_all_results.*\\.rds$", full.names = TRUE)
if (!length(models_rds)) {
  # Try alternative naming
  models_rds <- list.files(DIRS$out_paper1, pattern = "all_results.*\\.rds$", full.names = TRUE)
}
if (!length(models_rds)) stop("No saved models RDS found in ", DIRS$out_paper1, 
                               ". Re-run R09_models_event_study_v2.R with SAVE_MODELS_RDS=TRUE.")

# Use the most recent
models_path <- models_rds[which.max(file.info(models_rds)$mtime)]
log_msg("Loading models from: ", models_path)
all_results <- readRDS(models_path)

# ── FUNCTION: JOINT POST TEST ─────────────────────────────
joint_post_test <- function(m, model_label = "") {
  cn <- names(coef(m))
  # Find event_time terms
  et_terms <- cn[grepl("^event_time_es::", cn)]
  if (!length(et_terms)) return(NULL)
  
  # Parse event times
  k <- suppressWarnings(as.integer(sub("^event_time_es::\\(?(-?\\d+)\\)?.*$", "\\1", et_terms)))
  
  # Post-treatment terms (k >= 0)
  post_terms <- et_terms[!is.na(k) & k >= 0]
  pre_terms  <- et_terms[!is.na(k) & k < 0 & k != -1L]  # excluding ref = -1
  
  # Joint Wald test: all post-treatment = 0
  post_wald <- tryCatch(wald(m, post_terms), error = function(e) NULL)
  pre_wald  <- tryCatch(wald(m, pre_terms),  error = function(e) NULL)
  
  # Average post-treatment effect
  post_coefs <- coef(m)[post_terms]
  avg_post   <- mean(post_coefs)
  
  # Average pre-treatment effect (for comparison)
  pre_coefs <- coef(m)[pre_terms]
  avg_pre   <- mean(pre_coefs)
  
  data.table(
    model = model_label,
    n_post_periods = length(post_terms),
    n_pre_periods  = length(pre_terms),
    avg_post_coef  = round(avg_post, 4),
    avg_pre_coef   = round(avg_pre, 4),
    post_wald_F    = if (!is.null(post_wald)) round(as.numeric(post_wald$stat), 3) else NA_real_,
    post_wald_p    = if (!is.null(post_wald)) round(as.numeric(post_wald$p.value), 4) else NA_real_,
    pre_wald_F     = if (!is.null(pre_wald))  round(as.numeric(pre_wald$stat), 3) else NA_real_,
    pre_wald_p     = if (!is.null(pre_wald))  round(as.numeric(pre_wald$p.value), 4) else NA_real_
  )
}

# ── RUN ACROSS ALL ES MODELS ─────────────────────────────
results <- list()

for (panel_name in names(all_results)) {
  res <- all_results[[panel_name]]
  if (is.null(res$models_es) || !length(res$models_es)) next
  
  for (mname in names(res$models_es)) {
    m <- res$models_es[[mname]]
    label <- paste0(panel_name, " | ", mname)
    r <- joint_post_test(m, label)
    if (!is.null(r)) results[[length(results) + 1]] <- r
  }
}

if (length(results)) {
  out <- rbindlist(results, fill = TRUE)
  
  log_msg("")
  log_msg("============================================")
  log_msg("JOINT POST-TREATMENT TESTS")
  log_msg("============================================")
  print(out)
  
  write_audit_dt(out, "R09k_joint_post_tests")
  log_msg("")
  log_msg("Key: post_wald_p is the p-value for H0: all post-treatment coefficients = 0")
  log_msg("     If post_wald_p < 0.05, there is a statistically significant post-treatment effect")
  log_msg("     avg_post_coef is the simple average of post-treatment period estimates")
} else {
  log_msg("No event study models found in saved results.")
}

log_msg("Step 9k complete.")
