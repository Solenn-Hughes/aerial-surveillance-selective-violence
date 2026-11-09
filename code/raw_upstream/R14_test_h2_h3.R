# ============================================================
# R14_test_h2_h3.R
# Test H2 (territorial contestation) and H3 (rebel capability)
# via triple interactions in DiD and event-study specifications.
#
# H2: Extra-local intelligence → more selective violence in
#     areas of WEAK / CONTESTED territorial control
#     → expect POSITIVE triple interaction: post_eli × eli_feas × contestation_tv
#     (effect of ELI is stronger where contestation is higher)
#
# H3: Extra-local intelligence → stronger effect where rebels
#     have higher organisational capability (they're worth targeting)
#     → expect POSITIVE triple interaction: post_eli × eli_feas × rebel_cap_tv
#
# Input:  data/processed/panels/panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds (from R13)
# Output: audit tables and model objects
#
# Run AFTER R13.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 14: Test H2 (contestation) and H3 (rebel capability)")

# ── LOAD PANEL ───────────────────────────────────────────────
P_H2H3 <- file.path(DIRS$panels,
                     "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds")
if (!file.exists(P_H2H3)) stop("Missing: ", P_H2H3, " — run R13 first.")

dt <- readRDS(P_H2H3); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

# Outcome
dt[, log1p_y := log1p(events_type3)]

# Check H2/H3 variables exist
H2_VAR <- "contestation_tv"
H3_VAR <- "rebel_cap_tv"

assert_cols(dt, c(H2_VAR, H3_VAR, "eli_feas", "post_eli", "theatre_id", "pg_id"))
log_msg("H2 var: ", H2_VAR, "  H3 var: ", H3_VAR)

# Controls
CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

# ── HELPER FUNCTIONS ─────────────────────────────────────────
extract_terms <- function(m, pattern) {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  ct[grepl(pattern, term)]
}

# ══════════════════════════════════════════════════════════════
# PART 1: DIFFERENCE-IN-DIFFERENCES WITH TRIPLE INTERACTIONS
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 1: DiD TRIPLE INTERACTIONS\n")
cat("============================================================\n")

FE_SPECS <- list(
  global_month  = "pg_id + month_int",
  theatre_month = "pg_id + theatre_id^month_int"
)

did_results <- list()

for (fe_name in names(FE_SPECS)) {
  fe_part <- FE_SPECS[[fe_name]]
  log_msg("DiD with FE: ", fe_name)
  
  # ── M0: Baseline (no moderation) ──
  f0 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + ", ctrl_str,
    " | ", fe_part
  ))
  m0 <- feols(f0, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── M1: + H2 moderation ──
  f1 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:", H2_VAR,
    " + ", ctrl_str, " | ", fe_part
  ))
  m1 <- feols(f1, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── M2: + H3 moderation ──
  f2 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:", H3_VAR,
    " + ", ctrl_str, " | ", fe_part
  ))
  m2 <- feols(f2, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── M3: Both H2 + H3 ──
  f3 <- as.formula(paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:", H2_VAR,
    " + post_eli:eli_feas:", H3_VAR,
    " + ", ctrl_str, " | ", fe_part
  ))
  m3 <- feols(f3, data = dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── M4: With spatial contestation instead ──
  if ("contest_spatial_tv" %in% names(dt)) {
    f4 <- as.formula(paste0(
      "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas:contest_spatial_tv",
      " + post_eli:eli_feas:", H3_VAR,
      " + ", ctrl_str, " | ", fe_part
    ))
    m4 <- feols(f4, data = dt, cluster = ~theatre_id, notes = FALSE)
  }
  
  # Extract key coefficients
  get_row <- function(m, label, pattern) {
    hit <- extract_terms(m, pattern)
    if (nrow(hit) == 0) return(data.table(
      fe = fe_name, model = label, term = NA_character_,
      coef = NA_real_, se = NA_real_, p = NA_real_, nobs = nobs(m)
    ))
    data.table(
      fe = fe_name, model = label, term = hit$term[1],
      coef = round(hit$Estimate[1], 5),
      se   = round(hit$`Std. Error`[1], 5),
      p    = round(hit$`Pr(>|t|)`[1], 4),
      nobs = nobs(m)
    )
  }
  
  did_results[[length(did_results) + 1]] <- get_row(
    m0, "M0: Baseline", "post_eli.*eli_feas|eli_feas.*post_eli")
  did_results[[length(did_results) + 1]] <- get_row(
    m1, "M1: + H2 (contestation)", paste0("post_eli.*eli_feas.*", H2_VAR, "|", H2_VAR, ".*eli_feas.*post_eli"))
  did_results[[length(did_results) + 1]] <- get_row(
    m2, "M2: + H3 (rebel cap)", paste0("post_eli.*eli_feas.*", H3_VAR, "|", H3_VAR, ".*eli_feas.*post_eli"))
  did_results[[length(did_results) + 1]] <- get_row(
    m3, "M3: H2 + H3 (contestation)", paste0("post_eli.*eli_feas.*", H2_VAR, "|", H2_VAR, ".*eli_feas.*post_eli"))
  did_results[[length(did_results) + 1]] <- get_row(
    m3, "M3: H2 + H3 (rebel cap)", paste0("post_eli.*eli_feas.*", H3_VAR, "|", H3_VAR, ".*eli_feas.*post_eli"))
  
  if ("contest_spatial_tv" %in% names(dt)) {
    did_results[[length(did_results) + 1]] <- get_row(
      m4, "M4: Spatial H2 + H3 (spatial contest)",
      "post_eli.*eli_feas.*contest_spatial|contest_spatial.*eli_feas.*post_eli")
  }
  
  # Save model objects
  out_dir <- file.path(DIRS$out_paper1, "models")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(m0, file.path(out_dir, paste0("R14_DiD_M0_", fe_name, ".rds")))
  saveRDS(m1, file.path(out_dir, paste0("R14_DiD_M1_H2_", fe_name, ".rds")))
  saveRDS(m2, file.path(out_dir, paste0("R14_DiD_M2_H3_", fe_name, ".rds")))
  saveRDS(m3, file.path(out_dir, paste0("R14_DiD_M3_H2H3_", fe_name, ".rds")))
}

did_out <- rbindlist(did_results, fill = TRUE)
cat("\nDiD Results:\n")
print(did_out)
write_audit_dt(did_out, "R14_DiD_H2_H3")

# ══════════════════════════════════════════════════════════════
# PART 2: EVENT-STUDY WITH MODERATION
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 2: EVENT-STUDY WITH H2/H3 MODERATION\n")
cat("============================================================\n")

WINDOW <- 24L  # Use shorter window for triple interactions (power)

# Restrict to treated cells with event_time in window + all controls
es_dt <- copy(dt)
es_dt[, in_window := is.na(event_time) | (abs(event_time) <= WINDOW)]
es_dt <- es_dt[in_window == TRUE]

# Create event-time factor (set ref = -1)
es_dt[, et := fifelse(is.na(event_time), NA_integer_, event_time)]

ctrl_str_es <- ctrl_str

for (fe_name in names(FE_SPECS)) {
  fe_part <- FE_SPECS[[fe_name]]
  log_msg("Event-study with FE: ", fe_name, ", window=", WINDOW)
  
  # ── ES0: Base event study ──
  f_es0 <- as.formula(paste0(
    "log1p_y ~ i(et, eli_feas, ref = -1) + ", ctrl_str_es,
    " | ", fe_part
  ))
  m_es0 <- feols(f_es0, data = es_dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── ES1: + H2 moderated ──
  # Interact event-time dummies with eli_feas AND with eli_feas × H2
  f_es1 <- as.formula(paste0(
    "log1p_y ~ i(et, eli_feas, ref = -1) + ",
    "i(et, I(eli_feas * ", H2_VAR, "), ref = -1) + ",
    ctrl_str_es, " | ", fe_part
  ))
  m_es1 <- feols(f_es1, data = es_dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── ES2: + H3 moderated ──
  f_es2 <- as.formula(paste0(
    "log1p_y ~ i(et, eli_feas, ref = -1) + ",
    "i(et, I(eli_feas * ", H3_VAR, "), ref = -1) + ",
    ctrl_str_es, " | ", fe_part
  ))
  m_es2 <- feols(f_es2, data = es_dt, cluster = ~theatre_id, notes = FALSE)
  
  # ── Extract and test: joint significance of POST moderation terms ──
  # For H2: are the post-treatment interaction terms jointly != 0?
  
  extract_es_mod <- function(m, pattern, hyp_name) {
    ct <- as.data.table(coeftable(m), keep.rownames = "term")
    # Post-treatment terms only (event_time > 0)
    post_terms <- ct[grepl(pattern, term) & grepl("et::[1-9]|et::1[0-9]|et::2[0-4]", term)]
    
    if (nrow(post_terms) == 0) {
      log_msg("  No post-treatment moderation terms found for ", hyp_name)
      return(data.table(hypothesis = hyp_name, fe = fe_name,
                        n_post_terms = 0L, mean_coef = NA_real_,
                        share_positive = NA_real_, share_sig05 = NA_real_))
    }
    
    data.table(
      hypothesis = hyp_name, fe = fe_name,
      n_post_terms    = nrow(post_terms),
      mean_coef       = round(mean(post_terms$Estimate), 5),
      median_coef     = round(median(post_terms$Estimate), 5),
      share_positive  = round(mean(post_terms$Estimate > 0), 3),
      share_sig05     = round(mean(post_terms$`Pr(>|t|)` < 0.05), 3),
      min_p           = round(min(post_terms$`Pr(>|t|)`), 4)
    )
  }
  
  es1_h2 <- extract_es_mod(m_es1, H2_VAR, "H2: contestation")
  es2_h3 <- extract_es_mod(m_es2, H3_VAR, "H3: rebel capability")
  
  cat("\n  Event-study moderation (", fe_name, "):\n")
  print(rbind(es1_h2, es2_h3))
  
  # Save
  saveRDS(m_es0, file.path(out_dir, paste0("R14_ES0_base_", fe_name, ".rds")))
  saveRDS(m_es1, file.path(out_dir, paste0("R14_ES1_H2_", fe_name, ".rds")))
  saveRDS(m_es2, file.path(out_dir, paste0("R14_ES2_H3_", fe_name, ".rds")))
}

# ══════════════════════════════════════════════════════════════
# PART 3: SAMPLE SPLIT (ROBUSTNESS)
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("PART 3: SAMPLE SPLIT BY CONTESTATION LEVEL\n")
cat("============================================================\n")

# Split at median of contestation_tv (across treated cells only)
med_h2 <- median(dt[eli_feas > 0, contestation_tv], na.rm = TRUE)
dt[, h2_high := contestation_tv > med_h2]

split_results <- list()

for (fe_name in names(FE_SPECS)) {
  fe_part <- FE_SPECS[[fe_name]]
  
  for (grp in c("High contestation", "Low contestation")) {
    dt_sub <- if (grp == "High contestation") dt[h2_high == TRUE] else dt[h2_high == FALSE]
    
    f_split <- as.formula(paste0(
      "log1p_y ~ post_eli:eli_feas + ", ctrl_str,
      " | ", fe_part
    ))
    
    m_split <- tryCatch(
      feols(f_split, data = dt_sub, cluster = ~theatre_id, notes = FALSE),
      error = function(e) { log_msg("  Split failed for ", grp, ": ", e$message); NULL }
    )
    
    if (!is.null(m_split)) {
      hit <- extract_terms(m_split, "post_eli.*eli_feas|eli_feas.*post_eli")
      if (nrow(hit) > 0) {
        split_results[[length(split_results) + 1]] <- data.table(
          fe = fe_name, subsample = grp,
          coef = round(hit$Estimate[1], 5),
          se   = round(hit$`Std. Error`[1], 5),
          p    = round(hit$`Pr(>|t|)`[1], 4),
          nobs = nobs(m_split)
        )
      }
    }
  }
}

if (length(split_results)) {
  split_out <- rbindlist(split_results, fill = TRUE)
  cat("\nSample split results:\n")
  print(split_out)
  write_audit_dt(split_out, "R14_sample_split_H2")
}

# ══════════════════════════════════════════════════════════════
# PART 4: INTERPRETATION GUIDE
# ══════════════════════════════════════════════════════════════

cat("\n============================================================\n")
cat("INTERPRETATION GUIDE\n")
cat("============================================================\n")
cat("
H2 (Territorial contestation):
  Triple interaction post_eli × eli_feas × contestation_tv
  - POSITIVE coefficient: ELI has a STRONGER effect (more type 3
    violence against civilians) in contested areas. This would
    support H2: extra-local intelligence matters more where
    territorial control is weak.
  - NEGATIVE coefficient: ELI effect is weaker in contested areas.
  - NOT SIGNIFICANT: contestation does not moderate the effect.

H3 (Rebel capability):
  Triple interaction post_eli × eli_feas × rebel_cap_tv
  - POSITIVE coefficient: ELI has a STRONGER effect where rebels
    are more capable. Supports H3: information advantage matters
    more against organised adversaries.
  - NEGATIVE coefficient: effect is weaker against capable rebels.
  - NOT SIGNIFICANT: rebel capability does not moderate.

Sample split corroborates: if H2 is supported, the coefficient
in 'High contestation' should be larger than in 'Low contestation'.
")

log_msg("Step 14 complete.")
