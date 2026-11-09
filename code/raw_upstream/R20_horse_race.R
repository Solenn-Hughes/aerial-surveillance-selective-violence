# ============================================================
# R20_horse_race.R — Coefficient stability table
#
# Reports β from:
#   (a) Baseline: post_eli × eli_feas only
#   (b) + post_eli × population
#   (c) + post_eli × nightlights
#   (d) + post_eli × frontline proximity
#   (e) All three horse-race controls together
#
# If β barely moves across (a)–(e), cross-sectional confounders
# correlated with eli_feas are ruled out.
#
# Input:  panel_pg_month_ALLCW_paper1_shocks.rds
# Output: coefficient stability table
#
# Run AFTER R16. Fast (~5 min).
# ============================================================

rm(list = ls())
suppressMessages({ library(data.table); library(fixest) })

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 20: Horse-race coefficient stability table")

# ── LOAD ─────────────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
if (!file.exists(P_PANEL)) stop("Missing: ", P_PANEL)

dt <- readRDS(P_PANEL); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

FE <- "pg_id + theatre_id^month_int"

# Base controls (non-horse-race)
BASE_CTRL <- c("prec_gpcp", "temp", "water_ih", "forest_ih")
base_present <- intersect(BASE_CTRL, names(dt))
base_str <- paste(base_present, collapse = " + ")

# ── IDENTIFY HORSE-RACE VARIABLES ────────────────────────────
# Population
pop_col <- NULL
for (cc in c("pop_gpw_sum", "pop_gpw_max")) {
  if (cc %in% names(dt) && !all(is.na(dt[[cc]]))) { pop_col <- cc; break }
}

# Nightlights
nl_col <- NULL
for (cc in c("nlights_mean", "nlights_calib_mean")) {
  if (cc %in% names(dt) && !all(is.na(dt[[cc]]))) { nl_col <- cc; break }
}

# Frontline proximity
front_col <- NULL
for (cc in c("contestation_tv", "control_weak", "ttime_mean")) {
  if (cc %in% names(dt) && !all(is.na(dt[[cc]]))) { front_col <- cc; break }
}

log_msg("Horse-race vars: pop=", pop_col, " nl=", nl_col, " front=", front_col)

# ── BUILD SPECIFICATIONS ─────────────────────────────────────
specs <- list()

# (a) Baseline
specs[["(a) Baseline"]] <- paste0(
  "log1p_y ~ post_eli:eli_feas + ", base_str, " | ", FE
)

# (b) + post × population
if (!is.null(pop_col)) {
  specs[[paste0("(b) + post×", pop_col)]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:", pop_col, " + ", base_str, " | ", FE
  )
}

# (c) + post × nightlights
if (!is.null(nl_col)) {
  specs[[paste0("(c) + post×", nl_col)]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:", nl_col, " + ", base_str, " | ", FE
  )
}

# (d) + post × frontline
if (!is.null(front_col)) {
  specs[[paste0("(d) + post×", front_col)]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:", front_col, " + ", base_str, " | ", FE
  )
}

# (e) All three together
hr_terms <- c()
if (!is.null(pop_col))   hr_terms <- c(hr_terms, paste0("post_eli:", pop_col))
if (!is.null(nl_col))    hr_terms <- c(hr_terms, paste0("post_eli:", nl_col))
if (!is.null(front_col)) hr_terms <- c(hr_terms, paste0("post_eli:", front_col))

if (length(hr_terms) > 0) {
  specs[["(e) All horse-race"]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + ", paste(hr_terms, collapse = " + "),
    " + ", base_str, " | ", FE
  )
}

# ── RUN ALL SPECIFICATIONS ───────────────────────────────────
results <- list()

for (nm in names(specs)) {
  log_msg("Running: ", nm)
  f <- as.formula(specs[[nm]])
  m <- tryCatch(
    feols(f, data = dt, cluster = ~theatre_id, notes = FALSE),
    error = function(e) { log_msg("ERROR: ", nm, " — ", e$message); NULL }
  )
  
  if (!is.null(m)) {
    ct <- as.data.table(coeftable(m), keep.rownames = "term")
    hit <- ct[grepl("post_eli.*eli_feas|eli_feas.*post_eli", term)]
    if (nrow(hit) > 0) {
      results[[nm]] <- data.table(
        spec  = nm,
        beta  = round(hit$Estimate[1], 5),
        se    = round(hit$`Std. Error`[1], 5),
        p     = round(hit$`Pr(>|t|)`[1], 4),
        nobs  = nobs(m)
      )
    }
  }
}

out <- rbindlist(results)

cat("\n============================================================\n")
cat("COEFFICIENT STABILITY TABLE\n")
cat("============================================================\n")
print(out)

# Stability diagnostic
if (nrow(out) >= 2) {
  beta_range <- range(out$beta, na.rm = TRUE)
  cat(sprintf("\nβ range: [%.5f, %.5f]\n", beta_range[1], beta_range[2]))
  cat(sprintf("Max change from baseline: %.1f%%\n",
              100 * max(abs(out$beta - out$beta[1])) / abs(out$beta[1])))
  
  if (max(abs(out$beta - out$beta[1])) / abs(out$beta[1]) < 0.20) {
    cat("✓  β is stable (<20% change). Horse-race passed.\n")
  } else {
    cat("⚠️  β moved >20%. Investigate which control is absorbing variation.\n")
  }
}

write_audit_dt(out, "R20_horse_race_stability")

# Save all models
out_dir <- file.path(DIRS$out_paper1, "models")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
for (nm in names(specs)) {
  f <- as.formula(specs[[nm]])
  m <- tryCatch(feols(f, data = dt, cluster = ~theatre_id, notes = FALSE), error = function(e) NULL)
  if (!is.null(m)) saveRDS(m, file.path(out_dir, paste0("R20_", gsub("[^a-zA-Z0-9]", "_", nm), ".rds")))
}

log_msg("Step 20 complete.")
