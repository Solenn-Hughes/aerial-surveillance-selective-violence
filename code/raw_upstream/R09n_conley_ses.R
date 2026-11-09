# ============================================================
# R09n_conley_ses.R  (FIXED v2)
# FIX v2: replaced !!term with term_val variable (data.table
#         doesn't support !! bang-bang operator)
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 9n: Conley spatial HAC standard errors")

P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
if (!file.exists(P_ALLCW)) stop("Missing: ", P_ALLCW)

dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y := log1p(events_type3)]

# Coordinates from panel
if ("xcoord" %in% names(dt) && "ycoord" %in% names(dt)) {
  dt[, longitude := xcoord]
  dt[, latitude  := ycoord]
  log_msg("Using xcoord/ycoord from panel as longitude/latitude")
} else {
  stop("xcoord/ycoord not found in panel. Check your panel build.")
}

# Fit baseline model
CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

fml <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))

log_msg("Fitting baseline model...")
m <- feols(fml, data = dt, notes = FALSE)

cn <- names(coef(m))
term_val <- cn[grepl("post_eli.*eli_feas|eli_feas.*post_eli", cn)][1]
if (is.na(term_val)) stop("Interaction term not found.")

base_coef <- coef(m)[term_val]
log_msg("Baseline coefficient: ", round(base_coef, 5))

# ── VCOV COMPARISON ──────────────────────────────────────────
CUTOFFS_KM <- c(100, 200, 500)
results <- list()

# Helper to extract the interaction row from coeftable
get_hit <- function(m_summary) {
  ct <- as.data.table(coeftable(m_summary), keep.rownames = "term")
  ct[term == term_val]  # term_val is in parent scope
}

# 1) Theatre-clustered
se_theatre <- summary(m, cluster = ~theatre_id)
hit_th <- get_hit(se_theatre)
results[[1]] <- data.table(
  vcov_type = "Theatre-clustered (baseline)",
  cutoff_km = NA_integer_,
  coef  = round(hit_th$Estimate[1], 5),
  se    = round(hit_th$`Std. Error`[1], 5),
  t     = round(hit_th$`t value`[1], 3),
  p     = round(hit_th$`Pr(>|t|)`[1], 4)
)

# 2) Two-way clustered
se_twoway <- summary(m, cluster = ~pg_id + month_int)
hit_tw <- get_hit(se_twoway)
results[[2]] <- data.table(
  vcov_type = "Two-way (pg_id + month_int)",
  cutoff_km = NA_integer_,
  coef  = round(hit_tw$Estimate[1], 5),
  se    = round(hit_tw$`Std. Error`[1], 5),
  t     = round(hit_tw$`t value`[1], 3),
  p     = round(hit_tw$`Pr(>|t|)`[1], 4)
)

# 3) Conley spatial HAC
try_conley <- function(m, cutoff_km) {
  # Approach 1: auto-detect latitude/longitude columns
  res <- tryCatch(summary(m, vcov = conley(cutoff = cutoff_km)),
                  error = function(e) NULL)
  if (!is.null(res)) return(res)
  
  # Approach 2: vcov_conley
  res <- tryCatch(summary(m, vcov = vcov_conley(cutoff = cutoff_km)),
                  error = function(e) NULL)
  if (!is.null(res)) return(res)
  
  NULL
}

for (d_km in CUTOFFS_KM) {
  log_msg("Computing Conley SE with cutoff = ", d_km, " km...")
  se_conley <- try_conley(m, d_km)
  
  if (!is.null(se_conley)) {
    hit_con <- get_hit(se_conley)
    results[[length(results) + 1]] <- data.table(
      vcov_type = paste0("Conley (", d_km, " km)"),
      cutoff_km = d_km,
      coef  = round(hit_con$Estimate[1], 5),
      se    = round(hit_con$`Std. Error`[1], 5),
      t     = round(hit_con$`t value`[1], 3),
      p     = round(hit_con$`Pr(>|t|)`[1], 4)
    )
    log_msg("  Success.")
  } else {
    log_msg("  Conley failed at ", d_km, " km. fixest version: ",
            as.character(packageVersion("fixest")))
  }
}

out <- rbindlist(results)

cat("\n============================================================\n")
cat("STANDARD ERROR COMPARISON\n")
cat("============================================================\n")
print(out)

cat("\n── Interpretation ──\n")
cat("Theatre-clustered SE treats all cells within a theatre as perfectly correlated.\n")
cat("Two-way clustering allows for both spatial (pg_id) and temporal (month) correlation.\n")
cat("Conley allows spatial correlation to decay smoothly with distance.\n")
cat("If two-way SE is SMALLER → theatre-clustering is conservative (good).\n")

write_audit_dt(out, "R09n_conley_se_comparison")
log_msg("Step 9n complete.")