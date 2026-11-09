# ============================================================
# R09f_robustness_decay_radius.R  (FIXED — self-contained)
#
# Part 1: Decay-radius robustness — re-estimates eli_feas at
#   R = 200 / 400 / 800 km and checks β stability.
#
# Part 2: Component decomposition — separates eli_feas into
#   (a) distance-decay only, (b) terrain only, (c) composite.
#   Confirms β is not an artefact of index construction.
#
# FE:       pg_id + theatre_id^month_int
# Controls: nlights_calib_mean, nlights_mean, nlights_min,
#           prec_gpcp, temp, water_ih
# VCOV:     Conley 200km (primary) + theatre-clustered (secondary)
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_msg("R09f: Decay-radius robustness + component decomposition")

# ── SETTINGS (all defined here — no external dependencies) ───
OUTCOME      <- "events_type3"
CTRL_M0      <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
                  "prec_gpcp", "temp", "water_ih")
FE_SPEC      <- "pg_id + theatre_id^month_int"
R_VALS       <- c(200L, 400L, 800L)
R0_BASELINE  <- 400L   # radius used to build the panel's eli_feas

# ── HELPERS ──────────────────────────────────────────────────
assert_file <- function(p) {
  if (!file.exists(p)) stop("Missing file: ", p)
}
assert_cols2 <- function(dt, cols, nm = "dt") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(nm, " missing cols: ", paste(miss, collapse = ", "))
}

safe_conley <- function(m, cutoff_km = 200) {
  res <- tryCatch(summary(m, vcov = conley(cutoff = cutoff_km)), error = function(e) NULL)
  if (!is.null(res)) return(list(m = res, vcov_type = paste0("Conley_", cutoff_km, "km")))
  res <- tryCatch(summary(m, vcov = vcov_conley(cutoff = cutoff_km)), error = function(e) NULL)
  if (!is.null(res)) return(list(m = res, vcov_type = paste0("Conley_", cutoff_km, "km")))
  list(m = summary(m, cluster = ~theatre_id), vcov_type = "theatre_clustered_fallback")
}

extract_hit <- function(m_summary, pattern, label, component, vcov_type) {
  ct <- as.data.table(coeftable(m_summary), keep.rownames = "term")
  hit <- ct[grepl(pattern, term, perl = TRUE)][1]
  if (nrow(hit) == 0L || is.na(hit$Estimate)) {
    return(data.table(label = label, component = component, vcov_type = vcov_type,
                      coef = NA_real_, se = NA_real_, t_val = NA_real_, p = NA_real_,
                      nobs = NA_integer_))
  }
  data.table(
    label     = label,
    component = component,
    vcov_type = vcov_type,
    coef  = round(hit$Estimate[1],     5),
    se    = round(hit$`Std. Error`[1], 5),
    t_val = round(hit$`t value`[1],    3),
    p     = round(hit$`Pr(>|t|)`[1],   4),
    nobs  = nobs(m_summary)
  )
}

# ── 1) LOAD PANEL ─────────────────────────────────────────────
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
assert_file(P_PANEL)
dt <- readRDS(P_PANEL); setDT(dt)

if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id      := as.integer(pg_id)]
dt[, month      := if (exists("as_month_idate")) as_month_idate(month) else as.IDate(month)]
dt[, month_int  := year(month) * 12L + month(month)]
dt[, theatre_id := as.character(theatre_id)]

assert_cols2(dt, c("pg_id","theatre_id","month","month_int",
                   "post_eli","eli_feas","events_type3"), "panel")

# Coordinates for Conley
if ("xcoord" %in% names(dt)) {
  dt[, longitude := xcoord]; dt[, latitude := ycoord]
} else stop("xcoord/ycoord not found.")

# Check controls present
miss_ctrl <- setdiff(CTRL_M0, names(dt))
if (length(miss_ctrl)) stop("Missing controls: ", paste(miss_ctrl, collapse = ", "))
rhs_controls <- paste(CTRL_M0, collapse = " + ")

log_msg("Panel loaded: ", nrow(dt), " rows, ", uniqueN(dt$theatre_id), " theatres")

# ── 2) LOAD TERRAIN COMPONENT ─────────────────────────────────
cand_feas <- c(
  file.path(DIRS$eli, "eli_feasibility_pg.rds"),
  file.path("data", "processed", "eli", "eli_feasibility_pg.rds"),
  file.path("data", "processed", "eli", "eli_feasibility_pg_v2.rds")
)
P_FEAS <- cand_feas[file.exists(cand_feas)][1]
if (is.na(P_FEAS) || !nzchar(P_FEAS))
  stop("eli_feasibility_pg.rds not found. Check DIRS$eli.")

feas <- readRDS(P_FEAS); setDT(feas)
if (!("pg_id" %in% names(feas)) && "gid" %in% names(feas)) feas[, pg_id := gid]

terrain_col <- intersect(c("eli_feas_terrain","feas_terrain","eli_terrain"), names(feas))[1]
if (is.na(terrain_col))
  stop("Feasibility file has no terrain column. Expected eli_feas_terrain / feas_terrain / eli_terrain.")

feas_small <- feas[, .(pg_id = as.integer(pg_id),
                       eli_feas_terrain = as.numeric(get(terrain_col)))]
setkey(feas_small, pg_id); setkey(dt, pg_id)
dt <- feas_small[dt]
if (anyNA(dt$eli_feas_terrain)) stop("eli_feas_terrain has NA after merge.")

# ── 3) BACK OUT DISTANCE FROM BASELINE eli_feas ───────────────
# eli_feas = terrain × exp(-dist / R0)  =>  dist = -R0 × log(eli_feas / terrain)
dt[, range_base := eli_feas / eli_feas_terrain]
dt[!is.finite(range_base) | range_base <= 0, range_base := NA_real_]
dt[range_base > 1, range_base := 1]
dt[is.na(range_base),  dist_km := 1e6]          # terrain=0 → treat as very far
dt[!is.na(range_base), dist_km := -R0_BASELINE * log(range_base)]
if (any(!is.finite(dt$dist_km))) stop("Non-finite dist_km.")

log_msg("Distance back-calculation complete.")

# ── 4) PART 1: DECAY-RADIUS SWEEP ────────────────────────────
cat("\n====== PART 1: DECAY-RADIUS ROBUSTNESS ======\n")

ylog <- paste0("log1p_", OUTCOME)
dt[, (ylog) := log1p(get(OUTCOME))]

run_one_R <- function(R_km) {
  x <- copy(dt)
  x[, eli_feas_R := eli_feas_terrain * exp(-dist_km / R_km)]
  
  fml <- as.formula(paste0(
    ylog, " ~ post_eli:eli_feas_R + ", rhs_controls,
    " | ", FE_SPEC
  ))
  
  m_raw <- tryCatch(feols(fml, data = x, notes = FALSE),
                    error = function(e) { log_msg("ERROR R=", R_km, ": ", e$message); NULL })
  if (is.null(m_raw)) {
    return(data.table(R_km = as.integer(R_km), label = "decay_radius",
                      component = paste0("R=", R_km, "km"),
                      vcov_type = NA, coef = NA_real_, se = NA_real_,
                      t_val = NA_real_, p = NA_real_, nobs = NA_integer_))
  }
  
  # Primary: Conley 200km
  res_con <- safe_conley(m_raw, 200)
  hit_con <- extract_hit(res_con$m, "post_eli.*eli_feas_R|eli_feas_R.*post_eli",
                         "decay_radius", paste0("R=", R_km, "km"), res_con$vcov_type)
  hit_con[, R_km := as.integer(R_km)]
  
  # Secondary: theatre-clustered (for comparison only)
  m_th  <- tryCatch(summary(m_raw, cluster = ~theatre_id), error = function(e) NULL)
  if (!is.null(m_th)) {
    hit_th <- extract_hit(m_th, "post_eli.*eli_feas_R|eli_feas_R.*post_eli",
                          "decay_radius", paste0("R=", R_km, "km"), "theatre_clustered")
    hit_th[, R_km := as.integer(R_km)]
    return(rbind(hit_con, hit_th, fill = TRUE))
  }
  hit_con
}

decay_res <- rbindlist(lapply(R_VALS, run_one_R), fill = TRUE)

cat("\nDecay-radius results (Conley 200km primary):\n")
print(decay_res[grepl("Conley", vcov_type), .(R_km, component, coef, se, p, nobs)])
cat("\nTheatre-clustered (secondary / comparison only):\n")
print(decay_res[grepl("theatre", vcov_type), .(R_km, component, coef, se, p, nobs)])

out_decay <- file.path(DIRS$audit, paste0("R09f_decay_radius_", TS, ".csv"))
fwrite(decay_res, out_decay)
log_msg("Wrote: ", out_decay)

# ── 5) PART 2: COMPONENT DECOMPOSITION ───────────────────────
cat("\n====== PART 2: COMPONENT DECOMPOSITION ======\n")
cat("Tests whether β is driven by distance decay, terrain, or both.\n\n")

# Recompute range component at baseline R = 400km
dt[, eli_feas_range := exp(-dist_km / R0_BASELINE)]

# Normalise both components to [0,1] so coefficients are comparable
dt[, eli_feas_terrain_sc := eli_feas_terrain]   # already in [0,1]
dt[, eli_feas_range_sc   := eli_feas_range]      # already in [0,1]

decomp_rows <- list()

run_decomp <- function(feas_var, label) {
  fml <- as.formula(paste0(
    ylog, " ~ post_eli:", feas_var, " + ", rhs_controls,
    " | ", FE_SPEC
  ))
  m_raw <- tryCatch(feols(fml, data = dt, notes = FALSE),
                    error = function(e) { log_msg("ERROR decomp ", label, ": ", e$message); NULL })
  if (is.null(m_raw)) {
    return(data.table(label = label, component = feas_var,
                      vcov_type = NA, coef = NA_real_, se = NA_real_,
                      t_val = NA_real_, p = NA_real_, nobs = NA_integer_))
  }
  res_con <- safe_conley(m_raw, 200)
  extract_hit(res_con$m,
              paste0("post_eli.*", feas_var, "|", feas_var, ".*post_eli"),
              label, feas_var, res_con$vcov_type)
}

decomp_rows[[1]] <- run_decomp("eli_feas_range",   "Distance decay only")
decomp_rows[[2]] <- run_decomp("eli_feas_terrain",  "Terrain openness only")
decomp_rows[[3]] <- run_decomp("eli_feas",          "Full composite (range × terrain)")

decomp_dt <- rbindlist(decomp_rows, fill = TRUE)

cat("\nComponent decomposition results (Conley 200km):\n")
print(decomp_dt[, .(label, coef, se, p, nobs)])
cat("\nIf distance AND terrain are both positive → both components contribute.\n")
cat("If only one is significant → that component drives the index.\n")

out_decomp <- file.path(DIRS$audit, paste0("R09f_component_decomposition_", TS, ".csv"))
fwrite(decomp_dt, out_decomp)
log_msg("Wrote: ", out_decomp)

# ── 6) COMBINED SUMMARY ───────────────────────────────────────
cat("\n====== COMBINED ROBUSTNESS SUMMARY ======\n")
cat("PART 1 — Decay radius (Conley 200km):\n")
print(decay_res[grepl("Conley", vcov_type), .(R_km, coef, se, p)])
cat("\nPART 2 — Component decomposition (Conley 200km):\n")
print(decomp_dt[, .(label, coef, se, p)])

cat("\n✅ R09f complete.\n")
log_msg("R09f complete.")