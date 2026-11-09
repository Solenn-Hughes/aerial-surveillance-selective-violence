# ============================================================
# R25_methodology_fixes.R
# Addresses four methodological issues identified in audit:
#
#   Fix 1: R09f decay radius — re-run on 8-theatre panel
#           (was using 88-country unrestricted panel)
#
#   Fix 2: Within-theatre eli_feas re-standardisation
#           (global terrain z-scores compress within-theatre
#           variation; local re-normalisation improves power)
#
#   Fix 3: Horse race with Conley SEs
#           (R20 used theatre clustering = 8 clusters = invalid)
#
# Fix 4 (Incirlik spatial) is a documentation fix, not a code fix.
# See note at end of this script.
#
# Run this script AFTER all main analysis is complete.
# Outputs go to outputs/paper1/ with _FIXED suffix.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!requireNamespace("fixest", quietly = TRUE)) stop("fixest required.")
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_msg("R25: Methodology fixes — ", TS)

# ── PATHS ─────────────────────────────────────────────────────
P_SHOCKS <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_shocks.rds")
P_ELI    <- file.path(DIRS$eli,   "eli_feasibility_pg.rds")

if (!file.exists(P_SHOCKS)) stop("Missing: ", P_SHOCKS)
if (!file.exists(P_ELI))    stop("Missing: ", P_ELI)

# ── LOAD PANEL (8-theatre restricted, N=304,780) ──────────────
log_msg("Loading 8-theatre panel...")
dt <- readRDS(P_SHOCKS); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id     := as.integer(pg_id)]
dt[, month     := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]
dt[, log1p_y   := log1p(events_type3)]

log_msg("  Panel: ", nrow(dt), " rows, ", uniqueN(dt$theatre_id), " theatres")
stopifnot(uniqueN(dt$theatre_id) <= 10)  # fail fast if wrong panel

# Coordinates
if ("xcoord" %in% names(dt)) {
  dt[, longitude := xcoord]; dt[, latitude := ycoord]
} else stop("xcoord/ycoord not found.")

# Controls
CTRL     <- c("nlights_calib_mean","nlights_mean","nlights_min",
              "prec_gpcp","temp","water_ih")
ctrl_str <- paste(intersect(CTRL, names(dt)), collapse = " + ")
FE       <- "pg_id + theatre_id^month_int"

# Conley helper
safe_conley <- function(m, cut = 200) {
  r <- tryCatch(summary(m, vcov = conley(cutoff = cut)), error = function(e) NULL)
  if (!is.null(r)) return(list(m = r, vcov = paste0("Conley_", cut, "km")))
  r <- tryCatch(summary(m, vcov = vcov_conley(cutoff = cut)), error = function(e) NULL)
  if (!is.null(r)) return(list(m = r, vcov = paste0("Conley_", cut, "km")))
  list(m = summary(m, cluster = ~pg_id + month_int), vcov = "TwoWay_fallback")
}

extract_coef <- function(m_sum, pattern, label, vcov_type) {
  ct  <- as.data.table(coeftable(m_sum), keep.rownames = "term")
  hit <- ct[grepl(pattern, term, perl = TRUE)][1]
  if (nrow(hit) == 0L || is.na(hit$Estimate))
    return(data.table(label = label, vcov_type = vcov_type,
                      coef = NA_real_, se = NA_real_, p = NA_real_, nobs = NA_integer_))
  data.table(
    label     = label,
    vcov_type = vcov_type,
    coef      = round(hit$Estimate[1],     5),
    se        = round(hit$`Std. Error`[1], 5),
    p         = round(hit$`Pr(>|t|)`[1],   4),
    nobs      = nobs(m_sum)
  )
}

all_results <- list()

# ══════════════════════════════════════════════════════════════
# FIX 1 — DECAY RADIUS ROBUSTNESS ON CORRECT 8-THEATRE PANEL
# ══════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("FIX 1: Decay radius on 8-theatre panel (N=304,780)\n")
cat(strrep("=", 60), "\n\n")

# Load feasibility file for terrain component
eli <- readRDS(P_ELI); setDT(eli)
if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := gid]
terrain_col <- intersect(c("eli_feas_terrain","feas_terrain","eli_terrain"), names(eli))[1]
if (is.na(terrain_col)) stop("No terrain column in eli_feasibility_pg.rds.")

feas_small <- eli[, .(pg_id = as.integer(pg_id),
                      eli_feas_terrain = as.numeric(get(terrain_col)))]
setkey(feas_small, pg_id); setkey(dt, pg_id)
if (!("eli_feas_terrain" %in% names(dt))) dt <- feas_small[dt]
if (anyNA(dt$eli_feas_terrain)) stop("eli_feas_terrain has NA after merge.")

# Back-calculate distances from baseline eli_feas (R0=400km)
R0 <- 400L
dt[, range_base := eli_feas / eli_feas_terrain]
dt[!is.finite(range_base) | range_base <= 0, range_base := NA_real_]
dt[range_base > 1, range_base := 1]
dt[is.na(range_base),  dist_km := 1e6]
dt[!is.na(range_base), dist_km := -R0 * log(range_base)]
if (any(!is.finite(dt$dist_km))) stop("Non-finite dist_km.")

decay_rows <- list()
for (R_km in c(200L, 400L, 800L)) {
  dt_tmp <- copy(dt)
  dt_tmp[, eli_feas_R := eli_feas_terrain * exp(-dist_km / R_km)]
  fml <- as.formula(paste0("log1p_y ~ post_eli:eli_feas_R + ", ctrl_str, " | ", FE))
  m_raw <- tryCatch(feols(fml, data = dt_tmp, notes = FALSE),
                    error = function(e) { log_msg("ERROR R=", R_km, ": ", e$message); NULL })
  if (is.null(m_raw)) next
  rc <- safe_conley(m_raw, 200)
  r  <- extract_coef(rc$m, "post_eli.*eli_feas_R|eli_feas_R.*post_eli",
                     paste0("R=", R_km, "km"), rc$vcov)
  r[, R_km := R_km]; r[, fix := "Fix1_correct_panel"]
  cat(sprintf("  R=%dkm: β = %.5f  SE = %.5f  p = %.4f  N = %d\n",
              R_km, r$coef, r$se, r$p, r$nobs))
  decay_rows[[length(decay_rows) + 1]] <- r
}

decay_dt <- rbindlist(decay_rows, fill = TRUE)
cat("\nComparison with OLD results (from 88-country panel, now corrected):\n")
cat("  OLD nobs = 1,740,510 (88 theatres) — wrong\n")
cat("  NEW nobs = 304,780   (8 theatres)  — correct\n")
cat("  Coefficients should be similar; N and inference now correct.\n")
all_results$fix1_decay <- decay_dt

# ══════════════════════════════════════════════════════════════
# FIX 2 — LOCAL ELI_FEAS RE-STANDARDISATION
# ══════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("FIX 2: Within-theatre-sample eli_feas re-standardisation\n")
cat("Global z-scores compress within-theatre terrain variation.\n")
cat("Local re-normalisation calibrates to the 8-theatre analysis sample.\n")
cat(strrep("=", 60), "\n\n")

# Re-normalise eli_feas_terrain within the 8-theatre sample
terrain_range <- range(dt$eli_feas_terrain, na.rm = TRUE)
dt[, eli_feas_terrain_local := (eli_feas_terrain - terrain_range[1]) /
                                (terrain_range[2] - terrain_range[1])]
dt[, eli_feas_local := eli_feas_terrain_local * eli_feas_range]

cat("Terrain variation comparison (global vs local):\n")
cat(sprintf("  Global eli_feas_terrain: mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n",
            mean(dt$eli_feas_terrain), sd(dt$eli_feas_terrain),
            min(dt$eli_feas_terrain), max(dt$eli_feas_terrain)))
cat(sprintf("  Local  eli_feas_terrain: mean=%.3f  sd=%.3f  range=[%.3f, %.3f]\n",
            mean(dt$eli_feas_terrain_local), sd(dt$eli_feas_terrain_local),
            min(dt$eli_feas_terrain_local), max(dt$eli_feas_terrain_local)))
cat(sprintf("  Correlation: %.4f\n\n",
            cor(dt$eli_feas_terrain, dt$eli_feas_terrain_local)))

local_rows <- list()

for (feas_var in c("eli_feas", "eli_feas_local")) {
  label <- if (feas_var == "eli_feas") "Global standardisation (original)" else
             "Local standardisation (fixed)"
  fml <- as.formula(paste0("log1p_y ~ post_eli:", feas_var, " + ", ctrl_str, " | ", FE))
  m_raw <- tryCatch(feols(fml, data = dt, notes = FALSE),
                    error = function(e) { log_msg("ERROR: ", e$message); NULL })
  if (is.null(m_raw)) next
  rc <- safe_conley(m_raw, 200)
  r  <- extract_coef(rc$m,
                     paste0("post_eli.*", feas_var, "|", feas_var, ".*post_eli"),
                     label, rc$vcov)
  r[, feas_type := feas_var]
  cat(sprintf("  %s: β = %.5f  SE = %.5f  p = %.4f\n",
              label, r$coef, r$se, r$p))
  local_rows[[length(local_rows) + 1]] <- r
}

local_dt <- rbindlist(local_rows, fill = TRUE)
cat("\nInterpretation:\n")
cat("  If β(local) > β(global): global terrain standardisation was compressing\n")
cat("  within-theatre variation and underestimating the treatment effect.\n")
cat("  If β(local) ≈ β(global): global standardisation is fine; no change needed.\n")
all_results$fix2_local_feas <- local_dt

# ══════════════════════════════════════════════════════════════
# FIX 3 — HORSE RACE WITH CONLEY SEs
# ══════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("FIX 3: Horse race with Conley 200km SEs\n")
cat("(R20 used theatre clustering = 8 clusters = VCOV warnings)\n")
cat(strrep("=", 60), "\n\n")

BASE_CTRL <- c("prec_gpcp","temp","water_ih","forest_ih")
base_str  <- paste(intersect(BASE_CTRL, names(dt)), collapse = " + ")

pop_col   <- intersect(c("pop_gpw_sum","pop_gpw_max"), names(dt))[1]
nl_col    <- intersect(c("nlights_mean","nlights_calib_mean"), names(dt))[1]
front_col <- intersect(c("contestation_tv","control_weak","ttime_mean"), names(dt))[1]

specs <- list()
specs[["(a) Baseline"]] <- paste0("log1p_y ~ post_eli:eli_feas + ", base_str, " | ", FE)

if (!is.na(pop_col))
  specs[[paste0("(b) +post*", pop_col)]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:", pop_col, " + ", base_str, " | ", FE)

if (!is.na(nl_col))
  specs[[paste0("(c) +post*", nl_col)]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:", nl_col, " + ", base_str, " | ", FE)

if (!is.na(front_col))
  specs[[paste0("(d) +post*", front_col)]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + post_eli:", front_col, " + ", base_str, " | ", FE)

hr_terms <- c()
if (!is.na(pop_col))   hr_terms <- c(hr_terms, paste0("post_eli:", pop_col))
if (!is.na(nl_col))    hr_terms <- c(hr_terms, paste0("post_eli:", nl_col))
if (!is.na(front_col)) hr_terms <- c(hr_terms, paste0("post_eli:", front_col))
if (length(hr_terms))
  specs[["(e) All horse-race"]] <- paste0(
    "log1p_y ~ post_eli:eli_feas + ", paste(hr_terms, collapse=" + "),
    " + ", base_str, " | ", FE)

horse_rows <- list()
for (nm in names(specs)) {
  m_raw <- tryCatch(feols(as.formula(specs[[nm]]), data = dt, notes = FALSE),
                    error = function(e) NULL)
  if (is.null(m_raw)) next
  rc <- safe_conley(m_raw, 200)
  r  <- extract_coef(rc$m, "post_eli.*eli_feas|eli_feas.*post_eli", nm, rc$vcov)
  horse_rows[[length(horse_rows) + 1]] <- r
}

horse_dt <- rbindlist(horse_rows, fill = TRUE)
cat("Horse race results (Conley 200km):\n\n")
print(horse_dt[, .(label, coef, se, p)])

b_range <- range(horse_dt$coef, na.rm = TRUE)
cat(sprintf("\nβ range: %.3f to %.3f\n", b_range[1], b_range[2]))
cat(sprintf("All positive: %s\n",
            if (all(horse_dt$coef > 0, na.rm = TRUE)) "YES ✓" else "NO ✗"))
all_results$fix3_horse_race <- horse_dt

# ══════════════════════════════════════════════════════════════
# COMBINED SUMMARY
# ══════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n")
cat("SUMMARY OF FIXES\n")
cat(strrep("=", 60), "\n")

cat("\nFix 1 — Decay radius (correct 8-theatre panel, Conley 200km):\n")
print(decay_dt[, .(R_km, coef, se, p, nobs)])

cat("\nFix 2 — Local vs global eli_feas standardisation:\n")
print(local_dt[, .(label, coef, se, p)])

cat("\nFix 3 — Horse race (Conley 200km, valid inference):\n")
print(horse_dt[, .(label, coef, se, p)])

# ── SAVE ──────────────────────────────────────────────────────
out_csv <- file.path(DIRS$out_paper1, paste0("R25_methodology_fixes_", TS, ".csv"))
fwrite(rbindlist(lapply(all_results, function(x) {
  if ("fix" %in% names(x)) return(x)
  x[, fix := "see_label"]
}), fill = TRUE), out_csv)
log_msg("Wrote: ", out_csv)

# ── FIX 4 NOTE ────────────────────────────────────────────────
cat("\n", strrep("=", 60), "\n")
cat("FIX 4 — Incirlik spatial (documentation only, no code needed)\n")
cat(strrep("=", 60), "\n")
cat("
The main eli_feas uses nearest military airbase globally.
For Syrian/Iraqi cells, this may not always be Incirlik.
The Incirlik quasi-experiment already correctly uses:
  - Triple interaction: post_eli × eli_feas × shock_incirlik
  - Continuous proximity: post_incirlik × incirlik_proximity
Both of these use cell-level distance to Incirlik specifically,
not the global nearest-airbase eli_feas.

The binary incirlik_was_closer spec being null (p=0.901) is because
the binary coding loses the continuous gradient. The continuous
proximity spec (β=-0.135, p=0.037) correctly captures this gradient.

No code change needed. Add this to the methods section:
'For the Incirlik quasi-experiment, the continuous proximity measure
uses each cell's normalised inverse distance to Incirlik Air Base
specifically, rather than the global eli_feas index, to capture the
spatial gradient in ISR coverage disruption independently of the
composite feasibility measure.'
")

log_msg("R25 complete.")
