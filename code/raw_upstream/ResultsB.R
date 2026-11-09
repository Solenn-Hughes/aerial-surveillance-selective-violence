# ============================================================
# paper1_optionB_ALLCW_RUN_ALL.R
# RUNS EVERYTHING end-to-end to reproduce:
#  - CORE MODELS: post×exposure (PPML), post×exposure (logit), post×exposure (fat PPML)
#  - Binary robustness: post×near_asset
#  - Identification diagnostics (ctry_var + identifies count)
#  - Combined “Paper 1” model: post×exposure + post×ELI_geo_z + post×exposure×ELI_geo_z
#    + same model on identifying-country subset (DT_id)
#  - Macro H1 test: Dispersion ~ post×infra_mean_exposure + t + log1p(tot_events) | country FE
#
# IMPORTANT:
#  - This is Option B (assets proximity/intensity) + feasibility (ELI_geo_access).
#  - Uses ALL civil wars panel (no treated-country restriction).
#  - Designed to run once from a clean R session.
# ============================================================

rm(list = ls())

# -------------------------
# 0) CONFIG (EDIT THIS ONLY)
# -------------------------
PROJECT_ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
setwd(PROJECT_ROOT)

PANEL_RDS     <- "data/processed/panels/panel_pg_month_ALL_CIVIL_WARS_from_GED_FINAL.rds"
PG_STATIC_CSV <- "data/raw/priogrid/PRIO-GRID Static Variables - 2026-02-06.csv"

# Your Wikidata “template” (raw-ish pull) OR a curated list.
# If you already created eli_assets_used.csv and want to reuse it, set USE_EXISTING_ASSETS_USED=TRUE.
ASSETS_TEMPLATE_CSV <- "data/processed/eli/eli_assets_template.csv"
ASSETS_USED_CSV     <- "outputs/paper1_optionB_ALLCW/eli_assets_used.csv"

OUT_DIR <- "outputs/paper1_optionB_ALLCW"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Global “ELI era” cutoff
T0_POST <- as.Date("2011-01-01")

# Exposure definition
BUFFER_KM      <- 200
MAX_KM_TO_CW   <- 1500
CHUNK_SIZE     <- 500

# Filters
DROP_MOSTLY_WATER <- TRUE
WATER_THRESH <- 0.5

# Military-ish filter (Wikidata is noisy)
MIL_REGEX <- "(air\\s*base|airbase|air\\s*force\\s*base|\\bafb\\b|\\braf\\b|naval\\s*air\\s*station|\\bnas\\b|military\\s*air|army\\s*airfield|air\\s*station|\\bair\\s*base\\b)"

# Event-study window (not required for the results you posted, but kept if you want it)
ES_WIN <- 24

# Use existing precomputed assets_used file (if you want exact same list across reruns)
USE_EXISTING_ASSETS_USED <- FALSE

# Save artifacts
SAVE_CSV  <- TRUE
SAVE_TXT  <- TRUE

# -------------------------
# 1) Packages
# -------------------------
pkgs <- c("data.table","fixest","countrycode","stringi","geosphere")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

library(data.table)
library(fixest)
library(countrycode)
library(stringi)
library(geosphere)

# -------------------------
# 2) Helpers
# -------------------------
find_col <- function(nms, patterns){
  for (pat in patterns) {
    hit <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

ym_to_monthnum <- function(d){
  as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
}

get_exp_effect <- function(mod, term){
  b <- unname(coef(mod)[term])
  s <- unname(se(mod)[term])
  ci <- b + c(-1.96, 1.96) * s
  data.frame(
    term = term,
    beta = b,
    se = s,
    exp_effect = exp(b),
    ci_low = exp(ci[1]),
    ci_high = exp(ci[2])
  )
}

norm01 <- function(x){
  r <- range(x, na.rm = TRUE)
  if (!is.finite(r[1]) || !is.finite(r[2]) || r[1] == r[2]) return(rep(0, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

# Custom name fixes (only used if iso3_fe not already in panel)
CC_CUSTOM <- c(
  "Yemen (North Yemen)" = "Yemen"
)

# -------------------------
# 3) Load panel
# -------------------------
stopifnot(file.exists(PANEL_RDS))
DT <- as.data.table(readRDS(PANEL_RDS))

stopifnot(all(c("pg_id","month","events") %in% names(DT)))
DT[, month := as.Date(month)]
if (!"month_id" %in% names(DT)) DT[, month_id := ym_to_monthnum(month)]
if (!"any_event" %in% names(DT)) DT[, any_event := as.integer(events > 0)]
if (!"fat_best" %in% names(DT)) DT[, fat_best := 0]

# country key
if (!"ctry_key" %in% names(DT)) {
  if ("country_name" %in% names(DT)) DT[, ctry_key := as.character(country_name)]
  else stop("Panel must have country_name or ctry_key.", call.=FALSE)
}
DT[, ctry_key := stri_trim_both(as.character(ctry_key))]

# post
DT[, post := as.integer(month >= T0_POST)]

# Build iso3_fe ONLY if missing (so you reproduce exactly what you already had if it exists)
if (!"iso3_fe" %in% names(DT)) {
  DT[, iso3 := countrycode(ctry_key, "country.name", "iso3c", custom_match = CC_CUSTOM, warn = FALSE)]
  DT[, iso3_fe := fifelse(!is.na(iso3) & iso3 != "", iso3, ctry_key)]
}

cat("\n============================\n")
cat("PAPER 1 — OPTION B (ALL CW) RUN-ALL\n")
cat("============================\n\n")
cat("Rows:", nrow(DT), "\n")
cat("Cells:", uniqueN(DT$pg_id), "\n")
cat("Countries (iso3_fe):", uniqueN(DT$iso3_fe), "\n")
cat("Months:", uniqueN(DT$month), "\n")
cat("Month range:", as.character(min(DT$month)), "to", as.character(max(DT$month)), "\n")
cat("Share zero-event months:", mean(DT$events == 0), "\n\n")

# -------------------------
# 4) Merge PRIO static coords + mostly_water + (optional) ELI_geo_access
# -------------------------
stopifnot(file.exists(PG_STATIC_CSV))
pg <- fread(PG_STATIC_CSV)

if ("gid" %in% names(pg)) setnames(pg, "gid", "pg_id")
stopifnot("pg_id" %in% names(pg))

lon_col <- find_col(names(pg), c("^xcoord$","^lon$","longitude","^x$"))
lat_col <- find_col(names(pg), c("^ycoord$","^lat$","latitude","^y$"))
if (is.na(lon_col) || is.na(lat_col)) {
  stop("PRIO static must include lon/lat columns (e.g., xcoord/ycoord).", call.=FALSE)
}

# mostly_water
if (!"water_gc" %in% names(pg)) stop("PRIO static must include water_gc.", call.=FALSE)
pg[, mostly_water := as.integer(as.numeric(water_gc) >= WATER_THRESH)]

pg_coords <- unique(pg[, .(
  pg_id,
  lon = as.numeric(get(lon_col)),
  lat = as.numeric(get(lat_col)),
  mostly_water
)])
# Prevent suffix collisions BEFORE merge
for (v in c("lon","lat","mostly_water")) if (v %in% names(DT)) DT[, (v) := NULL]

# Merge coords
DT <- merge(DT, pg_coords, by = "pg_id", all.x = TRUE)
if (any(is.na(DT$lon) | is.na(DT$lat))) stop("Some pg_id missing coordinates after merge.", call.=FALSE)

# ---- FIX: resolve mostly_water suffix collisions robustly ----
mw_candidates <- c("mostly_water", "mostly_water.x", "mostly_water.y")

mw_found <- mw_candidates[mw_candidates %in% names(DT)]
if (length(mw_found) == 0) {
  stop("mostly_water not found after merge. Check PRIO static water_gc construction/merge.", call.=FALSE)
}

# Prefer non-missing if both exist; otherwise just rename the one that exists
if (all(c("mostly_water.x","mostly_water.y") %in% names(DT))) {
  DT[, mostly_water := fifelse(!is.na(mostly_water.x), mostly_water.x, mostly_water.y)]
  DT[, c("mostly_water.x","mostly_water.y") := NULL]
} else if ("mostly_water" %in% names(DT)) {
  # already fine
} else {
  # only one of the suffixed versions exists
  setnames(DT, mw_found[1], "mostly_water")
}

# ensure integer 0/1
DT[, mostly_water := as.integer(mostly_water)]
DT[is.na(mostly_water), mostly_water := 0L]


# Drop mostly-water if requested
if (DROP_MOSTLY_WATER) {
  DT <- DT[mostly_water == 0]
}

# ELI_geo_access: if missing, reconstruct using your earlier definition
if (!"ELI_geo_access" %in% names(DT)) {
  req_static <- c("forest_gc","mountains_mean","ttime_mean")
  miss <- setdiff(req_static, names(pg))
  if (length(miss) > 0) stop("Cannot reconstruct ELI_geo_access. Missing: ", paste(miss, collapse=", "), call.=FALSE)
  
  pg2 <- unique(pg[mostly_water == 0, .(pg_id, forest_gc, mountains_mean, ttime_mean)])
  pg2[, forest_n := norm01(as.numeric(forest_gc))]
  pg2[, mount_n  := norm01(as.numeric(mountains_mean))]
  pg2[, ttime_n  := norm01(as.numeric(ttime_mean))]
  pg2[, IMINT_UAV_access := (1 - forest_n) * (1 - mount_n)]
  pg2[, REMOTE_access    := (1 - ttime_n)]
  pg2[, ELI_geo_access   := IMINT_UAV_access * REMOTE_access]
  
  DT <- merge(DT, pg2[, .(pg_id, ELI_geo_access)], by = "pg_id", all.x = TRUE)
}

# Standardize feasibility at cell level
cell_f <- unique(DT[, .(pg_id, ELI_geo_access)])
cell_f[, ELI_geo_z := as.numeric(scale(ELI_geo_access))]
DT <- merge(DT, cell_f[, .(pg_id, ELI_geo_z)], by = "pg_id", all.x = TRUE)

# -------------------------
# 5) Build / load assets_used
# -------------------------
if (USE_EXISTING_ASSETS_USED) {
  stopifnot(file.exists(ASSETS_USED_CSV))
  A_use <- fread(ASSETS_USED_CSV)
} else {
  stopifnot(file.exists(ASSETS_TEMPLATE_CSV))
  A <- fread(ASSETS_TEMPLATE_CSV)
  
  stopifnot(all(c("asset_id","asset_name","country_iso3","lat","lon") %in% names(A)))
  
  A[, name_l := tolower(asset_name)]
  A[, wiki_l := tolower(if ("enwiki_url" %in% names(A)) enwiki_url else if ("enwiki" %in% names(A)) enwiki else "")]
  A[, mil_flag := as.integer(grepl(MIL_REGEX, name_l) | grepl(MIL_REGEX, wiki_l))]
  
  # civil-war country centroids (for proximity filter)
  cells0 <- unique(DT[, .(pg_id, lon, lat, iso3_fe)])
  # If iso3_fe mixes iso3 and names, the centroid filter still works using all iso3_fe groups:
  cw_cent <- cells0[, .(cent_lon = mean(lon), cent_lat = mean(lat)), by = iso3_fe]
  
  cent_mat <- as.matrix(cw_cent[, .(cent_lon, cent_lat)])
  A_mat <- as.matrix(A[, .(lon, lat)])
  
  min_km_to_cw <- rep(Inf, nrow(A))
  for (i in 1:nrow(cw_cent)) {
    d <- geosphere::distHaversine(A_mat, matrix(cent_mat[i,], nrow=nrow(A), ncol=2, byrow=TRUE)) / 1000
    min_km_to_cw <- pmin(min_km_to_cw, d)
  }
  A[, km_to_cw := min_km_to_cw]
  
  # Keep military-ish & close enough to any CW centroid
  A_use <- A[mil_flag == 1 & km_to_cw <= MAX_KM_TO_CW,
             .(asset_id, asset_name, country_iso3, lat, lon, enwiki_url = if ("enwiki_url" %in% names(A)) enwiki_url else wiki_l,
               km_to_cw, mil_flag)]
}

cat("Assets used rows:", nrow(A_use), "\n")
if (SAVE_CSV) fwrite(A_use, file.path(OUT_DIR, "eli_assets_used.csv"))

# -------------------------
# 6) Compute exposure (near_asset + exposure intensity)
# -------------------------
cells <- unique(DT[, .(pg_id, lon, lat, iso3_fe)])
cell_mat <- as.matrix(cells[, .(lon, lat)])
base_mat_all <- as.matrix(A_use[, .(lon, lat)])

n_cells <- nrow(cells)
n_bases <- nrow(A_use)

min_km   <- rep(Inf, n_cells)
n_within <- rep(0L, n_cells)

idx <- seq_len(n_bases)
chunks <- split(idx, ceiling(seq_along(idx) / CHUNK_SIZE))

cat("Computing exposure: ", n_cells, " cells × ", n_bases, " bases (", length(chunks), " chunks)\n")
for (ch in seq_along(chunks)) {
  idc <- chunks[[ch]]
  base_mat <- base_mat_all[idc, , drop = FALSE]
  
  Dkm <- geosphere::distm(cell_mat, base_mat, fun = geosphere::distHaversine) / 1000
  
  min_km <- pmin(min_km, apply(Dkm, 1, min))
  n_within <- n_within + rowSums(Dkm <= BUFFER_KM)
  
  if (ch %% 10 == 0 || ch == length(chunks)) cat("  chunk", ch, "/", length(chunks), "done\n")
}

cells[, min_km := min_km]
cells[, n_within := as.integer(n_within)]
cells[, near_asset := as.integer(min_km <= BUFFER_KM)]
cells[, exposure := log1p(n_within)]

if (SAVE_CSV) fwrite(cells[, .(pg_id, iso3_fe, lon, lat, min_km, n_within, near_asset, exposure)],
                     file.path(OUT_DIR, "cell_asset_exposure.csv"))

DT <- merge(DT, cells[, .(pg_id, near_asset, exposure)], by = "pg_id", all.x = TRUE)
DT[is.na(near_asset), near_asset := 0L]
DT[is.na(exposure), exposure := 0]

cat("Share near_asset==1 (cells):", mean(cells$near_asset == 1), "\n")
cat("Exposure summary:\n"); print(summary(cells$exposure)); cat("\n")

# -------------------------
# 7) CORE MODELS (the ones you pasted)
# -------------------------
cat("\n============================\n")
cat("CORE MODELS (cell FE + country×month FE)\n")
cat("============================\n")

m_events <- fepois(events ~ post:exposure | pg_id + iso3_fe^month_id,
                   data = DT, vcov = ~ pg_id)

m_any <- feglm(any_event ~ post:exposure | pg_id + iso3_fe^month_id,
               family = "binomial", data = DT, vcov = ~ pg_id)

m_fat <- fepois(fat_best ~ post:exposure | pg_id + iso3_fe^month_id,
                data = DT, vcov = ~ pg_id)

tab_core <- etable(m_events, m_any, m_fat)
print(tab_core)

cat("\nExponentiated effect (IRR/OR), post×exposure\n")
print(get_exp_effect(m_events, "post:exposure"))
print(get_exp_effect(m_any,    "post:exposure"))
print(get_exp_effect(m_fat,    "post:exposure"))

# Binary robustness (what you pasted)
m_events_bin <- fepois(events ~ post:near_asset | pg_id + iso3_fe^month_id,
                       data = DT, vcov = ~ pg_id)
cat("\nRobustness (binary): post×near_asset (PPML events)\n")
print(etable(m_events_bin))
print(get_exp_effect(m_events_bin, "post:near_asset"))

# -------------------------
# 8) Identification diagnostics (ctry_var + identifies count)
# -------------------------
cell_support <- unique(DT[, .(pg_id, iso3_fe, exposure, near_asset)])

ctry_var <- cell_support[, .(
  n_cells = .N,
  share_near = mean(near_asset == 1),
  has_near = any(near_asset == 1),
  has_far  = any(near_asset == 0),
  max_exposure = max(exposure),
  min_exposure = min(exposure)
), by = iso3_fe]

ctry_var[, identifies := (has_near & has_far & max_exposure > 0 & min_exposure == 0)]

cat("\nTOP 30 ctry_var\n")
print(ctry_var[order(-identifies, -share_near)][1:30])
cat("\n# countries that actually identify Option B:", sum(ctry_var$identifies), "\n\n")

if (SAVE_CSV) fwrite(ctry_var, file.path(OUT_DIR, "ctry_var_identification.csv"))

# -------------------------
# 9) Combined “Paper 1” model (what you pasted)
# -------------------------
m_comb <- fepois(
  events ~ post:exposure + post:ELI_geo_z + post:exposure:ELI_geo_z |
    pg_id + iso3_fe^month_id,
  data = DT, vcov = ~ pg_id
)

cat("\nCombined model (ALL CW)\n")
print(etable(m_comb))

id_ctry <- ctry_var[identifies == TRUE, iso3_fe]
DT_id <- DT[iso3_fe %in% id_ctry]

m_comb_id <- fepois(
  events ~ post:exposure + post:ELI_geo_z + post:exposure:ELI_geo_z |
    pg_id + iso3_fe^month_id,
  data = DT_id, vcov = ~ pg_id
)

cat("\nCombined model (identifying countries only)\n")
print(etable(m_comb, m_comb_id))

# -------------------------
# 10) Macro H1 test: Dispersion ~ post×infra_mean_exposure + trend + controls | country FE
# -------------------------
infra <- unique(DT[, .(pg_id, iso3_fe, exposure, near_asset)])
infra_ctry <- infra[, .(
  infra_mean_exposure = mean(exposure),
  infra_share_near    = mean(near_asset == 1)
), by = iso3_fe]

cm <- DT[, .(events), by = .(iso3_fe, month, month_id, pg_id, post)]
cm_tot <- cm[, .(tot = sum(events)), by = .(iso3_fe, month, month_id, post)]
cm <- merge(cm, cm_tot, by = c("iso3_fe","month","month_id","post"), all.x = TRUE)
cm <- cm[tot > 0]
cm[, share := events / tot]

out_cm <- cm[, .(
  tot_events = unique(tot),
  HHI = sum(share^2),
  Dispersion = 1 - sum(share^2)
), by = .(iso3_fe, month, month_id, post)][order(iso3_fe, month_id)]

out_cm <- merge(out_cm, infra_ctry, by = "iso3_fe", all.x = TRUE)
out_cm[, t := month_id - min(month_id)]

m_disp_H1 <- feols(
  Dispersion ~ post:infra_mean_exposure + t + log1p(tot_events) | iso3_fe,
  data = out_cm, vcov = "hetero"
)

cat("\nMacro H1 dispersion test\n")
print(summary(m_disp_H1))

if (SAVE_CSV) fwrite(out_cm, file.path(OUT_DIR, "country_month_dispersion.csv"))

# -------------------------
# 11) Save text snapshot
# -------------------------
if (SAVE_TXT) {
  capture.output({
    cat("PAPER 1 — OPTION B (ALL CW) RUN-ALL\n\n")
    cat("T0_POST:", as.character(T0_POST), "\n")
    cat("BUFFER_KM:", BUFFER_KM, " MAX_KM_TO_CW:", MAX_KM_TO_CW, "\n")
    cat("DROP_MOSTLY_WATER:", DROP_MOSTLY_WATER, " WATER_THRESH:", WATER_THRESH, "\n\n")
    
    cat("CORE MODELS\n")
    print(tab_core)
    cat("\nIRR/OR post×exposure\n")
    print(get_exp_effect(m_events, "post:exposure"))
    print(get_exp_effect(m_any,    "post:exposure"))
    print(get_exp_effect(m_fat,    "post:exposure"))
    
    cat("\nBinary robustness post×near_asset\n")
    print(get_exp_effect(m_events_bin, "post:near_asset"))
    
    cat("\nIdentification: # identifies\n")
    cat(sum(ctry_var$identifies), "\n\n")
    
    cat("\nCombined models\n")
    print(etable(m_comb, m_comb_id))
    
    cat("\nMacro dispersion H1\n")
    print(summary(m_disp_H1))
  }, file = file.path(OUT_DIR, "RUN_ALL_results_snapshot.txt"))
}

cat("\nDONE. Outputs in:", OUT_DIR, "\n")
