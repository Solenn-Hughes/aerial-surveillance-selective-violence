# ============================================================
# R15_modis_land_cover.R
# Build MODIS-based land openness measure and test as
# alternative terrain component for eli_feas.
#
# MODIS MCD12Q1 v6.1: yearly, 500m, IGBP classification
# IGBP classes reclassified into "open" vs "concealing":
#   Open (good IMINT):  Barren (16), Cropland (12,14),
#                       Grassland (10), Savanna (8,9),
#                       Urban (13), Snow/Ice (15)
#   Concealing (bad):   Forest (1-5), Shrubland (6,7),
#                       Wetland (11)
#   Water (0, 17):      excluded from denominator
#
# TWO APPROACHES:
#   A) AppEEARS (recommended — NASA web tool, no tile download)
#   B) Direct tile download + terra (if you prefer local processing)
#
# This script handles BOTH — it auto-detects which input you have.
#
# Run AFTER downloading MODIS data (see instructions below).
# ============================================================

# ╔══════════════════════════════════════════════════════════════╗
# ║  DOWNLOAD INSTRUCTIONS (DO THIS FIRST)                     ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                            ║
# ║  OPTION A: AppEEARS (RECOMMENDED — easiest)                ║
# ║                                                            ║
# ║  1. Register for NASA Earthdata:                           ║
# ║     https://urs.earthdata.nasa.gov/users/new               ║
# ║     (use institutional email, instant approval)            ║
# ║                                                            ║
# ║  2. Go to AppEEARS:                                        ║
# ║     https://appeears.earthdatacloud.nasa.gov/              ║
# ║     Log in with your Earthdata credentials                 ║
# ║                                                            ║
# ║  3. FIRST: Run Part 0 of this script (below) to generate  ║
# ║     your PRIO-GRID centroids CSV. This creates:            ║
# ║     data/raw/modis/prio_centroids_for_appeears.csv         ║
# ║                                                            ║
# ║  4. In AppEEARS, click "Extract" → "Point Sample"         ║
# ║     - Upload the CSV from step 3                           ║
# ║     - Under "Select the layers":                           ║
# ║       Search "MCD12Q1" → select "MCD12Q1.061"             ║
# ║       → tick "LC_Type1" (IGBP classification)              ║
# ║     - Start date: 01-01-2001                               ║
# ║     - End date: 12-31-2020                                 ║
# ║       (MODIS land cover is yearly; pick every few years    ║
# ║        or all years — it's fast to process)                ║
# ║     - Output: GeoTIFF & CSV (tick CSV)                     ║
# ║     - Click Submit                                         ║
# ║     - Wait for email notification (~10-30 min)             ║
# ║     - Download the CSV results                             ║
# ║                                                            ║
# ║  5. Save the downloaded CSV to:                            ║
# ║     data/raw/modis/appeears_MCD12Q1_output.csv             ║
# ║                                                            ║
# ║  6. Run the rest of this script.                           ║
# ║                                                            ║
# ║  NOTE: Point extraction gives you the land cover class     ║
# ║  at each centroid. For the FRACTION of open land within    ║
# ║  each 0.5° cell, use "Area Sample" instead of "Point      ║
# ║  Sample" — but this requires a GeoJSON polygon file.       ║
# ║  Part 0 below also generates the GeoJSON if you want      ║
# ║  area fractions (more accurate but slower to process).     ║
# ║                                                            ║
# ╚══════════════════════════════════════════════════════════════╝

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

log_msg("Step 15: MODIS land cover processing")

# ══════════════════════════════════════════════════════════════
# PART 0: Generate PRIO-GRID centroids for AppEEARS upload
# ══════════════════════════════════════════════════════════════
# Run this FIRST, BEFORE downloading from AppEEARS.
# It creates the CSV you upload to the AppEEARS web interface.

generate_appeears_input <- function() {
  log_msg("Generating PRIO-GRID centroids for AppEEARS...")
  
  # Get cell centroids from your panel
  P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
  if (!file.exists(P_ALLCW)) stop("Missing panel: ", P_ALLCW)
  
  dt <- readRDS(P_ALLCW); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  dt[, pg_id := as.integer(pg_id)]
  
  # Get unique cells with coordinates
  if ("xcoord" %in% names(dt) && "ycoord" %in% names(dt)) {
    cells <- unique(dt[, .(pg_id, Longitude = xcoord, Latitude = ycoord)])
  } else {
    # Reconstruct from pg_id
    cells <- unique(dt[, .(pg_id)])
    cells[, pg_row := ((pg_id - 1L) %/% 720L) + 1L]
    cells[, pg_col := ((pg_id - 1L) %% 720L) + 1L]
    # Centroid of each 0.5° cell
    cells[, Latitude  := (pg_row - 1L) * 0.5 - 90 + 0.25]
    cells[, Longitude := (pg_col - 1L) * 0.5 - 180 + 0.25]
    cells[, c("pg_row", "pg_col") := NULL]
  }
  
  # AppEEARS CSV format: ID, Latitude, Longitude
  appeears_csv <- cells[, .(ID = paste0("pg_", pg_id),
                             Latitude = Latitude,
                             Longitude = Longitude)]
  
  out_dir <- file.path(DIRS$data_raw, "modis")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  csv_path <- file.path(out_dir, "prio_centroids_for_appeears.csv")
  fwrite(appeears_csv, csv_path)
  log_msg("Wrote AppEEARS input: ", csv_path, " (", nrow(appeears_csv), " points)")
  
  # Also generate a GeoJSON for area extraction (if user wants fractions)
  if (requireNamespace("sf", quietly = TRUE)) {
    library(sf)
    
    # Create 0.5° × 0.5° polygons for each cell
    polys <- lapply(seq_len(nrow(cells)), function(i) {
      lon <- cells$Longitude[i]
      lat <- cells$Latitude[i]
      # Cell boundaries
      xmin <- lon - 0.25; xmax <- lon + 0.25
      ymin <- lat - 0.25; ymax <- lat + 0.25
      st_polygon(list(matrix(
        c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin),
        ncol = 2, byrow = TRUE
      )))
    })
    
    sf_cells <- st_sf(
      ID = paste0("pg_", cells$pg_id),
      pg_id = cells$pg_id,
      geometry = st_sfc(polys, crs = 4326)
    )
    
    geojson_path <- file.path(out_dir, "prio_cells_for_appeears.geojson")
    st_write(sf_cells, geojson_path, delete_dsn = TRUE, quiet = TRUE)
    log_msg("Wrote GeoJSON for area extraction: ", geojson_path)
    log_msg("Use this for 'Area Sample' in AppEEARS for fraction-of-open-land.")
  } else {
    log_msg("sf package not installed. Only CSV (point) input generated.")
    log_msg("Install sf if you want area-fraction extraction: install.packages('sf')")
  }
  
  cat("\n============================================================\n")
  cat("NEXT STEPS:\n")
  cat("1. Go to https://appeears.earthdatacloud.nasa.gov/\n")
  cat("2. Click 'Extract' → 'Point Sample'\n")
  cat("3. Upload: ", csv_path, "\n")
  cat("4. Search 'MCD12Q1' → select MCD12Q1.061 → tick LC_Type1\n")
  cat("5. Date range: 2001-01-01 to 2020-12-31\n")
  cat("6. Submit, wait for email, download CSV results\n")
  cat("7. Save to: data/raw/modis/appeears_MCD12Q1_output.csv\n")
  cat("8. Re-run this script to process.\n")
  cat("============================================================\n")
}

# ══════════════════════════════════════════════════════════════
# IGBP CLASSIFICATION LOOKUP
# ══════════════════════════════════════════════════════════════
# IGBP class → open (1) / concealing (0) / water (NA)
igbp_to_open <- function(class_val) {
  # Returns 1 for open terrain, 0 for concealing, NA for water/unclassified
  dplyr_case <- function(x) {
    out <- rep(NA_real_, length(x))
    # Open terrain (good for IMINT/surveillance):
    out[x == 8L]  <- 1  # Woody Savannas
    out[x == 9L]  <- 1  # Savannas
    out[x == 10L] <- 1  # Grasslands
    out[x == 12L] <- 1  # Croplands
    out[x == 13L] <- 1  # Urban and Built-Up
    out[x == 14L] <- 1  # Cropland/Natural Vegetation Mosaic
    out[x == 15L] <- 1  # Snow and Ice
    out[x == 16L] <- 1  # Barren or Sparsely Vegetated
    # Concealing terrain (bad for IMINT):
    out[x == 1L]  <- 0  # Evergreen Needleleaf Forest
    out[x == 2L]  <- 0  # Evergreen Broadleaf Forest
    out[x == 3L]  <- 0  # Deciduous Needleleaf Forest
    out[x == 4L]  <- 0  # Deciduous Broadleaf Forest
    out[x == 5L]  <- 0  # Mixed Forest
    out[x == 6L]  <- 0  # Closed Shrublands
    out[x == 7L]  <- 0  # Open Shrublands
    out[x == 11L] <- 0  # Permanent Wetlands
    # Water (excluded):
    out[x == 0L]  <- NA  # Water
    out[x == 17L] <- NA  # Water Bodies
    out[x == 255L] <- NA # Unclassified/Fill
    out
  }
  dplyr_case(as.integer(class_val))
}

# ══════════════════════════════════════════════════════════════
# PART 1: DETECT AND PROCESS INPUT
# ══════════════════════════════════════════════════════════════

MODIS_DIR <- file.path(DIRS$data_raw, "modis")

# Check for AppEEARS output
appeears_paths <- c(
  file.path(MODIS_DIR, "appeears_MCD12Q1_output.csv"),
  file.path(MODIS_DIR, "MCD12Q1-061-results.csv"),
  # AppEEARS sometimes names files differently
  list.files(MODIS_DIR, pattern = "MCD12Q1.*\\.csv$", full.names = TRUE)
)
appeears_paths <- unique(appeears_paths)
appeears_path <- NULL
for (pp in appeears_paths) {
  if (file.exists(pp)) { appeears_path <- pp; break }
}

if (is.null(appeears_path)) {
  log_msg("No AppEEARS output found. Generating upload files...")
  generate_appeears_input()
  log_msg("After downloading from AppEEARS, re-run this script.")
  stop("Waiting for AppEEARS download. See instructions above.")
}

log_msg("Found AppEEARS output: ", appeears_path)

# ── READ APPEEARS OUTPUT ─────────────────────────────────────
ap <- fread(appeears_path, showProgress = TRUE)
log_msg("AppEEARS rows: ", nrow(ap))

# Standardise column names
setnames(ap, tolower(names(ap)))

# AppEEARS CSV typically has columns like:
# ID, Latitude, Longitude, Date, MCD12Q1_061_LC_Type1, Category
# The exact names vary — detect them

# Find the land cover column
lc_cols <- names(ap)[grepl("lc_type1|lc.type1|igbp|land_cover", names(ap), ignore.case = TRUE)]
if (!length(lc_cols)) {
  # Try any column that looks like MODIS data
  lc_cols <- names(ap)[grepl("mcd12", names(ap), ignore.case = TRUE)]
}
if (!length(lc_cols)) {
  cat("Available columns: ", paste(names(ap), collapse = ", "), "\n")
  stop("Cannot find land cover column in AppEEARS output. ",
       "Check column names above.")
}
lc_col <- lc_cols[1]
log_msg("Using land cover column: ", lc_col)

# Find the ID column
id_col <- names(ap)[grepl("^id$|^point_id|^sample_id", names(ap), ignore.case = TRUE)]
if (!length(id_col)) id_col <- names(ap)[1]
id_col <- id_col[1]

# Find the date column
date_col <- names(ap)[grepl("^date$|^calendar_date|^start_date", names(ap), ignore.case = TRUE)]
if (!length(date_col)) date_col <- names(ap)[grepl("date", names(ap), ignore.case = TRUE)][1]

log_msg("ID column: ", id_col, "  Date column: ", date_col)

# Parse
ap[, lc_class := as.integer(get(lc_col))]
ap[, is_open := igbp_to_open(lc_class)]

# Extract pg_id from the ID column
# Our IDs were formatted as "pg_12345"
ap[, pg_id := as.integer(gsub("pg_", "", get(id_col)))]

# Parse year
if (!is.null(date_col) && date_col %in% names(ap)) {
  ap[, year := year(as.IDate(get(date_col)))]
} else {
  ap[, year := 2010L]  # default if no date
}

# ── AGGREGATE: fraction open per cell ────────────────────────
# For point extraction, each cell has one class per year.
# We'll take the modal class across years, or use the most recent year.

# If we have multiple years, compute the mode
lc_summary <- ap[!is.na(is_open), .(
  frac_open = mean(is_open, na.rm = TRUE),
  modal_class = as.integer(names(sort(table(lc_class), decreasing = TRUE))[1]),
  n_years = uniqueN(year)
), by = pg_id]

log_msg("Land cover computed for ", nrow(lc_summary), " cells")
log_msg("Mean openness: ", round(mean(lc_summary$frac_open), 3))
log_msg("SD openness: ", round(sd(lc_summary$frac_open), 3))

# Add text labels for audit
igbp_labels <- c(
  "0" = "Water", "1" = "Evergreen Needleleaf",
  "2" = "Evergreen Broadleaf", "3" = "Deciduous Needleleaf",
  "4" = "Deciduous Broadleaf", "5" = "Mixed Forest",
  "6" = "Closed Shrubland", "7" = "Open Shrubland",
  "8" = "Woody Savanna", "9" = "Savanna",
  "10" = "Grassland", "11" = "Wetland",
  "12" = "Cropland", "13" = "Urban",
  "14" = "Cropland/Veg Mosaic", "15" = "Snow/Ice",
  "16" = "Barren", "17" = "Water Bodies"
)

lc_summary[, modal_label := igbp_labels[as.character(modal_class)]]

# Audit: distribution of modal land cover classes
class_dist <- lc_summary[, .N, by = .(modal_class, modal_label)][order(-N)]
cat("\n============================================================\n")
cat("LAND COVER DISTRIBUTION (modal class across years)\n")
cat("============================================================\n")
print(class_dist)
write_audit_dt(class_dist, "R15_modis_class_distribution")

# ── MERGE INTO PANEL ─────────────────────────────────────────
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

n_before <- nrow(dt)
dt <- merge(dt, lc_summary[, .(pg_id, modis_open = frac_open)],
            by = "pg_id", all.x = TRUE)
stopifnot(nrow(dt) == n_before)

# Fill NAs (cells without MODIS match — shouldn't happen but be safe)
n_na <- sum(is.na(dt$modis_open))
if (n_na > 0) {
  log_msg("MODIS NAs: ", n_na, " — filling with median")
  dt[is.na(modis_open), modis_open := median(dt$modis_open, na.rm = TRUE)]
}

# Compare with existing forest_gc
if ("forest_gc" %in% names(dt)) {
  # forest_gc = fraction forest, so 1 - forest_gc ≈ openness
  r <- cor(dt$modis_open, 1 - dt$forest_gc, use = "complete.obs")
  log_msg("Cor(modis_open, 1 - forest_gc): ", round(r, 3))
}

# ── BUILD MODIS-ENHANCED ELI_FEAS_TERRAIN ────────────────────
# Replace the PRIO-GRID forest_gc with MODIS openness in the
# terrain component. Original formula:
#   terrain_raw = (-z(mountains) + -z(forest) + -z(water)) / 3
#   eli_feas_terrain = scale01(terrain_raw)
#
# New: use modis_open directly (already 0-1, higher = more open)
# instead of -z(forest_gc)

z_safe <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE); sdv <- sd(x, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[2] == rng[1]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

# Get mountains and water from panel or PRIO-GRID static
if (all(c("mountains_mean", "water_gc") %in% names(dt))) {
  # Columns already in panel — build enhanced terrain
  terrain_modis <- (
    -z_safe(dt$mountains_mean) +   # less mountainous = more open
     z_safe(dt$modis_open)     +   # more open = better IMINT (positive!)
    -z_safe(dt$water_gc)           # less water = easier access
  ) / 3
  
  dt[, eli_feas_terrain_modis := scale01(terrain_modis)]
  
  log_msg("MODIS-enhanced terrain component built.")
  cat("  Cor with original eli_feas_terrain:",
      round(cor(dt$eli_feas_terrain_modis,
                dt$eli_feas_terrain, use = "complete.obs"), 3), "\n")
} else {
  # Need to load from PRIO static
  P_PRIO <- file.path(DIRS$base, "priogrid_static.rds")
  if (file.exists(P_PRIO)) {
    pg <- readRDS(P_PRIO); setDT(pg)
    if ("gid" %in% names(pg)) pg[, pg_id := as.integer(gid)]
    
    pg_merge <- unique(pg[, .(pg_id, mountains_mean, water_gc)])
    dt <- merge(dt, pg_merge, by = "pg_id", all.x = TRUE,
                suffixes = c("", ".prio"))
    
    terrain_modis <- (
      -z_safe(dt$mountains_mean) +
       z_safe(dt$modis_open) +
      -z_safe(dt$water_gc)
    ) / 3
    
    dt[, eli_feas_terrain_modis := scale01(terrain_modis)]
    log_msg("MODIS terrain component built from PRIO-GRID static.")
  } else {
    log_msg("WARNING: Cannot find mountains/water data. Using modis_open directly.")
    dt[, eli_feas_terrain_modis := modis_open]
  }
}

# Build full MODIS-enhanced eli_feas
if ("eli_feas_range" %in% names(dt)) {
  dt[, eli_feas_modis := eli_feas_terrain_modis * eli_feas_range]
} else {
  P_ELI <- file.path(DIRS$eli, "eli_feasibility_pg.rds")
  eli <- readRDS(P_ELI); setDT(eli)
  if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := as.integer(gid)]
  n_before <- nrow(dt)
  dt <- merge(dt, unique(eli[, .(pg_id, eli_feas_range)]),
              by = "pg_id", all.x = TRUE, suffixes = c("", ".e2"))
  stopifnot(nrow(dt) == n_before)
  dt[, eli_feas_modis := eli_feas_terrain_modis * eli_feas_range]
}

# ── RUN COMPARISON MODELS ────────────────────────────────────
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

dt[, log1p_y := log1p(events_type3)]

CTRL <- c("nlights_calib_mean", "nlights_mean", "nlights_min",
          "prec_gpcp", "temp", "water_ih")
ctrl_present <- intersect(CTRL, names(dt))
ctrl_str <- paste(ctrl_present, collapse = " + ")

# M0: Original eli_feas
f0 <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m0 <- feols(f0, data = dt, cluster = ~theatre_id, notes = FALSE)

# M1: MODIS-enhanced eli_feas
f1 <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas_modis + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m1 <- feols(f1, data = dt, cluster = ~theatre_id, notes = FALSE)

# M2: Horse race — both
f2 <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + post_eli:eli_feas_modis + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m2 <- feols(f2, data = dt, cluster = ~theatre_id, notes = FALSE)

# M3: Just MODIS openness directly (no range component)
f3 <- as.formula(paste0(
  "log1p_y ~ post_eli:eli_feas + post_eli:modis_open + ", ctrl_str,
  " | pg_id + theatre_id^month_int"
))
m3 <- feols(f3, data = dt, cluster = ~theatre_id, notes = FALSE)

# ── EXTRACT AND REPORT ───────────────────────────────────────
extract_hit <- function(m, pattern) {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  ct[grepl(pattern, term)]
}

cat("\n============================================================\n")
cat("MODIS LAND COVER ROBUSTNESS\n")
cat("============================================================\n")

cat("\nM0: Baseline (original eli_feas)\n")
h0 <- extract_hit(m0, "post_eli.*eli_feas|eli_feas.*post_eli")
if (nrow(h0)) print(h0[, .(term, Est = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

cat("\nM1: MODIS-enhanced eli_feas\n")
h1 <- extract_hit(m1, "post_eli.*eli_feas_modis|eli_feas_modis.*post_eli")
if (nrow(h1)) print(h1[, .(term, Est = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

cat("\nM2: Horse race (original + MODIS)\n")
h2 <- extract_hit(m2, "post_eli|eli_feas")
if (nrow(h2)) print(h2[, .(term, Est = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

cat("\nM3: Original eli_feas + raw MODIS openness\n")
h3 <- extract_hit(m3, "post_eli|eli_feas|modis")
if (nrow(h3)) print(h3[, .(term, Est = round(Estimate, 5),
                           SE = round(`Std. Error`, 5), p = round(`Pr(>|t|)`, 4))])

results <- data.table(
  model = c("M0: Original eli_feas",
            "M1: MODIS-enhanced eli_feas",
            "M2: Horse race (original wins?)",
            "M3: Original + raw openness"),
  coef = c(h0$Estimate[1], h1$Estimate[1],
           h2$Estimate[1], h3$Estimate[1]),
  se = c(h0$`Std. Error`[1], h1$`Std. Error`[1],
         h2$`Std. Error`[1], h3$`Std. Error`[1]),
  p = c(h0$`Pr(>|t|)`[1], h1$`Pr(>|t|)`[1],
        h2$`Pr(>|t|)`[1], h3$`Pr(>|t|)`[1])
)

cat("\n── Summary ──\n")
print(results[, .(model, coef = round(coef, 5), se = round(se, 5), p = round(p, 4))])

cat("\n── Interpretation ──\n")
cat("If M1 ≈ M0: MODIS adds little beyond PRIO-GRID forest_gc. Expected.\n")
cat("If M1 > M0: MODIS-enhanced terrain is a better IMINT proxy.\n")
cat("M2 horse race shows which terrain measure dominates.\n")
cat("M3 tests raw MODIS openness as additional control.\n")

write_audit_dt(results, "R15_modis_robustness")

log_msg("Step 15 complete.")
