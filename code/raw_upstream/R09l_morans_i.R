# ============================================================
# R09l_morans_i.R  (FIXED)
# Spatial autocorrelation: Moran's I per theatre-month
#
# FIXES:
#   1) Fixed merge error: "y has duplicated column names"
#      — meta had theatre_id and month which duplicated moran_dt
#   2) Restricted to 8 treated theatres only (was computing
#      Moran's I for 80+ countries — slow and irrelevant)
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!requireNamespace("spdep", quietly = TRUE)) {
  stop("Package 'spdep' required. Install via: install.packages('spdep')")
}
suppressMessages(library(spdep))

if (!requireNamespace("fixest", quietly = TRUE)) stop("Need fixest.")
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

log_msg("Step 9l: Moran's I spatial autocorrelation analysis")

# ── LOAD PANEL + COORDINATES ─────────────────────────────────
P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
P_ELI   <- file.path(DIRS$eli, "eli_feasibility_pg.rds")

if (!file.exists(P_ALLCW)) stop("Missing: ", P_ALLCW)
if (!file.exists(P_ELI))   stop("Missing: ", P_ELI)

dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
dt[, pg_id := as.integer(pg_id)]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month) * 12L + month(month)]

# FIX: Restrict to 8 treated theatres only (not all 80+ countries)
dt <- dt[treated == 1L]
log_msg("Restricted to treated theatres: ", uniqueN(dt$theatre_id), " theatres, ", nrow(dt), " rows")

# Use xcoord/ycoord already in panel (avoid merge issues)
if ("xcoord" %in% names(dt) && "ycoord" %in% names(dt)) {
  dt[, lon := xcoord]
  dt[, lat := ycoord]
  log_msg("Using xcoord/ycoord from panel")
} else {
  eli <- readRDS(P_ELI); setDT(eli)
  if (!("pg_id" %in% names(eli)) && "gid" %in% names(eli)) eli[, pg_id := as.integer(gid)]
  coords <- unique(eli[, .(pg_id, lon, lat)])
  n_before <- nrow(dt)
  dt <- merge(dt, coords, by = "pg_id", all.x = TRUE)
  stopifnot(nrow(dt) == n_before)
}

n_no_coords <- sum(is.na(dt$lon) | is.na(dt$lat))
if (n_no_coords > 0) {
  log_msg("WARNING: ", n_no_coords, " rows missing coordinates. Dropping.")
  dt <- dt[!is.na(lon) & !is.na(lat)]
}

# ── BUILD SPATIAL WEIGHTS (per theatre) ──────────────────────
build_knn_weights <- function(lon, lat, k = 8) {
  coords_mat <- cbind(lon, lat)
  knn <- knearneigh(coords_mat, k = k)
  nb  <- knn2nb(knn)
  lw  <- nb2listw(nb, style = "W", zero.policy = TRUE)
  lw
}

theatres <- unique(dt$theatre_id)
theatre_weights <- list()
theatre_cells   <- list()

for (tid in theatres) {
  cells_t <- unique(dt[theatre_id == tid, .(pg_id, lon, lat)])
  setorder(cells_t, pg_id)
  
  if (nrow(cells_t) < 10) {
    log_msg("Skipping theatre ", tid, " — too few cells (", nrow(cells_t), ")")
    next
  }
  
  lw <- tryCatch(
    build_knn_weights(cells_t$lon, cells_t$lat, k = min(8, nrow(cells_t) - 1)),
    error = function(e) { log_msg("Weight build failed for ", tid, ": ", e$message); NULL }
  )
  
  if (!is.null(lw)) {
    theatre_weights[[tid]] <- lw
    theatre_cells[[tid]]   <- cells_t$pg_id
    log_msg("  Built weights for ", tid, ": ", nrow(cells_t), " cells")
  }
}

if (!length(theatre_weights)) stop("No theatres with valid spatial weights.")

# ── COMPUTE MORAN'S I PER THEATRE-MONTH ──────────────────────
log_msg("Computing Moran's I per theatre-month...")

tm_keys <- unique(dt[, .(theatre_id, month, month_int)])
setorder(tm_keys, theatre_id, month)

moran_rows <- list()
n_done <- 0L
n_total <- nrow(tm_keys)

for (i in seq_len(nrow(tm_keys))) {
  tid <- tm_keys$theatre_id[i]
  m_date <- tm_keys$month[i]
  
  if (!(tid %in% names(theatre_weights))) next
  
  lw       <- theatre_weights[[tid]]
  cell_ord <- theatre_cells[[tid]]
  
  sub <- dt[theatre_id == tid & month == m_date, .(pg_id, events_type3)]
  
  full <- data.table(pg_id = cell_ord)
  full <- merge(full, sub, by = "pg_id", all.x = TRUE)
  full[is.na(events_type3), events_type3 := 0L]
  setorder(full, pg_id)
  
  y <- full$events_type3
  
  if (sd(y) == 0) {
    moran_rows[[length(moran_rows) + 1]] <- data.table(
      theatre_id = tid, month = m_date, month_int = tm_keys$month_int[i],
      moran_I = NA_real_, moran_E = NA_real_, moran_p = NA_real_,
      sum_y = sum(y), n_cells = length(y), note = "zero variance"
    )
    next
  }
  
  mi <- tryCatch({
    moran.test(y, lw, alternative = "two.sided", zero.policy = TRUE)
  }, error = function(e) NULL)
  
  if (!is.null(mi)) {
    moran_rows[[length(moran_rows) + 1]] <- data.table(
      theatre_id = tid, month = m_date, month_int = tm_keys$month_int[i],
      moran_I = as.numeric(mi$estimate["Moran I statistic"]),
      moran_E = as.numeric(mi$estimate["Expectation"]),
      moran_p = as.numeric(mi$p.value),
      sum_y = sum(y), n_cells = length(y), note = "ok"
    )
  }
  
  n_done <- n_done + 1L
  if (n_done %% 200 == 0) log_msg("  Processed ", n_done, " / ", n_total, " theatre-months")
}

moran_dt <- rbindlist(moran_rows)
log_msg("Computed Moran's I for ", nrow(moran_dt), " theatre-months")
log_msg("  Of which defined (note='ok'): ", moran_dt[note == "ok", .N])

# ── MERGE POST_ELI AND RUN DiD ───────────────────────────────
# FIX: compute post_eli per theatre-month without creating duplicate columns
meta <- dt[, .(post_eli_meta = max(as.integer(post_eli), na.rm = TRUE)),
           by = .(theatre_id, month)]

moran_dt <- merge(moran_dt, meta, by = c("theatre_id", "month"), all.x = TRUE)
setnames(moran_dt, "post_eli_meta", "post_eli")

# Restrict to defined Moran's I
moran_def <- moran_dt[note == "ok" & !is.na(moran_I)]
log_msg("Theatre-months with defined Moran's I: ", nrow(moran_def))

# DiD: moran_I ~ post_eli | theatre_id + month_int
fml <- moran_I ~ post_eli | theatre_id + month_int
m_did <- feols(fml, data = moran_def, cluster = ~theatre_id,
               ssc = ssc(adj = TRUE, cluster.adj = TRUE, cluster.df = "min"),
               notes = FALSE)

cat("\n============================================\n")
cat("MORAN'S I DiD: Does spatial clustering change after ELI onset?\n")
cat("============================================\n")
print(summary(m_did))

ct <- as.data.table(coeftable(m_did), keep.rownames = "term")
hit <- ct[grepl("post_eli", term)]

result <- data.table(
  metric    = "Moran's I (global)",
  outcome   = "events_type3",
  coef      = if (nrow(hit)) round(hit$Estimate[1], 5) else NA,
  se        = if (nrow(hit)) round(hit$`Std. Error`[1], 5) else NA,
  t         = if (nrow(hit)) round(hit$`t value`[1], 3) else NA,
  p         = if (nrow(hit)) round(hit$`Pr(>|t|)`[1], 4) else NA,
  n         = nobs(m_did),
  fe        = "theatre + month",
  vcov      = "theatre-clustered + ssc"
)

cat("\n── Result ──\n")
print(result)

cat("\n── Interpretation ──\n")
if (!is.na(result$coef)) {
  if (result$coef < 0) {
    cat("Negative: Moran's I DECREASES after ELI onset → violence becomes\n")
    cat("LESS spatially clustered → consistent with H1 (dispersion).\n")
  } else {
    cat("Positive: Moran's I INCREASES after ELI onset → violence becomes\n")
    cat("MORE spatially clustered → inconsistent with simple H1 dispersion.\n")
    cat("This could indicate violence concentrating in newly accessible areas.\n")
  }
}

# Descriptive stats
cat("\n── Descriptive: Moran's I by pre/post ──\n")
desc <- moran_def[, .(
  mean_I  = round(mean(moran_I), 4),
  sd_I    = round(sd(moran_I), 4),
  median_I = round(median(moran_I), 4),
  n       = .N
), by = .(post_eli)]
print(desc)

write_audit_dt(result, "R09l_morans_i_did")
write_audit_dt(desc, "R09l_morans_i_descriptive")

out_path <- file.path(DIRS$out_paper1,
                      paste0("morans_i_theatre_month_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
fwrite(moran_dt, out_path)
log_msg("Saved Moran's I series: ", out_path)

log_msg("Step 9l complete.")