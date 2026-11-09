# ============================================================
# R09b_dispersion_metrics.R  (UPDATED — SUNAB PARSING FIXED)
# H1 tests using theatre×month dispersion/spread metrics.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found: R00_config.R or R/00_config.R")
source(cfg_path)

AUDIT_DIR  <- if (exists("DIRS") && "audit" %in% names(DIRS)) DIRS$audit else file.path("outputs","paper1","audit")
MODELS_DIR <- file.path("outputs","paper1","models")
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODELS_DIR, recursive = TRUE, showWarnings = FALSE)

write_audit_dt0 <- function(dt, stem) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  f <- file.path(AUDIT_DIR, paste0(stem, "_", ts, ".csv"))
  fwrite(as.data.table(dt), f)
  message("Wrote audit: ", f)
  invisible(f)
}

assert_file <- function(p) if (!file.exists(p)) stop("INTEGRITY FAIL: Missing file: ", p)
assert_cols2 <- function(dt, cols, nm="dt") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("INTEGRITY FAIL: ", nm, " missing cols: ", paste(miss, collapse=", "))
}

add_month_int <- function(dt) {
  if (!("month_int" %in% names(dt))) dt[, month_int := year(month) * 12L + month(month)]
  dt
}

# -------------------------
# Settings
# -------------------------
RUN_STRICT <- FALSE
OUTCOMES <- c("events_type3", "events_ged")
FE_THEATRE_MONTH <- "theatre_id + month"
REF_P <- -1L
LEAD_MIN <- -24L
LEAD_MAX <- -2L

VCOV_MAIN <- "theatre"   # "theatre" or "twoway"
USE_SSC   <- TRUE
SSC_SPEC  <- ssc(adj = TRUE, cluster.adj = TRUE, cluster.df = "min")

DO_M0BOOT_SANITYCHECK <- TRUE
R_SANITY_KM <- 200L
R0_BASE_KM  <- 400L

vcov_formula <- function(mode) {
  if (mode == "theatre") return(~theatre_id)
  if (mode == "twoway")  return(~theatre_id + month_int)
  stop("Unknown VCOV mode: ", mode)
}

# -------------------------
# Load panels
# -------------------------
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds")
assert_file(P_ALLCW)
if (RUN_STRICT) assert_file(P_STRICT)

load_panel <- function(path, label) {
  dt <- readRDS(path); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
  dt <- add_month_int(dt)
  dt[, theatre_id := as.character(theatre_id)]
  assert_cols2(dt, c("pg_id","month","month_int","theatre_id","treated","post_eli","event_time"), label)
  dt[, treated := as.integer(treated)]
  dt[, post_eli := as.integer(post_eli)]
  dt
}

PANELS <- list(ALLCW = load_panel(P_ALLCW, "ALLCW_posteli_v2"))
if (RUN_STRICT) PANELS$STRICT <- load_panel(P_STRICT, "STRICT_posteli_v2")

# -------------------------
# Build theatre×month metrics
# -------------------------
build_theatre_month_metrics <- function(dt, y, panel_label) {
  
  if (!(y %in% names(dt))) return(NULL)
  if (!is.numeric(dt[[y]])) stop("INTEGRITY FAIL: outcome not numeric: ", y)
  if (any(dt[[y]] < 0, na.rm=TRUE)) stop("INTEGRITY FAIL: outcome has negatives: ", y)
  
  x <- dt[, .(theatre_id, month, month_int, treated, post_eli, event_time, y = get(y))]
  
  tm <- x[, .(
    sum_y = sum(y, na.rm=TRUE),
    active_cells = sum(y > 0, na.rm=TRUE),
    any_events = as.integer(sum(y, na.rm=TRUE) > 0)
  ), by=.(theatre_id, month, month_int, treated, post_eli, event_time)]
  
  # HHI/dispersion
  x2 <- x[, .(theatre_id, month, y)]
  tm_sum <- x2[, .(sum_y = sum(y, na.rm=TRUE)), by=.(theatre_id, month)]
  setkey(tm_sum, theatre_id, month)
  setkey(x2, theatre_id, month)
  x2 <- tm_sum[x2]
  x2[, share := fifelse(sum_y > 0, y / sum_y, NA_real_)]
  
  hhi <- x2[sum_y > 0, .(HHI = sum(share^2, na.rm=TRUE)), by=.(theatre_id, month)]
  hhi[, dispersion := 1 - HHI]
  
  setkey(hhi, theatre_id, month)
  setkey(tm, theatre_id, month)
  tm <- hhi[tm]
  
  tm[, `:=`(log1p_active_cells = log1p(active_cells),
            sum_y0 = as.integer(sum_y == 0))]
  
  tm[, g_onset_int := fifelse(treated == 1L & !is.na(event_time),
                              as.integer(month_int - as.integer(event_time)),
                              0L)]
  
  badg <- tm[treated == 1L, .(ng = uniqueN(g_onset_int[!is.na(g_onset_int)])), by=theatre_id][ng > 1]
  if (nrow(badg)) stop("INTEGRITY FAIL: g_onset_int not constant within treated theatre (", panel_label, ", ", y, ").")
  
  tm[]
}

# -------------------------
# Extractors
# -------------------------
extract_post <- function(m, term_name="post_eli") {
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  row <- ct[term == term_name]
  if (!nrow(row)) return(data.table(coef=NA_real_, se=NA_real_, t=NA_real_, p=NA_real_))
  row <- row[1]
  data.table(coef=as.numeric(row$Estimate),
             se=as.numeric(row$`Std. Error`),
             t=as.numeric(row$`t value`),
             p=as.numeric(row$`Pr(>|t|)`))
}

# ---- FIXED: robust parsing of k from fixest sunab term names
extract_sunab <- function(m) {
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  sun <- ct[grepl("sunab", term)]
  if (!nrow(sun)) return(NULL)
  
  # Parse k (relative time). fixest formats can vary; try multiple patterns.
  sun[, k := NA_integer_]
  
  # Pattern A: last parentheses (...) at end, e.g. "...(-12)"
  kA <- suppressWarnings(as.integer(sub(".*\\((-?\\d+)\\)\\s*$", "\\1", sun$term)))
  sun[!is.na(kA), k := kA[!is.na(kA)]]
  
  # Pattern B: after final ::, e.g. "...::(-12)" or "...::(-12)"
  idx <- is.na(sun$k)
  if (any(idx)) {
    kB <- suppressWarnings(as.integer(sub(".*::\\(?(-?\\d+)\\)?\\s*$", "\\1", sun$term[idx])))
    sun[idx][!is.na(kB), k := kB[!is.na(kB)]]
  }
  
  # If still NA, leave NA (caller will handle)
  sun
}

sunab_wald_leads <- function(m, lead_min=LEAD_MIN, lead_max=LEAD_MAX) {
  sun <- extract_sunab(m)
  if (is.null(sun) || !nrow(sun)) return(data.table(n_leads=0L, wald_stat=NA_real_, wald_df=NA_real_, wald_p=NA_real_))
  
  leads <- sun[!is.na(k) & k >= lead_min & k <= lead_max]
  if (!nrow(leads)) return(data.table(n_leads=0L, wald_stat=NA_real_, wald_df=NA_real_, wald_p=NA_real_))
  
  wt <- try(fixest::wald(m, leads$term), silent=TRUE)
  if (inherits(wt, "try-error")) return(data.table(n_leads=nrow(leads), wald_stat=NA_real_, wald_df=nrow(leads), wald_p=NA_real_))
  data.table(n_leads=nrow(leads), wald_stat=as.numeric(wt$stat), wald_df=as.numeric(wt$df), wald_p=as.numeric(wt$p.value))
}

sunab_slope_leads <- function(m, lead_min=LEAD_MIN, lead_max=LEAD_MAX) {
  sun <- extract_sunab(m)
  if (is.null(sun) || !nrow(sun)) return(data.table(slope=NA_real_, slope_se=NA_real_, slope_p=NA_real_, n=0L))
  
  leads <- sun[!is.na(k) & k >= lead_min & k <= lead_max]
  leads <- leads[!is.na(`Std. Error`) & `Std. Error` > 0]
  if (nrow(leads) < 5) return(data.table(slope=NA_real_, slope_se=NA_real_, slope_p=NA_real_, n=nrow(leads)))
  
  w <- 1/(leads$`Std. Error`^2)
  fit <- lm(Estimate ~ k, data=leads, weights=w)
  cf <- summary(fit)$coefficients
  data.table(slope=as.numeric(cf["k","Estimate"]),
             slope_se=as.numeric(cf["k","Std. Error"]),
             slope_p=as.numeric(cf["k","Pr(>|t|)"]),
             n=nrow(leads))
}

# -------------------------
# Run models
# -------------------------
run_models_for_metric <- function(tm, panel_name, y, metric) {
  
  d <- copy(tm)
  if (metric == "dispersion") d <- d[sum_y > 0 & !is.na(dispersion)]
  
  vc <- vcov_formula(VCOV_MAIN)
  
  f_did <- as.formula(paste0(metric, " ~ post_eli | ", FE_THEATRE_MONTH))
  f_es  <- as.formula(paste0(metric, " ~ sunab(g_onset_int, month_int, ref.p=", REF_P, ") | ", FE_THEATRE_MONTH))
  
  m_did <- if (USE_SSC) feols(f_did, data=d, vcov=vc, ssc=SSC_SPEC) else feols(f_did, data=d, vcov=vc)
  m_es  <- if (USE_SSC) feols(f_es,  data=d, vcov=vc, ssc=SSC_SPEC) else feols(f_es,  data=d, vcov=vc)
  
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saveRDS(m_did, file.path(MODELS_DIR, paste0("R09b3_", panel_name, "_", y, "_", metric, "_DID_", stamp, ".rds")))
  saveRDS(m_es,  file.path(MODELS_DIR, paste0("R09b3_", panel_name, "_", y, "_", metric, "_SUNAB_", stamp, ".rds")))
  
  did_row <- cbind(
    data.table(panel=panel_name, outcome=y, metric=metric, model="DID", vcov=VCOV_MAIN, ssc=USE_SSC, nobs=nobs(m_did)),
    extract_post(m_did, "post_eli")
  )
  
  wald  <- sunab_wald_leads(m_es)
  slope <- sunab_slope_leads(m_es)
  
  diag_row <- cbind(
    data.table(panel=panel_name, outcome=y, metric=metric, model="SUNAB", vcov=VCOV_MAIN, ssc=USE_SSC, nobs=nobs(m_es)),
    wald, slope
  )
  
  coefs <- extract_sunab(m_es)
  if (!is.null(coefs) && nrow(coefs)) coefs[, `:=`(panel=panel_name, outcome=y, metric=metric)]
  
  list(did=did_row, diag=diag_row, coefs=coefs)
}

# -------------------------
# Optional sanity check (R=200)
# -------------------------
m0boot_sanitycheck <- function(panel_dt, panel_name) {
  if (!DO_M0BOOT_SANITYCHECK) return(NULL)
  
  assert_cols2(panel_dt, c("pg_id","theatre_id","month_int","post_eli","eli_feas","events_type3"), paste0(panel_name, "_panel"))
  
  CTRL_M0 <- c("nlights_calib_mean","nlights_mean","nlights_min","prec_gpcp","temp","water_ih")
  miss <- setdiff(CTRL_M0, names(panel_dt))
  if (length(miss)) stop("INTEGRITY FAIL: missing m0_boot controls: ", paste(miss, collapse=", "))
  
  cand_feas <- c(
    file.path("data","processed","eli","eli_feasibility_pg.rds"),
    file.path("data","processed","eli","eli_feasibility_pg_v2.rds"),
    file.path(DIRS$eli, "eli_feasibility_pg.rds")
  )
  P_FEAS <- cand_feas[file.exists(cand_feas)][1]
  assert_file(P_FEAS)
  
  feas <- readRDS(P_FEAS); setDT(feas)
  if (!("pg_id" %in% names(feas)) && "gid" %in% names(feas)) feas[, pg_id := gid]
  terrain_col <- c("eli_feas_terrain","feas_terrain","eli_terrain")
  terrain_col <- terrain_col[terrain_col %in% names(feas)][1]
  if (is.na(terrain_col)) stop("INTEGRITY FAIL: feasibility file lacks terrain component column.")
  
  feas_small <- feas[, .(pg_id=as.integer(pg_id), eli_feas_terrain=as.numeric(get(terrain_col)))]
  setkey(feas_small, pg_id)
  x <- copy(panel_dt)
  setkey(x, pg_id)
  x <- feas_small[x]
  if (anyNA(x$eli_feas_terrain)) stop("INTEGRITY FAIL: NA after merging eli_feas_terrain for sanity check")
  
  x[, range_base := eli_feas / eli_feas_terrain]
  x[!is.finite(range_base) | range_base <= 0, range_base := NA_real_]
  x[range_base > 1, range_base := 1]
  x[is.na(range_base), dist_km := 1e6]
  x[!is.na(range_base), dist_km := -R0_BASE_KM * log(range_base)]
  
  x[, eli_feas_R := eli_feas_terrain * exp(-dist_km / R_SANITY_KM)]
  
  fml <- as.formula(paste0(
    "log1p(events_type3) ~ post_eli*eli_feas_R + ",
    paste(CTRL_M0, collapse=" + "),
    " | pg_id + theatre_id^month_int"
  ))
  
  m <- if (USE_SSC) feols(fml, data=x, vcov=~theatre_id, ssc=SSC_SPEC) else feols(fml, data=x, vcov=~theatre_id)
  
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  row <- ct[term %in% c("post_eli:eli_feas_R","eli_feas_R:post_eli")][1]
  
  out <- data.table(panel=panel_name, R_km=R_SANITY_KM,
                    coef=as.numeric(row$Estimate),
                    se=as.numeric(row$`Std. Error`),
                    t=as.numeric(row$`t value`),
                    p=as.numeric(row$`Pr(>|t|)`),
                    nobs=as.integer(nobs(m)))
  
  write_audit_dt0(out, "R09b3_sanitycheck_m0boot_R200")
  out
}

# -------------------------
# Run
# -------------------------
all_did <- list()
all_diag <- list()
all_coefs <- list()
all_summ <- list()
all_sanity <- list()

for (pn in names(PANELS)) {
  base <- PANELS[[pn]]
  
  if (pn == "ALLCW") {
    sc <- m0boot_sanitycheck(base, pn)
    if (!is.null(sc)) all_sanity[[length(all_sanity)+1]] <- sc
  }
  
  for (y in OUTCOMES) {
    tm <- build_theatre_month_metrics(base, y, pn)
    if (is.null(tm)) next
    
    summ <- tm[, .(
      panel = pn,
      outcome = y,
      n_theatre_month = .N,
      theatres = uniqueN(theatre_id),
      treated_theatres = uniqueN(theatre_id[treated==1L]),
      share_sumy0 = mean(sum_y==0),
      share_disp_na = mean(is.na(dispersion))
    )]
    all_summ[[length(all_summ)+1]] <- summ
    
    for (mtr in c("log1p_active_cells","any_events","dispersion")) {
      res <- run_models_for_metric(tm, pn, y, mtr)
      all_did[[length(all_did)+1]] <- res$did
      all_diag[[length(all_diag)+1]] <- res$diag
      if (!is.null(res$coefs) && nrow(res$coefs)) all_coefs[[length(all_coefs)+1]] <- res$coefs
    }
  }
}

summ_tab <- rbindlist(all_summ, fill=TRUE)
did_tab  <- rbindlist(all_did, fill=TRUE)
diag_tab <- rbindlist(all_diag, fill=TRUE)
coef_tab <- if (length(all_coefs)) rbindlist(all_coefs, fill=TRUE) else data.table()

write_audit_dt0(summ_tab, "R09b3_H1_dispersion_summary")
write_audit_dt0(did_tab,  "R09b3_H1_dispersion_did")
write_audit_dt0(diag_tab, "R09b3_H1_dispersion_sunab_diag")
if (nrow(coef_tab)) write_audit_dt0(coef_tab, "R09b3_H1_dispersion_sunab_coefs")

cat("\n=== H1 DISPERSION METRICS: DiD results ===\n")
print(did_tab)

cat("\n=== H1 DISPERSION METRICS: SUNAB diagnostics ===\n")
print(diag_tab)

if (length(all_sanity)) {
  cat("\n=== SANITY CHECK: m0_boot R=200 ===\n")
  print(rbindlist(all_sanity, fill=TRUE))
}

cat("\n✅ Done.\n")