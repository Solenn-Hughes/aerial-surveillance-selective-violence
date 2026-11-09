# ============================================================
# R09d_export_for_chatgpt_paper1.R
# Exports the MINIMUM outputs needed to assess whether Paper 1
# answers the theory: DiD + ES pretrends + dispersion ES.
# Writes audits to outputs/paper1/audit/ and prints paste-ready tables.
# ============================================================

rm(list = ls())

suppressMessages(library(data.table))
if (!requireNamespace("fixest", quietly = TRUE)) stop("Need fixest installed.")
suppressMessages(library(fixest))

# -------------------------
# 0) Load config (DIRS, helpers)
# -------------------------
cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found at R00_config.R or R/00_config.R")
source(cfg_path)

# Fallbacks if some helpers aren't in config
if (!exists("latest_file")) {
  latest_file <- function(dir, pattern) {
    ff <- list.files(dir, pattern = pattern, full.names = TRUE)
    if (!length(ff)) return(NA_character_)
    ff[which.max(file.info(ff)$mtime)]
  }
}
if (!exists("write_audit_dt")) {
  write_audit_dt <- function(dt, stem) {
    out_dir <- if (exists("DIRS") && "audit" %in% names(DIRS)) DIRS$audit else file.path("outputs","paper1","audit")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ts_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
    f <- file.path(out_dir, paste0(stem, "_", ts_tag, ".csv"))
    fwrite(dt, f)
    message(sprintf("[%s] Wrote audit: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), f))
    invisible(f)
  }
}
assert_file <- function(path) { if (!file.exists(path)) stop("INTEGRITY FAIL: Missing file: ", path) }
assert_cols <- function(dt, cols, nm="dt") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("INTEGRITY FAIL: ", nm, " missing cols: ", paste(miss, collapse = ", "))
}

# -------------------------
# 1) Paths (v2 panels + controls + dispersion panels)
# -------------------------
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds")
P_CTRL   <- file.path(DIRS$base,   "controls_pg_month.rds")

# Dispersion panels (R09b outputs)
D_DISP_ALLCW_TYPE3 <- file.path(DIRS$panels, "dispersion_theatre_month_ALLCW_v2_events_type3.rds")
D_DISP_ALLCW_GED   <- file.path(DIRS$panels, "dispersion_theatre_month_ALLCW_v2_events_ged.rds")

assert_file(P_ALLCW)
assert_file(P_STRICT)
assert_file(P_CTRL)

ALLCW  <- readRDS(P_ALLCW);  setDT(ALLCW)
STRICT <- readRDS(P_STRICT); setDT(STRICT)

# Harmonise ids
if (!("pg_id" %in% names(ALLCW)) && "gid" %in% names(ALLCW)) ALLCW[, pg_id := gid]
if (!("pg_id" %in% names(STRICT)) && "gid" %in% names(STRICT)) STRICT[, pg_id := gid]

# Month handling (uses your config helper if present)
if (exists("as_month_idate")) {
  ALLCW[, month := as_month_idate(month)]
  STRICT[, month := as_month_idate(month)]
} else {
  # fall back: try IDate coercion
  ALLCW[, month := as.IDate(month)]
  STRICT[, month := as.IDate(month)]
}

ALLCW[,  month_int := year(month)*12L + month(month)]
STRICT[, month_int := year(month)*12L + month(month)]

# Required cols for core models
req_core <- c("pg_id","month","month_int","treated","event_time","eli_feas","post_eli")
assert_cols(ALLCW,  req_core, "ALLCW_posteli_v2")
assert_cols(STRICT, req_core, "STRICT_posteli_v2")

# theatre_id needed for theatre_month FE
has_theatre <- "theatre_id" %in% names(ALLCW)
if (has_theatre) {
  ALLCW[, theatre_id := as.character(theatre_id)]
  STRICT[, theatre_id := as.character(theatre_id)]
}

# -------------------------
# 2) Controls list (numeric controls from controls_pg_month.rds)
# -------------------------
ctrl <- readRDS(P_CTRL); setDT(ctrl)

DROP_NOT_CONTROLS <- c(
  "pg_id","month","gid","gid.x","gid.y","month_id","month_int",
  "country_id","country","lon","lat","xcoord","ycoord","row","col","year"
)
CONTROL_CAND <- setdiff(names(ctrl), DROP_NOT_CONTROLS)
CONTROL_CAND <- CONTROL_CAND[vapply(CONTROL_CAND, function(cc) is.numeric(ctrl[[cc]]), logical(1))]
if (!length(CONTROL_CAND)) stop("INTEGRITY FAIL: No numeric controls found in controls_pg_month.rds after drops.")

# Use only controls actually present in the panel file
ctrls_ALLCW  <- intersect(CONTROL_CAND, names(ALLCW))
ctrls_STRICT <- intersect(CONTROL_CAND, names(STRICT))

if (!length(ctrls_ALLCW))  stop("INTEGRITY FAIL: No controls intersected with ALLCW panel columns.")
if (!length(ctrls_STRICT)) stop("INTEGRITY FAIL: No controls intersected with STRICT panel columns.")

# -------------------------
# 3) Model builders
# -------------------------

# Extract DiD term from fixest model
extract_term <- function(m, term1="post_eli:eli_feas", term2="eli_feas:post_eli") {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  hit <- ct[term %in% c(term1, term2)]
  if (!nrow(hit)) return(data.table(term=NA_character_, est=NA_real_, se=NA_real_, p=NA_real_))
  data.table(term=hit$term[1], est=hit$Estimate[1], se=hit$`Std. Error`[1], p=hit$`Pr(>|t|)`[1])
}

# Build event-study variables keeping controls in sample (never-treated retained)
prep_es_vars <- function(dt, outcome, window=60L) {
  if (!(outcome %in% names(dt))) stop("INTEGRITY FAIL: Missing outcome: ", outcome)
  if (any(dt[[outcome]] < 0, na.rm=TRUE)) stop("INTEGRITY FAIL: Outcome has negatives: ", outcome)
  
  out <- copy(dt)
  out[, ylog := log1p(get(outcome))]
  
  # Keep never-treated: event_time_es set to 0, but eli_feas_t set to 0 so it contributes nothing
  out[, event_time_es := as.integer(event_time)]
  out[treated == 0L, event_time_es := 0L]
  
  # Window trimming applies ONLY to treated observations (outside window dropped)
  out[treated == 1L & (event_time_es < -24L | event_time_es > window), event_time_es := NA_integer_]
  
  # Interaction: 0 for never-treated
  out[, eli_feas_t := fifelse(treated == 1L, eli_feas, 0)]
  out
}

# Extract ES coefficients and compute joint Wald on leads
extract_es_and_wald <- function(m, lead_min=-24L, lead_max=-2L) {
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  # Match your naming convention: event_time_es::-24:eli_feas_t
  es <- ct[grepl("^event_time_es::\\-?\\d+:eli_feas_t$", term)]
  if (!nrow(es)) stop("INTEGRITY FAIL: No ES terms found with pattern event_time_es::<k>:eli_feas_t")
  
  es[, k := as.integer(sub("^event_time_es::(-?\\d+):eli_feas_t$", "\\1", term))]
  setorder(es, k)
  
  # Leads to test
  leads <- es[k >= lead_min & k <= lead_max]
  if (!nrow(leads)) stop("INTEGRITY FAIL: No lead terms found in requested window.")
  
  # Wald with numeric guard (don’t stop on PSD; return NA with note)
  terms <- leads$term
  b <- coef(m)[terms]
  V <- vcov(m)[terms, terms, drop=FALSE]
  
  min_eig <- tryCatch(min(eigen((V + t(V))/2, symmetric=TRUE, only.values=TRUE)$values),
                      error=function(e) NA_real_)
  
  w <- tryCatch({
    stat <- as.numeric(t(b) %*% solve(V) %*% b)
    df <- length(b)
    p <- 1 - pchisq(stat, df=df)
    list(stat=stat, df=df, p=p, status="OK", note="", min_eig=min_eig)
  }, error=function(e) {
    list(stat=NA_real_, df=length(b), p=NA_real_, status="WARN",
         note=paste0("Wald failed (likely PSD/singularity): ", e$message), min_eig=min_eig)
  })
  
  list(es=es, w=w)
}

fit_did <- function(dt, outcome, fe_mode=c("global_month","theatre_month"), ctrls) {
  fe_mode <- match.arg(fe_mode)
  dt2 <- copy(dt)
  dt2[, ylog := log1p(get(outcome))]
  
  fe_part <- if (fe_mode == "global_month") {
    "pg_id + month"
  } else {
    if (!("theatre_id" %in% names(dt2))) stop("INTEGRITY FAIL: theatre_id missing for theatre_month.")
    "pg_id + theatre_id^month_int"
  }
  
  rhs_controls <- paste(ctrls, collapse=" + ")
  fml <- as.formula(paste0(
    "ylog ~ post_eli * eli_feas + ", rhs_controls, " | ", fe_part
  ))
  
  feols(fml, data=dt2, vcov=~pg_id + month_int)
}

fit_es <- function(dt, outcome, fe_mode=c("global_month","theatre_month"), window=60L, ctrls) {
  fe_mode <- match.arg(fe_mode)
  dt2 <- prep_es_vars(dt, outcome, window=window)
  
  fe_part <- if (fe_mode == "global_month") {
    "pg_id + month"
  } else {
    if (!("theatre_id" %in% names(dt2))) stop("INTEGRITY FAIL: theatre_id missing for theatre_month.")
    "pg_id + theatre_id^month_int"
  }
  
  rhs_controls <- paste(ctrls, collapse=" + ")
  fml <- as.formula(paste0(
    "ylog ~ i(event_time_es, eli_feas_t, ref = -1) + ", rhs_controls, " | ", fe_part
  ))
  
  feols(fml, data=dt2, vcov=~pg_id + month_int)
}

# -------------------------
# 4) Run the MINIMUM set: focus on the pretrend-clean family you reported:
#    ALLCW + events_ged + w60, both FE modes.
#    (Also optionally run w24 for comparison if you want.)
# -------------------------
TARGET_PANEL_NAME <- "ALLCW"
TARGET_DT <- ALLCW
TARGET_CTRLS <- ctrls_ALLCW

OUTCOME_MAIN <- "events_ged"
WINDOW_MAIN  <- 60L
FE_MODES <- c("global_month", if (has_theatre) "theatre_month")

did_rows <- list()
es_rows  <- list()

for (fe_mode in FE_MODES) {
  # DiD
  m_did <- fit_did(TARGET_DT, OUTCOME_MAIN, fe_mode=fe_mode, ctrls=TARGET_CTRLS)
  d <- extract_term(m_did)
  did_rows[[length(did_rows)+1]] <- data.table(
    panel=TARGET_PANEL_NAME, outcome=OUTCOME_MAIN, window=WINDOW_MAIN, fe_mode=fe_mode,
    term=d$term, did_est=d$est, did_se=d$se, did_p=d$p, nobs=nobs(m_did)
  )
  
  # ES + Wald
  m_es <- fit_es(TARGET_DT, OUTCOME_MAIN, fe_mode=fe_mode, window=WINDOW_MAIN, ctrls=TARGET_CTRLS)
  tmp <- extract_es_and_wald(m_es, lead_min=-24L, lead_max=-2L)
  w <- tmp$w
  es_rows[[length(es_rows)+1]] <- data.table(
    panel=TARGET_PANEL_NAME, outcome=OUTCOME_MAIN, window=WINDOW_MAIN, fe_mode=fe_mode,
    n_leads=w$df, wald_stat=w$stat, wald_df=w$df, wald_p=w$p, wald_status=w$status,
    vcov_min_eig=w$min_eig, wald_note=w$note, nobs=nobs(m_es)
  )
  
  # Save ES coefficients (full) to audit
  es_coef <- tmp$es
  es_coef[, `:=`(panel=TARGET_PANEL_NAME, outcome=OUTCOME_MAIN, window=WINDOW_MAIN, fe_mode=fe_mode)]
  write_audit_dt(es_coef, paste0("R09d_ES_coefs_", TARGET_PANEL_NAME, "_", OUTCOME_MAIN, "_w", WINDOW_MAIN, "_", fe_mode))
}

did_tab <- rbindlist(did_rows)
es_tab  <- rbindlist(es_rows)

write_audit_dt(did_tab, "R09d_DiD_main_required_for_chatgpt")
write_audit_dt(es_tab,  "R09d_ES_pretrend_required_for_chatgpt")

cat("\n============================\n")
cat("PASTE THIS: Main DiD results\n")
cat("============================\n")
print(did_tab)

cat("\n=========================================\n")
cat("PASTE THIS: ES joint leads pretrend tests\n")
cat("=========================================\n")
print(es_tab)

# -------------------------
# 5) Dispersion ES (conditional on sum_y > 0)
# -------------------------
fit_dispersion_es <- function(path_rds, outcome_label, window=60L) {
  assert_file(path_rds)
  d <- readRDS(path_rds); setDT(d)
  
  # expected (best guess) cols:
  # theatre_id, month, month_int, treated, event_time, dispersion, sum_y
  if (exists("as_month_idate")) d[, month := as_month_idate(month)] else d[, month := as.IDate(month)]
  d[, month_int := year(month)*12L + month(month)]
  
  req <- c("theatre_id","month","month_int","treated","event_time","dispersion","sum_y")
  assert_cols(d, req, paste0("dispersion_panel_", outcome_label))
  
  # conditional sample
  d2 <- d[sum_y > 0 & !is.na(dispersion)]
  if (!nrow(d2)) stop("INTEGRITY FAIL: dispersion sample sum_y>0 is empty for ", outcome_label)
  
  # Need a theatre-level eli_feas; if missing, set to theatre mean of cell-level feasibility (fallback)
  if (!("eli_feas" %in% names(d2))) {
    # fallback from ALLCW grid panel (time-invariant)
    if (!("theatre_id" %in% names(ALLCW)) || !("eli_feas" %in% names(ALLCW))) {
      stop("INTEGRITY FAIL: dispersion panel lacks eli_feas and cannot build it from ALLCW.")
    }
    feas_map <- ALLCW[!is.na(theatre_id), .(eli_feas = mean(eli_feas, na.rm=TRUE)), by=theatre_id]
    d2 <- merge(d2, feas_map, by="theatre_id", all.x=TRUE)
    if (anyNA(d2$eli_feas)) stop("INTEGRITY FAIL: Could not impute eli_feas for all theatres in dispersion panel.")
  }
  
  # ES vars (theatre-level)
  d2[, event_time_es := as.integer(event_time)]
  d2[treated == 0L, event_time_es := 0L]
  d2[treated == 1L & (event_time_es < -24L | event_time_es > window), event_time_es := NA_integer_]
  d2[, eli_feas_t := fifelse(treated == 1L, eli_feas, 0)]
  
  # FE at theatre-month level: theatre + month
  m <- feols(dispersion ~ i(event_time_es, eli_feas_t, ref=-1) | theatre_id + month,
             data=d2, vcov=~theatre_id + month_int)
  
  tmp <- extract_es_and_wald(m, lead_min=-24L, lead_max=-2L)
  w <- tmp$w
  
  list(
    chk = data.table(
      panel="ALLCW_v2", outcome=outcome_label, window=window,
      sample="sum_y>0", nobs=nobs(m),
      n_leads=w$df, wald_stat=w$stat, wald_df=w$df, wald_p=w$p,
      wald_status=w$status, vcov_min_eig=w$min_eig, wald_note=w$note
    ),
    coefs = tmp$es
  )
}

disp_rows <- list()

if (file.exists(D_DISP_ALLCW_GED)) {
  res <- fit_dispersion_es(D_DISP_ALLCW_GED, "dispersion_from_events_ged", window=WINDOW_MAIN)
  write_audit_dt(res$chk,  "R09d_dispersion_ES_pretrend_required_for_chatgpt")
  write_audit_dt(res$coefs, "R09d_dispersion_ES_coefs_required_for_chatgpt")
  disp_rows[[length(disp_rows)+1]] <- res$chk
}

if (file.exists(D_DISP_ALLCW_TYPE3)) {
  res <- fit_dispersion_es(D_DISP_ALLCW_TYPE3, "dispersion_from_events_type3", window=WINDOW_MAIN)
  write_audit_dt(res$chk,  "R09d_dispersion_ES_pretrend_required_for_chatgpt")
  write_audit_dt(res$coefs, "R09d_dispersion_ES_coefs_required_for_chatgpt")
  disp_rows[[length(disp_rows)+1]] <- res$chk
}

if (length(disp_rows)) {
  disp_tab <- rbindlist(disp_rows, fill=TRUE)
  cat("\n========================================\n")
  cat("PASTE THIS: Dispersion ES pretrend tests\n")
  cat("========================================\n")
  print(disp_tab)
} else {
  cat("\nNOTE: Dispersion RDS files not found; skipping dispersion ES.\n")
}

cat("\n✅ Done. Paste the printed tables above here.\n")
