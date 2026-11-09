# ============================================================
# R09e_test_H2_H3.R
# Tests H2 & H3 via triple interactions in DiD and ES.
# Writes audits to outputs/paper1/audit/
# ============================================================

rm(list=ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

assert_file <- function(path) { if (!file.exists(path)) stop("INTEGRITY FAIL: Missing file: ", path) }
assert_cols <- function(dt, cols, nm="dt") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("INTEGRITY FAIL: ", nm, " missing cols: ", paste(miss, collapse=", "))
}

# ---- EDIT THESE ----
H2_VAR <- "control_weak"      # e.g., contested, weak_control, frontline_dist (numeric)
H3_VAR <- "rebel_capability"  # e.g., rebel_cap, cap_index, cap_high
# --------------------

P_ALLCW <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
P_CTRL  <- file.path(DIRS$base,   "controls_pg_month.rds")
assert_file(P_ALLCW); assert_file(P_CTRL)

dt <- readRDS(P_ALLCW); setDT(dt)
if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month)*12L + month(month)]

req <- c("pg_id","month","month_int","treated","event_time","eli_feas","post_eli","theatre_id")
assert_cols(dt, req, "ALLCW_posteli_v2")

# H2/H3 vars must exist and be numeric
assert_cols(dt, c(H2_VAR, H3_VAR), "ALLCW_posteli_v2")
if (!is.numeric(dt[[H2_VAR]])) stop("INTEGRITY FAIL: H2_VAR not numeric: ", H2_VAR)
if (!is.numeric(dt[[H3_VAR]])) stop("INTEGRITY FAIL: H3_VAR not numeric: ", H3_VAR)

# Controls list (same approach you used)
ctrl <- readRDS(P_CTRL); setDT(ctrl)
DROP_NOT_CONTROLS <- c("pg_id","month","gid","gid.x","gid.y","month_id","month_int",
                       "country_id","country","lon","lat","xcoord","ycoord","row","col","year")
CONTROL_CAND <- setdiff(names(ctrl), DROP_NOT_CONTROLS)
CONTROL_CAND <- CONTROL_CAND[vapply(CONTROL_CAND, function(cc) is.numeric(ctrl[[cc]]), logical(1))]
CTRL_USE <- intersect(CONTROL_CAND, names(dt))
if (!length(CTRL_USE)) stop("INTEGRITY FAIL: No usable numeric controls intersected with panel.")

# ---- helpers
write_audit <- function(x, stem) write_audit_dt(as.data.table(x), stem)

extract_term <- function(m, term_regex) {
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  out <- ct[grepl(term_regex, term)]
  out
}

# Event-study setup (keeps controls in sample, drops treated outside window)
prep_es <- function(dt, outcome, window=60L) {
  if (!(outcome %in% names(dt))) stop("Missing outcome: ", outcome)
  if (any(dt[[outcome]] < 0, na.rm=TRUE)) stop("Outcome has negatives: ", outcome)
  x <- copy(dt)
  x[, ylog := log1p(get(outcome))]
  x[, event_time_es := as.integer(event_time)]
  x[treated == 0L, event_time_es := 0L]
  x[treated == 1L & (event_time_es < -24L | event_time_es > window), event_time_es := NA_integer_]
  x[, eli_feas_t := fifelse(treated == 1L, eli_feas, 0)]
  # Moderators interacted with feasibility (0 for never-treated because eli_feas_t=0)
  x[, h2_feas_t := eli_feas_t * get(H2_VAR)]
  x[, h3_feas_t := eli_feas_t * get(H3_VAR)]
  x
}

# Joint Wald on ES leads for a given ES interaction variable name
joint_wald_es <- function(m, xvar, lead_min=-24L, lead_max=-2L) {
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  pat <- paste0("^event_time_es::\\-\\d+:", xvar, "$")
  es <- ct[grepl(pat, term)]
  if (!nrow(es)) stop("INTEGRITY FAIL: No ES terms matched: ", pat)
  es[, k := as.integer(sub(paste0("^event_time_es::(-?\\d+):", xvar, "$"), "\\1", term))]
  leads <- es[k >= lead_min & k <= lead_max]
  if (!nrow(leads)) stop("INTEGRITY FAIL: No lead terms for ", xvar)
  
  terms <- leads$term
  b <- coef(m)[terms]
  V <- vcov(m)[terms, terms, drop=FALSE]
  
  # WARN not stop on singularity
  tmp <- tryCatch({
    stat <- as.numeric(t(b) %*% solve(V) %*% b)
    df <- length(b)
    p <- 1 - pchisq(stat, df=df)
    list(stat=stat, df=df, p=p, n=df, status="OK")
  }, error=function(e) {
    list(stat=NA_real_, df=length(b), p=NA_real_, n=length(b), status="WARN", note=e$message)
  })
  tmp
}

# -------------------------
# Run baseline family: outcome/events_ged, w60, both FE modes
# -------------------------
OUTCOME <- "events_ged"
WINDOW  <- 60L
FE_MODES <- c("global_month","theatre_month")

rows <- list()

for (fe_mode in FE_MODES) {
  
  fe_part <- if (fe_mode=="global_month") {
    "pg_id + month"
  } else {
    "pg_id + theatre_id^month_int"
  }
  
  rhs_ctrl <- paste(CTRL_USE, collapse=" + ")
  
  # ---- DiD with triple interactions (H2 & H3)
  f_did <- as.formula(paste0(
    "log1p(", OUTCOME, ") ~ post_eli*eli_feas + post_eli*eli_feas*",
    H2_VAR, " + post_eli*eli_feas*", H3_VAR, " + ", rhs_ctrl, " | ", fe_part
  ))
  m_did <- feols(f_did, data=dt, vcov=~pg_id + month_int)
  
  # Extract the moderation terms you care about
  # H2 term: post_eli:eli_feas:H2_VAR (any ordering)
  h2_term <- extract_term(m_did, paste0("post_eli:eli_feas:", H2_VAR, "|post_eli:", H2_VAR, ":eli_feas|eli_feas:post_eli:", H2_VAR))
  h3_term <- extract_term(m_did, paste0("post_eli:eli_feas:", H3_VAR, "|post_eli:", H3_VAR, ":eli_feas|eli_feas:post_eli:", H3_VAR))
  
  # ---- ES with moderation components
  x <- prep_es(dt, OUTCOME, window=WINDOW)
  
  f_es <- as.formula(paste0(
    "ylog ~ i(event_time_es, eli_feas_t, ref=-1) + ",
    "i(event_time_es, h2_feas_t, ref=-1) + ",
    "i(event_time_es, h3_feas_t, ref=-1) + ",
    rhs_ctrl, " | ", fe_part
  ))
  m_es <- feols(f_es, data=x, vcov=~pg_id + month_int)
  
  w_base <- joint_wald_es(m_es, "eli_feas_t")
  w_h2   <- joint_wald_es(m_es, "h2_feas_t")
  w_h3   <- joint_wald_es(m_es, "h3_feas_t")
  
  rows[[length(rows)+1]] <- data.table(
    outcome=OUTCOME, window=WINDOW, fe_mode=fe_mode,
    did_post_feas_est = extract_term(m_did, "post_eli:eli_feas|eli_feas:post_eli")$Estimate[1],
    did_post_feas_p   = extract_term(m_did, "post_eli:eli_feas|eli_feas:post_eli")$`Pr(>|t|)`[1],
    
    H2_var = H2_VAR,
    did_h2_est = if (nrow(h2_term)) h2_term$Estimate[1] else NA_real_,
    did_h2_p   = if (nrow(h2_term)) h2_term$`Pr(>|t|)`[1] else NA_real_,
    es_h2_leads_p = w_h2$p,
    
    H3_var = H3_VAR,
    did_h3_est = if (nrow(h3_term)) h3_term$Estimate[1] else NA_real_,
    did_h3_p   = if (nrow(h3_term)) h3_term$`Pr(>|t|)`[1] else NA_real_,
    es_h3_leads_p = w_h3$p,
    
    es_base_leads_p = w_base$p,
    nobs_did = nobs(m_did),
    nobs_es  = nobs(m_es)
  )
  
  # Save models so R99 isolation doesn't matter
  out_dir <- DIRS$out_paper1
  dir.create(file.path(out_dir, "models"), showWarnings=FALSE, recursive=TRUE)
  saveRDS(m_did, file.path(out_dir, "models", paste0("DID_H2H3_", OUTCOME, "_w", WINDOW, "_", fe_mode, ".rds")))
  saveRDS(m_es,  file.path(out_dir, "models", paste0("ES_H2H3_",  OUTCOME, "_w", WINDOW, "_", fe_mode, ".rds")))
}

chk <- rbindlist(rows)
print(chk)
write_audit_dt(chk, "R09e_H2_H3_checkpoint")
