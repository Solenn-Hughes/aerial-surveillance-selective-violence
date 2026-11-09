# ============================================================
# R09c_checkpoint_fixed_ES_controls.R
# - Builds ES with controls kept in sample (no treated-only ES)
# - Produces a clean 4-row chk table
# ============================================================

suppressMessages(library(data.table))
if (!requireNamespace("fixest", quietly=TRUE)) stop("Need fixest installed.")
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

# -------- helpers
load_panel_posteli_v2 <- function(path, label){
  if (!file.exists(path)) stop("Missing: ", path)
  dt <- readRDS(path); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  assert_cols(dt, c("pg_id","month","treated","event_time","eli_feas"), label)
  dt[, pg_id := as.integer(pg_id)]
  dt[, month := as_month_idate(month)]
  dt[, month_int := year(month)*12L + month(month)]
  dt
}

# Controls list (same logic as your R09)
P_CTRL <- file.path(DIRS$base, "controls_pg_month.rds")
if (!file.exists(P_CTRL)) stop("Missing controls file: ", P_CTRL)
ctrl <- readRDS(P_CTRL); setDT(ctrl)
DROP_NOT_CONTROLS <- c("pg_id","month","gid","gid.x","gid.y","month_id","month_int",
                       "country_id","country","lon","lat","xcoord","ycoord","row","col","year")
CONTROL_CAND <- setdiff(names(ctrl), DROP_NOT_CONTROLS)
CONTROL_CAND <- CONTROL_CAND[vapply(CONTROL_CAND, function(cc) is.numeric(ctrl[[cc]]), logical(1))]
if (!length(CONTROL_CAND)) stop("No numeric controls found after drops.")

# Wald test for joint leads = 0 in [-24,-2]
joint_wald_leads <- function(m, lead_min=-24L, lead_max=-2L){
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  # Terms look like event_time_es::-24:eli_feas_t (after our construction)
  ct <- ct[grepl("^event_time_es::\\-\\d+:eli_feas_t$", term)]
  if (!nrow(ct)) stop("No ES terms found to test (unexpected).")
  
  ct[, k := as.integer(sub("^event_time_es::(-?\\d+):eli_feas_t$", "\\1", term))]
  leads <- ct[k >= lead_min & k <= lead_max]
  if (!nrow(leads)) stop("No lead terms found in requested window (unexpected).")
  
  # Manual Wald: b' V^{-1} b
  b <- coef(m)[leads$term]
  V <- vcov(m)[leads$term, leads$term, drop=FALSE]
  stat <- as.numeric(t(b) %*% solve(V) %*% b)
  df <- length(b)
  p <- 1 - pchisq(stat, df=df)
  list(stat=stat, df=df, p=p, n=df)
}

# -------- fit ES (with controls retained)
fit_es_with_controls <- function(dt, outcome, fe_mode=c("global_month","theatre_month"), window=24L){
  fe_mode <- match.arg(fe_mode)
  
  if (!(outcome %in% names(dt))) stop("Missing outcome: ", outcome)
  if (any(dt[[outcome]] < 0, na.rm=TRUE)) stop("Outcome has negatives: ", outcome)
  
  dt <- copy(dt)
  dt[, ylog := log1p(get(outcome))]
  
  # KEY FIX: keep never-treated in sample
  dt[, event_time_es := as.integer(event_time)]
  dt[treated == 0L, event_time_es := 0L]  # arbitrary, but safe
  
  # window only applies to treated observations
  dt[treated == 1L & (event_time_es < -window | event_time_es > window), event_time_es := NA_integer_]
  
  # interaction term: 0 for never-treated
  dt[, eli_feas_t := fifelse(treated == 1L, eli_feas, 0)]
  
  # FE
  fe_part <- if (fe_mode == "global_month") {
    "pg_id + month"
  } else {
    # theatre_month requires theatre_id
    if (!("theatre_id" %in% names(dt))) stop("Need theatre_id for theatre_month FE.")
    "pg_id + theatre_id^month"
  }
  
  # Controls present
  ctrls <- intersect(CONTROL_CAND, names(dt))
  if (!length(ctrls)) stop("No usable controls intersected with panel columns (unexpected).")
  rhs_controls <- paste(ctrls, collapse=" + ")
  
  fml <- as.formula(paste0(
    "ylog ~ i(event_time_es, eli_feas_t, ref = -1) + ", rhs_controls, " | ", fe_part
  ))
  
  # NOTE: rows with NA event_time_es (treated outside window) will drop; controls stay.
  m <- feols(fml, data=dt, vcov=~pg_id + month_int)
  
  m
}

# -------- run for both panels/outcomes
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds")

ALLCW  <- load_panel_posteli_v2(P_ALLCW,  "ALLCW_posteli_v2")
STRICT <- load_panel_posteli_v2(P_STRICT, "STRICT_posteli_v2")

PANELS <- list(ALLCW=ALLCW, STRICT=STRICT)
OUTCOMES <- c("events_type3","events_ged")

rows <- list()

for (pn in names(PANELS)) {
  dt <- PANELS[[pn]]
  for (y in OUTCOMES) {
    for (fe_mode in c("global_month","theatre_month")) {
      
      # theatre_month needs theatre_id; if missing, skip safely
      if (fe_mode=="theatre_month" && !("theatre_id" %in% names(dt))) next
      
      m <- fit_es_with_controls(dt, y, fe_mode=fe_mode, window=24L)
      w <- joint_wald_leads(m, -24L, -2L)
      
      rows[[length(rows)+1]] <- data.table(
        panel = pn,
        outcome = y,
        fe_mode = fe_mode,
        window = 24L,
        n_leads = w$n,
        wald_stat = w$stat,
        wald_df = w$df,
        wald_p = w$p
      )
    }
  }
}

chk2 <- rbindlist(rows)
print(chk2)

# write audit
write_audit_dt(chk2, "R09c_checkpoint_fixed_ES_controls")
