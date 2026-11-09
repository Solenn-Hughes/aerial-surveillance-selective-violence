# ============================================================
# R09_models_event_study.R — Step 9 (DiD + Event Study)
#
# Inputs (Step 8 outputs):
#   data/processed/panels/panel_pg_month_ALLCW_paper1_posteli.rds
#   data/processed/panels/panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds
#   data/processed/base/controls_pg_month.rds   (for canonical control list)
#
# Outputs:
#   outputs/paper1/tables/   (tables txt + tex)
#   outputs/paper1/figures/  (event-study plots)
#   outputs/paper1/          (event-study coef CSVs)
#   outputs/paper1/audit/    (audit CSV)
#
# FAIL-FAST only on true integrity failures:
#   - missing inputs
#   - duplicate keys pg_id×month
#   - missing key vars: outcomes, eli_feas, post_eli, event_time
#   - eli_feas out of [0,1]
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

# -------------------------
# 0) CONFIG + PACKAGES
# -------------------------
cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found (R00_config.R or R/00_config.R).")
source(cfg_path)

if (!requireNamespace("fixest", quietly = TRUE)) {
  stop("Package 'fixest' is required for Step 9. Install via install.packages('fixest').")
}
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

# -------------------------
# 1) SETTINGS
# -------------------------
RUN_ALLCW  <- TRUE
RUN_STRICT <- TRUE

# Your known outcomes
OUTCOMES <- c("events_type3", "events_ged")

# Main models
RUN_MAIN_FEOLS_LOG1P <- TRUE
RUN_MAIN_FEPOIS      <- TRUE

# Event study
RUN_EVENT_STUDY <- TRUE
EVENT_WINDOW <- c(-24L, 60L)
ES_REF <- -1L

# Two-way clustering
CLUSTER_TWOWAY <- TRUE

# -------------------------
# 2) PATHS
# -------------------------
P_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli.rds")
P_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds")
P_CTRL   <- file.path(DIRS$base,   "controls_pg_month.rds")

if (RUN_ALLCW  && !file.exists(P_ALLCW))  stop(sprintf("Missing required input: %s", P_ALLCW))
if (RUN_STRICT && !file.exists(P_STRICT)) stop(sprintf("Missing required input: %s", P_STRICT))
if (!file.exists(P_CTRL)) stop(sprintf("Missing required controls file: %s", P_CTRL))

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_AUDIT <- file.path(DIRS$audit, paste0("R09_models_event_study_", TS, ".csv"))

# -------------------------
# 3) HELPERS
# -------------------------
assert_has <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(sprintf("%s missing cols: %s", label, paste(miss, collapse = ", ")))
}

add_month_int <- function(dt) {
  dt[, month_int := (year(month) * 12L + month(month))]
  dt
}

assert_nonneg <- function(dt, col, label) {
  v <- dt[[col]]
  if (!is.numeric(v)) stop(sprintf("%s: %s is not numeric.", label, col))
  if (any(v < 0, na.rm = TRUE)) stop(sprintf("%s: %s has negative values.", label, col))
}

write_etable <- function(models, path_txt, path_tex = NULL, title = NULL) {
  txt <- capture.output(etable(models, headers = title))
  writeLines(txt, con = path_txt, useBytes = TRUE)
  if (!is.null(path_tex)) {
    tex <- capture.output(etable(models, tex = TRUE, headers = title))
    writeLines(tex, con = path_tex, useBytes = TRUE)
  }
  invisible(TRUE)
}

extract_es <- function(m, outcome, panel_name) {
  ct <- as.data.table(coeftable(m))
  ct[, term := rownames(coeftable(m))]
  ct <- ct[grepl("^event_time_es::", term)]
  ct[, event_time := as.integer(sub("^event_time_es::(-?\\d+).*", "\\1", term))]
  ct[, `:=`(panel = panel_name, outcome = outcome)]
  setcolorder(ct, c("panel","outcome","term","event_time","Estimate","Std. Error","t value","Pr(>|t|)"))
  ct[order(event_time)]
}

save_es_plot <- function(m, file_path, main_title) {
  png(file_path, width = 1200, height = 700, res = 150)
  on.exit(dev.off(), add = TRUE)
  iplot(m, ref.line = 0, main = main_title, xlab = "Event time (months)", ylab = "Effect (relative to ref)")
  invisible(TRUE)
}

# -------------------------
# 4) CONTROL LIST (robust, ignores merge junk)
# -------------------------
ctrl <- readRDS(P_CTRL); setDT(ctrl)
if (!("pg_id" %in% names(ctrl)) && "gid" %in% names(ctrl)) ctrl[, pg_id := gid]
assert_has(ctrl, c("pg_id","month"), "controls_pg_month")

# Drop obvious non-controls / merge artifacts
DROP_NOT_CONTROLS <- c(
  "pg_id","month","gid","gid.x","gid.y","month_id","month_int",
  "country_id","country","lon","lat","xcoord","ycoord",
  "row","col","year"   # <- add these
)

CONTROL_CAND <- setdiff(names(ctrl), DROP_NOT_CONTROLS)

# Keep only numeric controls (fixest will choke on stray character cols)
CONTROL_CAND <- CONTROL_CAND[vapply(CONTROL_CAND, function(cc) is.numeric(ctrl[[cc]]), logical(1))]

if (length(CONTROL_CAND) == 0) stop("No numeric control candidates found in controls_pg_month.rds after dropping junk.")

# -------------------------
# 5) LOAD PANELS + INTEGRITY (do NOT fail on missing junk controls)
# -------------------------
load_panel <- function(path, panel_name) {
  dt <- readRDS(path); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  dt[, pg_id := as.integer(pg_id)]
  dt[, month := as_month_idate(month)]
  assert_unique(dt, c("pg_id","month"), name = paste0(panel_name, " keys"))
  
  # Required core vars for modeling
  core <- c("eli_feas","post_eli","event_time")
  assert_has(dt, core, panel_name)
  
  # Bounds / types
  if (anyNA(dt$eli_feas)) stop(sprintf("%s: eli_feas has NA.", panel_name))
  if (any(dt$eli_feas < 0 | dt$eli_feas > 1)) stop(sprintf("%s: eli_feas out of [0,1].", panel_name))
  if (anyNA(dt$post_eli)) stop(sprintf("%s: post_eli has NA.", panel_name))
  if (any(!dt$post_eli %in% c(0L,1L))) stop(sprintf("%s: post_eli not binary.", panel_name))
  
  dt <- add_month_int(dt)
  dt
}

panels <- list()
if (RUN_ALLCW)  panels$ALLCW  <- load_panel(P_ALLCW,  "ALLCW_paper1_posteli")
if (RUN_STRICT) panels$STRICT <- load_panel(P_STRICT, "STRICT_paper1_posteli")

# Decide controls USED per panel: only those present in dt
controls_used <- lapply(panels, function(dt) intersect(CONTROL_CAND, names(dt)))
# Integrity: must use at least 1 control
if (any(vapply(controls_used, length, integer(1)) == 0L)) {
  stop("A panel ended up with zero usable controls after intersection (unexpected).")
}

# Controls missing (informational)
controls_missing <- lapply(panels, function(dt) setdiff(CONTROL_CAND, names(dt)))

# -------------------------
# 6) MODEL RUNNERS
# -------------------------
run_models_for_panel <- function(dt, panel_name, CONTROL_COLS) {
  results <- list(models_main = list(), models_es = list(), es_coefs = list(), audit = list())
  
  vcov_arg <- if (CLUSTER_TWOWAY) ~ pg_id + month_int else ~ pg_id
  
  rhs_controls <- paste(CONTROL_COLS, collapse = " + ")
  fe_part <- "pg_id + month"
  
  results$audit$panel_summary <- data.table(
    panel = panel_name,
    rows = nrow(dt),
    pg = uniqueN(dt$pg_id),
    months = uniqueN(dt$month),
    share_post = mean(dt$post_eli == 1L),
    eli_min = min(dt$eli_feas),
    eli_max = max(dt$eli_feas),
    n_controls_used = length(CONTROL_COLS)
  )
  
  for (y in OUTCOMES) {
    if (!(y %in% names(dt))) stop(sprintf("%s: missing outcome '%s'.", panel_name, y))
    assert_nonneg(dt, y, panel_name)
    
    if (RUN_MAIN_FEOLS_LOG1P) {
      ylog <- paste0("log1p_", y)
      if (!(ylog %in% names(dt))) dt[, (ylog) := log1p(get(y))]
      fml  <- as.formula(paste0(ylog, " ~ post_eli:eli_feas + ", rhs_controls, " | ", fe_part))
      results$models_main[[paste0(y, "_feols_log1p")]] <- feols(fml, data = dt, vcov = vcov_arg)
    }
    
    if (RUN_MAIN_FEPOIS) {
      fmlp <- as.formula(paste0(y,    " ~ post_eli:eli_feas + ", rhs_controls, " | ", fe_part))
      results$models_main[[paste0(y, "_fepois")]] <- fepois(fmlp, data = dt, vcov = vcov_arg)
    }
    
    if (RUN_EVENT_STUDY) {
      dt[, event_time_es := as.integer(event_time)]
      dt[event_time_es < EVENT_WINDOW[1] | event_time_es > EVENT_WINDOW[2], event_time_es := NA_integer_]
      dt_es <- dt[!is.na(event_time_es)]
      
      if (nrow(dt_es) == 0) stop(sprintf("%s: ES sample empty after windowing for %s.", panel_name, y))
      if (dt_es[event_time_es == ES_REF, .N] == 0) {
        stop(sprintf("%s: ES ref %d not present in window [%d,%d] for %s.",
                     panel_name, ES_REF, EVENT_WINDOW[1], EVENT_WINDOW[2], y))
      }
      
      ylog <- paste0("log1p_", y)
      if (!(ylog %in% names(dt_es))) dt_es[, (ylog) := log1p(get(y))]
      
      fml_es <- as.formula(paste0(
        ylog, " ~ i(event_time_es, eli_feas, ref = ", ES_REF, ") + ",
        rhs_controls, " | ", fe_part
      ))
      mes <- feols(fml_es, data = dt_es, vcov = vcov_arg)
      
      key <- paste0(y, "_eventstudy_feols_log1p")
      results$models_es[[key]] <- mes
      results$es_coefs[[key]] <- extract_es(mes, y, panel_name)
      
      fig_path <- file.path(DIRS$figures, paste0("event_study_", panel_name, "_", y, "_", TS, ".png"))
      save_es_plot(mes, fig_path, paste0("Event study: ", panel_name, " — ", y, " (log1p)"))
      
      results$audit[[paste0("es_", y)]] <- data.table(
        panel = panel_name,
        outcome = y,
        es_rows = nrow(dt_es),
        es_pg = uniqueN(dt_es$pg_id),
        es_months = uniqueN(dt_es$month),
        es_event_min = min(dt_es$event_time_es),
        es_event_max = max(dt_es$event_time_es),
        es_ref = ES_REF
      )
    }
  }
  
  results
}

all_results <- list()
for (nm in names(panels)) {
  all_results[[nm]] <- run_models_for_panel(panels[[nm]], nm, controls_used[[nm]])
}

# -------------------------
# 7) OUTPUTS + AUDIT
# -------------------------
# Tables + coef CSVs
for (panel_name in names(all_results)) {
  res <- all_results[[panel_name]]
  
  if (length(res$models_main) > 0) {
    txt_path <- file.path(DIRS$tables, paste0("R09_main_models_", panel_name, "_", TS, ".txt"))
    tex_path <- file.path(DIRS$tables, paste0("R09_main_models_", panel_name, "_", TS, ".tex"))
    write_etable(res$models_main, txt_path, tex_path, title = paste0("Step 9 — Main models (", panel_name, ")"))
    log_msg("Wrote: ", txt_path)
  }
  
  if (length(res$models_es) > 0) {
    txt_path <- file.path(DIRS$tables, paste0("R09_event_study_models_", panel_name, "_", TS, ".txt"))
    tex_path <- file.path(DIRS$tables, paste0("R09_event_study_models_", panel_name, "_", TS, ".tex"))
    write_etable(res$models_es, txt_path, tex_path, title = paste0("Step 9 — Event study models (", panel_name, ")"))
    log_msg("Wrote: ", txt_path)
  }
  
  if (length(res$es_coefs) > 0) {
    es_dt <- rbindlist(res$es_coefs, fill = TRUE)
    out_csv <- file.path(DIRS$out_paper1, paste0("event_study_coefs_", panel_name, "_", TS, ".csv"))
    fwrite(es_dt, out_csv)
    log_msg("Wrote: ", out_csv)
  }
}

# Build audit
audit_parts <- list()

# Control diagnostics
audit_parts$control_meta <- rbindlist(lapply(names(panels), function(nm) {
  data.table(
    section = "controls_used",
    panel = nm,
    n_controls_used = length(controls_used[[nm]]),
    controls_used = paste(controls_used[[nm]], collapse = ", ")
  )
}), fill = TRUE)

audit_parts$controls_missing <- rbindlist(lapply(names(panels), function(nm) {
  miss <- controls_missing[[nm]]
  data.table(
    section = "controls_missing_in_panel_ignored",
    panel = nm,
    n_missing = length(miss),
    missing = paste(head(miss, 50), collapse = ", ")
  )
}), fill = TRUE)

audit_parts$controls_drop_junk <- data.table(
  section = "controls_drop_junk",
  dropped = paste(intersect(DROP_NOT_CONTROLS, names(ctrl)), collapse = ", ")
)

for (panel_name in names(all_results)) {
  res <- all_results[[panel_name]]
  audit_parts[[paste0(panel_name, "_panel_summary")]] <- cbind(section = "panel_summary", res$audit$panel_summary)
  
  es_aud <- rbindlist(res$audit[grepl("^es_", names(res$audit))], fill = TRUE)
  if (nrow(es_aud) > 0) audit_parts[[paste0(panel_name, "_event_study_summary")]] <- cbind(section = "event_study_summary", es_aud)
  
  # Key coefficient snapshot
  if (length(res$models_main) > 0) {
    snap <- rbindlist(lapply(names(res$models_main), function(mn) {
      m <- res$models_main[[mn]]
      ct <- as.data.table(coeftable(m))
      ct[, term := rownames(coeftable(m))]
      keep <- ct[term %in% c("post_eli", "eli_feas", "post_eli:eli_feas")]
      if (nrow(keep) == 0) keep <- ct[grepl("post_eli|eli_feas", term)]
      keep[, `:=`(panel = panel_name, model = mn)]
      keep
    }), fill = TRUE)
    audit_parts[[paste0(panel_name, "_key_coefs_main")]] <- cbind(section = "key_coefs_main", snap)
  }
}

audit_dt <- rbindlist(audit_parts, fill = TRUE)
fwrite(audit_dt, OUT_AUDIT)
log_msg("Wrote audit: ", OUT_AUDIT)

log_msg("Step 9 complete.")
