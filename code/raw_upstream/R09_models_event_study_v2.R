# ============================================================
# R09_models_event_study_v2.R — Step 9 (DiD + Event Study) — v2
#
# Reads Step 8 v2 outputs if available (falls back to v1 if not).
# Saves models to RDS so you can inspect even after rm(list=ls()) scripts.
#
# Key upgrades vs v1:
#  - auto-picks *_posteli_v2.rds if present
#  - requires theatre_id for stronger FE option
#  - runs TWO FE variants:
#      (A) pg_id + month (global month FE)
#      (B) pg_id + theatre_id^month_int (theatre-by-month FE)  <-- recommended
#  - event study runs on treated units only (event_time defined only for treated)
#  - writes coefficient CSVs + audit with lead joint Wald test
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
  stop("Package 'fixest' is required. Install via install.packages('fixest').")
}
suppressMessages(library(fixest))
try(setFixest_nthreads(getDTthreads()), silent = TRUE)

# -------------------------
# 1) SETTINGS
# -------------------------
RUN_ALLCW  <- TRUE
RUN_STRICT <- TRUE

# Primary outcomes (event counts)
# Sensitivity outcomes (fatalities): fat_best_type3, civ_best_type3, fat_best_ged
# These are carried through from R04 via zero-fill in the base panel.
# We run all together; interpret fatality models as sensitivity checks.
OUTCOMES <- c("events_type3", "events_ged", "fat_best_type3", "civ_best_type3")

RUN_MAIN_FEOLS_LOG1P <- TRUE
RUN_MAIN_FEPOIS      <- TRUE

RUN_EVENT_STUDY <- TRUE
# Run two windows: tighter (balanced) + your original longer
EVENT_WINDOWS <- list(
  w24 = c(-24L, 24L),
  w60 = c(-24L, 60L)
)
ES_REF <- -1L

# Cluster
CLUSTER_TWOWAY <- TRUE

# Prefer theatre-month FE (recommended)
RUN_FE_SPECS <- c("global_month", "theatre_month")  # keep both for comparison

# Save models so all_results is inspectable later
SAVE_MODELS_RDS <- TRUE

# -------------------------
# 2) PATHS (auto-pick v2 if exists)
# -------------------------
pick_existing <- function(paths, label){
  for (p in paths) if (file.exists(p)) return(p)
  stop("No existing path found for ", label, ". Tried:\n", paste(paths, collapse = "\n"))
}

P_ALLCW <- pick_existing(
  c(
    file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds"),
    file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli.rds")
  ),
  "ALLCW posteli panel"
)

P_STRICT <- pick_existing(
  c(
    file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds"),
    file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds")
  ),
  "STRICT posteli panel"
)

P_CTRL <- file.path(DIRS$base, "controls_pg_month.rds")
if (!file.exists(P_CTRL)) stop("Missing required controls file: ", P_CTRL)

TS <- format(Sys.time(), "%Y%m%d_%H%M%S")

OUT_AUDIT <- file.path(DIRS$audit, paste0("R09_models_event_study_v2_", TS, ".csv"))
OUT_MODELS_RDS <- file.path(DIRS$out_paper1, paste0("R09_models_event_study_v2_models_", TS, ".rds"))

# -------------------------
# 3) HELPERS
# -------------------------
assert_has <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(sprintf("%s missing cols: %s", label, paste(miss, collapse = ", ")))
}

add_month_int <- function(dt) {
  # month is IDate (YYYY-MM-01)
  dt[, month_int := (data.table::year(month) * 12L + data.table::month(month))]
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

extract_es <- function(m, outcome, panel_name, fe_spec, window_tag) {
  ct <- as.data.table(coeftable(m))
  ct[, term := rownames(coeftable(m))]
  # fixest i() terms look like: event_time_es::(-24):eli_feas  (or similar)
  ct <- ct[grepl("^event_time_es::", term)]
  ct[, event_time := as.integer(sub("^event_time_es::\\(?(-?\\d+)\\)?.*$", "\\1", term))]
  ct[, `:=`(panel = panel_name, outcome = outcome, fe_spec = fe_spec, window = window_tag)]
  setcolorder(ct, c("panel","outcome","fe_spec","window","term","event_time","Estimate","Std. Error","t value","Pr(>|t|)"))
  ct[order(event_time)]
}

save_es_plot <- function(m, file_path, main_title) {
  png(file_path, width = 1200, height = 700, res = 150)
  on.exit(dev.off(), add = TRUE)
  iplot(m, ref.line = 0, main = main_title, xlab = "Event time (months)", ylab = "Effect (relative to ref)")
  invisible(TRUE)
}

# Pre-trend joint test on leads: k in [lead_min, lead_max]
wald_leads <- function(m, lead_min=-24L, lead_max=-2L) {
  cn <- names(coef(m))
  if (is.null(cn) || !length(cn)) return(list(n_leads=0L, p=NA_real_, stat=NA_real_, note="no coefficients"))
  
  # keep only event_time_es i() terms
  et_terms <- cn[grepl("^event_time_es::", cn)]
  if (!length(et_terms)) return(list(n_leads=0L, p=NA_real_, stat=NA_real_, note="no event_time terms"))
  
  # parse k
  k <- suppressWarnings(as.integer(sub("^event_time_es::\\(?(-?\\d+)\\)?.*$", "\\1", et_terms)))
  keep <- et_terms[!is.na(k) & k >= lead_min & k <= lead_max]
  if (!length(keep)) return(list(n_leads=0L, p=NA_real_, stat=NA_real_, note="no leads in window"))
  
  # Drop any that got aliased out (shouldn’t happen, but be safe)
  keep <- intersect(keep, names(coef(m)))
  if (!length(keep)) return(list(n_leads=0L, p=NA_real_, stat=NA_real_, note="lead terms collinear/dropped"))
  
  w <- tryCatch(wald(m, keep), error = function(e) e)
  if (inherits(w, "error")) return(list(n_leads=length(keep), p=NA_real_, stat=NA_real_, note=paste0("wald error: ", w$message)))
  
  # fixest::wald returns a list with p.value / stat (naming can vary slightly across versions)
  p <- if (!is.null(w$p.value)) w$p.value else if (!is.null(w$pval)) w$pval else NA_real_
  s <- if (!is.null(w$stat)) w$stat else if (!is.null(w$statistic)) w$statistic else NA_real_
  
  list(n_leads=length(keep), p=as.numeric(p), stat=as.numeric(s), note="ok")
}

# -------------------------
# 4) CONTROL LIST (robust)
# -------------------------
ctrl <- readRDS(P_CTRL); setDT(ctrl)
if (!("pg_id" %in% names(ctrl)) && "gid" %in% names(ctrl)) ctrl[, pg_id := gid]
assert_has(ctrl, c("pg_id","month"), "controls_pg_month")

DROP_NOT_CONTROLS <- c(
  "pg_id","month","gid","gid.x","gid.y","month_id","i.month_id","month_int",
  "country_id","country","lon","lat","xcoord","ycoord",
  "row","col","year",
  # posteli fields that might appear
  "theatre_id","treated","post_eli","event_time","eli_onset_month","event_time_es"
)

CONTROL_CAND <- setdiff(names(ctrl), DROP_NOT_CONTROLS)
CONTROL_CAND <- CONTROL_CAND[vapply(CONTROL_CAND, function(cc) is.numeric(ctrl[[cc]]), logical(1))]

if (!length(CONTROL_CAND)) stop("No numeric controls found after dropping junk columns.")

# -------------------------
# 5) LOAD PANELS + INTEGRITY
# -------------------------
load_panel <- function(path, panel_name) {
  dt <- readRDS(path); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  dt[, pg_id := as.integer(pg_id)]
  dt[, month := as_month_idate(month)]
  
  assert_unique(dt, c("pg_id","month"), name = paste0(panel_name, " keys"))
  
  # Required vars
  core <- c("eli_feas","post_eli","event_time")
  assert_has(dt, core, panel_name)
  
  # theatre_id strongly preferred for the improved FE spec
  if (!("theatre_id" %in% names(dt))) {
    # attempt to build from country if present
    if ("country" %in% names(dt)) {
      dt[, theatre_id := as.character(country)]
    } else if ("country_id" %in% names(dt)) {
      dt[, theatre_id := paste0("country_id_", as.integer(country_id))]
    } else {
      stop(panel_name, ": no theatre_id and no country/country_id to construct it.")
    }
  }
  
  # Bounds/types
  if (anyNA(dt$eli_feas)) stop(panel_name, ": eli_feas has NA.")
  if (any(dt$eli_feas < 0 | dt$eli_feas > 1)) stop(panel_name, ": eli_feas out of [0,1].")
  if (anyNA(dt$post_eli)) stop(panel_name, ": post_eli has NA.")
  if (any(!dt$post_eli %in% c(0L,1L))) stop(panel_name, ": post_eli not binary.")
  
  dt <- add_month_int(dt)
  dt
}

panels <- list()
if (RUN_ALLCW)  panels$ALLCW  <- load_panel(P_ALLCW,  "ALLCW_posteli")
if (RUN_STRICT) panels$STRICT <- load_panel(P_STRICT, "STRICT_posteli")

controls_used <- lapply(panels, function(dt) intersect(CONTROL_CAND, names(dt)))
if (any(vapply(controls_used, length, integer(1)) == 0L)) {
  stop("A panel ended up with zero usable controls after intersection (unexpected).")
}
controls_missing <- lapply(panels, function(dt) setdiff(CONTROL_CAND, names(dt)))

# -------------------------
# 6) MODEL RUNNERS
# -------------------------
run_models_for_panel <- function(dt, panel_name, CONTROL_COLS) {
  results <- list(models_main = list(), models_es = list(), es_coefs = list(), audit = list())
  
  vcov_arg <- if (CLUSTER_TWOWAY) ~ pg_id + month_int else ~ pg_id
  
  rhs_controls <- paste(CONTROL_COLS, collapse = " + ")
  
  results$audit$panel_summary <- data.table(
    panel = panel_name,
    rows = nrow(dt),
    pg = uniqueN(dt$pg_id),
    months = uniqueN(dt$month),
    theatres = uniqueN(dt$theatre_id),
    treated_theatres = dt[, uniqueN(theatre_id[post_eli %in% c(0L,1L) & !is.na(event_time)])],
    share_post = mean(dt$post_eli == 1L),
    eli_min = min(dt$eli_feas),
    eli_max = max(dt$eli_feas),
    n_controls_used = length(CONTROL_COLS)
  )
  
  # FE specs
  fe_map <- list(
    global_month  = "pg_id + month",
    theatre_month = "pg_id + theatre_id^month_int"
  )
  
  for (y in OUTCOMES) {
    if (!(y %in% names(dt))) {
      log_msg("SKIP outcome '", y, "' — not found in panel '", panel_name, "'. Check R04 zero_cols.")
      next
    }
    assert_nonneg(dt, y, panel_name)
    
    for (fe_name in RUN_FE_SPECS) {
      fe_part <- fe_map[[fe_name]]
      if (is.null(fe_part)) stop("Unknown FE spec: ", fe_name)
      
      # MAIN models
      if (RUN_MAIN_FEOLS_LOG1P) {
        ylog <- paste0("log1p_", y)
        if (!(ylog %in% names(dt))) dt[, (ylog) := log1p(get(y))]
        fml  <- as.formula(paste0(ylog, " ~ post_eli:eli_feas + ", rhs_controls, " | ", fe_part))
        key  <- paste0(y, "_feols_log1p__", fe_name)
        results$models_main[[key]] <- feols(fml, data = dt, vcov = vcov_arg)
      }
      
      if (RUN_MAIN_FEPOIS) {
        fmlp <- as.formula(paste0(y, " ~ post_eli:eli_feas + ", rhs_controls, " | ", fe_part))
        keyp <- paste0(y, "_fepois__", fe_name)
        results$models_main[[keyp]] <- fepois(fmlp, data = dt, vcov = vcov_arg)
      }
      
      # EVENT STUDY (treated-only; event_time defined there)
      if (RUN_EVENT_STUDY) {
        for (wtag in names(EVENT_WINDOWS)) {
          w <- EVENT_WINDOWS[[wtag]]
          
          dt[, event_time_es := as.integer(event_time)]
          dt[event_time_es < w[1] | event_time_es > w[2], event_time_es := NA_integer_]
          
          dt_es <- dt[!is.na(event_time_es)]
          
          # require ref period exists
          if (dt_es[event_time_es == ES_REF, .N] == 0) {
            # Don’t hard fail: just log and skip this ES run
            results$audit[[paste0("es_skip__", y, "__", fe_name, "__", wtag)]] <- data.table(
              panel = panel_name, outcome = y, fe_spec = fe_name, window = wtag,
              note = paste0("SKIP: ref ", ES_REF, " not present in window [", w[1], ",", w[2], "]")
            )
            next
          }
          
          ylog <- paste0("log1p_", y)
          if (!(ylog %in% names(dt_es))) dt_es[, (ylog) := log1p(get(y))]
          
          fml_es <- as.formula(paste0(
            ylog, " ~ i(event_time_es, eli_feas, ref = ", ES_REF, ") + ",
            rhs_controls, " | ", fe_part
          ))
          
          mes <- feols(fml_es, data = dt_es, vcov = vcov_arg)
          
          key_es <- paste0(y, "_eventstudy_feols_log1p__", fe_name, "__", wtag)
          results$models_es[[key_es]] <- mes
          results$es_coefs[[key_es]]  <- extract_es(mes, y, panel_name, fe_name, wtag)
          
          # Pre-trend joint test (leads)
          wl <- wald_leads(mes, lead_min = w[1], lead_max = -2L)
          results$audit[[paste0("wald_leads__", y, "__", fe_name, "__", wtag)]] <- data.table(
            panel = panel_name, outcome = y, fe_spec = fe_name, window = wtag,
            n_leads = wl$n_leads, wald_stat = wl$stat, wald_p = wl$p, note = wl$note
          )
          
          # Plot
          fig_path <- file.path(
            DIRS$figures,
            paste0("event_study_v2_", panel_name, "_", y, "__", fe_name, "__", wtag, "_", TS, ".png")
          )
          save_es_plot(mes, fig_path, paste0("ES v2: ", panel_name, " — ", y, " — ", fe_name, " — ", wtag))
        }
      }
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
audit_parts <- list()

# Controls diagnostics
audit_parts$control_used <- rbindlist(lapply(names(panels), function(nm) {
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

audit_parts$panel_summary <- rbindlist(lapply(names(all_results), function(panel_name) {
  cbind(section = "panel_summary", all_results[[panel_name]]$audit$panel_summary)
}), fill = TRUE)

# Wald lead tests + ES skips
audit_parts$wald_and_skips <- rbindlist(lapply(names(all_results), function(panel_name) {
  aud <- all_results[[panel_name]]$audit
  out <- rbindlist(aud[grepl("^(wald_leads__|es_skip__)", names(aud))], fill = TRUE)
  if (nrow(out)) cbind(section = "wald_and_skips", out) else NULL
}), fill = TRUE)

# Write tables + coef CSVs per panel
for (panel_name in names(all_results)) {
  res <- all_results[[panel_name]]
  
  if (length(res$models_main) > 0) {
    txt_path <- file.path(DIRS$tables, paste0("R09v2_main_models_", panel_name, "_", TS, ".txt"))
    tex_path <- file.path(DIRS$tables, paste0("R09v2_main_models_", panel_name, "_", TS, ".tex"))
    write_etable(res$models_main, txt_path, tex_path, title = paste0("Step 9 v2 — Main models (", panel_name, ")"))
    log_msg("Wrote: ", txt_path)
  }
  
  if (length(res$models_es) > 0) {
    txt_path <- file.path(DIRS$tables, paste0("R09v2_event_study_models_", panel_name, "_", TS, ".txt"))
    tex_path <- file.path(DIRS$tables, paste0("R09v2_event_study_models_", panel_name, "_", TS, ".tex"))
    write_etable(res$models_es, txt_path, tex_path, title = paste0("Step 9 v2 — Event study models (", panel_name, ")"))
    log_msg("Wrote: ", txt_path)
  }
  
  if (length(res$es_coefs) > 0) {
    es_dt <- rbindlist(res$es_coefs, fill = TRUE)
    out_csv <- file.path(DIRS$out_paper1, paste0("event_study_coefs_v2_", panel_name, "_", TS, ".csv"))
    fwrite(es_dt, out_csv)
    log_msg("Wrote: ", out_csv)
  }
}

audit_dt <- rbindlist(audit_parts, fill = TRUE)
fwrite(audit_dt, OUT_AUDIT)
log_msg("Wrote audit: ", OUT_AUDIT)

if (SAVE_MODELS_RDS) {
  saveRDS(all_results, OUT_MODELS_RDS)
  log_msg("Saved models RDS: ", OUT_MODELS_RDS)
}

log_msg("Step 9 v2 complete.")