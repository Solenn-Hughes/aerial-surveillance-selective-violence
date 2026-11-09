# ============================================================
# R99_run_all_v2.R — One-button run (Paper 1 pipeline, v2)
#
# Runs:
#  Step 1..7 as before
#  Step 8: theatre-specific onsets (v2)
#  Step 9: models + ES (v2)
#  Step 9b: dispersion metrics + ES
#  Step 10: audits
#
# Robust to step scripts that call rm(list=ls()) because each
# step is sourced into its own environment.
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

# -------------------------
# 0) SETTINGS
# -------------------------
FORCE_ALL     <- FALSE  # TRUE = rerun everything
FORCE_MODEL   <- TRUE   # TRUE = always rerun Step 9 / 9b / 10
ECHO_SOURCE   <- FALSE  # TRUE = echo source() lines

# -------------------------
# 1) SCRIPT LOCATOR
# -------------------------
locate_script <- function(candidates) {
  for (p in candidates) if (file.exists(p)) return(p)
  stop("Could not locate script. Tried:\n", paste(candidates, collapse = "\n"))
}

# -------------------------
# 2) SOURCE CONFIG (Step 0) IN GLOBAL (sets DIRS, helpers, etc.)
# -------------------------
cfg_path <- locate_script(c("R00_config.R", "R/00_config.R", "00_config.R"))
source(cfg_path, echo = ECHO_SOURCE)

TS <- if (exists(".now_stamp", inherits = TRUE)) .now_stamp() else format(Sys.time(), "%Y%m%d_%H%M%S")
if (!dir.exists(DIRS$audit)) dir.create(DIRS$audit, recursive = TRUE, showWarnings = FALSE)

runlog <- data.table(
  step = character(),
  script = character(),
  status = character(),
  note = character(),
  started = character(),
  finished = character()
)

# -------------------------
# 3) SAFE SOURCE (child script cannot wipe runner state)
# -------------------------
safe_source <- function(script_path) {
  e <- new.env(parent = .GlobalEnv)
  sys.source(script_path, envir = e, keep.source = FALSE)
  invisible(TRUE)
}

# -------------------------
# 4) STEP RUNNER
# -------------------------
run_step <- function(step_name, script_candidates, outputs = character(0), force = FALSE) {
  started <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  script_path <- locate_script(script_candidates)
  
  need_run <- isTRUE(force)
  if (!need_run) {
    if (length(outputs) > 0) need_run <- any(!file.exists(outputs)) else need_run <- FALSE
  }
  
  if (!need_run) {
    runlog <<- rbind(runlog, data.table(
      step = step_name,
      script = script_path,
      status = "SKIP",
      note = "All expected outputs exist (or step is timestamped).",
      started = started,
      finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    ))
    log_msg(step_name, " — skip.")
    return(invisible(TRUE))
  }
  
  log_msg(step_name, " — running: ", script_path)
  
  ok <- TRUE
  err <- NULL
  tryCatch(
    safe_source(script_path),
    error = function(e) { ok <<- FALSE; err <<- conditionMessage(e) }
  )
  
  finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  if (!ok) {
    runlog <<- rbind(runlog, data.table(
      step = step_name,
      script = script_path,
      status = "FAIL",
      note = err,
      started = started,
      finished = finished
    ))
    write_audit_dt(runlog, paste0("R99_run_all_v2_RUNLOG_FAIL_", TS), subdir = DIRS$audit)
    stop(step_name, " failed:\n", err)
  }
  
  if (length(outputs) > 0) {
    missing <- outputs[!file.exists(outputs)]
    if (length(missing)) {
      runlog <<- rbind(runlog, data.table(
        step = step_name,
        script = script_path,
        status = "FAIL",
        note = paste0("Step ran but outputs missing: ", paste(missing, collapse = " | ")),
        started = started,
        finished = finished
      ))
      write_audit_dt(runlog, paste0("R99_run_all_v2_RUNLOG_FAIL_", TS), subdir = DIRS$audit)
      stop(step_name, " ran but expected output(s) not found:\n", paste(missing, collapse = "\n"))
    }
  }
  
  runlog <<- rbind(runlog, data.table(
    step = step_name,
    script = script_path,
    status = "OK",
    note = "",
    started = started,
    finished = finished
  ))
  
  log_msg(step_name, " — OK.")
  invisible(TRUE)
}

# -------------------------
# 5) DEFINE STEPS + OUTPUTS
# -------------------------
steps <- list(
  list(
    name = "Step 1 — ingest PRIO-GRID",
    scripts = c("R01_ingest_priogrid.R", "R/01_ingest_priogrid.R"),
    outs = c(file.path(DIRS$base, "priogrid_static.rds"),
             file.path(DIRS$base, "priogrid_yearly.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Step 2 — ingest airbases",
    scripts = c("R02_ingest_airbases.R", "R/02_ingest_airbases.R"),
    outs = c(file.path(DIRS$base, "airbases_points.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Step 3 — ingest GED + aggregate",
    scripts = c("R03_ingest_ged.R", "R/03_ingest_ged.R",
                "R03_ingest_events_ged.R", "R/03_ingest_events_ged.R"),
    outs = c(file.path(DIRS$base, "ged_clean.rds"),
             file.path(DIRS$base, "ged_pg_month.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Step 4 — build base panel",
    scripts = c("R04_build_base_panel_pg_month.R", "R/04_build_base_panel_pg_month.R"),
    outs = c(file.path(DIRS$base, "panel_pg_month_base.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Extra — build civil war panels",
    scripts = c("R05_build_cw_panels.R", "R/05_build_cw_panels.R"),
    outs = c(file.path(DIRS$panels, "panel_pg_month_ALLCW.rds"),
             file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Step 6 — build ELI feasibility",
    scripts = c("R06_build_eli_feasibility.R", "R/06_build_eli_feasibility.R"),
    outs = c(file.path(DIRS$eli, "eli_feasibility_pg.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Step 5 — build controls",
    scripts = c("R05_build_controls_pg_month.R", "R/05_build_controls_pg_month.R"),
    outs = c(file.path(DIRS$base, "controls_pg_month.rds")),
    force = FORCE_ALL
  ),
  list(
    name = "Step 7 — merge into paper1 panels",
    scripts = c("R07_build_panel_pg_month.R", "R/07_build_panel_pg_month.R"),
    outs = c(file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1.rds"),
             file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1.rds")),
    force = FORCE_ALL
  ),
  
  # --- v2 onset definition
  list(
    name = "Step 8 — define post-ELI + event_time (v2)",
    scripts = c("R08_define_post_eli_v2.R", "R/08_define_post_eli_v2.R"),
    outs = c(file.path(DIRS$eli, "eli_onsets_v2.csv"),
             file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds"),
             file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds")),
    force = FORCE_ALL
  ),
  
  # --- v2 models
  list(
    name = "Step 9 — models + event study (v2)",
    scripts = c("R09_models_event_study_v2.R", "R/09_models_event_study_v2.R"),
    outs = character(0),
    force = (FORCE_ALL || FORCE_MODEL)
  ),
  
  # --- dispersion step
  list(
    name = "Step 9b — dispersion metrics + ES",
    scripts = c("R09b_dispersion_metrics.R", "R/09b_dispersion_metrics.R"),
    outs = character(0),                 # timestamped outputs; rely on failure-on-error
    force = (FORCE_ALL || FORCE_MODEL)
  ),
  
  list(
    name = "Step 10 — automated audits",
    scripts = c("R10_audits.R", "R/10_audits.R"),
    outs = character(0),
    force = (FORCE_ALL || FORCE_MODEL)
  )
)

# -------------------------
# 6) RUN ALL STEPS
# -------------------------
log_msg("=== R99_run_all_v2 starting ===")
log_msg("FORCE_ALL=", FORCE_ALL, " FORCE_MODEL=", FORCE_MODEL)

invisible(lapply(steps, function(st) {
  run_step(st$name, st$scripts, st$outs, force = isTRUE(st$force))
}))

# -------------------------
# 7) FINAL “MUST EXIST” CHECKS (v2)
# -------------------------
final_must <- c(
  file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds"),
  file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds"),
  file.path(DIRS$eli,    "eli_feasibility_pg.rds"),
  file.path(DIRS$base,   "controls_pg_month.rds")
)

missing_final <- final_must[!file.exists(final_must)]
if (length(missing_final)) {
  runlog <- rbind(runlog, data.table(
    step = "FINAL_CHECK",
    script = "R99_run_all_v2.R",
    status = "FAIL",
    note = paste0("Missing final outputs: ", paste(missing_final, collapse = " | ")),
    started = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ))
  write_audit_dt(runlog, paste0("R99_run_all_v2_RUNLOG_FAIL_", TS), subdir = DIRS$audit)
  stop("R99 finished but final required outputs are missing:\n", paste(missing_final, collapse = "\n"))
}

write_audit_dt(runlog, paste0("R99_run_all_v2_RUNLOG_OK_", TS), subdir = DIRS$audit)
log_msg("=== R99_run_all_v2 complete (OK) ===")
