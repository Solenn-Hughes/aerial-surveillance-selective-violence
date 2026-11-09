# From-raw upstream builder. Uses only raw source data in isolated workspace.
# This original upstream branch does NOT itself construct the published FINAL30 sample.
PROJECT <- 'C:/Users/USERNAME/Downloads/JPR_REPLICATION_RELEASE_FROM_ACTUAL_CODE'
if (!dir.exists(file.path(PROJECT,'data','raw'))) stop('Run INSTALL_LOCAL.R first')
ROOT <- file.path(PROJECT,'_fresh_raw_upstream')
if (dir.exists(ROOT)) unlink(ROOT,recursive=TRUE,force=TRUE)
dir.create(ROOT,recursive=TRUE)
copy_tree <- function(from,to) {
  f <- list.files(from,recursive=TRUE,full.names=TRUE,include.dirs=FALSE)
  f <- f[file.exists(f)&!dir.exists(f)]
  from <- normalizePath(from,winslash='/',mustWork=TRUE)
  for (p in f) {
    rel <- substring(gsub('\\\\','/',p),nchar(from)+2L)
    dest <- file.path(to,rel)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(p,dest,overwrite=TRUE)) stop('Copy failed: ',p)
  }
}
copy_tree(file.path(PROJECT,'data','raw'),file.path(ROOT,'data','raw'))
up <- file.path(PROJECT,'code','raw_upstream')
files <- c('R00_config.R','R01_ingest_priogrid.R','R02_ingest_airbases.R',
 'R03_ingest_ged.R','R04_build_base_panel_pg_month.R','R05_build_cw_panels.R',
 'R06_build_eli_feasibility.R','R05_build_controls_pg_month.R',
 'R07_build_panel_pg_month.R','R08_define_post_eli_v2.R',
 'R09_models_event_study_v2.R','R09b_dispersion_metrics.R','R10_audits.R','R99_run_all.R')
if (!all(file.exists(file.path(up,files)))) stop('Missing original upstream script: ',
 paste(files[!file.exists(file.path(up,files))],collapse=', '))
for (name in files) {
  lines <- readLines(file.path(up,name),warn=FALSE,encoding='UTF-8')
  lines <- gsub('C:/Users/USERNAME/OneDrive/Desktop/THESIS',
    gsub('\\\\','/',ROOT),lines,fixed=TRUE)
  writeLines(lines,file.path(ROOT,name),useBytes=TRUE)
}
cat('Executing original from-raw R99 pipeline in isolated directory...\n')
old <- getwd(); on.exit(setwd(old),add=TRUE)
setwd(ROOT)
raw_env <- new.env(parent=.GlobalEnv)
sys.source(file.path(ROOT,'R99_run_all.R'),envir=raw_env)
cat('Raw upstream pipeline finished. Verifying FINAL30 stage...\n')
final <- file.path(ROOT,'outputs','final_airfield_screen_20260809_FINAL30','11_final_analysis_panel.rds')
if (!file.exists(final)) stop(paste(
 'The original upstream R99 code ran, but the separate FINAL30 selection/construction step',
 'was not included in the uploaded code. This is NOT a full publication replication.',
 'The stage must be recovered and linked; do not substitute a saved final panel.'))
cat('FINAL30 panel found in fresh build:',final,'\n')
