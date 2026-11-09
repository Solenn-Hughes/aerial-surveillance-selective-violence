# Local-only staging of original datasets and saved model panels.
# Run in repository root: source('PREPARE_LOCAL_DATA.R')
PROJECT <- normalizePath(getwd(), winslash='/', mustWork=TRUE)
ARCHIVE <- file.path(Sys.getenv('USERPROFILE'), 'Downloads', 'JPR_COMPLETE_RESEARCH_ARCHIVE')
if (!dir.exists(ARCHIVE)) stop('Local research archive missing. Set ARCHIVE in PREPARE_LOCAL_DATA.R to its location.')
copy_tree <- function(from, to) {
  if (!dir.exists(from)) stop('Required source directory not found: ', from)
  fs <- list.files(from, recursive=TRUE, full.names=TRUE)
  fs <- fs[file.exists(fs) & !dir.exists(fs)]
  root <- normalizePath(from, winslash='/', mustWork=TRUE)
  for (f in fs) {
    rel <- substring(gsub('\\\\','/',f), nchar(root)+2L)
    dest <- file.path(to, rel)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if (!file.exists(dest) && !file.copy(f,dest)) stop('Copy failed: ',f)
  }
  length(fs)
}
cat('Copying raw data locally (excluded by .gitignore)...\n')
copy_tree(file.path(ARCHIVE,'data','raw'),file.path(PROJECT,'data','raw'))
copy_tree(file.path(ARCHIVE,'data','external_LWJ'),file.path(PROJECT,'data','external_LWJ'))
copy_tree(file.path(ARCHIVE,'outputs','final_airfield_screen_20260809_FINAL30'),file.path(PROJECT,'reference','historical_outputs','final_airfield_screen_20260809_FINAL30'))
copy_tree(file.path(ARCHIVE,'outputs','paper1_geo_afghanistan_executable_acled','panels'),file.path(PROJECT,'reference','historical_outputs','paper1_geo_afghanistan_executable_acled','panels'))
copy_tree(file.path(ARCHIVE,'outputs','paper1_geo_afghanistan_executable_acled','diagnostics_h2_spatial_adm2'),file.path(PROJECT,'reference','historical_outputs','diagnostics_h2_spatial_adm2'))
p <- file.path(PROJECT,'reference','historical_outputs','final_airfield_screen_20260809_FINAL30','11_final_analysis_panel.rds')
if (!file.exists(p)) stop('Saved FINAL30 panel not found')
dest <- file.path(PROJECT,'outputs','final_airfield_screen_20260809_FINAL30','11_final_analysis_panel.rds')
dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
if (!file.exists(dest) && !file.copy(p,dest)) stop('Cannot stage final panel')
cat('Local data ready. Saved-panel check: source("RUN_VERIFIED_14.R")\n')
