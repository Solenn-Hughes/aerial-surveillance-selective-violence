# Re-run original published H2 supplementary analysis scripts without changing estimators.
# Requires completed INSTALL_LOCAL.R and installed data.table, fixest packages.
PROJECT <- 'C:/Users/USERNAME/Downloads/JPR_REPLICATION_RELEASE_FROM_ACTUAL_CODE'
if (!requireNamespace('data.table',quietly=TRUE) || !requireNamespace('fixest',quietly=TRUE))
 stop('Install packages: install.packages(c("data.table","fixest"))')
SRC <- file.path(PROJECT,'code','published_H2')
PATCHED <- file.path(PROJECT,'_patched_original_H2')
dir.create(PATCHED,recursive=TRUE,showWarnings=FALSE)
orig <- c('JPR_H2_LWJ_00_COMMON.R',
 'JPR_H2_FINAL_17_PPML_SAMPLE_RECONCILIATION.R',
 'JPR_H2_FINAL_18B_FROZEN_VALIDATION_OLS_FIX.R')
if (!all(file.exists(file.path(SRC,orig)))) stop('Original H2 scripts missing')
for(n in orig) {
 s <- readLines(file.path(SRC,n),warn=FALSE,encoding='UTF-8')
 s <- gsub('C:/Users/USERNAME/OneDrive/Desktop/THESIS',PROJECT,s,fixed=TRUE)
 s <- gsub('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R',
           file.path(PATCHED,'JPR_H2_LWJ_00_COMMON.R'),s,fixed=TRUE)
 # ROOT and DOWNLOADS data dependencies are redirected to the replicated package.
 s <- gsub('C:/Users/USERNAME/Downloads',file.path(PROJECT,'data','external_LWJ'),s,fixed=TRUE)
 writeLines(s,file.path(PATCHED,n),useBytes=TRUE)
}
# Verify that the original source still defines the expected reference data.
if (!file.exists(file.path(PROJECT,'data','external_LWJ','LWJ_Dudek_district_periods.csv')) ||
    !file.exists(file.path(PROJECT,'data','external_LWJ','LWJ_FINAL30_VERIFIED_DISTRICT_CROSSWALK.csv')))
 stop('External LWJ data missing. Run INSTALL_LOCAL.R and verify archive inputs')
if (!file.exists(file.path(PROJECT,'outputs','paper1_geo_afghanistan_executable_acled',
    'diagnostics_h2_spatial_adm2','district_month_contestedness_spatial_adm2.csv')))
 stop('Historical ACLED district-month input missing. Cannot calculate frozen appendix models')
cat('Re-running original final H2 event-support and frozen-timing scripts...\n')
source(file.path(PATCHED,'JPR_H2_FINAL_17_PPML_SAMPLE_RECONCILIATION.R'),local=FALSE)
source(file.path(PATCHED,'JPR_H2_FINAL_18B_FROZEN_VALIDATION_OLS_FIX.R'),local=FALSE)
cat('Appendix model runs completed; inspect outputs/JPR_H2_LWJ_REMAINING_DIAGNOSTICS_20261106.\n')
