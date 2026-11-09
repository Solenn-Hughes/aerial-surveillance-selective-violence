# Re-estimate publication model results from the saved FINAL30 panel.
# This tier does NOT assert construction from raw source data.
PROJECT <- normalizePath(getOption('jpr.project_root', getwd()), winslash='/', mustWork=TRUE)
if (!requireNamespace('data.table',quietly=TRUE) || !requireNamespace('fixest',quietly=TRUE))
  stop('Install packages first: install.packages(c("data.table","fixest"))')
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(fixest))
PANEL <- file.path(PROJECT,'outputs','final_airfield_screen_20260809_FINAL30','11_final_analysis_panel.rds')
if (!file.exists(PANEL)) stop('Run INSTALL_LOCAL.R first. Published panel missing.')
DEST <- file.path(PROJECT,'outputs','publication_model_replication')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
D <- data.table::as.data.table(readRDS(PANEL))
need <- c('cell_id','month','events_type3','post_closure','eli_feas','eli_final','prec_gpcp','temp','nlights_calib_mean','district_id_for_h2')
miss <- setdiff(need,names(D))
if (length(miss)) stop('Published panel missing: ',paste(miss,collapse=', '))
D[,month:=data.table::as.IDate(month)]
D[,cell_id:=factor(cell_id)]
if (nrow(D)!=104112L || data.table::uniqueN(D$cell_id)!=241L || anyDuplicated(D,by=c('cell_id','month')))
  stop('Published FINAL30 panel size or cell-month keys mismatch. Do not estimate.')

results <- list()
add_result <- function(label,X,fml,term,vcov_spec,target_beta,target_se,target_n,tol=0.005) {
  mod <- fixest::fepois(fml,data=X,notes=FALSE,warn=FALSE)
  ct <- fixest::coeftable(mod,vcov=vcov_spec)
  if (!(term %in% rownames(ct))) stop(label,': requested coefficient not identified: ',term)
  beta <- unname(ct[term,1]); se <- unname(ct[term,2]); p <- unname(ct[term,4]); n <- as.integer(stats::nobs(mod))
  ok <- is.finite(beta) && is.finite(se) && abs(beta-target_beta)<tol && abs(se-target_se)<tol && n==target_n
  results[[length(results)+1L]] <<- data.frame(result=label,term=term,beta=beta,se=se,p=p,n=n,
    target_beta=target_beta,target_se=target_se,target_n=target_n,within_tolerance=ok)
  cat(label,': beta=',sprintf('%.5f',beta),' SE=',sprintf('%.5f',se),' N=',n,' match=',ok,'\n',sep='')
  invisible(mod)
}
FULL <- D
DROP <- D[month < as.IDate('2020-03-01') | month >= as.IDate('2021-07-01')]
controls <- 'prec_gpcp + temp + nlights_calib_mean'
H1 <- stats::as.formula(paste('events_type3 ~ post_closure:eli_final +',controls,'| cell_id + month'))
add_result('H1_full',FULL,H1,'post_closure:eli_final',~cell_id+month,-2.770786,1.089019,54312L)
# Full H1 exclusion estimated with the same specification and comparison benchmark.
add_result('H1_drop_anticipation',DROP,H1,'post_closure:eli_final',~cell_id+month,-2.719,1.124,49956L,tol=0.005)

# Appendix A1: deliberately restricted original ACLED interactions.
mods <- c('pre_contested_share_spatial','pre_taliban_activity_share_spatial','pre_contested_intensity_spatial')
for(i in seq_along(mods)) {
  m <- mods[i]
  form <- stats::as.formula(paste('events_type3 ~ post_closure:eli_final +',
    paste0('post_closure:eli_final:',m),'+',controls,'| cell_id + month'))
  targets_full <- c(-16.102380,-17.666752,-3.120677)
  se_full <- c(4.273193,4.603180,0.880878)
  targets_drop <- c(-17.408208,-19.070708,-3.388981)
  se_drop <- c(4.867394,5.268160,0.984560)
  add_result(paste0('A1_',m,'_full'),FULL,form,paste0('post_closure:eli_final:',m),
    ~cell_id+month,targets_full[i],se_full[i],54312L)
  add_result(paste0('A1_',m,'_drop'),DROP,form,paste0('post_closure:eli_final:',m),
    ~cell_id+month,targets_drop[i],se_drop[i],49956L)
}

# Preferred H2: joint ACLED-LWJ, with required lower-order post interactions.
LWJ <- file.path(PROJECT,'data','external_LWJ','LWJ_Dudek_district_periods.csv')
CROSS <- file.path(PROJECT,'data','external_LWJ','LWJ_FINAL30_VERIFIED_DISTRICT_CROSSWALK.csv')
if (!file.exists(LWJ) || !file.exists(CROSS)) stop('Required LWJ/crosswalk files missing; see data/external_LWJ')
L <- data.table::fread(LWJ)
CW <- data.table::fread(CROSS,colClasses=list(character=c('adm2_id','lwj_object_id')))
L[,lwj_object_id:=as.character(lwj_object_id)]
CW[,lwj_object_id:=as.character(lwj_object_id)]
if (anyDuplicated(CW$adm2_id) || anyDuplicated(CW$lwj_object_id)) stop('Crosswalk not unique')
C <- unique(D[,.(cell_id,district_id_for_h2)])
if (anyDuplicated(C$cell_id)) stop('Cell has multiple H2 districts')
C[,lwj_object_id:=CW$lwj_object_id[match(as.character(district_id_for_h2),CW$adm2_id)]]
if (anyNA(C$lwj_object_id)) stop('Unmapped historical districts')
H <- L[period %in% 1:5,.(lwj_contested_share=mean(control=='Contested')),
        by=lwj_object_id]
D[,lwj_object_id:=C$lwj_object_id[match(cell_id,C$cell_id)]]
D[,lwj_contested_share:=H$lwj_contested_share[match(lwj_object_id,H$lwj_object_id)]]
if (anyNA(D$lwj_contested_share)) stop('Missing historical moderator')
D[,district_id_for_h2:=factor(district_id_for_h2)]
F2 <- events_type3 ~ post_closure:eli_final +
 post_closure:pre_contested_share_spatial + post_closure:lwj_contested_share +
 post_closure:eli_final:pre_contested_share_spatial +
 post_closure:eli_final:lwj_contested_share | cell_id + month
H2FULL <- D[month>=as.IDate('2015-01-01') & month<=as.IDate('2024-12-01')]
H2DROP <- H2FULL[month<as.IDate('2020-03-01') | month>=as.IDate('2021-07-01')]
add_result('H2_joint_full',H2FULL,F2,'post_closure:eli_final:lwj_contested_share',
 ~district_id_for_h2+month,-4.997379,2.441607,19920L)
add_result('H2_joint_drop',H2DROP,F2,'post_closure:eli_final:lwj_contested_share',
 ~district_id_for_h2+month,-6.309611,2.512772,16536L)

# H3: match the source panel and exact original estimator, not FINAL30 H1/H2 panel.
# Historical source: original Afghanistan executable ACLED branch.
# The following candidates permit an installed clone; local data are not included in Git.
H3_SOURCE_CANDIDATES <- c(
 file.path(PROJECT,'reference','historical_outputs','paper1_geo_afghanistan_executable_acled','panels','panel_core_cell_month.rds'),
 file.path('C:/Users/USERNAME/Downloads/JPR_COMPLETE_RESEARCH_ARCHIVE',
   'outputs','paper1_geo_afghanistan_executable_acled','panels','panel_core_cell_month.rds'),
 file.path('C:/Users/USERNAME/OneDrive/Desktop/THESIS',
   'outputs','paper1_geo_afghanistan_executable_acled','panels','panel_core_cell_month.rds'),
 file.path(PROJECT,'outputs','paper1_geo_afghanistan_executable_acled','panels','panel_core_cell_month.rds')
)
H3_SOURCE <- H3_SOURCE_CANDIDATES[file.exists(H3_SOURCE_CANDIDATES)][1L]
if (is.na(H3_SOURCE)) stop('Original H3 core panel is missing from previously known locations')
cat('H3 original core panel: ',H3_SOURCE,'\n',sep='')
H3DATA <- data.table::as.data.table(readRDS(H3_SOURCE))
H3DATA[,month:=data.table::as.IDate(month)]
H3DATA[,cell_id:=factor(cell_id)]
# The publication H3 models used the final cell-level geography from the
# FINAL30 publication output, merged onto the original month-level panel.
# The H3 variables embedded in the historical panel are NOT authoritative.
H3_GEO_CANDIDATES <- c(
  file.path(PROJECT, 'reference','historical_outputs','final_airfield_screen_20260809_FINAL30',
            'final_replacement_values','05_H3_FINAL_CELL_GEOGRAPHY.csv'),
  file.path(PROJECT, 'outputs','final_airfield_screen_20260809_FINAL30',
            'final_replacement_values','05_H3_FINAL_CELL_GEOGRAPHY.csv'),
  file.path('C:/Users/USERNAME/OneDrive/Desktop/THESIS','outputs',
            'final_airfield_screen_20260809_FINAL30','final_replacement_values',
            '05_H3_FINAL_CELL_GEOGRAPHY.csv'),
  file.path('C:/Users/USERNAME/Downloads/JPR_COMPLETE_RESEARCH_ARCHIVE','outputs',
            'final_airfield_screen_20260809_FINAL30','final_replacement_values',
            '05_H3_FINAL_CELL_GEOGRAPHY.csv')
)
H3_GEO <- H3_GEO_CANDIDATES[file.exists(H3_GEO_CANDIDATES)][1L]
if (is.na(H3_GEO)) stop('Missing publication H3 geography: 05_H3_FINAL_CELL_GEOGRAPHY.csv')
GEO <- data.table::fread(H3_GEO)
if ('cell_id_num' %in% names(GEO) && !'cell_id' %in% names(GEO))
  data.table::setnames(GEO,'cell_id_num','cell_id')
if (!'cell_id' %in% names(GEO) || anyDuplicated(as.character(GEO$cell_id)))
  stop('Invalid publication H3 geography: cell identifiers missing or duplicated')
GEO_COLS <- c('low_eli_cell','high_eli_cell','source_cluster_id',
              'adjacent_0_50','fringe_50_100','h3_local_sample',
              'refuge_ring_0_50','refuge_ring_50_100','refuge_ring_100_150')
if (length(setdiff(GEO_COLS,names(GEO))))
  stop('Publication H3 geography missing columns: ',
       paste(setdiff(GEO_COLS,names(GEO)),collapse=', '))
geo_index <- match(as.character(H3DATA$cell_id),as.character(GEO$cell_id))
if (anyNA(geo_index) || data.table::uniqueN(H3DATA$cell_id) != 241L)
  stop('Publication H3 geography does not cover all 241 analysis cells')
for (col in GEO_COLS) H3DATA[, (col) := GEO[[col]][geo_index]]
cat('H3 final publication geography: ',H3_GEO,' (',
    data.table::uniqueN(H3DATA$cell_id),' matched cells)\n',sep='')
H3_REQUIRED <- c('events_type3','post_closure','low_eli_cell','h3_local_sample','adjacent_0_50',
                 'fringe_50_100','source_cluster_id','refuge_ring_0_50','refuge_ring_50_100',
                 'refuge_ring_100_150','prec_gpcp','temp','nlights_calib_mean')
if (length(setdiff(H3_REQUIRED,names(H3DATA))))
  stop('Original core panel lacks columns: ',paste(setdiff(H3_REQUIRED,names(H3DATA)),collapse=', '))
# Canonical file accompanies the actual source archive, not freshly computed outputs.
CANONICAL <- file.path(PROJECT,'reference','FINAL30_CANONICAL_all_model_results.csv')
if (!file.exists(CANONICAL)) stop('Missing original canonical H3 output: ',CANONICAL)
# Read canonical results as an ordinary data.frame to avoid data.table i-column
# evaluation ambiguity in the publication-reference lookup.
CAN <- utils::read.csv(CANONICAL, stringsAsFactors=FALSE, check.names=FALSE,
                       fileEncoding="UTF-8-BOM")
stopifnot(all(c("result_block","term","estimate","std_error","nobs") %in% names(CAN)))
CAN$result_block <- trimws(as.character(CAN$result_block))
CAN$term <- trimws(as.character(CAN$term))
get_ref <- function(block, target_term) {
  key1 <- trimws(as.character(block))
  key2 <- trimws(as.character(target_term))
  matching <- which(!is.na(CAN$result_block) & !is.na(CAN$term) &
                    CAN$result_block == key1 & CAN$term == key2)
  if (length(matching) != 1L) {
    alternatives <- CAN$term[CAN$result_block == key1]
    stop(sprintf("Canonical reference lookup found %d matching rows for %s / %s. Available terms for block: %s. File: %s",
                 length(matching), key1, key2, paste(alternatives,collapse=" | "), CANONICAL),
         call.=FALSE)
  }
  CAN[matching, , drop=FALSE]
}
# Fail before fitting costly models if any expected canonical row is absent.
check_h3_refs <- list(
  c("afghanistan_h3_local_cluster_ppml","post_closure:adjacent_0_50"),
  c("afghanistan_h3_global_rings_ppml","post_closure:refuge_ring_0_50"),
  c("afghanistan_h3_global_rings_ppml","post_closure:refuge_ring_50_100"),
  c("afghanistan_h3_global_rings_ppml","post_closure:refuge_ring_100_150")
)
for (key in check_h3_refs) get_ref(key[[1]],key[[2]])
cat("H3 canonical reference lookup: 4 of 4 keys located.\n")
LOCAL <- H3DATA[low_eli_cell==1L & h3_local_sample==1L &
               (adjacent_0_50==1L | fringe_50_100==1L)]
LOCAL[,source_cluster_id:=factor(source_cluster_id)]
LOCAL[,source_cluster_month:=interaction(source_cluster_id,month,drop=TRUE)]
H3F <- stats::as.formula(paste('events_type3 ~ post_closure:adjacent_0_50 +', controls,
                             '| cell_id + source_cluster_month'))
fit3 <- fixest::fepois(H3F,data=LOCAL,notes=FALSE,warn=FALSE)
ct3 <- fixest::coeftable(fit3,vcov=~cell_id+source_cluster_id)
lt <- 'post_closure:adjacent_0_50'
if (!(lt %in% rownames(ct3))) stop('Original H3 local term not identified')
local_ref <- get_ref('afghanistan_h3_local_cluster_ppml',lt)
local3 <- data.frame(result='H3_local',term=lt,beta=unname(ct3[lt,1]),
  se=unname(ct3[lt,2]),n=as.integer(stats::nobs(fit3)),
  target_beta=local_ref$estimate,target_se=local_ref$std_error,
  target_n=as.integer(local_ref$nobs),source_panel=H3_SOURCE)
local3$within_tolerance <- abs(local3$beta-local3$target_beta)<0.005 &
 abs(local3$se-local3$target_se)<0.005 & local3$n==local3$target_n
write.csv(local3,file.path(DEST,'H3_LOCAL.csv'),row.names=FALSE)
cat('H3_local: beta=',round(local3$beta,5),' SE=',round(local3$se,5),
 ' N=',local3$n,' match=',local3$within_tolerance,'\n',sep='')
G <- H3DATA[low_eli_cell==1L]
terms <- paste0('post_closure:',c('refuge_ring_0_50','refuge_ring_50_100','refuge_ring_100_150'))
H3G <- stats::as.formula(paste('events_type3 ~',paste(terms,collapse=' + '),'+',controls,
                             '| cell_id + month'))
fitg <- fixest::fepois(H3G,data=G,notes=FALSE,warn=FALSE)
ctg <- fixest::coeftable(fitg,vcov=~cell_id+month)
if (!all(terms %in% rownames(ctg))) stop('Original H3 global terms not identified')
refs <- lapply(terms,function(t)get_ref('afghanistan_h3_global_rings_ppml',t))
glob <- data.frame(result='H3_global_ring',term=terms,beta=unname(ctg[terms,1]),
  se=unname(ctg[terms,2]),n=as.integer(stats::nobs(fitg)),
  target_beta=vapply(refs,function(r)r$estimate[1],numeric(1)),
  target_se=vapply(refs,function(r)r$std_error[1],numeric(1)),
  target_n=vapply(refs,function(r)as.integer(r$nobs[1]),integer(1)),source_panel=H3_SOURCE)
glob$within_tolerance <- abs(glob$beta-glob$target_beta)<0.005 &
 abs(glob$se-glob$target_se)<0.005 & glob$n==glob$target_n
write.csv(glob,file.path(DEST,'H3_GLOBAL_RING.csv'),row.names=FALSE)
for(i in seq_len(nrow(glob))) cat('H3_global_',i,': beta=',round(glob$beta[i],5),
  ' SE=',round(glob$se[i],5),' N=',glob$n[i],
  ' match=',glob$within_tolerance[i],'\n',sep='')
comparison <- data.table::rbindlist(results,fill=TRUE)
data.table::fwrite(comparison,file.path(DEST,'PUBLISHED_COEFFICIENT_COMPARISON.csv'))
if (!all(comparison$within_tolerance) || !all(glob$within_tolerance) || !all(local3$within_tolerance))
 stop('Some published model results did not match; inspect PUBLISHED_COEFFICIENT_COMPARISON.csv and H3_GLOBAL_RING.csv')
# Save the complete checked specification table in one authoritative report.
all_checks <- data.table::rbindlist(list(comparison,local3,glob),fill=TRUE)
if (nrow(all_checks) != 14L || !all(all_checks$within_tolerance))
 stop('Internal check failed: expected all 14 model comparisons to pass')
data.table::fwrite(all_checks,file.path(DEST,'ALL_14_PUBLISHED_MODEL_CHECKS.csv'))
cat('\nPASS: all 14 independently re-estimated saved-panel model checks match.\n')
cat('NOT a from-raw verification. Remaining A2-A4 and figure-generation checks still need execution.\n')
