# JPR H2 — historical Long War Journal territorial-control measurement
# Inputs: published FINAL30, Dudek LWJ periods, documented district crosswalk.
# Does not overwrite FINAL30. No model selection based on statistical significance.
# Runs six published H2 benchmark replications, then matched LWJ H2 tests.

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
DOWNLOADS <- 'C:/Users/USERNAME/Downloads'
PANEL <- file.path(ROOT, 'outputs', 'final_airfield_screen_20260809_FINAL30',
                   '11_final_analysis_panel.rds')
LWJ <- file.path(DOWNLOADS, 'LWJ_Dudek_district_periods.csv')
CROSSWALK <- file.path(DOWNLOADS, 'LWJ_FINAL30_VERIFIED_DISTRICT_CROSSWALK.csv')
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_LWJ_REMAINING_DIAGNOSTICS_20261106')
for (p in c(PANEL, LWJ, CROSSWALK)) if (!file.exists(p)) stop('Input missing: ', p)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

L <- fread(LWJ, encoding='UTF-8')
CW <- fread(CROSSWALK, encoding='UTF-8', colClasses=list(character=c('adm2_id','lwj_object_id')))
D <- as.data.table(readRDS(PANEL))
D[, `:=`(cell_id=as.character(cell_id),
         month=as.IDate(month),
         district_id_for_h2=as.character(district_id_for_h2))]
L[, lwj_object_id := as.character(lwj_object_id)]
CW[, lwj_object_id := as.character(lwj_object_id)]

required <- c('cell_id','month','district_id_for_h2','events_type3', 'post_closure',
              'eli_final','prec_gpcp','temp','nlights_calib_mean',
              'pre_contested_share_spatial','pre_taliban_activity_share_spatial',
              'pre_contested_intensity_spatial')
if (length(missing <- setdiff(required,names(D)))) stop('FINAL30 columns missing: ',paste(missing,collapse=', '))
if (length(missing <- setdiff(c('adm2_id','lwj_object_id','match_basis'),names(CW))))
  stop('Crosswalk columns missing: ',paste(missing,collapse=', '))
if (length(missing <- setdiff(c('lwj_object_id','period','control'),names(L))))
  stop('LWJ columns missing: ',paste(missing,collapse=', '))
if (nrow(D)!=104112L || uniqueN(D$cell_id)!=241L || anyDuplicated(D,by=c('cell_id','month')))
  stop('FINAL30 integrity mismatch')
if (nrow(L)!=3663L || uniqueN(L$lwj_object_id)!=407L ||
    anyDuplicated(L,by=c('lwj_object_id','period')) || !setequal(L$period,1:9))
  stop('LWJ source structure mismatch')
if (!all(L$control %in% c('GOA','Taliban','Contested','Unconfirmed')))
  stop('Unexpected LWJ control value')
if (anyDuplicated(CW$adm2_id) || anyDuplicated(CW$lwj_object_id))
  stop('Crosswalk is not one-to-one: investigate before estimating')
if (!all(CW$lwj_object_id %in% L$lwj_object_id))
  stop('Crosswalk contains LWJ IDs not in source')
cell <- unique(D[,.(cell_id,district_id_for_h2)])
if (anyDuplicated(cell$cell_id)) stop('Panel cell assigned to multiple districts')
cell[, lwj_object_id := CW$lwj_object_id[match(district_id_for_h2, CW$adm2_id)]]
cat('\nCROSSWALK: ',uniqueN(cell$district_id_for_h2),' panel districts, ',
    nrow(cell),' cells; matched cells ',sum(!is.na(cell$lwj_object_id)),'/',nrow(cell),'\n',sep='')
fwrite(cell,file.path(OUT,'01_cell_join_audit.csv'))
if (anyNA(cell$lwj_object_id)) stop('Unmatched panel cells: see 01_cell_join_audit.csv')
if (uniqueN(cell$district_id_for_h2)!=168L) stop('Expected 168 panel districts')
if (anyDuplicated(cell$lwj_object_id)) {
  # Several PRIO-GRID cells can share one historical district: normal.
  if (uniqueN(cell[,.(district_id_for_h2,lwj_object_id)]) != uniqueN(cell$district_id_for_h2))
    stop('Distinct panel districts share source ID')
}
if (uniqueN(cell[,.(district_id_for_h2,lwj_object_id)]) != uniqueN(cell$district_id_for_h2) ||
    anyDuplicated(unique(cell[,.(district_id_for_h2,lwj_object_id)]),by='lwj_object_id'))
  stop('Panel districts are not in one-to-one relation with LWJ districts')

# Only periods 1–5, pre-Doha. Do NOT use outcome-dependent control categories.
P <- L[period %in% 1:5,.(
  lwj_contested_share=mean(control=='Contested'),
  lwj_taliban_control_share=mean(control=='Taliban'),
  lwj_goa_control_share=mean(control=='GOA'),
  lwj_unconfirmed_share=mean(control=='Unconfirmed'),
  lwj_confirmed_contested_share=if (sum(control!='Unconfirmed')>0)
       sum(control=='Contested')/sum(control!='Unconfirmed') else NA_real_,
  n_pre_doha_periods=.N
),by=lwj_object_id]
if (any(P$n_pre_doha_periods!=5L)) stop('Historical periods incomplete')
if (any(abs(P$lwj_contested_share+P$lwj_taliban_control_share+
             P$lwj_goa_control_share+P$lwj_unconfirmed_share-1)>1e-10))
  stop('LWJ statuses fail composition check')
D[, lwj_object_id := cell$lwj_object_id[match(cell_id, cell$cell_id)]]
for (v in c('lwj_contested_share','lwj_taliban_control_share',
            'lwj_goa_control_share','lwj_unconfirmed_share','lwj_confirmed_contested_share'))
  D[, (v) := P[[v]][match(lwj_object_id,P$lwj_object_id)]]
if (anyNA(D$lwj_contested_share)) stop('Missing mapped LWJ district histories')
if (nrow(D)!=104112L) stop('Panel rows unexpectedly changed')
D[,cell_id:=factor(cell_id)]
D[,month:=as.IDate(month)]
fwrite(unique(D[,.(cell_id,district_id_for_h2,lwj_object_id,
                    lwj_contested_share,lwj_taliban_control_share,
                    lwj_goa_control_share,lwj_unconfirmed_share,
                    lwj_confirmed_contested_share)]),
       file.path(OUT,'02_cell_lwj_moderators.csv'))
cat('LWJ pre-Doha status shares prepared, no missing unconditional measures.\n')

fit_h2 <- function(X,moderator,expanded=FALSE){
  rhs <- c('post_closure:eli_final',paste0('post_closure:eli_final:',moderator),
           'prec_gpcp','temp','nlights_calib_mean')
  if (expanded) rhs <- c(rhs,paste0('post_closure:',moderator))
  form <- as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),'| cell_id + month'))
  fit <- fepois(form,data=X,notes=FALSE,warn=FALSE)
  term <- paste0('post_closure:eli_final:',moderator)
  z <- coeftable(fit,vcov=~cell_id+month)
  if (!(term %in% rownames(z))) stop('Three-way interaction dropped: ',term)
  list(fit=fit,beta=unname(z[term,1]),se=unname(z[term,2]),
       p=unname(z[term,4]),n=as.integer(nobs(fit)))
}


# Common helpers for supplementary audits (not an estimator selection procedure).
DIAG_VARS <- c('lwj_contested_share','lwj_taliban_control_share','lwj_unconfirmed_share')
DIAG_TRIPLES <- paste0('post_closure:eli_final:',DIAG_VARS)
DIAG_LOWER <- paste0('post_closure:',DIAG_VARS)
DIAG_CONTROLS <- c('prec_gpcp','temp','nlights_calib_mean')
DIAG_FORM <- as.formula(paste('events_type3 ~', paste(c('post_closure:eli_final',DIAG_TRIPLES,DIAG_LOWER, DIAG_CONTROLS),collapse=' + '),'| cell_id + month'))
DIAG_WINDOWS <- list(full=D, drop_anticipation=D[month<as.IDate('2020-03-01') | month>=as.IDate('2021-07-01')])
DIAG_TERM <- 'post_closure:eli_final:lwj_contested_share'
diag_fit <- function(dat,formula=DIAG_FORM) {
  f <- fixest::fepois(formula,data=dat,notes=FALSE,warn=FALSE)
  zz <- fixest::coeftable(f,vcov=~cell_id+month)
  if (!DIAG_TERM %in% rownames(zz)) stop('Key interaction dropped: ',DIAG_TERM)
  list(fit=f,beta=unname(zz[DIAG_TERM,1]),se=unname(zz[DIAG_TERM,2]),p=unname(zz[DIAG_TERM,4]),n=as.integer(nobs(f)))
}
cat('Common LWJ data ready: ',uniqueN(D$cell_id),' cells, ',uniqueN(D$district_id_for_h2),' assigned districts.\n',sep='')
