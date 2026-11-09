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
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_LWJ_MATCHED_20261106')
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

bench <- data.table(
  window=rep(c('full','drop_anticipation'),each=3),
  moderator=rep(c('pre_contested_share_spatial',
                  'pre_taliban_activity_share_spatial',
                  'pre_contested_intensity_spatial'),times=2),
  beta_expected=c(-16.1023802,-17.6667523,-3.120677,
                  -17.40820781,-19.07071,-3.38898),
  se_expected=c(4.273193509,4.60318,0.88088,
                4.867393903,5.26816,0.98456),
  n_expected=rep(c(54312L,49956L),each=3))
windows <- list(full=D,
  drop_anticipation=D[month < as.IDate('2020-03-01') |
                      month >= as.IDate('2021-07-01')])
cat('\n=== REQUIRED ORIGINAL FINAL30 REPRODUCTION ===\n')
reps <- list()
for (w in names(windows)) for (m in unique(bench$moderator)) {
  f <- fit_h2(windows[[w]],m,FALSE)
  b <- bench[window==w & moderator==m]
  if (nrow(b)!=1L) stop('Benchmark lookup ambiguous')
  good <- abs(f$beta-b$beta_expected)<0.002 &&
          abs(f$se-b$se_expected)<0.002 && f$n==b$n_expected
  reps[[length(reps)+1L]] <- data.table(window=w,moderator=m,
       beta=f$beta,se=f$se,n=f$n,replicated=good)
  cat(w,m,'beta=',round(f$beta,4),'SE=',round(f$se,4),
      'N=',f$n,'reproduced=',good,'\n')
  if (!good) stop('Published model not reproduced; LWJ results NOT estimated')
}
fwrite(rbindlist(reps),file.path(OUT,'03_original_H2_replication.csv'))

# PRIMARY: independent LWJ pre-Doha territorial-contestation measure.
# SECONDARY: Taliban-control and GOA-control shares, and exclusion of unconfirmed periods.
mods <- c('lwj_contested_share','lwj_taliban_control_share',
          'lwj_goa_control_share','lwj_confirmed_contested_share')
results <- list()
cat('\n=== LWJ H2: RESTRICTED vs FULLY SPECIFIED ===\n')
for (w in names(windows)) for (m in mods) {
  X <- windows[[w]]
  if (anyNA(X[[m]])) {
    cat('SKIP: ',w,' / ',m,' contains unconfirmed-only districts; no listwise deletion concealed.\n')
    next
  }
  a <- tryCatch(fit_h2(X,m,FALSE),error=function(e)e)
  b <- tryCatch(fit_h2(X,m,TRUE),error=function(e)e)
  if (inherits(a,'error') || inherits(b,'error')) {
    results[[length(results)+1L]] <- data.table(window=w,moderator=m,
       specification=c('restricted','expanded'),beta=NA_real_,se=NA_real_,
       p=NA_real_,n=NA_integer_,error=if(inherits(a,'error'))conditionMessage(a) else conditionMessage(b))
    next
  }
  aa <- tryCatch(as.integer(fixest::obs(a$fit)),error=function(e)integer(0))
  bb <- tryCatch(as.integer(fixest::obs(b$fit)),error=function(e)integer(0))
  same <- length(aa)>0 && identical(aa,bb)
  if (!same || a$n!=b$n) stop('Different model samples for ',w,' / ',m)
  for (spec in c('restricted','expanded')) {
    z <- if(spec=='restricted') a else b
    results[[length(results)+1L]] <- data.table(window=w,moderator=m,
       specification=spec,beta=z$beta,se=z$se,p=z$p,n=z$n,
       same_estimation_rows=same,error='')
    cat(w,' / ',m,' / ',spec,' beta=',round(z$beta,4),
        ' SE=',round(z$se,4),' p=',signif(z$p,4),' N=',z$n,'\n',sep='')
  }
}
R <- rbindlist(results,fill=TRUE)
fwrite(R,file.path(OUT,'04_LWJ_H2_models.csv'))
print(R[,.(window,moderator,specification,beta,se,p,n,same_estimation_rows)])
cat('\nDONE. Outputs at: ',OUT,'\n',sep='')
cat('PRIMARY INTERPRETATION: expanded lwj_contested_share result; compare its sign, SE and p.\n')
cat('LWJ status proxies territorial control, not directly observed access to informants.\n')
