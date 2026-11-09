# JPR H2 / 15 — joint ACLED-LWJ timing diagnostic with frozen 2017 LWJ classification.
# Calendar periods separate pre-closure, anticipation, and post-closure.
# This is a timing / pre-trend diagnostic, NOT proof of an ISR-specific causal shock.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST<-file.path(OUT,'15_joint_timing_event_study')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
FIRST<-L[period==1L,.(lwj_object_id,early_status=control)]
if(anyDuplicated(FIRST$lwj_object_id))stop('Duplicated first-period source district')
X<-copy(D[month>=as.IDate('2018-01-01') & month<=as.IDate('2024-12-01')])
X[,early_status:=FIRST$early_status[match(lwj_object_id,FIRST$lwj_object_id)]]
if(anyNA(X$early_status))stop('Missing 2017 LWJ classification')
X[,`:=`(early_contested=as.integer(early_status=='Contested'),
        early_taliban=as.integer(early_status=='Taliban'),
        early_unconfirmed=as.integer(early_status=='Unconfirmed'))]
# The ACLED moderator is inherited from the original panel and uses information through 2021.
# As such, estimates of pre-trends conditional on it are diagnostic, NOT genuinely prospective.
X[,`:=`(eli_acled=eli_final*pre_contested_share_spatial,
        eli_lwj=eli_final*early_contested,
        eli_taliban=eli_final*early_taliban,
        eli_unconfirmed=eli_final*early_unconfirmed)]
X[,period:=fifelse(month<as.IDate('2019-01-01'),'2018',
              fifelse(month<as.IDate('2020-01-01'),'2019',
              fifelse(month<as.IDate('2020-03-01'),'2020_JanFeb',
              fifelse(month<as.IDate('2021-07-01'),'anticipation',
              fifelse(month<as.IDate('2022-01-01'),'2021_H2',
              fifelse(month<as.IDate('2023-01-01'),'2022',
              fifelse(month<as.IDate('2024-01-01'),'2023','2024')))))))]
levels_period<-c('2018','2019','2020_JanFeb','anticipation','2021_H2','2022','2023','2024')
X[,period:=factor(period,levels=levels_period)]
if(anyNA(X$period))stop('Unclassified calendar period')

# Lower-order period x raw ACLED and LWJ status and period x ELI effects included.
# Pre-shock 2019 is the reference period. Month FE absorb common month changes.
terms<-c('eli_final','pre_contested_share_spatial','early_contested',
         'early_taliban','early_unconfirmed','eli_acled','eli_lwj',
         'eli_taliban','eli_unconfirmed')
F<-as.formula(paste('events_type3 ~',
  paste(sprintf('i(period, %s, ref="2019")',terms),collapse=' + '),
  '| cell_id + month'))
results<-list(); wald_results<-list(); sample_results<-list()
for(w in c('full','exclude_anticipation')) {
  Y<-if(w=='full')copy(X) else copy(X[period!='anticipation'])
  model<-tryCatch(fixest::fepois(F,data=Y,notes=FALSE,warn=TRUE),error=function(e)e)
  if(inherits(model,'error')) {
    results[[w]]<-data.table(window=w,status=conditionMessage(model))
    next
  }
  vc<-tryCatch(fixest::vcov(model,vcov=~district_id_for_h2+month),error=function(e)e)
  if(inherits(vc,'error')) {results[[w]]<-data.table(window=w,status=conditionMessage(vc));next}
  b<-coef(model)
  labs<-names(b)
  # Built-in i() names commonly use period::2018:eli_lwj, but save all names for verification.
  cc<-tryCatch(fixest::coeftable(model,vcov=~district_id_for_h2+month),error=function(e)e)
  if(inherits(cc,'error'))stop('coeftable failed: ',conditionMessage(cc))
  Z<-data.table(window=w,term=rownames(cc),beta=cc[,1],se=cc[,2],p=cc[,4],
                n=nobs(model))
  results[[w]]<-Z
  saveRDS(model,file.path(DEST,paste0(w,'_fit.rds')))
  sample_results[[w]]<-data.table(window=w,input_rows=nrow(Y),estimated_n=nobs(model),
     retained_cells=uniqueN(Y[fixest::obs(model),cell_id]),
     retained_districts=uniqueN(Y[fixest::obs(model),district_id_for_h2]))
  # Exact 2018 vs 2019 pre-period comparison, both before closure.
  # 2020 Jan-Feb comparison is potentially informative but has only 2 months.
  for(pr in c('2018','2020_JanFeb')) {
    hits<-labs[grepl(paste0('period::',pr,':eli_lwj'),labs,fixed=TRUE)]
    if(length(hits)!=1L) {
      wald_results[[length(wald_results)+1L]]<-data.table(window=w,pre_period=pr,
        status='not_estimable_or_label_mismatch',term=paste(hits,collapse=';'))
      next
    }
    nm<-hits[1]; se<-sqrt(vc[nm,nm]); z<-b[nm]/se
    wald_results[[length(wald_results)+1L]]<-data.table(window=w,pre_period=pr,
      status='ok',term=nm,beta=unname(b[nm]),se=se,p=2*pnorm(-abs(z)))
  }
  cat('\n=== ',w,' JOINT EVENT-STUDY: LWJ X ELI DIFFERENCES VS 2019 ===\n',sep='')
  print(Z[grepl(':eli_lwj',term,fixed=TRUE),.(term,beta,se,p,n)])
}
RES<-rbindlist(results,fill=TRUE)
WD<-rbindlist(wald_results,fill=TRUE)
SS<-rbindlist(sample_results,fill=TRUE)
fwrite(RES,file.path(DEST,'01_all_period_interactions.csv'))
fwrite(WD,file.path(DEST,'02_pre_2018_and_2020_comparisons.csv'))
fwrite(SS,file.path(DEST,'03_event_study_sample_sizes.csv'))
cat('\n=== PRE-CLOSURE CHECKS ===\n');print(WD)
cat('\n=== SAMPLE SIZES ===\n');print(SS)
cat('CAUTION: 2017 status snapshot was observed in late 2017 and may be misclassified.\n')
cat('CAUTION: ACLED moderator uses data through 2021, including later than 2018-19.\n')
cat('CAUTION: No evidence about ISR-specific exclusion restrictions is supplied by this test.\n')
cat('DONE: ',DEST,'\n')
