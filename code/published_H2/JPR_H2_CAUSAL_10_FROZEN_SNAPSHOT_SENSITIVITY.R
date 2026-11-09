# JPR H2 / 10: Frozen October 2017 control vs 2017-19 averages, holding sample constant.
# The 2017 baseline uses one noisy snapshot; five-snapshot baseline is more precise but
# later-dated relative to some placebo periods. Not a randomised treatment.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST<-file.path(OUT,'10_frozen_snapshot_sensitivity');dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
FIRST<-L[period==1L,.(lwj_object_id,
    early_contested=as.integer(control=='Contested'),
    early_taliban=as.integer(control=='Taliban'),
    early_unconfirmed=as.integer(control=='Unconfirmed'))]
X<-copy(D[month>=as.IDate('2018-01-01') & month<=as.IDate('2024-12-01')])
for(v in c('early_contested','early_taliban','early_unconfirmed')) {
 X[,(v):=FIRST[[v]][match(lwj_object_id,FIRST$lwj_object_id)]]
 if(anyNA(X[[v]]))stop('Missing baseline classification: ',v)
}
models<-list(
  frozen_2017=c('early_contested','early_taliban','early_unconfirmed'),
  average_2017_2019=c('lwj_contested_share','lwj_taliban_control_share','lwj_unconfirmed_share'))
R<-list()
for(w in c('full','drop_anticipation')) {
 Y<-if(w=='full')X else X[month<as.IDate('2020-03-01')|month>=as.IDate('2021-07-01')]
 for(spec in names(models)){
  vv<-models[[spec]]
  rhs<-c('post_closure:eli_final',paste0('post_closure:eli_final:',vv),paste0('post_closure:',vv))
  fm<-as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),'|cell_id + month'))
  fit<-tryCatch(fixest::fepois(fm,data=Y,notes=FALSE,warn=FALSE),error=function(e)e)
  if(inherits(fit,'error')){R[[length(R)+1L]]<-data.table(window=w,spec=spec,status=conditionMessage(fit));next}
  ct<-tryCatch(fixest::coeftable(fit,vcov=~district_id_for_h2+month),error=function(e)e)
  if(inherits(ct,'error')){R[[length(R)+1L]]<-data.table(window=w,spec=spec,status=conditionMessage(ct));next}
  term<-paste0('post_closure:eli_final:',vv[1]);z<-if(term%in%rownames(ct))unname(ct[term,]) else rep(NA_real_,4)
  R[[length(R)+1L]]<-data.table(window=w,spec=spec,beta=z[1],se=z[2],p=z[4],
  n=as.integer(nobs(fit)),input_rows=nrow(Y),term_estimable=term%in%rownames(ct),status='ok')
  saveRDS(fit,file.path(DEST,paste0(w,'_',spec,'.rds')))
 }
}
RES<-rbindlist(R,fill=TRUE);print(RES);fwrite(RES,file.path(DEST,'01_frozen_vs_average.csv'))
cat('DONE:',DEST,'\n')
