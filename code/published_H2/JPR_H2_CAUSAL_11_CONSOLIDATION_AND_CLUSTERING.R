# JPR H2 / 11: Province-by-month adjustments and inference at the moderator district level.
# Province-month FE absorb province-wide changes; district-month model is stricter (script 08).
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST<-file.path(OUT,'11_province_month_and_clustering');dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
# Province names come from existing documented crosswalk, not generated from ADM2 code.
if(!'shapefile_province'%in%names(CW))stop('No province identifier in crosswalk')
X<-copy(D[month>=as.IDate('2015-01-01')&month<=as.IDate('2024-12-01')])
X[,province:=CW$shapefile_province[match(district_id_for_h2,CW$adm2_id)]]
if(anyNA(X$province))stop('Missing province matches')
X[,province_month:=interaction(province,as.character(month),drop=TRUE)]
# Adjust for status-dependent post changes. province-month absorbs national month FE.
rhs<-c('post_closure:eli_final',DIAG_TRIPLES,DIAG_LOWER)
FORM<-as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),'|cell_id + province_month'))
baseForm<-as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),'|cell_id + month'))
R<-list()
for(w in c('full','drop_anticipation')){
 Y<-if(w=='full')X else X[month<as.IDate('2020-03-01') |month>=as.IDate('2021-07-01')]
 for(s in c('month_FE','province_by_month_FE')){
  fm<-if(s=='month_FE')baseForm else FORM
  fit<-tryCatch(fixest::fepois(fm,data=Y,notes=FALSE,warn=FALSE),error=function(e)e)
  if(inherits(fit,'error')){R[[length(R)+1L]]<-data.table(window=w,spec=s,status=conditionMessage(fit));next}
  for(vv in c('cell_by_month','district_by_month')){
   ct<-tryCatch(fixest::coeftable(fit,vcov=if(vv=='cell_by_month') ~cell_id+month else ~district_id_for_h2+month),error=function(e)e)
   if(inherits(ct,'error')){R[[length(R)+1L]]<-data.table(window=w,spec=s,cluster=vv,status=conditionMessage(ct));next}
   zz<-if(DIAG_TERM%in%rownames(ct))unname(ct[DIAG_TERM,]) else rep(NA_real_,4)
   R[[length(R)+1L]]<-data.table(window=w,spec=s,cluster=vv,beta=zz[1],se=zz[2],p=zz[4],
     n=as.integer(nobs(fit)),input_rows=nrow(Y),status='ok')
  }
  saveRDS(fit,file.path(DEST,paste0(w,'_',s,'.rds')))
 }
}
R<-rbindlist(R,fill=TRUE);print(R);fwrite(R,file.path(DEST,'01_province_month_cluster.csv'))
cat('DONE:',DEST,'\n')
