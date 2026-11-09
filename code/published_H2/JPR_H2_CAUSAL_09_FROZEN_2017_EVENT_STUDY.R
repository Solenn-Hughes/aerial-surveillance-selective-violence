# JPR H2 / 09: Prospective LWJ moderator frozen at first observed 2017 snapshot.
# Differential annual ELI slopes relative to 2019; no future status for 2018 placebo.
# 2021 combines pre- and post-closure months: DO NOT read as a clean post-treatment year.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST<-file.path(OUT,'09_frozen_2017_event_study');dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
FIRST<-L[period==1L,.(lwj_object_id,early_control=control)]
if(anyDuplicated(FIRST$lwj_object_id))stop('Duplicate period-one district')
X<-copy(D[month >= as.IDate('2018-01-01') & month <= as.IDate('2024-12-01')])
X[,early_control:=FIRST$early_control[match(lwj_object_id,FIRST$lwj_object_id)]]
if(anyNA(X$early_control))stop('Some first-period statuses missing')
for (v in c('Contested','Taliban','Unconfirmed')) {
 nm<-paste0('early_',tolower(v)); X[,(nm):=as.integer(early_control==v)]
 X[,(paste0('eli_',tolower(v))):=eli_final*get(nm)]
}
X[,year:=as.integer(format(month,'%Y'))]
# By-year lower-order status effects absorb flexible changes in status-related baselines.
# By-year ELI effects absorb country-wide time-varying ELI gradients.
xs<-c('eli_final','early_contested','early_taliban','early_unconfirmed',
       'eli_contested','eli_taliban','eli_unconfirmed')
ints<-paste0('i(year,',xs,',ref=2019)')
F<-as.formula(paste('events_type3 ~',paste(ints,collapse=' + '),'| cell_id + month'))
res<-list()
for(w in c('full','exclude_anticipation')) {
 Y<-if(w=='full')X else X[month < as.IDate('2020-03-01') | month >= as.IDate('2021-07-01')]
 fit<-tryCatch(fixest::fepois(F,data=Y,notes=FALSE,warn=FALSE),error=function(e)e)
 if(inherits(fit,'error')) {res[[length(res)+1L]]<-data.table(window=w,status=conditionMessage(fit));next}
 ct<-tryCatch(fixest::coeftable(fit,vcov=~district_id_for_h2+month),error=function(e)e)
 if(inherits(ct,'error')) {res[[length(res)+1L]]<-data.table(window=w,status=conditionMessage(ct));next}
 terms<-grep('eli_contested',rownames(ct),value=TRUE,fixed=TRUE)
 # Term labels depend on installed fixest version: save all model terms for audit.
 fwrite(data.table(term=rownames(ct),beta=unname(ct[,1]),se=unname(ct[,2]),p=unname(ct[,4])),
        file.path(DEST,paste0(w,'_all_coefficients.csv')))
 for(term in terms){v<-unname(ct[term,]);res[[length(res)+1L]]<-data.table(window=w,term=term,
  beta=v[1],se=v[2],p=v[4],n=as.integer(nobs(fit)),status='ok')}
 if(!length(terms))res[[length(res)+1L]]<-data.table(window=w,status='eli_contested interactions not estimated')
 saveRDS(fit,file.path(DEST,paste0(w,'_model.rds')))
 cat('\n',w,'estimated N=',nobs(fit),' contested-year terms=',length(terms),'\n',sep='')
}
R<-rbindlist(res,fill=TRUE);print(R);fwrite(R,file.path(DEST,'01_contested_year_interactions.csv'))
cat('First LWJ period is October 2017; verify district coding dates before claiming a prospective test.\n')
cat('2018 is the pre-treatment check, relative to 2019; zero difference does not prove parallel trends.\n')
cat('DONE:',DEST,'\n')
