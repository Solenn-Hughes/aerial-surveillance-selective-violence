# JPR H2 / 12: Do historical battle exposure and territorial control separately
# moderate Post x ELI once both are allowed independent post-change effects?
# These are correlated proxies, NOT exogenous instruments. Collinearity is reported.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST<-file.path(OUT,'12_joint_acled_lwj');dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
# Focus on one pre-existing ACLED moderator and one LWJ territorial measure.
# Do not choose among ACLED measures based on significance.
A<-'pre_contested_share_spatial'; B<-'lwj_contested_share'
BASE<-c('post_closure:eli_final',paste0('post_closure:',c(A,B)),
         paste0('post_closure:eli_final:',c(A,B)))
F<-as.formula(paste('events_type3 ~',paste(BASE,collapse=' + '),'|cell_id + month'))
R<-list()
for(w in c('full','drop_anticipation')){
 X<-copy(D[month>=as.IDate('2015-01-01') & month<=as.IDate('2024-12-01')])
 if(w=='drop_anticipation')X<-X[month<as.IDate('2020-03-01')|month>=as.IDate('2021-07-01')]
 fit<-tryCatch(fixest::fepois(F,data=X,notes=FALSE,warn=FALSE),error=function(e)e)
 if(inherits(fit,'error')){R[[length(R)+1L]]<-data.table(window=w,status=conditionMessage(fit));next}
 ct<-tryCatch(fixest::coeftable(fit,vcov=~district_id_for_h2+month),error=function(e)e)
 if(inherits(ct,'error')){R[[length(R)+1L]]<-data.table(window=w,status=conditionMessage(ct));next}
 for(t in paste0('post_closure:eli_final:',c(A,B))) {
  z<-if(t%in%rownames(ct))unname(ct[t,]) else rep(NA_real_,4)
  R[[length(R)+1L]]<-data.table(window=w,term=t,beta=z[1],se=z[2],p=z[4],
    n=as.integer(nobs(fit)),input_rows=nrow(X),status='ok')
 }
 saveRDS(fit,file.path(DEST,paste0(w,'_joint.rds')))
}
R<-rbindlist(R,fill=TRUE);print(R)
fwrite(R,file.path(DEST,'01_joint_acled_lwj.csv'))
C<-unique(D[,.(cell_id,pre_contested_share_spatial,lwj_contested_share)])
correlation<-cor(C$pre_contested_share_spatial,C$lwj_contested_share,use='complete.obs')
cat('Cell-level correlation between moderators:',round(correlation,4),'\n')
cat('WARNING: Comparison of two noisy correlated proxies is not causal identification.\n')
cat('DONE:',DEST,'\n')
