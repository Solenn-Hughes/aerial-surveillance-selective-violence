# 03: targeted district influence analysis; not a coefficient-significance search.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
OUT3 <- file.path(OUT,'03_influence')
dir.create(OUT3,recursive=TRUE,showWarnings=FALSE)
# Examine all unconfirmed-status districts plus biggest outcome-contributing districts;
# fixed selection uses levels/outcomes, NOT regression p-values.
K <- unique(D[,.(cell_id,district_id_for_h2,lwj_unconfirmed_share)])
C <- D[,.(events=sum(events_type3,na.rm=TRUE)),by=district_id_for_h2]
U <- unique(K[lwj_unconfirmed_share>0,district_id_for_h2])
target <- unique(c(U,C[order(-events),head(district_id_for_h2,15L)]))
cat('Districts in leave-one-out audit: ',length(target),' (unconfirmed or top 15 by outcome counts)\n',sep='')
fwrite(data.table(district_id_for_h2=target),file.path(OUT3,'01_target_districts.csv'))
R<-list()
for(w in names(DIAG_WINDOWS)) {
  X <- DIAG_WINDOWS[[w]]
  baseline <- diag_fit(X)
  R[[length(R)+1L]]<-data.table(window=w,omitted_district='NONE',beta=baseline$beta,se=baseline$se,p=baseline$p,n=baseline$n,cells=uniqueN(X$cell_id),status='ok')
  for(id in target){
    Y <- X[district_id_for_h2!=id]
    rr <- tryCatch(diag_fit(Y),error=function(e){warning(w,' missing ',id,': ',conditionMessage(e));NULL})
    if(is.null(rr)){
      R[[length(R)+1L]] <- data.table(window=w,omitted_district=id,beta=NA_real_,se=NA_real_,p=NA_real_,n=nrow(Y),cells=uniqueN(Y$cell_id),status='fit_failed')
    } else {
      R[[length(R)+1L]] <- data.table(window=w,omitted_district=id,beta=rr$beta,se=rr$se,p=rr$p,n=rr$n,cells=uniqueN(Y$cell_id),status='ok')
    }
  }
  cat('Completed ',w,': ',length(target),' leave-one-district-out attempts.\n',sep='')
}
T<-rbindlist(R,fill=TRUE)
fwrite(T,file.path(OUT3,'02_district_influence.csv'))
print(T[,.(window,omitted_district,beta=round(beta,3),se=round(se,3),p=round(p,4),status)])
cat('Influence output: ',OUT3,'\n')
