# JPR H2 FINAL 18 -- SINGLE FROZEN, PROSPECTIVE H2 VALIDATION
# Does not overwrite prior models. Reports results irrespective of significance.
# Uses original district-month ACLED source, ending 2017-12-31, and first LWJ period.
# 2018 and Jan-Feb 2020 vs 2019 test pre-Doha contrasts, not an omnibus proof of trends.
# Note that LWJ first period is late 2017; ACLED pre-2018 is backward-looking.

source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
SRC <- file.path(ROOT,'outputs','paper1_geo_afghanistan_executable_acled',
                 'diagnostics_h2_spatial_adm2','district_month_contestedness_spatial_adm2.csv')
DEST <- file.path(OUT,'18_FINAL_frozen_prospective_validation')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
if (!file.exists(SRC)) stop('Missing original ACLED district-month file: ',SRC)
DM <- fread(SRC)
need <- c('district_id','month','contested_any')
if (length(bad <- setdiff(need,names(DM)))) stop('Missing ACLED fields: ',paste(bad,collapse=', '))
DM[, `:=`(district_id=as.character(district_id),month=as.IDate(month))]
if (anyDuplicated(DM,by=c('district_id','month'))) stop('Duplicate ACLED district-month rows')
if (anyNA(DM$contested_any)) stop('Missing contested_any; inspect data source')
if (any(!DM$contested_any %in% c(0,1))) stop('contested_any must be binary')
PRE <- DM[month < as.IDate('2018-01-01'),.(acled_frozen=mean(contested_any),
                                           months_acled=.N,
                                           first_month=min(month),
                                           last_month=max(month)),by=district_id]
if (!nrow(PRE)) stop('No ACLED pre-2018 history found')
FIRST <- L[period==1L,.(lwj_object_id,first_status=control)]
if (anyDuplicated(FIRST$lwj_object_id)) stop('Duplicated LWJ period-one district')
X <- copy(D[month>=as.IDate('2018-01-01') & month<=as.IDate('2024-12-01')])
X[,acled_frozen:=PRE$acled_frozen[match(as.character(district_id_for_h2),PRE$district_id)]]
X[,months_acled:=PRE$months_acled[match(as.character(district_id_for_h2),PRE$district_id)]]
X[,first_status:=FIRST$first_status[match(lwj_object_id,FIRST$lwj_object_id)]]
if (anyNA(X$acled_frozen) || anyNA(X$first_status)) stop('Unmatched early ACLED or LWJ district; do not fill missing as zero')
X[,`:=`(lwj_contested=as.numeric(first_status=='Contested'),
        lwj_taliban=as.numeric(first_status=='Taliban'),
        lwj_unconfirmed=as.numeric(first_status=='Unconfirmed'))]
X[,`:=`(eli_acled=eli_final*acled_frozen,
        eli_lwj=eli_final*lwj_contested,
        eli_taliban=eli_final*lwj_taliban,
        eli_unconfirmed=eli_final*lwj_unconfirmed,
        log1p_events=log1p(events_type3))]
if (anyNA(X[,.(eli_acled,eli_lwj,eli_taliban,eli_unconfirmed,log1p_events)])) stop('Missing frozen model variables')
X[,phase:=fifelse(month<as.IDate('2019-01-01'),'2018',
          fifelse(month<as.IDate('2020-01-01'),'2019',
          fifelse(month<as.IDate('2020-03-01'),'2020_JanFeb',
          fifelse(month<as.IDate('2021-07-01'),'anticipation',
          fifelse(month<as.IDate('2022-01-01'),'2021_H2',
          fifelse(month<as.IDate('2023-01-01'),'2022',
          fifelse(month<as.IDate('2024-01-01'),'2023','2024')))))))]
X[,phase:=factor(phase,levels=c('2018','2019','2020_JanFeb','anticipation','2021_H2','2022','2023','2024'))]
if(anyNA(X$phase)) stop('Unexpected phase')
CELL <- unique(X[,.(cell_id,district_id_for_h2,lwj_object_id,
                     first_status,acled_frozen,months_acled)])
if(anyDuplicated(CELL$cell_id))stop('Nonconstant cell mapping')
fwrite(CELL,file.path(DEST,'01_FROZEN_MODERATORS_AND_COVERAGE.csv'))
cat('\n=== FROZEN EXPOSURE COVERAGE ===\n')
print(CELL[,.(cells=.N,districts=uniqueN(district_id_for_h2),
               min_acled_months=min(months_acled),median_acled_months=median(months_acled),
               max_acled_months=max(months_acled),cells_under_12_months=sum(months_acled<12),
               correlation=cor(acled_frozen,as.numeric(first_status=='Contested')))])

# Match the prior phase-diagnostic lower order terms; no model chosen by p-value.
SLOPES <- c('eli_final','acled_frozen','lwj_contested','lwj_taliban',
            'lwj_unconfirmed','eli_acled','eli_lwj','eli_taliban','eli_unconfirmed')
PRE_TERMS <- list(
 LWJ_2df=c('phase::2018:eli_lwj','phase::2020_JanFeb:eli_lwj'),
 JOINT_4df=c('phase::2018:eli_lwj','phase::2020_JanFeb:eli_lwj',
             'phase::2018:eli_acled','phase::2020_JanFeb:eli_acled'))
P_FORM <- as.formula(paste('events_type3 ~',
     paste(sprintf('i(phase,%s,ref="2019")',SLOPES),collapse=' + '),
     '| cell_id + month'))
O_FORM <- update(P_FORM, log1p_events ~ .)
# Frozen main model includes independent post changes for ACLED and ALL status groups;
# GOA is the category reference. Both proxies are measured before 2018 outcomes.
MAIN_RHS <- c('post_closure:eli_final',
 'post_closure:acled_frozen','post_closure:lwj_contested',
 'post_closure:lwj_taliban','post_closure:lwj_unconfirmed',
 'post_closure:eli_final:acled_frozen',
 'post_closure:eli_final:lwj_contested',
 'post_closure:eli_final:lwj_taliban',
 'post_closure:eli_final:lwj_unconfirmed')
M_FORM <- as.formula(paste('events_type3 ~',paste(MAIN_RHS,collapse=' + '),'| cell_id + month'))

# Confidence intervals and joint Wald tests use stats::vcov; carry warnings forward.
# Two-way district + month is primary to match prior tests; district-only is a
# pre-labelled sensitivity if two-way has covariance repair warnings.
getV <- function(f, spec) {
 w <- character()
 v <- tryCatch(withCallingHandlers(stats::vcov(f,vcov=spec),
  warning=function(e){w<<-unique(c(w,conditionMessage(e)));invokeRestart('muffleWarning')}),
  error=function(e)e)
 list(v=v,w=w)
}
wald <- function(f,V,terms) {
 if(inherits(V,'error')) return(list(status='VCOV_ERROR',stat=NA_real_,p=NA_real_,min_eig=NA_real_))
 if(any(!terms %in% names(coef(f))))return(list(status='TERMS_DROPPED',stat=NA_real_,p=NA_real_,min_eig=NA_real_))
 C <- as.matrix(V[terms,terms,drop=FALSE]);C <- (C+t(C))/2
 eig <- eigen(C,symmetric=TRUE,only.values=TRUE)$values
 mn <- min(eig)
 if(any(!is.finite(eig)) || mn<=1e-10)return(list(status='INVALID_SUBMATRIX',stat=NA_real_,p=NA_real_,min_eig=mn))
 b <- as.numeric(coef(f)[terms]); stat <- as.numeric(crossprod(b,solve(C,b)))
 list(status='computed',stat=stat,p=pchisq(stat,length(terms),lower.tail=FALSE),min_eig=mn)
}
OUTPUT <- list(); COEFS <- list(); SAMPLES <- list(); k <- 0L;j <- 0L;s <- 0L
for (win in c('full','exclude_anticipation')) {
 Y <- if(win=='full') copy(X) else copy(X[phase!='anticipation'])
 for (analysis in c('MAIN_PPML','PRE_PPML','PRE_OLS')) {
  form <- switch(analysis, MAIN_PPML=M_FORM,PRE_PPML=P_FORM,PRE_OLS=O_FORM)
  fw <- character()
  f <- tryCatch(withCallingHandlers(
   if (analysis=='PRE_OLS') feols(form,data=Y,notes=FALSE,warn=TRUE)
   else fepois(form,data=Y,notes=FALSE,warn=TRUE),
   warning=function(e){fw<<-unique(c(fw,conditionMessage(e)));invokeRestart('muffleWarning')}),
   error=function(e)e)
  if(inherits(f,'error'))stop('MODEL FAILURE ',win,' ',analysis,': ',conditionMessage(f))
  o <- as.integer(fixest::obs(f))
  s<-s+1L
  SAMPLES[[s]] <- data.table(window=win,analysis=analysis,input_rows=nrow(Y),
    retained_rows=nobs(f),input_cells=uniqueN(Y$cell_id),
    retained_cells=uniqueN(Y$cell_id[o]),
    input_events=sum(Y$events_type3),retained_events=sum(Y$events_type3[o]),
    model_warnings=paste(fw,collapse=' | '))
  for (vcname in c('district_month','district_only')) {
   vc <- if(vcname=='district_month') ~district_id_for_h2+month else ~district_id_for_h2
   z <- getV(f,vc)
   if(inherits(z$v,'error'))stop('VCOV extraction failed ',win,' ',analysis,' ',vcname,': ',conditionMessage(z$v))
   two_way_repair <- any(grepl('positive definite|fixed|adjust',z$w,ignore.case=TRUE))
   targets <- if(analysis=='MAIN_PPML')list(MAIN_LWJ='post_closure:eli_final:lwj_contested') else PRE_TERMS
   for(test in names(targets)) {
    terms <- targets[[test]]
    if(analysis=='MAIN_PPML') {
     if(!terms %in% names(coef(f))) {
      ww<-list(status='TERMS_DROPPED',stat=NA_real_,p=NA_real_,min_eig=NA_real_)
     } else {
      C<-z$v[terms,terms,drop=FALSE]
      ww<-wald(f,z$v,terms)
     }
    } else ww<-wald(f,z$v,terms)
    k<-k+1L
    OUTPUT[[k]]<-data.table(window=win,analysis=analysis,vcov=vcname,test=test,
      status=ww$status,chi_square=ww$stat,df=length(terms),p=ww$p,
      min_eigenvalue=ww$min_eig,vcov_repaired=two_way_repair,
      vcov_warnings=paste(z$w,collapse=' | '),n=nobs(f))
   }
   tt <- if(analysis=='MAIN_PPML')c('post_closure:eli_final:lwj_contested')
         else unlist(PRE_TERMS,use.names=FALSE)
   tt<-unique(tt)
   for (t in tt) {
    if(!t %in% names(coef(f)))next
    j<-j+1L
    se<-sqrt(z$v[t,t]);b<-unname(coef(f)[t])
    COEFS[[j]]<-data.table(window=win,analysis=analysis,vcov=vcname,
      term=t,beta=b,se=se,z=b/se,p=2*pnorm(-abs(b/se)),n=nobs(f),vcov_repaired=two_way_repair)
   }
  }
 }
}
RES<-rbindlist(OUTPUT,fill=TRUE); TT<-rbindlist(COEFS,fill=TRUE); SS<-rbindlist(SAMPLES,fill=TRUE)
fwrite(RES,file.path(DEST,'02_JOINT_TESTS_AND_COVARIANCE.csv'))
fwrite(TT,file.path(DEST,'03_COEFFICIENTS.csv'))
fwrite(SS,file.path(DEST,'04_PPML_AND_OLS_SAMPLES.csv'))
cat('\n=== AUTHORITATIVE FROZEN MAIN H2 AND PRE-TREND JOINT TESTS ===\n')
print(RES[,.(window,analysis,vcov,test,status,chi_square,df,p,min_eigenvalue,vcov_repaired,n)])
cat('\n=== FROZEN LWJ COEFFICIENTS AND PRE-TREND COEFFICIENTS ===\n')
print(TT)
cat('\n=== SAMPLE AUDIT ===\n')
print(SS)
cat('\nInterpretation: the original main joint H2 (-4.997 / -6.310) uses pre-Doha LWJ shares,\n')
cat('whereas this diagnostic uses late-2017 LWJ status and pre-2018 ACLED histories.\n')
cat('Neither failure to reject pre-trends nor a covariance sensitivity establishes a uniquely ISR-causal effect.\n')
cat('DONE: ',DEST,'\n',sep='')
