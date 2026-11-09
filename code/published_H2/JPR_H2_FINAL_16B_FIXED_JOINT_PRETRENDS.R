# JPR H2 FINAL CHECK 16 — joint, pre-specified pre-treatment Wald tests
# Uses frozen first LWJ assessment (late 2017) and existing FINAL30 inputs.
# Pre-periods: calendar 2018 and Jan-Feb 2020, both compared with calendar 2019.
# NOTE: ACLED contestedness was assembled through June 2021; it is NOT
# prospectively predetermined for the 2018/2020 comparisons.
# This is a diagnostic for differential trends, NOT proof of parallel trends.
# Keep this file alongside JPR_H2_LWJ_00_COMMON.R and its two LWJ CSV inputs.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT,'16_FINAL_joint_pretrend_Wald')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)

FIRST <- L[period == 1L, .(lwj_object_id, first_status=control)]
stopifnot(!anyDuplicated(FIRST$lwj_object_id))
X <- copy(D[month >= as.IDate('2018-01-01') & month <= as.IDate('2024-12-01')])
X[, first_status := FIRST$first_status[match(lwj_object_id,FIRST$lwj_object_id)]]
if(anyNA(X$first_status)) stop('Frozen LWJ status missing: check district join.')
X[, `:=`(first_contested=as.numeric(first_status=='Contested'),
         first_taliban=as.numeric(first_status=='Taliban'),
         first_unconfirmed=as.numeric(first_status=='Unconfirmed'))]
X[, `:=`(eli_acled=eli_final*pre_contested_share_spatial,
         eli_lwj=eli_final*first_contested,
         eli_taliban=eli_final*first_taliban,
         eli_unconfirmed=eli_final*first_unconfirmed)]
X[, phase := data.table::fifelse(month < as.IDate('2019-01-01'),'2018',
     data.table::fifelse(month < as.IDate('2020-01-01'),'2019',
     data.table::fifelse(month < as.IDate('2020-03-01'),'2020_JanFeb',
     data.table::fifelse(month < as.IDate('2021-07-01'),'anticipation',
     data.table::fifelse(month < as.IDate('2022-01-01'),'2021_H2',
     data.table::fifelse(month < as.IDate('2023-01-01'),'2022',
     data.table::fifelse(month < as.IDate('2024-01-01'),'2023','2024')))))))]
X[,phase:=factor(phase,levels=c('2018','2019','2020_JanFeb','anticipation','2021_H2','2022','2023','2024'))]
if(anyNA(X$phase)) stop('Phase construction failed.')
X[,log1p_events:=log1p(events_type3)]
# Parallel to Test 15B, with lower-order interactions included for each phase.
SLOPES <- c('eli_final','pre_contested_share_spatial',
            'first_contested','first_taliban','first_unconfirmed',
            'eli_acled','eli_lwj','eli_taliban','eli_unconfirmed')
make_formula <- function(y) as.formula(paste0(y,' ~ ',
  paste(sprintf('i(phase, %s, ref="2019")',SLOPES),collapse=' + '),
  ' | cell_id + month'))
F_P <- make_formula('events_type3')
F_O <- make_formula('log1p_events')

# Joint Wald test using the SAME covariance matrix used for individual terms.
# Records non-positive-definite or automatically repaired covariance warnings.
# Correction: stats::vcov dispatches to fixest's model method; fixest::vcov is not exported.
wald_terms <- function(fit,terms,vc,scenario,family,test_label,warning_text) {
  vv_warn <- character()
  V <- tryCatch(withCallingHandlers(
    stats::vcov(fit,vcov=vc),
    warning=function(w){vv_warn <<- c(vv_warn,conditionMessage(w));
                        invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(V,'error')) return(data.table(scenario=scenario,family=family,
    test=test_label,status='VCOV_ERROR',error=conditionMessage(V),
    chi_square=NA_real_,df=length(terms),p=NA_real_,estimated_rows=nobs(fit),covariance_warning=TRUE,
    warnings=paste(c(warning_text,vv_warn),collapse=' | ')))
  missing <- setdiff(terms, names(coef(fit)))
  if(length(missing)) return(data.table(scenario=scenario,family=family,
    test=test_label,status='TERMS_NOT_ESTIMABLE',error=paste(missing,collapse='; '),
    chi_square=NA_real_,df=length(terms),p=NA_real_,estimated_rows=nobs(fit),covariance_warning=NA,
    warnings=paste(c(warning_text,vv_warn),collapse=' | ')))
  b <- as.numeric(coef(fit)[terms]); C <- as.matrix(V[terms,terms,drop=FALSE])
  ev <- eigen((C+t(C))/2,symmetric=TRUE,only.values=TRUE)$values
  if(any(!is.finite(ev)) || min(ev)<=1e-10) return(data.table(
    scenario=scenario,family=family,test=test_label,status='SINGULAR_VCOV',
    error='Covariance submatrix is not positive definite',minimum_eigenvalue=min(ev),
    chi_square=NA_real_,df=length(terms),p=NA_real_,estimated_rows=nobs(fit),covariance_warning=TRUE,
    warnings=paste(c(warning_text,vv_warn),collapse=' | ')))
  W <- as.numeric(t(b)%*%solve(C,b)); p <- pchisq(W,df=length(terms),lower.tail=FALSE)
  data.table(scenario=scenario,family=family,test=test_label,status='ok',
    chi_square=W,df=length(terms),p=p,input_rows=NA_integer_,
    estimated_rows=nobs(fit),minimum_eigenvalue=min(ev),
    covariance_warning=any(grepl('positive definite|fixed|adjust',c(warning_text,vv_warn),ignore.case=TRUE)),
    warnings=paste(unique(c(warning_text,vv_warn)),collapse=' | '),error='')
}

TESTS <- list(
  'LWJ_H2_preperiods_2df'=c('phase::2018:eli_lwj','phase::2020_JanFeb:eli_lwj'),
  'ACLED_and_LWJ_H2_preperiods_4df'=c('phase::2018:eli_lwj','phase::2020_JanFeb:eli_lwj',
                                     'phase::2018:eli_acled','phase::2020_JanFeb:eli_acled'))
all_tests <- list(); all_terms <- list(); sample <- list(); errors <- list()
for(win in c('full','exclude_anticipation')) {
  Y <- if(win=='full') copy(X) else copy(X[phase!='anticipation'])
  for(family in c('PPML','OLS_log1p')) {
    key <- paste(win,family,sep='__'); warnings_fit <- character()
    fit <- tryCatch(withCallingHandlers(
      if(family=='PPML') fepois(F_P,data=Y,notes=FALSE,warn=TRUE)
      else feols(F_O,data=Y,notes=FALSE,warn=TRUE),
      warning=function(w){warnings_fit <<- c(warnings_fit,conditionMessage(w));
                          invokeRestart('muffleWarning')}),error=function(e)e)
    if(inherits(fit,'error')) {
      msg <- conditionMessage(fit); cat('FAILED ',key,': ',msg,'\n',sep='')
      errors[[length(errors)+1L]] <- data.table(run=key,error=msg)
      next
    }
    if(nobs(fit)<1) stop('No observations retained in ',key)
    vc <- ~district_id_for_h2 + month
    ct_warn<-character()
    ct<-tryCatch(withCallingHandlers(fixest::coeftable(fit,vcov=vc),
      warning=function(w){ct_warn<<-c(ct_warn,conditionMessage(w));
                          invokeRestart('muffleWarning')}),error=function(e)e)
    if(inherits(ct,'error')) stop('Covariance error in ',key,': ',conditionMessage(ct))
    z <- data.table(term=rownames(ct),beta=as.numeric(ct[,1]),se=as.numeric(ct[,2]),
                    p=as.numeric(ct[,4]))
    z <- z[grepl('phase::(2018|2020_JanFeb):eli_(lwj|acled)$',term)]
    z[,`:=`(scenario=win,family=family)]
    all_terms[[key]]<-z
    sample[[key]]<-data.table(scenario=win,family=family,input_rows=nrow(Y),
       estimated_rows=nobs(fit),retained_cells=uniqueN(Y$cell_id[fixest::obs(fit)]),
       retained_events=sum(Y$events_type3[fixest::obs(fit)]),
       fit_warnings=paste(unique(warnings_fit),collapse=' | '),
       vcov_warnings=paste(unique(ct_warn),collapse=' | '))
    for(test_name in names(TESTS)){
      a <- wald_terms(fit,TESTS[[test_name]],vc,win,family,test_name,
                      paste(unique(c(warnings_fit,ct_warn)),collapse=' | '))
      a[,input_rows:=nrow(Y)]
      all_tests[[paste(key,test_name,sep='__')]]<-a
    }
    cat('\n',key,': input=',nrow(Y),' retained=',nobs(fit),'\n',sep='')
    print(z[,.(term,beta,se,p)])
  }
}
if(!length(all_tests)) stop('All models failed; see printed errors.')
RES<-rbindlist(all_tests,fill=TRUE); TERMS<-rbindlist(all_terms,fill=TRUE)
SAMPLE<-rbindlist(sample,fill=TRUE)
fwrite(RES,file.path(DEST,'01_JOINT_WALD_PRETRENDS.csv'))
fwrite(TERMS,file.path(DEST,'02_INDIVIDUAL_PRE_COEFFICIENTS.csv'))
fwrite(SAMPLE,file.path(DEST,'03_SAMPLE_AND_VCOV_AUDIT.csv'))
fwrite(if(length(errors))rbindlist(errors) else data.table(run=character(),error=character()),
       file.path(DEST,'04_MODEL_ERRORS.csv'))
cat('\n=== FINAL JOINT PRE-TREND TESTS ===\n')
print(RES[,.(scenario,family,test,status,chi_square,df,p,estimated_rows,minimum_eigenvalue,covariance_warning,warnings,error)])
cat('\nCAUTION: non-rejection is not proof of parallel trends.\n')
cat('Frozen LWJ assessment dates from late 2017; ACLED measure includes information through 2021.\n')
cat('Any covariance warnings must be assessed before relying on p-values.\n')
cat('These tests are NOT identical to the dissertation Doha H1 Wald (3 df, p=0.062).\n')
cat('Output: ',DEST,'\n',sep='')
