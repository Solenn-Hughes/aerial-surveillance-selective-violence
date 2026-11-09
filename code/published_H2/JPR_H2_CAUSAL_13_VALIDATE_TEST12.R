# JPR H2 — TEST 13: Validate the joint ACLED + LWJ result, NOT a causal proof.
# Requires: JPR_H2_LWJ_00_COMMON.R, LWJ CSV, documented crosswalk, FINAL30 panel.
# Place this script and those existing inputs directly in Windows Downloads.
# Run: source('C:/Users/USERNAME/Downloads/JPR_H2_CAUSAL_13_VALIDATE_TEST12.R')

source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT, '13_validate_joint_test12')
dir.create(DEST, recursive=TRUE, showWarnings=FALSE)

A <- 'pre_contested_share_spatial'
B <- 'lwj_contested_share'
T_A <- paste0('post_closure:eli_final:', A)
T_B <- paste0('post_closure:eli_final:', B)
# All models include the independent Post x moderator terms.
form_acled <- as.formula(paste0('events_type3 ~ post_closure:eli_final + post_closure:', A,
  ' + ', T_A, ' | cell_id + month'))
form_lwj <- as.formula(paste0('events_type3 ~ post_closure:eli_final + post_closure:', B,
  ' + ', T_B, ' | cell_id + month'))
form_joint <- as.formula(paste0('events_type3 ~ post_closure:eli_final + post_closure:', A,
  ' + post_closure:', B, ' + ', T_A, ' + ', T_B, ' | cell_id + month'))

# The controls were perfectly collinear in the prior 2015-2024 tests and
# made no difference to estimates or retained observations. Omit explicitly.
# This is NOT a claim that time-varying confounding is absent.
models <- list(ACLED_only=form_acled, LWJ_only=form_lwj, JOINT=form_joint)

runfit <- function(form, dat) {
  tryCatch(fixest::fepois(form, data=dat, notes=FALSE, warn=FALSE),
           error=function(e) e)
}
extract <- function(fit, vc, scenario, modelname, n_input, used_expected=NA_integer_) {
  if (inherits(fit,'error')) return(data.table(scenario=scenario, model=modelname,
       cluster=vc, term=NA_character_, status=conditionMessage(fit), n_input=n_input))
  warned <- character()
  cc <- tryCatch(withCallingHandlers(fixest::coeftable(fit, vcov=as.formula(vc)),
             warning=function(w) {warned <<- c(warned,conditionMessage(w)); invokeRestart('muffleWarning')}),
             error=function(e)e)
  if (inherits(cc,'error')) return(data.table(scenario=scenario, model=modelname,
       cluster=vc, term=NA_character_, status=conditionMessage(cc), n_input=n_input,
       n=as.integer(nobs(fit))))
  out <- lapply(c(T_A,T_B), function(term) {
    z <- if (term %in% rownames(cc)) as.numeric(cc[term,]) else rep(NA_real_,4)
    data.table(scenario=scenario,model=modelname,cluster=vc,term=term,
      beta=z[1],se=z[2],p=z[4],n=as.integer(nobs(fit)),n_input=n_input,
      kept_same_n=if(is.na(used_expected)) NA else as.integer(nobs(fit))==used_expected,
      warning=if(length(warned)) paste(unique(warned),collapse=' | ') else '',
      status=if(term %in% rownames(cc)) 'ok' else 'not_estimable')
  })
  rbindlist(out,fill=TRUE)
}

ALL <- list(); SAMPLE <- list(); DROP <- list()
for (win in c('full_2015_2024','drop_anticipation_2015_2024')) {
  X <- copy(D[month >= as.IDate('2015-01-01') & month <= as.IDate('2024-12-01')])
  if (win=='drop_anticipation_2015_2024')
    X <- X[month < as.IDate('2020-03-01') | month >= as.IDate('2021-07-01')]
  X[, rid13 := .I]
  joint0 <- runfit(form_joint,X)
  if (inherits(joint0,'error')) stop('Joint baseline failed: ',conditionMessage(joint0))
  retained <- fixest::obs(joint0)
  if(length(retained)!=nobs(joint0)) stop('fixest::obs indexing inconsistent with N')
  Y <- copy(X[retained])
  fwrite(X[,.(rid13,cell_id,month,district_id_for_h2,
             retained=rid13 %in% retained)],
         file.path(DEST,paste0(win,'_observation_audit.csv')))
  SAMPLE[[length(SAMPLE)+1L]] <- data.table(window=win,input_rows=nrow(X),
       retained_joint=nrow(Y),dropped_joint=nrow(X)-nrow(Y),
       cells_retained=uniqueN(Y$cell_id),districts_retained=uniqueN(Y$district_id_for_h2),
       cells_with_positive_outcome=uniqueN(Y[events_type3>0]$cell_id),
       retained_total_events=sum(Y$events_type3))
  cat('\nWINDOW ',win,': input=',nrow(X),' joint retained=',nrow(Y),
       ' cells=',uniqueN(Y$cell_id),' districts=',uniqueN(Y$district_id_for_h2),'\n',sep='')
  for (nam in names(models)) {
    fit <- runfit(models[[nam]],Y)
    if (!inherits(fit,'error')) {
      for (vc in c('~cell_id+month','~district_id_for_h2+month'))
        ALL[[length(ALL)+1L]] <- extract(fit,vc,paste0(win,'_COMMON_SAMPLE'),nam,nrow(X),nrow(Y))
      # Actual row identity must match, not merely N.
      ids <- fixest::obs(fit)
      DROP[[length(DROP)+1L]] <- data.table(window=win,model=nam,
        same_rows=identical(sort(ids),seq_len(nrow(Y))),
        n=as.integer(nobs(fit)),expected=nrow(Y))
    } else ALL[[length(ALL)+1L]] <- extract(fit,'~cell_id+month',paste0(win,'_COMMON_SAMPLE'),nam,nrow(X))
  }
  # Reproduce prior joint-model result exactly before adding sensitivities.
  for (vc in c('~cell_id+month','~district_id_for_h2+month'))
    ALL[[length(ALL)+1L]] <- extract(joint0,vc,paste0(win,'_ORIGINAL_SAMPLE'),
                                     'JOINT',nrow(X),nobs(joint0))

  # Geographic sensitivity: exact-name subset is different geography, NOT
  # definitive evidence of a mislinked district. Compare and label accordingly.
  exact_ids <- CW[match_basis=='exact_name_province',adm2_id]
  EXACT <- X[district_id_for_h2 %in% exact_ids]
  if (nrow(EXACT)>0L) {
    fe <- runfit(form_joint,EXACT)
    ALL[[length(ALL)+1L]] <- extract(fe,'~district_id_for_h2+month',
                          paste0(win,'_EXACT_NAME_SUBPOPULATION'), 'JOINT',nrow(EXACT))
  }
  # Temporal concentration: post-year truncations retain unchanged pre-period.
  for (year_end in c(2021L,2022L,2023L)) {
    sub <- X[month <= as.IDate(sprintf('%d-12-01',year_end))]
    fs <- runfit(form_joint,sub)
    ALL[[length(ALL)+1L]] <- extract(fs,'~district_id_for_h2+month',
               paste0(win,'_ENDING_',year_end),'JOINT',nrow(sub))
  }
}
RESULTS <- rbindlist(ALL,fill=TRUE)
SAMPLES <- rbindlist(SAMPLE,fill=TRUE)
IDENTICAL <- rbindlist(DROP,fill=TRUE)
fwrite(RESULTS,file.path(DEST,'01_all_estimates.csv'))
fwrite(SAMPLES,file.path(DEST,'02_sample_audit.csv'))
fwrite(IDENTICAL,file.path(DEST,'03_common_sample_identity.csv'))
cat('\n=== BASELINE AND FIXED-SAMPLE TEST 12 COMPARISON ===\n')
print(RESULTS[grepl('COMMON_SAMPLE|ORIGINAL_SAMPLE',scenario) &
    term==T_B,.(scenario,model,cluster,beta,se,p,n,status)])
cat('\n=== IS THE COMMON RETAINED SAMPLE TRULY IDENTICAL? ===\n')
print(IDENTICAL)
cat('\n=== EXACT-MATCH SUBPOPULATION AND TIMING SENSITIVITIES ===\n')
print(RESULTS[grepl('EXACT_NAME|ENDING_',scenario) & term==T_B,
   .(scenario,cluster,beta,se,p,n,status)])
cat('\n=== SAMPLE LOSS ===\n');print(SAMPLES)
if (any(!IDENTICAL$same_rows)) cat('\nWARNING: Models do not use identical retained rows; DO NOT treat as same-sample evidence.\n')
if (any(grepl('VCOV|positive definite|not positive',RESULTS$warning,ignore.case=TRUE)))
 cat('\nWARNING: At least one covariance matrix required adjustment; inspect 01_all_estimates.csv.\n')
cat('\nINTERPRETATION: Significant coefficients are not proof of a surveillance mechanism.\n')
cat('No surveillance-intensity instrument or independent surveillance shock is in these data.\n')
cat('FILES:',DEST,'\n')
