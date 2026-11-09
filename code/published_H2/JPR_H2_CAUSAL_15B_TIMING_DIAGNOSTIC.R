# H2 Test 15B: auditable, joint ACLED-LWJ calendar-phase diagnostics.
# Unlike original Test 15, EVERY failed model writes and prints its actual error.
# This is not, by itself, evidence of an ISR-specific causal shock.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT, '15B_joint_timing_repaired')
dir.create(DEST, recursive=TRUE, showWarnings=FALSE)

FIRST <- L[period == 1L, .(lwj_object_id, early_status = control)]
stopifnot(!anyDuplicated(FIRST$lwj_object_id))
X <- copy(D[month >= as.IDate('2018-01-01') & month <= as.IDate('2024-12-01')])
X[, early_status := FIRST$early_status[match(lwj_object_id, FIRST$lwj_object_id)]]
if (anyNA(X$early_status)) stop('Some early LWJ statuses are missing')
X[, `:=`(early_contested = as.integer(early_status == 'Contested'),
         early_taliban = as.integer(early_status == 'Taliban'),
         early_unconfirmed = as.integer(early_status == 'Unconfirmed'))]
X[, `:=`(eli_acled = eli_final * pre_contested_share_spatial,
         eli_lwj = eli_final * early_contested,
         eli_taliban = eli_final * early_taliban,
         eli_unconfirmed = eli_final * early_unconfirmed)]
X[, phase := fifelse(month < as.IDate('2019-01-01'), '2018',
             fifelse(month < as.IDate('2020-01-01'), '2019',
             fifelse(month < as.IDate('2020-03-01'), '2020_JanFeb',
             fifelse(month < as.IDate('2021-07-01'), 'anticipation',
             fifelse(month < as.IDate('2022-01-01'), '2021_H2',
             fifelse(month < as.IDate('2023-01-01'), '2022',
             fifelse(month < as.IDate('2024-01-01'), '2023','2024')))))))]
X[, phase := factor(phase, levels=c('2018','2019','2020_JanFeb','anticipation','2021_H2','2022','2023','2024'))]
stopifnot(!anyNA(X$phase))
X[, log1p_events := log1p(events_type3)]

# Phases are defined before estimation. 2019 is the pre-closure reference.
# Full model allows the ACLED, LWJ, ELI, and early-control categories to have
# distinct phase effects. Simpler model is explicitly labelled diagnostic.
slopes_full <- c('eli_final', 'pre_contested_share_spatial',
  'early_contested','early_taliban','early_unconfirmed',
  'eli_acled','eli_lwj','eli_taliban','eli_unconfirmed')
slopes_simple <- c('eli_final','pre_contested_share_spatial',
  'early_contested','eli_acled','eli_lwj')
make_f <- function(outcome, slopes) {
  rhs <- paste(sprintf('i(phase, %s, ref="2019")', slopes), collapse=' + ')
  as.formula(paste(outcome, '~',rhs,'| cell_id + month'))
}
records <- list(); fits <- list(); k <- 0L
for (win in c('full','exclude_anticipation')) {
  Y <- if(win == 'full') copy(X) else copy(X[phase != 'anticipation'])
  for (spec in c('full_controls','simpler_diagnostic')) {
    slopes <- if(spec == 'full_controls') slopes_full else slopes_simple
    for (family in c('PPML','OLS_log1p')) {
      key <- paste(win,spec,family,sep='__')
      cat('\n--- ',key,' ---\n',sep='')
      f <- make_f(if(family=='PPML')'events_type3' else 'log1p_events', slopes)
      warnings <- character()
      model <- tryCatch(withCallingHandlers({
        if(family=='PPML') fixest::fepois(f,data=Y,notes=FALSE,warn=TRUE)
        else fixest::feols(f,data=Y,notes=FALSE,warn=TRUE)
      }, warning=function(w){
        warnings <<- unique(c(warnings,conditionMessage(w)))
        invokeRestart('muffleWarning')
      }), error=function(e)e)
      k<-k+1L
      if(inherits(model,'error')) {
        msg <- conditionMessage(model)
        cat('MODEL ERROR: ',msg,'\n',sep='')
        records[[k]] <- data.table(window=win,specification=spec,family=family,
          status='ERROR',error=msg,input_rows=nrow(Y),estimated_n=NA_integer_,
          period=NA_character_,term=NA_character_,beta=NA_real_,se=NA_real_,p=NA_real_,
          warnings=paste(warnings,collapse=' | '))
        next
      }
      fits[[key]] <- model
      saveRDS(model,file.path(DEST,paste0(key,'.rds')))
      ct <- tryCatch(fixest::coeftable(model,vcov=~district_id_for_h2+month),error=function(e)e)
      if(inherits(ct,'error')) {
        msg<-conditionMessage(ct);cat('VCOV ERROR: ',msg,'\n',sep='')
        records[[k]]<-data.table(window=win,specification=spec,family=family,
          status='VCOV_ERROR',error=msg,input_rows=nrow(Y),estimated_n=nobs(model),
          period=NA_character_,term=NA_character_,beta=NA_real_,se=NA_real_,p=NA_real_,
          warnings=paste(warnings,collapse=' | '))
        next
      }
      T <- data.table(term=rownames(ct),beta=unname(ct[,1]),se=unname(ct[,2]),
                      p=unname(ct[,4]))
      # fixest labels i(phase,...): phase::2018:eli_lwj; retain raw labels for audit.
      T <- T[grepl('eli_lwj',term,fixed=TRUE)]
      if(!nrow(T)) cat('WARNING: eli_lwj terms absent; see saved fit and coefficient names.\n')
      T[, `:=`(window=win,specification=spec,family=family,
               status=if(nrow(T)) 'ok' else 'NO_LWJ_TERMS',error='',
               input_rows=nrow(Y),estimated_n=nobs(model),
               period=as.character(sub('.*phase::([^:]+):eli_lwj.*','\\1',term)),
               warnings=paste(warnings,collapse=' | '))]
      records[[k]] <- T
      cat('Input:',nrow(Y),'Estimated:',nobs(model),'ELI x LWJ phase terms:',nrow(T),'\n')
      if(nrow(T))print(T[,.(period,beta,se,p)])
      if(length(warnings)) cat('WARNINGS: ',paste(warnings,collapse=' | '),'\n')
    }
  }
}
R <- rbindlist(records,fill=TRUE)
fwrite(R,file.path(DEST,'01_all_attempts_including_errors.csv'))
PRE <- R[status=='ok' & period %in% c('2018','2020_JanFeb')]
fwrite(PRE,file.path(DEST,'02_pre_period_comparisons.csv'))
cat('\n=== ERRORS (IF ANY) ===\n')
print(R[status!='ok',.(window,specification,family,status,error)])
cat('\n=== PRE-PERIOD LWJ INTERACTIONS, RELATIVE TO 2019 ===\n')
print(PRE[,.(window,specification,family,period,beta,se,p,estimated_n)])
cat('\nNOTE: Late-2017 LWJ measurement is fixed; ACLED moderator includes later information.\n')
cat('NOTE: PPML separation can change the effective sample. Compare PPML and all-cell OLS descriptively, not as identical estimands.\n')
cat('NOTE: An insignificant pre-period test is not proof of parallel trends.\n')
cat('OUTPUT: ',DEST,'\n',sep='')
