# JPR H2 / 07: Explain sample loss and dropped controls in Test 06.
# Descriptive audit; DOES NOT establish causal identification.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT, '07_sample_audit'); dir.create(DEST, recursive=TRUE, showWarnings=FALSE)
X <- copy(D[month >= as.IDate('2015-01-01') & month <= as.IDate('2024-12-01')])
COV <- c('prec_gpcp','temp','nlights_calib_mean')
X[, pre_period := ifelse(month < as.IDate('2020-03-01'), 'pre_doha',
                        ifelse(month < as.IDate('2021-07-01'), 'anticipation', 'post_closure'))]
MISS <- X[, lapply(.SD, function(z) sum(is.na(z))), .SDcols=c(COV,'events_type3'), by=pre_period]
print(MISS); fwrite(MISS, file.path(DEST,'01_missing_by_period.csv'))
STATS <- rbindlist(lapply(COV,function(v) {
 z <- X[[v]]; data.table(variable=v,non_missing=sum(!is.na(z)),
 unique_nonmissing=uniqueN(z[!is.na(z)]), minimum=if(all(is.na(z))) NA_real_ else min(z,na.rm=TRUE),
 maximum=if(all(is.na(z))) NA_real_ else max(z,na.rm=TRUE))
}))
print(STATS); fwrite(STATS, file.path(DEST,'02_control_variation.csv'))
F0 <- as.formula(paste('events_type3 ~',paste(c('post_closure:eli_final',DIAG_TRIPLES,DIAG_LOWER),collapse=' + '),'| cell_id + month'))
F1 <- DIAG_FORM
res <- list(); omitted <- list()
for (w in c('full','drop_anticipation')) {
 dat <- if(w=='full') X else X[pre_period!='anticipation']
 for (s in c('no_controls','original_controls')) {
  fm <- if(s=='no_controls') F0 else F1
  fit <- tryCatch(fixest::fepois(fm,data=dat,notes=FALSE,warn=FALSE),error=function(e)e)
  if(inherits(fit,'error')) {
   res[[length(res)+1L]] <- data.table(window=w,model=s,input_rows=nrow(dat),status=conditionMessage(fit)); next
  }
  ct <- tryCatch(fixest::coeftable(fit,vcov=~cell_id+month),error=function(e)e)
  key <- if(!inherits(ct,'error') && DIAG_TERM %in% rownames(ct)) unname(ct[DIAG_TERM,]) else rep(NA_real_,4)
  res[[length(res)+1L]] <- data.table(window=w,model=s,input_rows=nrow(dat),input_cells=uniqueN(dat$cell_id),
    estimated_rows=as.integer(nobs(fit)),excluded_rows=nrow(dat)-as.integer(nobs(fit)),
    dropped_controls=paste(intersect(COV,fit$collin.var),collapse=';'),
    coefficient=key[1],se=key[2],p=key[4],status='ok')
  saveRDS(fit,file.path(DEST,paste0(w,'_',s,'.rds')))
  # fixest stores successive exclusions in obs_selection; count without assuming its encoding.
  ss <- fit$obs_selection
  omitted[[length(omitted)+1L]] <- data.table(window=w,model=s,
    selection_stages=if(is.null(ss)) 'none' else paste(names(ss),collapse=';'),
    selection_lengths=if(is.null(ss)) '' else paste(vapply(ss,length,integer(1)),collapse=';'))
 }
}
R <- rbindlist(res,fill=TRUE); print(R)
fwrite(R,file.path(DEST,'03_model_sample_comparison.csv'))
fwrite(rbindlist(omitted,fill=TRUE),file.path(DEST,'04_fixest_obs_selection.csv'))
cat('\nAUDIT DONE:',DEST,'\n')
cat('Do not claim control collinearity explains sample loss unless these diagnostics establish it.\n')
