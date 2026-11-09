# 04: negative-time placebo with pre-2020 outcomes. Descriptive falsification only.
# NOTE: LWJ moderator uses information through 2019 (periods 1:5), which can be
# FUTURE INFORMATION relative to a 2017/2018 pseudo-shock. Placebo cannot
# identify pre-trends causally; interpret as a warning diagnostic only.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
OUT4 <- file.path(OUT,'04_placebo')
dir.create(OUT4,recursive=TRUE,showWarnings=FALSE)
# Hold pre-closure sample fixed and impose two early placebo shifts. Never fit on post-Doha.
X0 <- copy(D[month>=as.IDate('2014-01-01') & month<as.IDate('2020-03-01')])
if(!nrow(X0))stop('Pre-treatment sample empty')
placebos<-c('2017-01-01','2019-01-01')
R<-list()
for(date_txt in placebos) {
  date <- as.IDate(date_txt)
  X<-copy(X0)
  X[,placebo_post:=as.integer(month>=date)]
  fm<-as.formula(paste('events_type3 ~',paste(c('placebo_post:eli_final',paste0('placebo_post:eli_final:',DIAG_VARS),paste0('placebo_post:',DIAG_VARS),DIAG_CONTROLS),collapse=' + '),'| cell_id + month'))
  fit<-tryCatch(fepois(fm,data=X,notes=FALSE,warn=FALSE),error=function(e){warning(as.character(date),': ',conditionMessage(e));NULL})
  if(is.null(fit))next
  z<-coeftable(fit,vcov=~cell_id+month)
  term<-'placebo_post:eli_final:lwj_contested_share'
  if(!term %in% rownames(z)){warning('Placebo coefficient absent ',date);next}
  R[[length(R)+1L]]<-data.table(fake_shock=as.character(date),beta=unname(z[term,1]),se=unname(z[term,2]),p=unname(z[term,4]),n=nobs(fit),
     limitation='Moderator defined with 2017-2019 information, potentially after placebo date')
}
if(length(R)){
 print(rbindlist(R));fwrite(rbindlist(R),file.path(OUT4,'01_pre_doha_fake_shocks.csv'))
} else warning('No placebo estimates succeeded')
cat('Placebo output: ',OUT4,'\n')
