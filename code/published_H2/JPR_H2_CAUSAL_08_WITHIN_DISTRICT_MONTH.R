# JPR H2 / 08: Within-district-by-month test (most demanding consolidation adjustment).
# This is NOT automatically causal: cell-specific ELI-correlated shocks remain possible.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT, '08_within_district_month'); dir.create(DEST, recursive=TRUE, showWarnings=FALSE)
X <- copy(D[month >= as.IDate('2015-01-01') & month <= as.IDate('2024-12-01')])
# Only districts with at least two cells and nonzero within-district ELI dispersion can identify ELI slopes.
Z <- unique(X[,.(district_id_for_h2,cell_id,eli_final,lwj_contested_share)])
SUM <- Z[,.(n_cells=.N,eli_range=max(eli_final)-min(eli_final),contested_share=first(lwj_contested_share)),by=district_id_for_h2]
SUM[,informative := n_cells>=2L & eli_range>1e-6]
print(SUM[,.(districts=.N,with_multi_cells=sum(n_cells>=2),with_eli_variation=sum(informative),
 cells_in_informative=sum(n_cells[informative]))])
fwrite(SUM,file.path(DEST,'01_district_identification_support.csv'))
X[,district_month:=interaction(district_id_for_h2,as.character(month),drop=TRUE)]
# The district x month FE absorb all district-wide changes and the Post x status main effects.
# Include Post x ELI and differential ELI by status; omit government status as the reference.
trip <- paste(DIAG_TRIPLES,collapse=' + ')
noctrl <- as.formula(paste('events_type3 ~ post_closure:eli_final +',trip,'| cell_id + district_month'))
withctrl <- as.formula(paste('events_type3 ~ post_closure:eli_final +',trip,'+ prec_gpcp + temp + nlights_calib_mean | cell_id + district_month'))
res <- list()
for (w in c('full','drop_anticipation')) {
 Y <- if(w=='full') X else X[month < as.IDate('2020-03-01') | month >= as.IDate('2021-07-01')]
 for (spec in c('no_controls','with_controls')) {
  fit <- tryCatch(fixest::fepois(if(spec=='no_controls') noctrl else withctrl,
                      data=Y,notes=FALSE,warn=FALSE),error=function(e)e)
  if(inherits(fit,'error')) {
   msg<-conditionMessage(fit);cat('FAILED',w,spec,msg,'\n')
   res[[length(res)+1L]]<-data.table(window=w,specification=spec,status=msg);next
  }
  ct<-tryCatch(fixest::coeftable(fit,vcov=~district_id_for_h2+month),error=function(e)e)
  if(inherits(ct,'error')) {res[[length(res)+1L]]<-data.table(window=w,specification=spec,status=conditionMessage(ct));next}
  for(term in DIAG_TRIPLES) {
   v<-if(term %in% rownames(ct)) unname(ct[term,]) else rep(NA_real_,4)
   res[[length(res)+1L]]<-data.table(window=w,specification=spec,term=term,
       beta=v[1],se=v[2],p=v[4],n=as.integer(nobs(fit)),
       n_input=nrow(Y),term_estimable=term %in% rownames(ct),
       dropped=paste(fit$collin.var,collapse=';'),status='ok')
  }
  saveRDS(fit,file.path(DEST,paste0(w,'_',spec,'.rds')))
 }
}
R<-rbindlist(res,fill=TRUE);print(R)
fwrite(R,file.path(DEST,'02_within_district_month_H2.csv'))
cat('\nDONE:',DEST,'\n')
cat('Interpret only with district support, estimated N, and coefficient estimability.\n')
