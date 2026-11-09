# 02: sensitivity to the historical-name crosswalk, without constructing alternative invented matches.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
OUT2 <- file.path(OUT,'02_crosswalk_sensitivity')
dir.create(OUT2,recursive=TRUE,showWarnings=FALSE)
stopifnot(all(c('adm2_id','match_basis','documented_reason') %in% names(CW)))
cat('\nMatch methods used by cells:\n')
J <- unique(D[,.(cell_id,district_id_for_h2)])
J[, match_basis:=CW$match_basis[match(district_id_for_h2,CW$adm2_id)]]
print(J[,.(cells=.N,districts=uniqueN(district_id_for_h2)),by=match_basis][order(-cells)])
fwrite(J[,.(cells=.N,districts=uniqueN(district_id_for_h2)),by=match_basis],file.path(OUT2,'01_crosswalk_method_support.csv'))
# Main result retained on all cells. Then use only exact normalized-name+province matches.
# Different samples answer different questions; this is a sensitivity check, NOT a superior estimator.
exact_basis <- c('exact_name_province','normalized_province_name','exact_after_province_normalization')
cat('Observed crosswalk match_basis labels: ',paste(sort(unique(CW$match_basis)),collapse=', '),'\n')
exact_ids <- CW[match_basis %in% exact_basis,adm2_id]
runset <- list(all_cells=unique(J$cell_id))
if(length(exact_ids)>0L) runset$exact_name_cells <- J[district_id_for_h2 %in% exact_ids,cell_id]
if(length(runset)<2L) warning('No exact-match subset found. Check basis labels in 01_crosswalk_method_support.csv')
R<-list()
for(w in names(DIAG_WINDOWS)) for(s in names(runset)) {
  dx <- DIAG_WINDOWS[[w]][cell_id %in% runset[[s]]]
  if(uniqueN(dx$cell_id)<30L || sum(dx$events_type3,na.rm=TRUE)==0L) {cat('Skipping small sample: ',w,' ',s,'\n');next}
  rr<-tryCatch(diag_fit(dx),error=function(e){warning(w,' ',s,': ',conditionMessage(e));NULL})
  if(is.null(rr))next
  R[[length(R)+1L]]<-data.table(window=w,sample=s,cells=uniqueN(dx$cell_id),districts=uniqueN(dx$district_id_for_h2),beta=rr$beta,se=rr$se,p=rr$p,n=rr$n)
  cat(w,s,'cells=',uniqueN(dx$cell_id),'beta=',round(rr$beta,3),'SE=',round(rr$se,3),'p=',round(rr$p,4),'\n')
}
if(length(R))fwrite(rbindlist(R),file.path(OUT2,'02_crosswalk_results.csv'))
fwrite(CW[,.(adm2_id,shapefile_province,shapefile_district,lwj_object_id,lwj_province,lwj_district,match_basis,documented_reason)],file.path(OUT2,'03_crosswalk_audit.csv'))
cat('Crosswalk sensitivity outputs: ',OUT2,'\n')
