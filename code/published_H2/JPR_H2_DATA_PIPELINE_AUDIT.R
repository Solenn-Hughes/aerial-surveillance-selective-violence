# ============================================================
# JPR H2 — DATA AND MEASUREMENT AUDIT (read-only)
# Purpose: check for data / timing / mapping mistakes BEFORE any
# specification search. Does not replace published estimates.
# Required packages: data.table, fixest
# ============================================================
suppressPackageStartupMessages({library(data.table);library(fixest)})
ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
PANEL <- file.path(ROOT, 'outputs', 'final_airfield_screen_20260809_FINAL30', '11_final_analysis_panel.rds')
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_DATA_PIPELINE_AUDIT')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
if(!file.exists(PANEL)) stop('FINAL30 panel not found: ', PANEL)
D <- as.data.table(readRDS(PANEL))
D[,month:=as.IDate(month)]
mods <- c('pre_contested_share_spatial','pre_taliban_activity_share_spatial','pre_contested_intensity_spatial')
need <- c('cell_id','month','post_closure','events_type3','eli_final',mods,
          'prec_gpcp','temp','nlights_calib_mean')
if(length(x <- setdiff(need,names(D)))) stop('Missing: ',paste(x,collapse=', '))
D[,cell_id:=as.character(cell_id)]
report <- list(); flag <- function(check,value,note='') {
  report[[length(report)+1L]] <<- data.table(check=check,value=as.character(value),note=note)
  cat(check,': ',value,if(nzchar(note))paste0(' | ',note) else '', '\n',sep='')
}
cat('\n===== 1. BASIC PANEL INTEGRITY =====\n')
flag('Rows',nrow(D));flag('Distinct cells',uniqueN(D$cell_id))
flag('Distinct months',uniqueN(D$month))
flag('Months span',paste(min(D$month),max(D$month)))
flag('Duplicate cell-months',anyDuplicated(D,by=c('cell_id','month')))
flag('Missing outcome',sum(is.na(D$events_type3)))
flag('Negative/noninteger outcome',sum(is.finite(D$events_type3) & (D$events_type3<0 | D$events_type3!=floor(D$events_type3))))
flag('Post code inconsistent with July 2021',sum(!is.na(D$post_closure) & D$post_closure != as.integer(D$month>=as.IDate('2021-07-01'))))
flag('Full panel outcome events',sum(D$events_type3,na.rm=TRUE))
flag('Rows in anticipation period',nrow(D[month>=as.IDate('2020-03-01') & month<as.IDate('2021-07-01')]))
for(v in c('eli_final',mods)){
  n_vary <- D[,.(vary=uniqueN(get(v),na.rm=TRUE)>1L),by=cell_id][vary==TRUE,.N]
  flag(paste('Time variation in',v,'(cells)'),n_vary)
  flag(paste('Missing',v),sum(is.na(D[[v]])))
}
for(v in mods[1:2]) flag(paste('Outside 0-1:',v),sum(!is.na(D[[v]]) & (D[[v]]<0 | D[[v]]>1)))

cat('\n===== 2. DISTINCT CELLS, MODERATORS, AND OUTCOMES =====\n')
cell <- unique(D[,c('cell_id','eli_final',mods,intersect(c('district_id_for_h2','pre_contested_total_spatial','pre_months_observed_spatial'),names(D))),with=FALSE])
if(anyDuplicated(cell$cell_id)) flag('Multiple moderator/ELI records per cell',nrow(cell)-uniqueN(cell$cell_id))
else flag('Multiple moderator/ELI records per cell',0)
for(v in c('eli_final',mods)){
  s<-cell[[v]]
  flag(paste('Cell summary',v),paste('unique',uniqueN(s,na.rm=TRUE),'zero',sum(s==0,na.rm=TRUE),
      'p10',round(quantile(s,.10,na.rm=TRUE),4),'median',round(median(s,na.rm=TRUE),4),
      'p90',round(quantile(s,.90,na.rm=TRUE),4)))
}
if('district_id_for_h2' %in% names(D)){
  distcheck <- unique(D[,.(cell_id,district_id_for_h2)])[,.(n_districts=uniqueN(district_id_for_h2,na.rm=TRUE)),by=cell_id]
  flag('Cells assigned multiple district IDs',nrow(distcheck[n_districts>1]))
  flag('Cells lacking district IDs',nrow(distcheck[n_districts==0]))
  flag('Distinct assigned districts',uniqueN(D$district_id_for_h2,na.rm=TRUE))
  district_mods<-unique(D[,c('district_id_for_h2',mods),with=FALSE])
  flag('Districts with multiple moderator records',nrow(district_mods)-uniqueN(district_mods$district_id_for_h2))
}
period <- D[,.(rows=.N,events=sum(events_type3,na.rm=TRUE),
              positive_cell_months=sum(events_type3>0,na.rm=TRUE),
              cells_with_any_events=uniqueN(cell_id[events_type3>0])),
           by=.(period=fifelse(month<as.IDate('2020-03-01'),'pre_Doha',
                 fifelse(month<as.IDate('2021-07-01'),'anticipation','post_closure')))]
fwrite(period,file.path(OUT,'02_period_outcomes.csv'));print(period)
fwrite(cell,file.path(OUT,'02_cell_measurements.csv'))

cat('\n===== 3. ORIGINAL H2 REPRODUCTION =====\n')
bench <- data.table(moderator=mods,expected=c(-16.1023802013192,-17.6667522797993,-3.12067660696125),
                    expected_se=c(4.27319350902624,4.60318002832200,0.880878277343997))
fits <- list()
for(m in mods){
  tri<-paste('post_closure','eli_final',m,sep=':')
  f<-as.formula(paste('events_type3 ~',paste(c('post_closure:eli_final',tri,
       'prec_gpcp','temp','nlights_calib_mean'),collapse=' + '),'| cell_id + month'))
  fit <- fepois(f,data=D,notes=FALSE,warn=FALSE)
  zz <- coeftable(fit,vcov=~cell_id+month)
  if(!(tri %in% rownames(zz))) stop('Model dropped term: ',tri)
  b<-as.numeric(zz[tri,1]);se<-as.numeric(zz[tri,2]);n<-nobs(fit)
  tar<-bench[moderator==m]
  pass<-abs(b-tar$expected[[1]])<1e-3 && abs(se-tar$expected_se[[1]])<1e-3 && n==54312
  cat(m,' | beta=',round(b,5),' SE=',round(se,5),' N=',n,' | reproduced=',pass,'\n',sep='')
  fits[[length(fits)+1L]]<-data.table(moderator=m,beta=b,se=se,n=n,reproduced=pass)
  if(!pass) stop('Stop: original published model did not reproduce. No further inference.')
}
fwrite(rbindlist(fits),file.path(OUT,'03_original_replication.csv'))

cat('\n===== 4. AUDIT ORIGINAL MODERATOR CONSTRUCTION =====\n')
# IMPORTANT SOURCE-CODE FACT: H2 IMPROVEMENT.R uses month < 2021-07-01
# for construction of pre-shock shares and uses log1p(mean(battles_GT))
# for the intensity variable, not log1p(total events).
flag('Original construction window','All observed district-months before July 2021',
     'Includes March 2020–June 2021; dropping outcome rows does not rebuild moderators')
flag('Original intensity transformation','log1p(mean(monthly contested battle count))',
     'Verify manuscript describes this rather than log number of battles')
if('pre_months_observed_spatial' %in% names(cell)){
  flag('Number of distinct pre-month denominators',uniqueN(cell$pre_months_observed_spatial,na.rm=TRUE))
  print(cell[,.(cells=.N),by=pre_months_observed_spatial][order(pre_months_observed_spatial)])
}
if(all(c('pre_contested_total_spatial','pre_months_observed_spatial') %in% names(cell))){
  cell[,intensity_from_source:=log1p(pre_contested_total_spatial/pre_months_observed_spatial)]
  cell[,intensity_from_total:=log1p(pre_contested_total_spatial)]
  flag('Intensity mismatches vs log1p(total/months)',sum(abs(cell$intensity_from_source-cell$pre_contested_intensity_spatial)>1e-6,na.rm=TRUE))
  flag('Intensity differs from log1p(total) in cells',sum(abs(cell$intensity_from_total-cell$pre_contested_intensity_spatial)>1e-6,na.rm=TRUE))
  fwrite(cell,file.path(OUT,'04_intensity_consistency.csv'))
}

cat('\n===== 5. OPTIONAL DISTRICT-MONTH RECONSTRUCTION =====\n')
DM_PATH <- file.path(ROOT,'outputs','paper1_geo_afghanistan_executable_acled',
                     'diagnostics_h2_spatial_adm2','district_month_contestedness_spatial_adm2.csv')
if(!file.exists(DM_PATH)){
  flag('District-month source file','NOT FOUND',DM_PATH)
  cat('Cannot independently audit the ACLED district-month aggregation or rebuild pre-Doha moderators.\n')
} else if(!('district_id_for_h2' %in% names(D))){
  flag('District-month reconstruction','SKIPPED','Panel has no district_id_for_h2')
} else {
  dm<-fread(DM_PATH)
  dm[,month:=as.IDate(month)]
  for(v in c('district_id','contested_any','taliban_active_any','battles_GT'))
    if(!v %in% names(dm)) stop('District-month file missing ',v)
  flag('Duplicate district-months',anyDuplicated(dm,by=c('district_id','month')))
  flag('Source district-month date range',paste(min(dm$month),max(dm$month)))
  calc<-function(end){dm[month<end,.(share=mean(contested_any),tal_share=mean(taliban_active_any),
                 intensity=log1p(mean(battles_GT)),n_months=.N),by=district_id]}
  src<-calc(as.IDate('2021-07-01'))
  predoha<-calc(as.IDate('2020-03-01'))
  setnames(src,c('share','tal_share','intensity'),mods)
  setnames(predoha,c('share','tal_share','intensity'),paste0('preDoha_',mods))
  x<-merge(cell,src,by.x='district_id_for_h2',by.y='district_id',all.x=TRUE,suffixes=c('_panel','_rebuilt'))
  for(v in mods){
    a<-paste0(v,'_panel');b<-paste0(v,'_rebuilt')
    flag(paste('Moderator reconstruction mismatches:',v),sum(!is.na(x[[a]]) & !is.na(x[[b]]) & abs(x[[a]]-x[[b]])>1e-7))
  }
  fwrite(merge(cell,predoha,by.x='district_id_for_h2',by.y='district_id',all.x=TRUE),
         file.path(OUT,'05_preDoha_moderator_candidates.csv'))
  fwrite(x,file.path(OUT,'05_moderator_reconstruction.csv'))
  flag('Pre-Doha moderator candidates','SAVED','This is a measurement sensitivity, not a discovered error')
}

cat('\n===== 6. PROVENANCE AND SPATIAL ASSIGNMENT =====\n')
MAPAUDIT <- file.path(ROOT,'outputs','paper1_geo_afghanistan_executable_acled',
                      'diagnostics_h2_spatial_adm2','spatial_acled_adm2_match_summary.csv')
if(file.exists(MAPAUDIT)){
  cat('Spatial join summary from original pipeline:\n'); print(fread(MAPAUDIT))
  file.copy(MAPAUDIT,file.path(OUT,'06_original_spatial_join_summary.csv'),overwrite=TRUE)
} else cat('Spatial join summary not found; cannot verify within-polygon vs nearest-within-25km matching.\n')

fwrite(rbindlist(report),file.path(OUT,'01_audit_summary.csv'))
cat('\nDONE. Read-only audit saved to: ',OUT,'\n',sep='')
cat('Start with 01_audit_summary.csv and 03_original_replication.csv.\n')
cat('No statistical conclusion about H2 should be drawn from measurement discrepancies until the source data are verified.\n')
