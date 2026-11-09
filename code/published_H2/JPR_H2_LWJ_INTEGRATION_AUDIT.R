# =============================================================
# JPR H2 — LWJ / Dudek district-control integration and audit
# Verified workbook: Civilian Victimization Data (Final).xlsx
# Main sheet: 407 source districts x 9 six-month periods.
# Uses periods 1-5 (late 2017 through late 2019) as pre-Doha.
# Does NOT treat source-specific LWJ object IDs as Afghan ADM2 IDs.
# Stage 1: audited join; Stage 2: sensitivity only if join passes.
# Never overwrites the FINAL30 panel or original estimates.
# =============================================================
suppressPackageStartupMessages({library(data.table); library(sf); library(fixest)})
ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
PANEL <- file.path(ROOT,'outputs','final_airfield_screen_20260809_FINAL30','11_final_analysis_panel.rds')
SHAPE <- file.path(ROOT,'data','raw','afg_admin','afg_admin2.shp')
# Download the accompanying LWJ_Dudek_district_periods.csv to Downloads.
LWJ <- 'C:/Users/USERNAME/Downloads/LWJ_Dudek_district_periods.csv'
OUT <- file.path(ROOT,'outputs','JPR_H2_LWJ_AUDIT')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
for(x in c(PANEL,SHAPE,LWJ)) if(!file.exists(x)) stop('Missing file: ',x)

pick <- function(x,choices,required=TRUE){
 ix<-match(tolower(choices),tolower(x)); ix<-ix[!is.na(ix)]
 if(!length(ix)){if(required)stop('Missing column. Available: ',paste(x,collapse=', '));return(NA_character_)}
 x[ix[1]]
}
norm <- function(x){
 x<-tolower(iconv(as.character(x),to='ASCII//TRANSLIT',sub=''))
 x<-gsub('&',' and ',x,fixed=TRUE)
 x<-gsub('[^a-z0-9]','',x)
 x
}
L<-fread(LWJ,encoding='UTF-8')
stopifnot(nrow(L)==407L*9L,uniqueN(L$lwj_object_id)==407L)
if(anyDuplicated(L,by=c('lwj_object_id','period')))stop('Duplicate source district-period')
if(!setequal(unique(L$period),1:9))stop('Unexpected period numbers')
if(!all(L$control %in% c('GOA','Taliban','Contested','Unconfirmed')))stop('Unexpected control status')
L[,lwj_object_id:=as.character(lwj_object_id)]
# IMPORTANT: only periods 1-5 are unambiguously before Doha.
LP<-L[period %in% 1:5,.(
 lwj_contested_share=mean(control=='Contested'),
 lwj_goa_share=mean(control=='GOA'),
 lwj_taliban_share=mean(control=='Taliban'),
 lwj_unconfirmed_share=mean(control=='Unconfirmed'),
 lwj_n_periods=.N),by=.(lwj_object_id,province,district)]
if(any(LP$lwj_n_periods!=5L))stop('Incomplete historical source')
LP[,`:=`(province_key=norm(province),district_key=norm(district))]
fwrite(LP,file.path(OUT,'01_lwj_pre_doha_measurements.csv'))

S<-st_read(SHAPE,quiet=TRUE)
name_col<-pick(names(S),c('adm2_en','adm2_name','district','district_name','name_2'))
prov_col<-pick(names(S),c('adm1_en','adm1_name','province','province_name','name_1'))
id_col<-pick(names(S),c('adm2_pcode','adm2_code','dist_id','district_id','adm2_code_en'),FALSE)
if(is.na(id_col))stop('Shapefile lacks stable ADM2 ID; cannot link confidently to panel')
SM<-unique(data.table(adm2_id=as.character(S[[id_col]]),
 adm2_name=as.character(S[[name_col]]),province_name=as.character(S[[prov_col]])))
SM[,`:=`(province_key=norm(province_name),district_key=norm(adm2_name))]
if(anyDuplicated(SM$adm2_id))stop('Shapefile ADM2 IDs not unique')

# Exact normalized province + district name matching only.
# Do not guess near matches. Save unresolved source districts for manual review.
LX<-LP[,.(lwj_object_id,province,district,province_key,district_key)]
SF<-SM[,.(adm2_id,adm2_name,province_name,province_key,district_key)]
SX<-merge(SF,LX,by=c('province_key','district_key'),all.x=TRUE,allow.cartesian=TRUE)
SX[,match_type:=fifelse(is.na(lwj_object_id),'unmatched','exact_province_district')]
fwrite(SX,file.path(OUT,'02_shapefile_to_lwj_review.csv'))
if(SX[!is.na(lwj_object_id),anyDuplicated(adm2_id)] > 0 ||
   SX[!is.na(lwj_object_id),anyDuplicated(lwj_object_id)] > 0)
  stop('Ambiguous exact district matches; review 02_shapefile_to_lwj_review.csv')

D<-as.data.table(readRDS(PANEL))
D[,`:=`(cell_id=as.character(cell_id),month=as.IDate(month),
        district_id_for_h2=as.character(district_id_for_h2))]
needed<-c('cell_id','month','district_id_for_h2','post_closure','eli_final','events_type3',
          'prec_gpcp','temp','nlights_calib_mean')
if(length(miss<-setdiff(needed,names(D))))stop('Missing FINAL30 fields: ',paste(miss,collapse=', '))
cell<-unique(D[,.(cell_id,district_id_for_h2)])
if(anyDuplicated(cell$cell_id))stop('Ambiguous panel district assignment')
cat('\nSource rows:',nrow(L),'districts:',nrow(LP),'\n')
cat('Shapefile districts:',nrow(SM),'\n')
cat('Exact province/district matches:',SX[!is.na(lwj_object_id),.N],'\n')
cat('Sample panel ADM2 IDs:',paste(head(unique(cell$district_id_for_h2),6),collapse=', '),'\n')
cat('Sample shapefile ADM2 IDs:',paste(head(SM$adm2_id,6),collapse=', '),'\n')

# The original H2 script matches the panel district_id_for_h2 to
# a shapefile ADM2 ID. Verify that this remains true here.
match_id_fraction<-mean(cell$district_id_for_h2 %in% SM$adm2_id)
cat('Panel district IDs found in shapefile:',round(match_id_fraction*100,1),'% of cells\n')
if(match_id_fraction < .99){
 fwrite(cell,file.path(OUT,'03_panel_ids_need_crosswalk.csv'))
 stop('Panel ADM2 IDs do not align with shapefile IDs. Wrote crosswalk review files; STOPPING rather than creating false matches.')
}
valid_map<-SX[!is.na(lwj_object_id),.(district_id_for_h2=adm2_id,lwj_object_id)]
cell<-merge(cell,valid_map,by='district_id_for_h2',all.x=TRUE,sort=FALSE)
match_cells<-sum(!is.na(cell$lwj_object_id))
cat('FINAL30 cells matched by exact name+province:',match_cells,'/',nrow(cell),'\n')
fwrite(cell,file.path(OUT,'03_cell_match_audit.csv'))
fwrite(cell[is.na(lwj_object_id)],file.path(OUT,'03_unmatched_cells.csv'))
# SOURCE names can differ by spelling. Do not proceed until matched cells audited.
if(match_cells < nrow(cell)){
 cat('\nSTOP: Some panel cells are unmatched. Review 02_shapefile_to_lwj_review.csv and 03_unmatched_cells.csv.\n')
 cat('Once the spelling crosswalk is documented, re-run with approved matches. No models estimated.\n')
 stop('Matching incomplete: no H2 models run. Review the saved crosswalk files.')
}

D<-merge(D,cell[,.(cell_id,lwj_object_id)],by='cell_id',all.x=TRUE,sort=FALSE)
D<-merge(D,LP[,.(lwj_object_id,lwj_contested_share,lwj_goa_share,lwj_taliban_share,
                  lwj_unconfirmed_share)],by='lwj_object_id',all.x=TRUE,sort=FALSE)
if(anyNA(D$lwj_contested_share)||nrow(D)!=104112L)stop('Panel change or missing LWJ moderators')
cat('Matched FULL FINAL30 panel, with zero unmatched cells. Starting comparison.\n')

mods<-c('lwj_contested_share','lwj_taliban_share')
windows<-list(full=D,drop_anticipation=D[month<'2020-03-01'|month>='2021-07-01'])
rows<-list()
for(w in names(windows))for(m in mods)for(spec in c('restricted','expanded')){
 X<-windows[[w]]
 terms<-c('post_closure:eli_final',paste0('post_closure:eli_final:',m))
 if(spec=='expanded')terms<-c(terms,paste0('post_closure:',m))
 form<-as.formula(paste('events_type3 ~',paste(c(terms,'prec_gpcp','temp','nlights_calib_mean'),collapse=' + '),'| cell_id + month'))
 fit<-fepois(form,data=X,notes=FALSE,warn=FALSE)
 tri<-paste0('post_closure:eli_final:',m)
 z<-coeftable(fit,vcov=~cell_id+month)
 if(!(tri %in% rownames(z)))stop('Missing three-way coefficient ',tri)
 rows[[length(rows)+1L]]<-data.table(window=w,moderator=m,specification=spec,
   beta=unname(z[tri,1]),se=unname(z[tri,2]),p=unname(z[tri,4]),n=nobs(fit))
 cat(w,m,spec,'beta=',round(z[tri,1],4),'se=',round(z[tri,2],4),'p=',round(z[tri,4],4),'N=',nobs(fit),'\n')
}
fwrite(rbindlist(rows),file.path(OUT,'04_lwj_h2_models.csv'))
cat('DONE — matching verified; see',OUT,'\n')
