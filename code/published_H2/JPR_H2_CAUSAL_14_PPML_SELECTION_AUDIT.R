# JPR H2 / 14 — PPML selection, support, and all-cell outcome sensitivity.
# Does not modify the original panel. Depends on 00_COMMON.R and LWJ CSVs in Downloads.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT, '14_ppml_selection_audit')
dir.create(DEST, recursive=TRUE, showWarnings=FALSE)

F <- events_type3 ~ post_closure:eli_final +
  post_closure:pre_contested_share_spatial + post_closure:lwj_contested_share +
  post_closure:eli_final:pre_contested_share_spatial +
  post_closure:eli_final:lwj_contested_share | cell_id + month
KEY <- 'post_closure:eli_final:lwj_contested_share'

fit_safe <- function(form, x, poisson=TRUE) {
  tryCatch(if (poisson) fixest::fepois(form, data=x, notes=FALSE, warn=TRUE)
           else fixest::feols(form, data=x, notes=FALSE, warn=TRUE),
           error=function(e)e)
}
get_result <- function(m, label, input_n, vc=~district_id_for_h2+month) {
  if(inherits(m,'error')) return(data.table(model=label, status=conditionMessage(m), input_n=input_n))
  zz <- tryCatch(fixest::coeftable(m,vcov=vc),error=function(e)e)
  if(inherits(zz,'error')) return(data.table(model=label,status=conditionMessage(zz),input_n=input_n,n=nobs(m)))
  if(!(KEY %in% rownames(zz))) return(data.table(model=label,status='KEY_DROPPED',input_n=input_n,n=nobs(m)))
  data.table(model=label,status='ok',input_n=input_n,n=nobs(m),beta=zz[KEY,1],se=zz[KEY,2],p=zz[KEY,4])
}

all_models <- list(); all_groups <- list(); all_cells <- list(); all_months <- list()
for(w in c('full','exclude_anticipation')) {
  X <- copy(D[month>=as.IDate('2015-01-01') & month<=as.IDate('2024-12-01')])
  if(w=='exclude_anticipation') X <- X[month<as.IDate('2020-03-01') | month>=as.IDate('2021-07-01')]
  X[, rid14 := .I]
  m <- fit_safe(F,X)
  if(inherits(m,'error')) stop('PPML joint fit failed in ',w,': ',conditionMessage(m))
  kept <- fixest::obs(m)
  if(length(kept)!=nobs(m) || anyNA(kept) || any(kept<1L | kept>nrow(X))) stop('Unexpected obs() indices')
  X[, retained := rid14 %in% kept]
  if(sum(X$retained)!=nobs(m)) stop('Retained observation audit failed')
  # Verify if entire cells, rather than individual observations, are excluded.
  C <- X[,.(cell_rows=.N,retained_rows=sum(retained),total_events=sum(events_type3),
           post_events=sum(events_type3[post_closure==1L]),
           eli=eli_final[1],lwj_contested=lwj_contested_share[1],
           acled_contested=pre_contested_share_spatial[1],
           district=district_id_for_h2[1],
           status_match=CW$match_basis[match(district_id_for_h2[1],CW$adm2_id)]),by=cell_id]
  C[,`:=`(group=fifelse(retained_rows==cell_rows,'all_retained',
                       fifelse(retained_rows==0L,'all_dropped','partly_retained')),
          zero_events=total_events==0)]
  C[,window:=w]; all_cells[[w]] <- C
  G <- C[,.(cells=.N, districts=uniqueN(district),zero_outcome_cells=sum(zero_events),
            events=sum(total_events),avg_eli=mean(eli),median_eli=median(eli),
            avg_lwj_contested=mean(lwj_contested),avg_acled_contested=mean(acled_contested),
            non_exact_cells=sum(status_match!='exact_name_province')),by=group]
  G[,window:=w]; all_groups[[w]] <- G
  M <- X[,.(input_rows=.N,retained_rows=sum(retained),total_events=sum(events_type3),
           retained_events=sum(events_type3[retained])),by=month]
  M[,window:=w]; all_months[[w]]<-M
  all_models[[length(all_models)+1L]] <- get_result(m,paste0(w,'_PPML_retained'),nrow(X))
  # A linear fixed-effects sensitivity keeps cells with no violence.
  # This changes the estimand and is NOT a causal identification strategy.
  X[,log1p_events:=log1p(events_type3)]
  f_ols <- log1p_events ~ post_closure:eli_final +
    post_closure:pre_contested_share_spatial + post_closure:lwj_contested_share +
    post_closure:eli_final:pre_contested_share_spatial +
    post_closure:eli_final:lwj_contested_share | cell_id + month
  ols <- fit_safe(f_ols,X,poisson=FALSE)
  all_models[[length(all_models)+1L]] <- get_result(ols,paste0(w,'_OLS_log1p_all_cells'),nrow(X))
  if(!inherits(ols,'error')) saveRDS(ols,file.path(DEST,paste0(w,'_OLS_all_cells.rds')))
  saveRDS(m,file.path(DEST,paste0(w,'_PPML.rds')))
  fwrite(X[,.(rid14,cell_id,month,district_id_for_h2,retained,events_type3)],
         file.path(DEST,paste0(w,'_observation_retention.csv')))
  cat('\n',w,': input ',nrow(X),' retained ',sum(X$retained),' events retained ',sum(X[retained]$events_type3),
      ' of ',sum(X$events_type3),'\n',sep='')
  print(G)
}
fwrite(rbindlist(all_cells,fill=TRUE),file.path(DEST,'01_cell_level_retention.csv'))
fwrite(rbindlist(all_groups,fill=TRUE),file.path(DEST,'02_retained_vs_dropped_groups.csv'))
fwrite(rbindlist(all_months,fill=TRUE),file.path(DEST,'03_month_level_retention.csv'))
R <- rbindlist(all_models,fill=TRUE)
fwrite(R,file.path(DEST,'04_PPML_vs_all_cell_OLS.csv'))
cat('\n=== SAME-H2-FORMULA OUTCOME SENSITIVITY ===\n');print(R)
cat('\nIMPORTANT: All-cell log(1+y) OLS and selected-sample PPML estimate different parameters.\n')
cat('Zero-outcome cell exclusion is compatible with PPML FE separation, not proof that selection is ignorable.\n')
cat('DONE: ',DEST,'\n')
