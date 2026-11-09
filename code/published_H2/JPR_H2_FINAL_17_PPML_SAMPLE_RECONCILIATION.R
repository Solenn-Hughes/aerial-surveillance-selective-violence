# JPR H2 FINAL CHECK 17 — definitive joint PPML sample/event reconciliation
# Uses exact same expanded joint H2 as Test 13, and BOTH 2015-2024 windows.
# This is a sample audit, not another specification search.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
DEST <- file.path(OUT,'17_FINAL_PPML_reconciliation')
dir.create(DEST,recursive=TRUE,showWarnings=FALSE)
F <- events_type3 ~ post_closure:eli_final +
  post_closure:pre_contested_share_spatial + post_closure:lwj_contested_share +
  post_closure:eli_final:pre_contested_share_spatial +
  post_closure:eli_final:lwj_contested_share | cell_id + month
key <- 'post_closure:eli_final:lwj_contested_share'
res <- list(); cells_all <- list(); months_all <- list()
for(win in c('full_2015_2024','exclude_anticipation_2015_2024')){
  X <- copy(D[month>=as.IDate('2015-01-01') & month<=as.IDate('2024-12-01')])
  if(win=='exclude_anticipation_2015_2024')
    X <- X[month<as.IDate('2020-03-01') | month>=as.IDate('2021-07-01')]
  X[,audit_id:=seq_len(.N)]
  if(anyNA(X$events_type3) || any(X$events_type3<0)) stop('Outcome invalid in ',win)
  fit <- fixest::fepois(F,data=X,notes=FALSE,warn=TRUE)
  obs <- fixest::obs(fit)
  if(anyNA(obs) || length(obs)!=nobs(fit) || anyDuplicated(obs) ||
     any(obs<1L | obs>nrow(X))) stop('obs() invalid in ',win)
  X[,retained:=audit_id %in% obs]
  if(sum(X$retained)!=nobs(fit)) stop('Retained n mismatch')
  # Explicit predicate (rather than X[retained]) prevents ambiguous evaluation.
  input_total <- X[,sum(events_type3)]
  retained_total <- X[retained==TRUE,sum(events_type3)]
  removed_total <- X[retained==FALSE,sum(events_type3)]
  if(input_total != retained_total + removed_total) stop('Event conservation failed')
  bycell<-X[,.(rows=.N,kept_rows=sum(retained),events=sum(events_type3),
       kept_events=sum(events_type3[retained==TRUE]),
       removed_events=sum(events_type3[retained==FALSE]),
       eli=eli_final[1],lwj_contested=lwj_contested_share[1],
       acled_contested=pre_contested_share_spatial[1],
       district=district_id_for_h2[1]),by=cell_id]
  bycell[,group:=fifelse(kept_rows==0,'dropped',
                  fifelse(kept_rows==rows,'retained','partially_retained'))]
  bycell[,window:=win];cells_all[[win]]<-bycell
  monthdt<-X[,.(input_rows=.N,retained_rows=sum(retained),
                 input_events=sum(events_type3),
                 retained_events=sum(events_type3[retained==TRUE])),by=month]
  monthdt[,window:=win];months_all[[win]]<-monthdt
  if(sum(bycell$events)!=input_total || sum(bycell$kept_events)!=retained_total ||
     sum(monthdt$input_events)!=input_total || sum(monthdt$retained_events)!=retained_total)
    stop('Cell/month reconciliation failed: ',win)
  coeff <- fixest::coeftable(fit,vcov=~district_id_for_h2+month)
  if(!key %in% rownames(coeff)) stop('Main LWJ interaction dropped')
  res[[win]]<-data.table(window=win,input_rows=nrow(X),estimated_rows=nobs(fit),
    input_cells=uniqueN(X$cell_id),retained_cells=uniqueN(X[retained==TRUE,cell_id]),
    fully_dropped_cells=sum(bycell$group=='dropped'),
    partly_dropped_cells=sum(bycell$group=='partially_retained'),
    total_input_events=input_total,total_retained_events=retained_total,
    total_dropped_events=removed_total,
    dropped_cells_with_events=sum(bycell$group=='dropped' & bycell$events>0),
    coefficient=unname(coeff[key,1]),se=unname(coeff[key,2]),p=unname(coeff[key,4]))
  cat('\n=== ',win,' ===\n',sep='')
  print(res[[win]])
  print(bycell[,.(cells=.N,events=sum(events),avg_eli=mean(eli),
                 avg_lwj_contested=mean(lwj_contested),
                 avg_acled_contested=mean(acled_contested)),by=group])
  fwrite(X[,.(cell_id,month,district_id_for_h2,events_type3,retained)],
         file.path(DEST,paste0(win,'_OBSERVATION_RETENTION.csv')))
}
R<-rbindlist(res); C<-rbindlist(cells_all); M<-rbindlist(months_all)
fwrite(R,file.path(DEST,'01_VERIFIED_SAMPLE_TOTALS.csv'))
fwrite(C,file.path(DEST,'02_CELL_RETENTION.csv'))
fwrite(M,file.path(DEST,'03_MONTH_RETENTION.csv'))
cat('\n=== AUTHORITATIVE TOTALS (do not reuse inconsistent Test 14 console line) ===\n')
print(R)
cat('Zero-only cell exclusion is compatible with PPML separation. It does NOT prove representativeness.\n')
cat('Output: ',DEST,'\n',sep='')
