# ============================================================================
# JPR H2: COMMON SUPPORT, PRE-TRENDS/PLACEBOS, AND MARGINAL EFFECTS
# Verified FINAL30 source; no search for significant specifications.
# R packages: data.table, fixest (install.packages(c('data.table','fixest')))
# Run: source('C:/Users/USERNAME/Downloads/JPR_H2_FINAL_DIAGNOSTICS.R')
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(fixest) })
ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
PANEL <- file.path(ROOT, 'outputs', 'final_airfield_screen_20260809_FINAL30',
                   '11_final_analysis_panel.rds')
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_FINAL_DIAGNOSTICS')
if (!file.exists(PANEL)) stop('Missing FINAL30 file: ', PANEL)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
D <- as.data.table(readRDS(PANEL))
D[, month := as.IDate(month)]
D[, cell_id := as.factor(cell_id)]
mods <- c('pre_contested_share_spatial', 'pre_taliban_activity_share_spatial',
          'pre_contested_intensity_spatial')
ctrl <- c('prec_gpcp', 'temp', 'nlights_calib_mean')
need <- c('events_type3','post_closure','eli_final','cell_id','month',mods,ctrl)
if (length(missing <- setdiff(need,names(D))))
  stop('Missing required columns: ',paste(missing,collapse=', '))
if (anyDuplicated(D[,.(cell_id,month)])) stop('Repeated cell-month observations')
if (!all(na.omit(unique(D$post_closure)) %in% c(0,1))) stop('Post is not binary')

# Locked original H2 benchmarks: coefficient, clustered SE, sample size.
bench <- data.table(
 window=rep(c('full','drop_anticipation'),each=3L), moderator=rep(mods,2L),
 beta=c(-16.1023802013192,-17.6667522797993,-3.12067660696125,
        -17.4082078077671,-19.0707077670960,-3.38898087940946),
 se=c(4.27319350902624,4.60318002832200,0.880878277343997,
      4.86739390345531,5.26815958334537,0.984560054320593),
 n=rep(c(54312L,49956L), each=3L))
formula_h2 <- function(m, extra=character()) {
  x <- c('post_closure:eli_final', paste('post_closure',m,sep=':'),
         paste('post_closure','eli_final',m,sep=':'), extra, ctrl)
  as.formula(paste('events_type3 ~',paste(x,collapse=' + '),'| cell_id + month'))
}
formula_original <- function(m) as.formula(paste('events_type3 ~',
  paste(c('post_closure:eli_final',paste('post_closure','eli_final',m,sep=':'),ctrl),
        collapse=' + '),'| cell_id + month'))
coefrow <- function(fit,term){
  tab <- coeftable(fit,vcov=~cell_id+month)
  if (!term %in% rownames(tab)) stop('Missing model coefficient: ',term)
  as.numeric(tab[term,c(1,2,4)])
}
fitfun <- function(f,d) fepois(f,data=d,notes=FALSE,warn=FALSE)
# Cell-level variables are expected to be fixed. Reject silent time variation.
for (m in c('eli_final',mods)) {
  n_per_cell <- D[, .(n_values=uniqueN(get(m),na.rm=TRUE)), by=cell_id]
  if (any(n_per_cell$n_values>1L))
    stop('Not time-invariant within cell: ',m)
}

win <- list(full=D,
  drop_anticipation=D[month < as.IDate('2020-03-01') |
                      month >= as.IDate('2021-07-01')])
reproduction <- list(); common <- list(); bins <- list(); marg <- list(); placebo <- list()

# Group labels use distinct estimation cells and moderator values, never outcomes.
make_groups <- function(cell,m){
  positive <- cell[get(m)>0 & is.finite(get(m)),get(m)]
  if (length(positive)<2L) stop('Too few positive moderator cells for ',m)
  cut <- median(positive)
  cell[,group:=fifelse(get(m)==0,'zero',
     fifelse(get(m)<=cut,'positive_lower','positive_upper'))]
  cell[,group:=factor(group,levels=c('zero','positive_lower','positive_upper'))]
  cell
}
for(w in names(win)) for(m in mods){
  cat('\n===== ',w,' / ',m,' =====\n',sep='')
  data <- win[[w]]
  tri <- paste('post_closure','eli_final',m,sep=':')
  term_eli <- 'post_closure:eli_final'
  original <- fitfun(formula_original(m),data)
  b <- coefrow(original,tri)
  ww<-w; mm<-m
  target <- bench[window==ww & moderator==mm]
  if(nrow(target)!=1L) stop('Benchmark selection error')
  verified <- nobs(original)==target$n[[1]] &&
    abs(b[1]-target$beta[[1]])<0.001 && abs(b[2]-target$se[[1]])<0.001
  cat('Published H2 exact benchmark reproduced: ',verified,
      ' | beta=',round(b[1],4),' SE=',round(b[2],4),' N=',nobs(original),'\n',sep='')
  if(!verified) stop('Benchmark mismatch: do NOT interpret following results')
  expanded <- fitfun(formula_h2(m),data)
  if(!identical(sort(as.integer(obs(original))),sort(as.integer(obs(expanded)))))
    stop('Expanded model has different estimation rows: ',w,' / ',m)
  ce <- coefrow(expanded,tri)
  reproduction[[length(reproduction)+1L]] <- data.table(window=w,moderator=m,
     original_beta=b[1],original_se=b[2],original_n=nobs(original),
     expanded_beta=ce[1],expanded_se=ce[2],expanded_p=ce[3],
     expanded_n=nobs(expanded),matching_rows=TRUE)

  # 1) COMMON SUPPORT.  Exactly one row per estimation cell, not 432 repeated months.
  cells <- unique(data[as.integer(obs(expanded)),
                  .(cell_id,eli_final,m_value=get(m))])
  if(anyDuplicated(cells$cell_id)) stop('Cell exposure/moderator not constant')
  setnames(cells,'m_value',m)
  cells <- make_groups(cells,m)
  # Within each group: 10-90 ELI band; common 10-90 overlap across observed groups.
  grp <- cells[!is.na(group) & is.finite(eli_final), .(
     n_cells=.N,eli_min=min(eli_final),eli_p10=as.numeric(quantile(eli_final,.1)),
     eli_median=median(eli_final),eli_p90=as.numeric(quantile(eli_final,.9)),
     eli_max=max(eli_final),moderator_median=median(get(m))),by=group]
  group_levels <- as.character(grp$group[grp$n_cells>=5L])
  if(length(group_levels)>=2L){
    low <- max(grp[group %in% group_levels,eli_p10]); high <- min(grp[group %in% group_levels,eli_p90])
    has_overlap <- is.finite(low) && is.finite(high) && low<high
  } else {low<-NA_real_; high<-NA_real_; has_overlap<-FALSE}
  frac <- if(has_overlap) mean(cells$eli_final>=low & cells$eli_final<=high,na.rm=TRUE) else NA_real_
  common[[length(common)+1L]] <- data.table(window=w, moderator=m,
    n_distinct_cells=nrow(cells),n_groups_with_5plus_cells=length(group_levels),
    common_band_low=low,common_band_high=high,nonempty_common_band=has_overlap,
    fraction_cells_in_common_band=frac,
    raw_corr_eli_m=cor(cells$eli_final,cells[[m]],use='complete.obs'))
  grp[,`:=`(window=w,moderator=m,common_low=low,common_high=high)]
  bins[[length(bins)+1L]] <- grp
  cat('Common 10-90 ELI band across groups: ',
      if(has_overlap) paste(round(low,3),'to',round(high,3)) else 'NONE',
      ' ; fraction of cells in band: ',round(frac,3),'\n',sep='')

  # 2) MARGINAL EFFECTS from *fully specified expanded* model.
  # Slope of post/pre log-rate contrast in ELI at M=m0: beta_postELI + m0*beta_triple.
  # The moderation contrast is (Mhigh-Mlow)*beta_triple; test explicitly.
  V <- vcov(expanded,vcov=~cell_id+month)
  cf <- coef(expanded)
  if(!all(c(term_eli,tri)%in% names(cf))) stop('Model dropped required term')
  mvals <- c(q25=as.numeric(quantile(cells[[m]],.25,na.rm=TRUE)),
             q75=as.numeric(quantile(cells[[m]],.75,na.rm=TRUE)))
  # If the 25th percentile equals 75th, retain 0 and positive-median instead.
  if(mvals[1]==mvals[2]){
     mvals <- c(zero=0,positive_median=median(cells[[m]][cells[[m]]>0],na.rm=TRUE))
  }
  effect <- function(a,b,contrast){
    est <- as.numeric(a*cf[[term_eli]]+b*cf[[tri]])
    se <- sqrt(max(0,a*a*V[term_eli,term_eli]+b*b*V[tri,tri]+2*a*b*V[term_eli,tri]))
    data.table(estimate_log_rate=est,se_log_rate=se,
      rate_ratio=exp(est),ci95_low=exp(est-1.96*se),ci95_high=exp(est+1.96*se),
      p_normal_approx=2*pnorm(-abs(est/se)),meaning=contrast)
  }
  nm<-names(mvals)
  for(j in seq_along(mvals)) marg[[length(marg)+1L]] <- cbind(
     data.table(window=w,moderator=m,moderator_level=nm[j],moderator_value=mvals[j]),
     effect(1,mvals[j],'Post-change ELI slope at specified M (per 1 ELI unit)'))
  delta<-unname(diff(mvals))
  marg[[length(marg)+1L]] <- cbind(
     data.table(window=w,moderator=m,moderator_level='HIGH_MINUS_LOW',moderator_value=delta),
     effect(0,delta,'Difference between high and low M ELI slopes (H2 contrast)'))
  cat('Expanded H2 p=',signif(ce[3],4),
      ' ; marginal contrast p=',signif(2*pnorm(-abs(ce[1]/ce[2])),4),
      ' (normal approximation)\n',sep='')

  # 3) PRE-WITHDRAWAL PLACEBO SPLITS. Use only PRE-DOHA months (no anticipation).
  # A genuine pre-trend diagnostic cannot prove parallel trends. These two
  # falsification tests ask whether *similar* triple-interaction shifts already
  # happened before the actual shock, using the complete lower-order structure.
  # In placebo subsamples, post_fake=1 after the indicated cut-off.
  if(w=='full'){
    pre <- copy(D[month < as.IDate('2020-03-01')])
    for(date in c('2017-01-01','2019-01-01')){
      pre[,placebo_post:=as.integer(month>=as.IDate(date))]
      trp<-paste('placebo_post','eli_final',m,sep=':')
      ff<-as.formula(paste('events_type3 ~',paste(c(
          'placebo_post:eli_final',paste('placebo_post',m,sep=':'),trp,ctrl),
          collapse=' + '),'| cell_id + month'))
      tryCatch({
        pf<-fitfun(ff,pre)
        rr<-coefrow(pf,trp)
        placebo[[length(placebo)+1L]] <- data.table(
          placebo_start=date,moderator=m,beta=rr[1],se=rr[2],p=rr[3],
          n=nobs(pf),error='')
        cat('Pre-Doha placebo ',date,': beta=',round(rr[1],3),
            ', p=',signif(rr[3],3),'\n',sep='')
      }, error=function(e){
        placebo[[length(placebo)+1L]] <<- data.table(placebo_start=date,
          moderator=m,beta=NA_real_,se=NA_real_,p=NA_real_,
          n=NA_integer_,error=conditionMessage(e))
        cat('Placebo failed for ',m,' / ',date,': ',conditionMessage(e),'\n',sep='')
      })
    }
  }
}
res1<-rbindlist(reproduction,fill=TRUE)
res2<-rbindlist(common,fill=TRUE)
res3<-rbindlist(bins,fill=TRUE)
res4<-rbindlist(marg,fill=TRUE)
res5<-rbindlist(placebo,fill=TRUE)
fwrite(res1,file.path(OUT,'01_REPRODUCTION_AND_EXPANDED_H2.csv'))
fwrite(res2,file.path(OUT,'02_COMMON_SUPPORT_SUMMARY.csv'))
fwrite(res3,file.path(OUT,'03_ELI_BY_CONTESTEDNESS_GROUP.csv'))
fwrite(res4,file.path(OUT,'04_MARGINAL_EFFECTS_AND_H2_CONTRAST.csv'))
fwrite(res5,file.path(OUT,'05_PRE_DOHA_PLACEBO_TESTS.csv'))
cat('\n========== SUMMARY =========\n')
cat('Reproduced benchmarks:',nrow(res1),'of 6.\n')
cat('All paired samples identical:',all(res1$matching_rows),'\n')
cat('\nCommon support:\n'); print(res2[,.(window,moderator,n_distinct_cells,
  nonempty_common_band,fraction_cells_in_common_band)],digits=3)
cat('\nExpanded-model H2 marginal contrasts (NOT individual slopes):\n')
print(res4[moderator_level=='HIGH_MINUS_LOW',
  .(window,moderator,estimate_log_rate,se_log_rate,p_normal_approx)],digits=3)
cat('\nPlacebo tests (diagnostic, not proof of parallel trends):\n')
print(res5[,.(placebo_start,moderator,beta,se,p,n,error)],digits=3)
cat('\nSaved five CSVs to ',OUT,'\n',sep='')
cat('Important: significant marginal ELI slopes within one contestedness level\n',
    'do NOT establish H2. The HIGH_MINUS_LOW difference tests moderation.\n',
    'Two placebo dates are correlated checks and should not be treated as\n',
    'independent confirmation. Empty common support limits identification.\n',sep='')
