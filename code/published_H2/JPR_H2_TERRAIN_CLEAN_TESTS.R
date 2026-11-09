# ============================================================================
# JPR H2: PRE-SPECIFIED TERRAIN-BASED TESTS (FINAL30 ONLY)
# Goal: test whether terrain observability moderates the post-closure change
# differently in more-contested areas, while ALWAYS allowing Post x M.
# This script does NOT select specifications based on significance.
# Requires: install.packages(c('data.table','fixest'))
# Run: source('C:/Users/USERNAME/Downloads/JPR_H2_TERRAIN_CLEAN_TESTS.R')
# ============================================================================
suppressPackageStartupMessages({library(data.table); library(fixest)})

ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
PANEL <- file.path(ROOT, 'outputs', 'final_airfield_screen_20260809_FINAL30',
                   '11_final_analysis_panel.rds')
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_TERRAIN_CLEAN_TESTS')
if (!file.exists(PANEL)) stop('FINAL30 panel not found: ', PANEL)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
D <- as.data.table(readRDS(PANEL))
mods <- c('pre_contested_share_spatial','pre_taliban_activity_share_spatial',
          'pre_contested_intensity_spatial')
controls <- c('prec_gpcp','temp','nlights_calib_mean')
terrain <- 'eli_feas_terrain'  # Existing source variable; NOT composite ELI.
need <- c('events_type3','post_closure','eli_final',terrain,'cell_id','month',mods,controls)
missing <- setdiff(need,names(D))
if (length(missing)) stop('Required FINAL30 column(s) missing: ',paste(missing,collapse=', '),
                         '\nAvailable columns: ',paste(names(D),collapse=', '))
D[,month:=as.IDate(month)]
D[,cell_id:=as.factor(cell_id)]
if (!all(na.omit(unique(D$post_closure)) %in% c(0,1)))
  stop('post_closure must be 0/1')

# Confirm the terrain variable is pre-shock/time-invariant at the CELL level.
# Scale and cutoff are determined ONCE, across unique cells, not per model.
terrain_audit <- D[,.(n_terrain=uniqueN(na.omit(get(terrain))),
                     min_terrain=if(all(is.na(get(terrain)))) NA_real_ else min(get(terrain),na.rm=TRUE),
                     max_terrain=if(all(is.na(get(terrain)))) NA_real_ else max(get(terrain),na.rm=TRUE)),by=cell_id]
if (any(terrain_audit$n_terrain > 1L))
  stop('Terrain changes within some cells; verify its measurement before proceeding.')
cells <- unique(D[,.(cell_id,terrain_raw=get(terrain),eli_final)])
if (anyDuplicated(cells$cell_id)) stop('ELI or terrain changes within cells. Audit first.')
mu_t <- mean(cells$terrain_raw,na.rm=TRUE)
sd_t <- sd(cells$terrain_raw,na.rm=TRUE)
med_t <- median(cells$terrain_raw,na.rm=TRUE)
if (!is.finite(sd_t) || sd_t <= 0) stop('Terrain variable has no variation.')
D[,terrain_z := (get(terrain)-mu_t)/sd_t]
D[,terrain_high := as.integer(get(terrain) >= med_t)]
cat('Terrain measure:', terrain, '\n')
cat('Cell-level mean:',signif(mu_t,6),'SD:',signif(sd_t,6),
    'median:',signif(med_t,6),'cells:',nrow(cells),'\n')

# Source-verified exact benchmark for ORIGINAL published H2; fails closed.
expected <- data.table(
  window=rep(c('full','drop_anticipation'),each=3),
  moderator=rep(mods,2),
  beta=c(-16.1023802013192,-17.6667522797993,-3.12067660696125,
         -17.4082078077671,-19.0707077670960,-3.38898087940946),
  se=c(4.27319350902624,4.60318002832200,0.880878277343997,
       4.86739390345531,5.26815958334537,0.984560054320593),
  n=c(54312L,54312L,54312L,49956L,49956L,49956L))

extract <- function(model,term){
  ct <- coeftable(model,vcov=~cell_id+month)
  if (!term %in% rownames(ct)) return(c(beta=NA_real_,se=NA_real_,p=NA_real_))
  setNames(as.numeric(ct[term,c(1,2,4)]),c('beta','se','p'))
}
make_fml <- function(terms){
  as.formula(paste('events_type3 ~',paste(c(terms,controls),collapse=' + '),
                   '| cell_id + month'))
}
fit <- function(formula,data,label){
  cat('  Estimating:',label,'\n')
  tryCatch(fepois(formula,data=data,notes=FALSE,warn=FALSE),
           error=function(e){
             cat('  MODEL ERROR:',conditionMessage(e),'\n'); NULL
           })
}
record <- function(model,window,mod,spec,term,reference=NULL){
  if (is.null(model)) return(data.table(window=window,moderator=mod,
    specification=spec,term=term,beta=NA_real_,se=NA_real_,p=NA_real_,
    n=NA_integer_,sample_matches_reference=NA))
  z <- extract(model,term)
  same <- if (is.null(reference)) NA else
    identical(sort(as.integer(obs(model))),sort(as.integer(obs(reference))))
  data.table(window=window,moderator=mod,specification=spec,term=term,
             beta=unname(z['beta']),se=unname(z['se']),p=unname(z['p']),
             n=as.integer(nobs(model)),sample_matches_reference=same)
}

windows <- list(full=D,
  drop_anticipation=D[month < as.IDate('2020-03-01') |
                      month >= as.IDate('2021-07-01')])
results <- list(); summaries <- list(); k<-0L
for(w in names(windows)){
  for(m in mods){
    cat('\n=== ',w,' / ',m,' ===\n',sep='')
    dat <- windows[[w]]
    oldterm <- paste('post_closure','eli_final',m,sep=':')
    old <- fit(make_fml(c('post_closure:eli_final',oldterm)),dat,'published benchmark')
    if(is.null(old)) stop('Published benchmark model failed. Aborting.')
    cw <- w; cm <- m
    bench <- expected[window==cw & moderator==cm]
    if(nrow(bench)!=1L) stop('Benchmark lookup failed.')
    ob <- extract(old,oldterm)
    ok <- nobs(old)==bench$n[[1L]] &&
      abs(ob['beta']-bench$beta[[1L]]) < 0.001 &&
      abs(ob['se']-bench$se[[1L]]) < 0.001
    cat('Published benchmark reproduced:',ok,
        '; beta:',round(ob['beta'],4),'; SE:',round(ob['se'],4),
        '; N:',nobs(old),'\n')
    if(!ok) stop('Benchmark differs; no new H2 results should be interpreted.')
    k <- k+1L
    results[[length(results)+1L]] <- record(old,w,m,'published_eli_restricted',oldterm)

    # PRIMARY: fully specified continuous terrain interaction.
    # Post x M is always included. Cell FE absorb M, terrain and M x terrain;
    # month FE absorb Post. This tests H2 using the terrain component only.
    tri <- paste('post_closure','terrain_z',m,sep=':')
    pm <- paste('post_closure',m,sep=':')
    pt <- 'post_closure:terrain_z'
    primary <- fit(make_fml(c(pt,pm,tri)),dat,'PRIMARY continuous terrain + Post x M')
    results[[length(results)+1L]] <- record(primary,w,m,'PRIMARY_terrain_continuous',tri,old)
    results[[length(results)+1L]] <- record(primary,w,m,'PRIMARY_post_x_m',pm,old)
    results[[length(results)+1L]] <- record(primary,w,m,'PRIMARY_post_x_terrain',pt,old)

    # SECONDARY: same full hierarchy, but median-split terrain (a coarse
    # functional-form check, NOT an independent identification strategy).
    trihi <- paste('post_closure','terrain_high',m,sep=':')
    pthi <- 'post_closure:terrain_high'
    binary <- fit(make_fml(c(pthi,pm,trihi)),dat,
                  'SECONDARY binary terrain + Post x M')
    results[[length(results)+1L]] <- record(binary,w,m,'SECONDARY_terrain_median_split',trihi,old)

    # DIAGNOSTICS: distinct-cell terrain/moderator correlation and support.
    # Not hypothesis tests; useful for checking whether moderation is
    # identified across the actual range, rather than extrapolation.
    S <- unique(dat[,.(cell_id,terrain_z,terrain_high,m_value=get(m),eli_final)])
    if (anyDuplicated(S$cell_id)) stop('Moderator varies within cell: ',m)
    row <- data.table(window=w,moderator=m,original_n=nobs(old),
      primary_n=if(is.null(primary)) NA_integer_ else nobs(primary),
      binary_n=if(is.null(binary)) NA_integer_ else nobs(binary),
      primary_same_rows=if(is.null(primary)) NA else
        identical(sort(as.integer(obs(primary))),sort(as.integer(obs(old)))),
      binary_same_rows=if(is.null(binary)) NA else
        identical(sort(as.integer(obs(binary))),sort(as.integer(obs(old)))),
      n_cells=nrow(S),n_nonzero_m_cells=sum(S$m_value>0,na.rm=TRUE),
      corr_terrain_m=if(sd(S$m_value,na.rm=TRUE)>0)
        cor(S$terrain_z,S$m_value,use='complete.obs') else NA_real_,
      terrain_p25=as.numeric(quantile(S$terrain_z,.25,na.rm=TRUE)),
      terrain_p75=as.numeric(quantile(S$terrain_z,.75,na.rm=TRUE)),
      m_p25=as.numeric(quantile(S$m_value,.25,na.rm=TRUE)),
      m_p75=as.numeric(quantile(S$m_value,.75,na.rm=TRUE)))
    summaries[[length(summaries)+1L]] <- row
    if (!is.null(primary)){
      r <- extract(primary,tri)
      cat('PRIMARY three-way beta:',round(r['beta'],4),
          'SE:',round(r['se'],4),'p:',signif(r['p'],4),
          '; N:',nobs(primary),'; same rows:',row$primary_same_rows,'\n')
    }
    if(!is.null(binary)){
      r <- extract(binary,trihi)
      cat('SECONDARY median-split beta:',round(r['beta'],4),
          'SE:',round(r['se'],4),'p:',signif(r['p'],4),
          '; N:',nobs(binary),'; same rows:',row$binary_same_rows,'\n')
    }
  }
}
R <- rbindlist(results,fill=TRUE)
S <- rbindlist(summaries,fill=TRUE)
fwrite(R,file.path(OUT,'01_ALL_H2_TERRAIN_COEFFICIENTS.csv'))
fwrite(S,file.path(OUT,'02_SAMPLE_AND_SUPPORT_CHECKS.csv'))
cat('\n=== PRIMARY TERRAIN THREE-WAY ESTIMATES ===\n')
print(R[specification=='PRIMARY_terrain_continuous',
        .(window,moderator,beta,se,p,n,sample_matches_reference)],digits=4)
cat('\n=== SECONDARY MEDIAN-SPLIT THREE-WAY ESTIMATES ===\n')
print(R[specification=='SECONDARY_terrain_median_split',
        .(window,moderator,beta,se,p,n,sample_matches_reference)],digits=4)
cat('\nFiles saved to:',OUT,'\n')
cat('Interpretation: a terrain-based interaction consistent with H2 is not\n',
    'by itself evidence that the mechanism is uniquely causal. The two tests\n',
    'are not independent confirmation; inspect all moderators/windows.\n',sep='')
