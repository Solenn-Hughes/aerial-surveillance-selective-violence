# =====================================================================
# JPR H2 PRE-DOHA MODERATOR SENSITIVITY — READ ONLY
#
# 1. Independently rebuilds all 3 moderators from the original ADM2
#    district-month source, ending BEFORE 2020-03-01.
# 2. Verifies the original pre-July-2021 moderator reconstruction.
# 3. Reproduces all SIX published FINAL30 H2 benchmarks.
# 4. On identical observation windows, estimates:
#       A original restricted, original moderator (benchmark)
#       B original expanded, original moderator (comparison)
#       C PRE-DOHA restricted, rebuilt moderator
#       D PRE-DOHA expanded, rebuilt moderator (key sensitivity)
# Nothing overwrites the FINAL30 panel or original files.
# Packages: data.table, fixest. Run with source(...).
# =====================================================================

suppressPackageStartupMessages({library(data.table); library(fixest)})
ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
PANEL <- file.path(ROOT, 'outputs', 'final_airfield_screen_20260809_FINAL30',
                   '11_final_analysis_panel.rds')
DISTRICT_MONTH <- file.path(ROOT, 'outputs',
  'paper1_geo_afghanistan_executable_acled', 'diagnostics_h2_spatial_adm2',
  'district_month_contestedness_spatial_adm2.csv')
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_PREDOHA_MODERATOR_TEST')
if (!file.exists(PANEL)) stop('Missing FINAL30 panel: ', PANEL)
if (!file.exists(DISTRICT_MONTH)) stop('Missing original district-month source: ', DISTRICT_MONTH)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

MODS <- c('pre_contested_share_spatial', 'pre_taliban_activity_share_spatial',
          'pre_contested_intensity_spatial')
CTRL <- c('prec_gpcp', 'temp', 'nlights_calib_mean')
D <- as.data.table(readRDS(PANEL))
D[, month := as.IDate(month)]
D[, cell_id := as.factor(cell_id)]
need <- c('cell_id','district_id_for_h2','month','events_type3','post_closure',
          'eli_final', MODS, CTRL)
if (length(missing <- setdiff(need, names(D))))
  stop('FINAL30 missing variables: ', paste(missing, collapse=', '))
if (anyDuplicated(D, by=c('cell_id','month')))
  stop('Duplicate FINAL30 cell-month records')
if (any(is.na(D$district_id_for_h2))) stop('Some cells lack ADM2 assignment')

DM <- fread(DISTRICT_MONTH)
if (length(missing <- setdiff(c('district_id','month','contested_any',
                             'taliban_active_any','battles_GT'),names(DM))))
  stop('District-month source missing: ',paste(missing,collapse=', '))
DM[, month := as.IDate(month)]
if (anyDuplicated(DM, by=c('district_id','month')))
  stop('District-month source has duplicate district-months')
for (v in c('contested_any','taliban_active_any','battles_GT')) {
  if (anyNA(DM[[v]])) stop('Source has missing ', v, ' — stopping rather than silently dropping')
}
# Match original audit definition exactly: mean of the binary flags and
# log1p(mean MONTHLY contested battle counts), not log1p(total battles).
aggregate_mod <- function(cutoff) {
  DM[month < as.IDate(cutoff), .(
    contested = mean(contested_any),
    taliban = mean(taliban_active_any),
    intensity = log1p(mean(battles_GT)),
    months = .N), by=district_id]
}
ORIG <- aggregate_mod('2021-07-01')
PRE  <- aggregate_mod('2020-03-01')
if (anyNA(ORIG) || anyNA(PRE)) stop('Missing reconstructed district moderator')
if (any(PRE$months > ORIG$months)) stop('Pre-Doha denominators exceed original')
if (uniqueN(ORIG$district_id) != uniqueN(PRE$district_id))
  stop('ADM2 coverage differs across windows; investigate before interpreting')

# Exact character matching guards against factor/integer district ID issues.
idx_orig <- match(as.character(D$district_id_for_h2),as.character(ORIG$district_id))
idx_pre  <- match(as.character(D$district_id_for_h2),as.character(PRE$district_id))
if (anyNA(idx_orig) || anyNA(idx_pre))
  stop('Unmatched district IDs in one of the reconstructed datasets')

for (j in seq_along(MODS)) {
  m <- MODS[[j]]
  src_var <- c('contested','taliban','intensity')[[j]]
  reconstructed_orig <- ORIG[[src_var]][idx_orig]
  delta <- max(abs(D[[m]] - reconstructed_orig),na.rm=TRUE)
  if (!is.finite(delta) || delta > 1e-7)
    stop('Original moderator does NOT reconstruct: ',m,
         ' ; maximum absolute difference=',format(delta,digits=12))
  D[, (paste0('predoha_',m)) := PRE[[src_var]][idx_pre]]
}
cat('Original district-month moderator construction verified for all three variables.\n')
cat('Pre-Doha moderators rebuilt using months before March 2020.\n')

# Per-cell comparison, no monthly duplication in summaries.
CELL <- unique(D[,c('cell_id','district_id_for_h2',MODS,paste0('predoha_',MODS)),with=FALSE])
if (anyDuplicated(CELL$cell_id)) stop('Cell-level moderation variables are not constant')
change <- rbindlist(lapply(MODS,function(m) {
  a <- CELL[[m]]; b <- CELL[[paste0('predoha_',m)]]
  data.table(moderator=m, n_cells=length(a), n_changed=sum(abs(a-b)>1e-10),
    share_changed=mean(abs(a-b)>1e-10),
    correlation=cor(a,b), median_original=median(a),median_predoha=median(b),
    median_absolute_change=median(abs(a-b)),
    max_absolute_change=max(abs(a-b)))
}))
fwrite(change,file.path(OUT,'01_MODERATOR_MEASUREMENT_CHANGES.csv'))
cat('\n===== MEASUREMENT CHANGES (one row per cell) =====\n')
print(change)

# Source-verified published FINAL30 coefficient/SE/N benchmarks.
EXPECT <- data.table(
  window=rep(c('full','drop_anticipation'),each=3L),
  moderator=rep(MODS,2L),
  beta=c(-16.1023802013192,-17.6667522797993,-3.12067660696125,
         -17.4082078077671,-19.0707077670960,-3.38898087940946),
  se=c(4.27319350902624,4.60318002832200,0.880878277343997,
       4.86739390345531,5.26815958334537,0.984560054320593),
  n=rep(c(54312L,49956L),each=3L))

fit_one <- function(df, m, modname, expanded, window, specification) {
  tri <- paste('post_closure','eli_final',modname,sep=':')
  low <- paste('post_closure',modname,sep=':')
  rhs <- c('post_closure:eli_final', if(expanded) low,tri,CTRL)
  f <- as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),
                         '| cell_id + month'))
  fit <- fepois(f,data=df,notes=FALSE,warn=FALSE)
  z <- coeftable(fit,vcov=~cell_id+month)
  if (!(tri %in% rownames(z))) stop('Three-way term dropped: ',tri)
  if (expanded && !(low %in% rownames(z))) stop('Post x moderator term dropped: ',low)
  rec <- data.table(window=window,moderator=m,specification=specification,
    moderator_measure=ifelse(grepl('predoha',specification),'pre_Doha','pre_July2021'),
    expanded=expanded,beta=as.numeric(z[tri,1]),
    se=as.numeric(z[tri,2]),p=as.numeric(z[tri,4]),
    post_m_beta=if(expanded)as.numeric(z[low,1]) else NA_real_,
    post_m_se=if(expanded)as.numeric(z[low,2]) else NA_real_,
    post_m_p=if(expanded)as.numeric(z[low,4]) else NA_real_,
    n=as.integer(nobs(fit)),
    formula=paste(deparse(f),collapse=''))
  list(record=rec,fit=fit)
}

# Both models within a pair use the identical dataset rows and same model
# terms except the one added post x moderator interaction.
# Note fixest::obs() returns row positions within the supplied df.
extract_obs <- function(fit) as.integer(fixest::obs(fit))
windows <- list(
  full = D,
  drop_anticipation = D[month < as.IDate('2020-03-01') |
                         month >= as.IDate('2021-07-01')]
)
all_res <- list(); all_checks <- list(); k <- 0L
for (w in names(windows)) for (m in MODS) {
  cat('\n===== ',w,' / ',m,' =====\n',sep='')
  dat <- windows[[w]]
  pm <- paste0('predoha_',m)
  model_defs <- list(
    original_restricted=list(var=m,expanded=FALSE),
    original_expanded=list(var=m,expanded=TRUE),
    predoha_restricted=list(var=pm,expanded=FALSE),
    predoha_expanded=list(var=pm,expanded=TRUE)
  )
  fitted <- list()
  for (nm in names(model_defs)) {
    spec <- model_defs[[nm]]
    fitted[[nm]] <- fit_one(dat,m,spec$var,spec$expanded,w,nm)
    k <- k+1L; all_res[[k]] <- fitted[[nm]]$record
  }
  # Explicitly scoped benchmark lookup — avoids data.table scoping bug.
  target <- EXPECT[window == w & moderator == m]
  if(nrow(target) != 1L) stop('Benchmark lookup returned ',nrow(target),' rows')
  actual <- fitted$original_restricted$record
  reproduced <- abs(actual$beta-target$beta)<1e-3 &&
    abs(actual$se-target$se)<1e-3 && actual$n==target$n
  cat('Published model reproduced: ',reproduced,
      ' ; beta=',round(actual$beta,5),' ; SE=',round(actual$se,5),
      ' ; N=',actual$n,'\n',sep='')
  if(!reproduced) stop('Published coefficient did NOT reproduce; stopping without interpreting new results.')
  refrows <- sort(extract_obs(fitted$original_restricted$fit))
  rowcompare <- vapply(fitted,function(item) {
    identical(refrows,sort(extract_obs(item$fit)))
  },logical(1))
  all_checks[[length(all_checks)+1L]] <- data.table(
    window=w,moderator=m,published_reproduced=reproduced,
    original_expanded_same_rows=rowcompare[['original_expanded']],
    predoha_restricted_same_rows=rowcompare[['predoha_restricted']],
    predoha_expanded_same_rows=rowcompare[['predoha_expanded']],
    n_original=actual$n,n_predoha_expanded=fitted$predoha_expanded$record$n)
  cat('Same estimation rows (all 4 models): ',all(rowcompare), '\n',sep='')
  if(!all(rowcompare)) cat('WARNING: samples differ. Comparisons need qualification.\n')
  pre_r <- fitted$predoha_restricted$record
  pre_e <- fitted$predoha_expanded$record
  cat('Pre-Doha restricted triple: beta=',round(pre_r$beta,4),
      ' SE=',round(pre_r$se,4),' p=',signif(pre_r$p,4),'\n',sep='')
  cat('Pre-Doha EXPANDED triple: beta=',round(pre_e$beta,4),
      ' SE=',round(pre_e$se,4),' p=',signif(pre_e$p,4),'\n',sep='')
}
RES <- rbindlist(all_res,fill=TRUE)
CHECKS <- rbindlist(all_checks,fill=TRUE)
fwrite(RES,file.path(OUT,'02_ALL_MODEL_RESULTS.csv'))
fwrite(CHECKS,file.path(OUT,'03_REPRODUCTION_AND_SAMPLE_CHECKS.csv'))
WIDE <- dcast(RES,window+moderator~specification,
              value.var=c('beta','se','p','n'))
fwrite(WIDE,file.path(OUT,'04_SIDE_BY_SIDE_COMPARISON.csv'))
cat('\n===== PRIMARY RESULT: PRE-DOHA EXPANDED MODELS =====\n')
print(RES[specification=='predoha_expanded',
          .(window,moderator,beta=round(beta,4),
            se=round(se,4),p=signif(p,4),n)])
cat('\n===== AUDIT CHECKS =====\n');print(CHECKS)
cat('\nOutputs saved to: ',OUT,'\n',sep='')
cat('Interpretation rule: compare expanded pre-Doha to expanded original,\n',
    'and verify same samples; do not choose a moderator definition by its p-value.\n',sep='')
