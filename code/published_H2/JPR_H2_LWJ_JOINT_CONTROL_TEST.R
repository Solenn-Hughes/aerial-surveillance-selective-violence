# JPR H2 — historical Long War Journal territorial-control measurement
# Inputs: published FINAL30, Dudek LWJ periods, documented district crosswalk.
# Does not overwrite FINAL30. No model selection based on statistical significance.
# Runs six published H2 benchmark replications, then matched LWJ H2 tests.

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

ROOT <- 'C:/Users/USERNAME/OneDrive/Desktop/THESIS'
DOWNLOADS <- 'C:/Users/USERNAME/Downloads'
PANEL <- file.path(ROOT, 'outputs', 'final_airfield_screen_20260809_FINAL30',
                   '11_final_analysis_panel.rds')
LWJ <- file.path(DOWNLOADS, 'LWJ_Dudek_district_periods.csv')
CROSSWALK <- file.path(DOWNLOADS, 'LWJ_FINAL30_VERIFIED_DISTRICT_CROSSWALK.csv')
OUT <- file.path(ROOT, 'outputs', 'JPR_H2_LWJ_JOINT_CONTROL_20261106')
for (p in c(PANEL, LWJ, CROSSWALK)) if (!file.exists(p)) stop('Input missing: ', p)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

L <- fread(LWJ, encoding='UTF-8')
CW <- fread(CROSSWALK, encoding='UTF-8', colClasses=list(character=c('adm2_id','lwj_object_id')))
D <- as.data.table(readRDS(PANEL))
D[, `:=`(cell_id=as.character(cell_id),
         month=as.IDate(month),
         district_id_for_h2=as.character(district_id_for_h2))]
L[, lwj_object_id := as.character(lwj_object_id)]
CW[, lwj_object_id := as.character(lwj_object_id)]

required <- c('cell_id','month','district_id_for_h2','events_type3', 'post_closure',
              'eli_final','prec_gpcp','temp','nlights_calib_mean',
              'pre_contested_share_spatial','pre_taliban_activity_share_spatial',
              'pre_contested_intensity_spatial')
if (length(missing <- setdiff(required,names(D)))) stop('FINAL30 columns missing: ',paste(missing,collapse=', '))
if (length(missing <- setdiff(c('adm2_id','lwj_object_id','match_basis'),names(CW))))
  stop('Crosswalk columns missing: ',paste(missing,collapse=', '))
if (length(missing <- setdiff(c('lwj_object_id','period','control'),names(L))))
  stop('LWJ columns missing: ',paste(missing,collapse=', '))
if (nrow(D)!=104112L || uniqueN(D$cell_id)!=241L || anyDuplicated(D,by=c('cell_id','month')))
  stop('FINAL30 integrity mismatch')
if (nrow(L)!=3663L || uniqueN(L$lwj_object_id)!=407L ||
    anyDuplicated(L,by=c('lwj_object_id','period')) || !setequal(L$period,1:9))
  stop('LWJ source structure mismatch')
if (!all(L$control %in% c('GOA','Taliban','Contested','Unconfirmed')))
  stop('Unexpected LWJ control value')
if (anyDuplicated(CW$adm2_id) || anyDuplicated(CW$lwj_object_id))
  stop('Crosswalk is not one-to-one: investigate before estimating')
if (!all(CW$lwj_object_id %in% L$lwj_object_id))
  stop('Crosswalk contains LWJ IDs not in source')
cell <- unique(D[,.(cell_id,district_id_for_h2)])
if (anyDuplicated(cell$cell_id)) stop('Panel cell assigned to multiple districts')
cell[, lwj_object_id := CW$lwj_object_id[match(district_id_for_h2, CW$adm2_id)]]
cat('\nCROSSWALK: ',uniqueN(cell$district_id_for_h2),' panel districts, ',
    nrow(cell),' cells; matched cells ',sum(!is.na(cell$lwj_object_id)),'/',nrow(cell),'\n',sep='')
fwrite(cell,file.path(OUT,'01_cell_join_audit.csv'))
if (anyNA(cell$lwj_object_id)) stop('Unmatched panel cells: see 01_cell_join_audit.csv')
if (uniqueN(cell$district_id_for_h2)!=168L) stop('Expected 168 panel districts')
if (anyDuplicated(cell$lwj_object_id)) {
  # Several PRIO-GRID cells can share one historical district: normal.
  if (uniqueN(cell[,.(district_id_for_h2,lwj_object_id)]) != uniqueN(cell$district_id_for_h2))
    stop('Distinct panel districts share source ID')
}
if (uniqueN(cell[,.(district_id_for_h2,lwj_object_id)]) != uniqueN(cell$district_id_for_h2) ||
    anyDuplicated(unique(cell[,.(district_id_for_h2,lwj_object_id)]),by='lwj_object_id'))
  stop('Panel districts are not in one-to-one relation with LWJ districts')

# Only periods 1–5, pre-Doha. Do NOT use outcome-dependent control categories.
P <- L[period %in% 1:5,.(
  lwj_contested_share=mean(control=='Contested'),
  lwj_taliban_control_share=mean(control=='Taliban'),
  lwj_goa_control_share=mean(control=='GOA'),
  lwj_unconfirmed_share=mean(control=='Unconfirmed'),
  lwj_confirmed_contested_share=if (sum(control!='Unconfirmed')>0)
       sum(control=='Contested')/sum(control!='Unconfirmed') else NA_real_,
  n_pre_doha_periods=.N
),by=lwj_object_id]
if (any(P$n_pre_doha_periods!=5L)) stop('Historical periods incomplete')
if (any(abs(P$lwj_contested_share+P$lwj_taliban_control_share+
             P$lwj_goa_control_share+P$lwj_unconfirmed_share-1)>1e-10))
  stop('LWJ statuses fail composition check')
D[, lwj_object_id := cell$lwj_object_id[match(cell_id, cell$cell_id)]]
for (v in c('lwj_contested_share','lwj_taliban_control_share',
            'lwj_goa_control_share','lwj_unconfirmed_share','lwj_confirmed_contested_share'))
  D[, (v) := P[[v]][match(lwj_object_id,P$lwj_object_id)]]
if (anyNA(D$lwj_contested_share)) stop('Missing mapped LWJ district histories')
if (nrow(D)!=104112L) stop('Panel rows unexpectedly changed')
D[,cell_id:=factor(cell_id)]
D[,month:=as.IDate(month)]
fwrite(unique(D[,.(cell_id,district_id_for_h2,lwj_object_id,
                    lwj_contested_share,lwj_taliban_control_share,
                    lwj_goa_control_share,lwj_unconfirmed_share,
                    lwj_confirmed_contested_share)]),
       file.path(OUT,'02_cell_lwj_moderators.csv'))
cat('LWJ pre-Doha status shares prepared, no missing unconditional measures.\n')

fit_h2 <- function(X,moderator,expanded=FALSE){
  rhs <- c('post_closure:eli_final',paste0('post_closure:eli_final:',moderator),
           'prec_gpcp','temp','nlights_calib_mean')
  if (expanded) rhs <- c(rhs,paste0('post_closure:',moderator))
  form <- as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),'| cell_id + month'))
  fit <- fepois(form,data=X,notes=FALSE,warn=FALSE)
  term <- paste0('post_closure:eli_final:',moderator)
  z <- coeftable(fit,vcov=~cell_id+month)
  if (!(term %in% rownames(z))) stop('Three-way interaction dropped: ',term)
  list(fit=fit,beta=unname(z[term,1]),se=unname(z[term,2]),
       p=unname(z[term,4]),n=as.integer(nobs(fit)))
}

bench <- data.table(
  window=rep(c('full','drop_anticipation'),each=3),
  moderator=rep(c('pre_contested_share_spatial',
                  'pre_taliban_activity_share_spatial',
                  'pre_contested_intensity_spatial'),times=2),
  beta_expected=c(-16.1023802,-17.6667523,-3.120677,
                  -17.40820781,-19.07071,-3.38898),
  se_expected=c(4.273193509,4.60318,0.88088,
                4.867393903,5.26816,0.98456),
  n_expected=rep(c(54312L,49956L),each=3))
windows <- list(full=D,
  drop_anticipation=D[month < as.IDate('2020-03-01') |
                      month >= as.IDate('2021-07-01')])
cat('\n=== REQUIRED ORIGINAL FINAL30 REPRODUCTION ===\n')
reps <- list()
for (w in names(windows)) for (m in unique(bench$moderator)) {
  f <- fit_h2(windows[[w]],m,FALSE)
  b <- bench[window==w & moderator==m]
  if (nrow(b)!=1L) stop('Benchmark lookup ambiguous')
  good <- abs(f$beta-b$beta_expected)<0.002 &&
          abs(f$se-b$se_expected)<0.002 && f$n==b$n_expected
  reps[[length(reps)+1L]] <- data.table(window=w,moderator=m,
       beta=f$beta,se=f$se,n=f$n,replicated=good)
  cat(w,m,'beta=',round(f$beta,4),'SE=',round(f$se,4),
      'N=',f$n,'reproduced=',good,'\n')
  if (!good) stop('Published model not reproduced; LWJ results NOT estimated')
}
fwrite(rbindlist(reps),file.path(OUT,'03_original_H2_replication.csv'))


# ===============================================================
# JOINT LWJ TERRITORIAL-CONTROL MODERATION
# GOA share is omitted reference; Contested, Taliban, Unconfirmed
# status shares all included, hence no perfect compositional collinearity.
# Both models include Post x ELI. The FULL model also includes
# Post x each included status share, allowing independent post shifts.
# ===============================================================
status_vars <- c('lwj_contested_share','lwj_taliban_control_share','lwj_unconfirmed_share')
status_labels <- c('Contested vs GOA','Taliban vs GOA','Unconfirmed vs GOA')
status_terms <- paste0('post_closure:eli_final:',status_vars)
post_terms <- paste0('post_closure:',status_vars)

# In a five-period share, e.g. 0.6 means 3 of 5 intervals in that status.
# A triple coefficient represents the change in the ELI slope when moving
# from all GOA periods to all periods of that other status (on log-rate scale).
cell_status <- unique(D[,c('cell_id','district_id_for_h2','lwj_object_id',
                           'lwj_contested_share','lwj_taliban_control_share',
                           'lwj_goa_control_share','lwj_unconfirmed_share'),with=FALSE])
if (nrow(cell_status)!=241L || anyNA(cell_status[,..status_vars]))
  stop('Missing/duplicated cell-level territorial-status shares')
if (any(abs(cell_status$lwj_contested_share+cell_status$lwj_taliban_control_share+
            cell_status$lwj_goa_control_share+cell_status$lwj_unconfirmed_share - 1)>1e-9))
  stop('LWJ status shares do not sum to 1')

cat('\n=== STATUS COMPOSITION AND SUPPORT (ONE ROW PER CELL) ===\n')
summary_status <- rbindlist(lapply(c('lwj_goa_control_share',status_vars),function(m) {
  v <- cell_status[[m]]
  data.table(status=m,cells_with_positive=sum(v>0),cells_all_five=sum(v==1),
             min=min(v),median=median(v),max=max(v),distinct_values=uniqueN(v))
}))
print(summary_status)
fwrite(summary_status,file.path(OUT,'04_joint_status_support.csv'))

# Strictly same underlying sample and status measurement for all models.
# Compare the full model to a restricted specification only as context.
run_joint <- function(X, expanded=TRUE) {
  rhs <- c('post_closure:eli_final',
           status_terms,
           if(expanded) post_terms,
           'prec_gpcp','temp','nlights_calib_mean')
  formula <- as.formula(paste('events_type3 ~',paste(rhs,collapse=' + '),
                              '| cell_id + month'))
  fit <- fepois(formula,data=X,notes=FALSE,warn=FALSE)
  zz <- coeftable(fit,vcov=~cell_id+month)
  needed <- c('post_closure:eli_final', status_terms, if(expanded) post_terms)
  missing_terms <- setdiff(needed,rownames(zz))
  if (length(missing_terms)) stop('Term omitted/renamed: ',paste(missing_terms,collapse=', '),
                                  '\nActual terms: ',paste(rownames(zz),collapse=', '))
  list(fit=fit,tab=zz,formula=deparse(formula,width.cutoff=500L))
}

# Contrasts use the same TWO-WAY CLUSTER-ROBUST covariance matrix as the model.
# Wald normal approximation; includes correlation between estimated coefficients.
contrast <- function(fit, weights, label, window, model) {
  b <- coef(fit)
  V <- vcov(fit,vcov=~cell_id+month)
  if (!all(names(weights) %in% names(b))) stop('Contrast term absent from coefficient vector')
  if (!all(names(weights) %in% rownames(V))) stop('Contrast term absent from covariance matrix')
  w <- setNames(rep(0,length(b)),names(b)); w[names(weights)] <- weights
  est <- sum(w*b)
  v <- as.numeric(t(w)%*%V[names(b),names(b),drop=FALSE]%*%w)
  if (!is.finite(v) || v<=0) stop('Contrast covariance invalid: ',label)
  se <- sqrt(v); z <- est/se
  data.table(window=window,model=model,contrast=label,
             beta=est,se=se,z=z,p=2*pnorm(-abs(z)),
             ci_low=est-qnorm(.975)*se,ci_high=est+qnorm(.975)*se,
             n=as.integer(nobs(fit)))
}

joint_results <- list()
contrast_results <- list()
meta <- list()
for (w in names(windows)) {
  X <- windows[[w]]
  cat('\n=== JOINT LWJ H2 / ',w,' ===\n',sep='')
  restricted <- run_joint(X,FALSE)
  full <- run_joint(X,TRUE)
  a <- as.integer(fixest::obs(restricted$fit))
  b <- as.integer(fixest::obs(full$fit))
  if (!identical(a,b) || length(a)==0 || nobs(restricted$fit)!=nobs(full$fit))
    stop('Different estimation rows in restricted versus expanded joint models for ',w)
  cat('Matched estimation observations: ',nobs(full$fit),
      '; identical models samples: TRUE\n',sep='')
  for (spec in c('restricted','expanded')) {
    f <- if (spec=='restricted') restricted else full
    for (term in c('post_closure:eli_final',status_terms,if(spec=='expanded')post_terms)) {
      z <- f$tab[term,]
      joint_results[[length(joint_results)+1L]] <-
        data.table(window=w,model=spec,term=term,beta=unname(z[1]),
                   se=unname(z[2]),p=unname(z[4]),n=as.integer(nobs(f$fit)))
    }
    # These are the conditional post-closure ELI slopes for PURE control histories.
    # They should not be confused with the differences between these slopes.
    contrast_results[[length(contrast_results)+1L]] <- contrast(f$fit,
      c('post_closure:eli_final'=1),'ELI slope: GOA only',w,spec)
    for (i in seq_along(status_vars)) {
      wts <- c('post_closure:eli_final'=1)
      wts[status_terms[i]] <- 1
      contrast_results[[length(contrast_results)+1L]] <- contrast(f$fit,wts,
        paste0('ELI slope: ',status_labels[i]),w,spec)
    }
    for (i in seq_along(status_vars)) {
      wts <- setNames(1,status_terms[i])
      contrast_results[[length(contrast_results)+1L]] <- contrast(f$fit,wts,
        paste0('ELI moderation: ',status_labels[i]),w,spec)
    }
    wts <- setNames(c(1,-1),status_terms[c(1,2)])
    contrast_results[[length(contrast_results)+1L]] <- contrast(f$fit,wts,
      'ELI moderation difference: Contested minus Taliban',w,spec)
  }
  meta[[length(meta)+1L]] <- data.table(window=w,
      restricted_formula=restricted$formula,
      expanded_formula=full$formula,
      n=as.integer(nobs(full$fit)),same_rows=TRUE)
  cat('Expanded joint triple estimates (relative to GOA):\n')
  print(full$tab[status_terms,,drop=FALSE])
  cat('Expanded joint Post x status terms:\n')
  print(full$tab[post_terms,,drop=FALSE])
}
R <- rbindlist(joint_results,fill=TRUE)
C <- rbindlist(contrast_results,fill=TRUE)
M <- rbindlist(meta,fill=TRUE)
fwrite(R,file.path(OUT,'05_joint_model_coefficients.csv'))
fwrite(C,file.path(OUT,'06_joint_model_contrasts.csv'))
fwrite(M,file.path(OUT,'07_model_formulas_and_sample_checks.csv'))
cat('\n=== MAIN EXPANDED MODEL CONTRASTS (p-values are for differences) ===\n')
print(C[model=='expanded' & grepl('moderation',contrast),
        .(window,contrast,beta=round(beta,4),se=round(se,4),p=signif(p,4),n)])
cat('\nOutput directory: ',OUT,'\n',sep='')
cat('IMPORTANT: LWJ join is a documented name/history crosswalk, not polygon-verified;\n',
    'these estimates remain conditional on mapping accuracy.\n',sep='')
cat('Interpret expanded joint tests before restricted results.\n')
cat('No R execution was performed in preparation of this file.\n')
