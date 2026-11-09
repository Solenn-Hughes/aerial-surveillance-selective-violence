# JPR H2 / diagnostic 06: separate territorial-status time shocks from ELI moderation.
# Needs JPR_H2_LWJ_00_COMMON.R + LWJ input CSVs in Windows Downloads.
# This is a sensitivity diagnostic, NOT independent causal identification.
# NEVER automatically change specifications in response to significance.

COMMON <- 'C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R'
if (!file.exists(COMMON)) stop('Missing required common file: ', COMMON)
source(COMMON, local = TRUE)

DIR06 <- file.path(OUT, '06_status_specific_month_shocks')
dir.create(DIR06, recursive = TRUE, showWarnings = FALSE)

# Use a post-2015 comparison era, avoiding four decades of unrelated conflicts.
# LWJ pre-Doha status is frozen at periods 1-5 (2017-2019).
START <- as.IDate('2015-01-01')
END <- as.IDate('2024-12-01')
BASE <- copy(D[month >= START & month <= END])
if (uniqueN(BASE$cell_id) < 200L) stop('Unexpected loss of cells in analysis period')
if (!all(c('lwj_contested_share','lwj_taliban_control_share',
           'lwj_unconfirmed_share') %in% names(BASE))) stop('Missing LWJ shares')

# The LWJ shares sum to one. GOA is omitted from status time-varying slopes.
# For each calendar month, model A absorbs *different violence baselines* as a
# flexible function of pre-Doha contested / Taliban / unconfirmed shares.
# In these specifications, Post x status is absorbed by month-specific status
# slopes; do not include it again.

TRIP <- c('post_closure:eli_final:lwj_contested_share',
          'post_closure:eli_final:lwj_taliban_control_share',
          'post_closure:eli_final:lwj_unconfirmed_share')
CTRL <- c('prec_gpcp','temp','nlights_calib_mean')

# A: Status-specific month shocks, common Post x ELI retained.
FORM_A <- as.formula(paste(
 'events_type3 ~',
 paste(c('post_closure:eli_final', TRIP, CTRL), collapse=' + '),
 '| cell_id + month[lwj_contested_share,lwj_taliban_control_share,lwj_unconfirmed_share]'))

# B: Additional ELI-specific month shocks, allowing *all* districts with
# higher ELI to experience their own arbitrary monthly changes. The common
# Post x ELI is now absorbed, but differential ELI x status (the H2 contrasts)
# remain as explicit coefficients. This is demanding and may be imprecise.
FORM_B <- as.formula(paste(
 'events_type3 ~', paste(c(TRIP,CTRL), collapse=' + '),
 '| cell_id + month[lwj_contested_share,lwj_taliban_control_share,lwj_unconfirmed_share,eli_final]'))

# C: Comparison, conventional expanded joint model on the EXACT same period.
FORM_C <- DIAG_FORM

WINDOWS <- list(
  'full_2015_2024'=BASE,
  'exclude_anticipation_2015_2024'=BASE[month < as.IDate('2020-03-01') |
                                           month >= as.IDate('2021-07-01')]
)
MODELS <- list('C_expanded_reference'=FORM_C,
               'A_month_by_status'=FORM_A,
               'B_month_by_status_and_ELI'=FORM_B)

get_n <- function(fit) as.integer(stats::nobs(fit))
results <- list(); checks <- list(); failures <- list(); k <- 0L
for (w in names(WINDOWS)) {
  dat <- WINDOWS[[w]]
  cat('\n=== WINDOW:', w, 'cells:', uniqueN(dat$cell_id), 'rows:', nrow(dat), '===\n')
  for (s in names(MODELS)) {
    cat('\nEstimating', s, '...\n')
    fit <- tryCatch(
      fixest::fepois(MODELS[[s]], data=dat, notes=FALSE, warn=FALSE),
      error=function(e) e)
    if (inherits(fit,'error')) {
      msg <- conditionMessage(fit)
      warning(w, '/', s, ' failed: ', msg)
      failures[[length(failures)+1L]] <- data.table(window=w,specification=s,error=msg)
      next
    }
    ct <- tryCatch(fixest::coeftable(fit, vcov=~cell_id+month),
                   error=function(e) e)
    if (inherits(ct,'error')) {
      failures[[length(failures)+1L]] <- data.table(window=w,specification=s,
                                                    error=conditionMessage(ct))
      next
    }
    for (term in TRIP) {
      v <- if (term %in% rownames(ct)) ct[term,] else rep(NA_real_,4L)
      k <- k+1L
      results[[k]] <- data.table(window=w,specification=s,term=term,
          beta=unname(v[1]),se=unname(v[2]),p=unname(v[4]),
          n=get_n(fit),n_cells=uniqueN(dat$cell_id),
          term_estimable=term %in% rownames(ct))
      cat(sprintf('  %s: beta=%.4f SE=%.4f p=%.4g N=%d\n',
                  term, v[1], v[2], v[4], get_n(fit)))
    }
    # Explicit check for dropped regressors (esp. due to collinearity).
    cv <- fit$collin.var
    if (is.null(cv)) cv <- character(0)
    checks[[length(checks)+1L]] <- data.table(window=w,specification=s,
        n=get_n(fit), dropped_terms=paste(cv,collapse='; '))
    saveRDS(fit, file.path(DIR06,paste0(w,'__',s,'.rds')))
  }
}

if (length(results)) {
  R <- rbindlist(results,fill=TRUE)
  fwrite(R,file.path(DIR06,'06_H2_month_specific_results.csv'))
  print(R[grepl('lwj_contested_share$',term),
          .(window,specification,beta,se,p,n,term_estimable)])
  if (any(!R$term_estimable)) warning('Some H2 contrasts not estimable; do not interpret missing terms')
  nn <- unique(R[,.(window,specification,n)])
  for (w in unique(nn$window)) {
    if (uniqueN(nn[window==w]$n)>1L) warning('Different estimation samples across models in window ',w)
  }
} else warning('No models estimated; see failures file')
if (length(checks)) fwrite(rbindlist(checks),file.path(DIR06,'06_collin_and_N_audit.csv'))
if (length(failures)) fwrite(rbindlist(failures),file.path(DIR06,'06_errors.csv'))
cat('\nDONE. Outputs:',DIR06,'\n')
cat('CAUSAL CAUTION: this absorbs group-specific monthly shocks but does not isolate surveillance\n')
cat('from unmeasured shocks whose effects vary jointly with territorial status and ELI.\n')
