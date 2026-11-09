# JPR H2: why does adding Post x Moderator reduce precision?
# Uses ONLY the verified published FINAL30 panel and eli_final.
# Does not alter the original analyses or select a new preferred specification.
# Run: source("C:/Users/USERNAME/Downloads/JPR_H2_PRECISION_DIAGNOSTIC_CORRECTED.R")

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
PANEL <- file.path(ROOT, "outputs", "final_airfield_screen_20260809_FINAL30",
                   "11_final_analysis_panel.rds")
OUT <- file.path(ROOT, "outputs", "JPR_H2_PRECISION_DIAGNOSTIC")
if (!file.exists(PANEL)) stop("FINAL30 panel not found: ", PANEL)
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

DT <- as.data.table(readRDS(PANEL))
mods <- c("pre_contested_share_spatial", "pre_taliban_activity_share_spatial",
          "pre_contested_intensity_spatial")
controls <- c("prec_gpcp", "temp", "nlights_calib_mean")
required <- c("events_type3", "post_closure", "eli_final", "cell_id", "month", mods, controls)
missing_cols <- setdiff(required, names(DT))
if (length(missing_cols)) stop("Missing expected FINAL30 columns: ", paste(missing_cols,collapse=", "))
DT[, month := as.IDate(month)]
DT[, cell_id := as.factor(cell_id)]

# Source-verified benchmark values: original three-way coefficients, SE and N.
expected <- data.table(
  window=rep(c("full","drop_anticipation"), each=3),
  moderator=rep(mods,2),
  coefficient=c(-16.1023802013192,-17.6667522797993,-3.12067660696125,
                -17.4082078077671,-19.0707077670960,-3.38898087940946),
  se=c(4.27319350902624,4.60318002832200,0.880878277343997,
       4.86739390345531,5.26815958334537,0.984560054320593),
  n=rep(c(54312L,49956L), each=3))

coef_record <- function(fit, term) {
  z <- coeftable(fit, vcov=~cell_id+month)
  if (!term %in% rownames(z)) stop("Term not estimated: ",term,
                                  "; available: ",paste(rownames(z),collapse=", "))
  c(estimate=unname(z[term,1]),se=unname(z[term,2]),p=unname(z[term,4]))
}

all_main <- list()
all_overlap <- list()
all_ranges <- list()
k <- 0L

windows <- list(
  full=DT,
  drop_anticipation=DT[month < as.IDate("2020-03-01") |
                        month >= as.IDate("2021-07-01")]
)

for (window in names(windows)) for (mod in mods) {
  cat("\n=== ",window," / ",mod," ===\n",sep="")
  D <- copy(windows[[window]])
  # Explicit generated regressors prevent parsing ambiguities.
  D[, x_post_eli := post_closure * eli_final]
  D[, x_post_m := post_closure * get(mod)]
  D[, x_triple := post_closure * eli_final * get(mod)]
  D[, x_m := get(mod)]

  # Match the EXACT formula syntax used by JPR_H2_PROPER_RETEST.R,
  # which reproduced FINAL30 estimates in all six baseline models.
  term_eli <- paste("post_closure", "eli_final", sep=":")
  term_m <- paste("post_closure", mod, sep=":")
  term_triple <- paste("post_closure", "eli_final", mod, sep=":")
  rhs_base <- paste(c(term_eli,term_triple,controls), collapse=" + ")
  rhs_full <- paste(c(term_eli,term_m,term_triple,controls),collapse=" + ")
  f_base <- as.formula(paste("events_type3 ~",rhs_base,"| cell_id + month"))
  f_full <- as.formula(paste("events_type3 ~",rhs_full,"| cell_id + month"))
  base <- fepois(f_base,data=D,notes=FALSE)
  full <- fepois(f_full,data=D,notes=FALSE)
  id_base <- as.integer(obs(base))
  id_full <- as.integer(obs(full))
  same_sample <- identical(sort(id_base),sort(id_full))
  if (!same_sample) stop("Different estimation observations for ",mod," / ",window)
  b0 <- coef_record(base,term_triple)
  b1 <- coef_record(full,term_triple)
  gamma <- coef_record(full,term_m)
  current_window <- window
  current_mod <- mod
  benchmark <- expected[window == current_window & moderator == current_mod]
  verified <- nrow(benchmark)==1L && nobs(base)==benchmark$n[[1L]] &&
              abs(b0["estimate"]-benchmark$coefficient)<0.001 &&
              abs(b0["se"]-benchmark$se)<0.001
  if (!verified) stop(
    "Published H2 NOT reproduced for ",window," / ",mod,".\n",
    "  Estimated beta=",signif(b0["estimate"],10),
    ", SE=",signif(b0["se"],10),", N=",nobs(base),"\n",
    "  Expected beta=",signif(benchmark$coefficient,10),
    ", SE=",signif(benchmark$se,10),", N=",benchmark$n,"\n",
    "  Formula: ",paste(deparse(f_base),collapse=""),"\n",
    "  No diagnostic estimates should be interpreted until reproduced.")

  # Weight/residualization diagnostics on EXACT expanded-model rows. The
  # weights are the fitted PPML conditional means, reflecting the information
  # weighting in the Poisson model. These are DESCRIPTIVE collinearity checks,
  # not separate causal models or hypothesis tests.
  S <- copy(D[id_full])
  S[, ppml_weight := as.numeric(predict(full, type="response"))]
  if (anyNA(S$ppml_weight) || any(S$ppml_weight <= 0))
    stop("Nonpositive/missing fitted PPML weights")

  residual_formula <- function(target) as.formula(paste(
    target, "~", paste(c("x_post_eli",controls),collapse=" + "), "| cell_id + month"))
  residualize <- function(target) {
    m <- feols(residual_formula(target), data=S, weights=~ppml_weight,
               notes=FALSE, warn=FALSE)
    as.numeric(resid(m))
  }
  r_triple <- residualize("x_triple")
  r_post_m <- residualize("x_post_m")
  weighted_cor <- function(x,y,w) {
    mx <- weighted.mean(x,w); my <- weighted.mean(y,w)
    num <- sum(w*(x-mx)*(y-my))
    den <- sqrt(sum(w*(x-mx)^2)*sum(w*(y-my)^2))
    if (den==0) return(NA_real_)
    num/den
  }
  corr_adjusted <- weighted_cor(r_triple,r_post_m,S$ppml_weight)
  # When the correlation is nearly +/-1 the independent three-way signal
  # is tiny, so the variance inflation proxy becomes large.
  vif_proxy <- if (is.finite(corr_adjusted) && abs(corr_adjusted)<1)
    1/(1-corr_adjusted^2) else NA_real_
  corr_raw <- weighted_cor(S$x_triple,S$x_post_m,S$ppml_weight)

  # The robust covariance is useful because it reports correlation in the
  # estimated two-parameter sampling uncertainty after all adjustments.
  V <- vcov(full, vcov=~cell_id+month)
  idx <- c(term_m,term_triple)
  corr_coef <- V[idx[1],idx[2]]/sqrt(V[idx[1],idx[1]]*V[idx[2],idx[2]])

  k <- k+1L
  all_main[[k]] <- data.table(
    window=window,moderator=mod,n=nobs(full),published_reproduced=verified,
    original_beta=b0["estimate"],original_se=b0["se"],original_p=b0["p"],
    expanded_beta=b1["estimate"],expanded_se=b1["se"],expanded_p=b1["p"],
    expanded_post_m=gamma["estimate"],expanded_post_m_se=gamma["se"],
    expanded_post_m_p=gamma["p"],
    se_inflation=unname(b1["se"]/b0["se"]),
    weighted_raw_corr=unname(corr_raw),
    weighted_partial_corr=unname(corr_adjusted),
    vif_proxy_from_partial_corr=unname(vif_proxy),
    coefficient_error_corr=unname(corr_coef),
    same_estimation_rows=same_sample)
  all_overlap[[k]] <- data.table(window=window,moderator=mod,
    metric=c("weighted_raw_corr_post_m_vs_triple",
             "weighted_partial_corr_after_post_eli_controls_cell_month_FE",
             "VIF_proxy_1_over_1_minus_partial_corr_squared",
             "cluster_robust_coefficient_error_corr"),
    value=c(corr_raw,corr_adjusted,vif_proxy,corr_coef))

  # Actual overlap between ELI and moderator in estimation sample. Restrict
  # to distinct cells: time-invariant exposure should not count 200 times.
  cells <- unique(S[,.(cell_id,eli_final,x_m)])
  e <- as.numeric(quantile(cells$eli_final,c(.05,.25,.5,.75,.95),na.rm=TRUE))
  mm <- as.numeric(quantile(cells$x_m,c(.05,.25,.5,.75,.95),na.rm=TRUE))
  all_ranges[[k]] <- data.table(window=window,moderator=mod,
    n_unique_cells=nrow(cells),
    eli_p05=e[1],eli_p25=e[2],eli_median=e[3],eli_p75=e[4],eli_p95=e[5],
    m_p05=mm[1],m_p25=mm[2],m_median=mm[3],m_p75=mm[4],m_p95=mm[5],
    corr_eli_m=if (sd(cells$eli_final)>0 && sd(cells$x_m)>0)
      cor(cells$eli_final,cells$x_m) else NA_real_)

  cat("Published reproduced:",verified,"; N:",nobs(full),"\n")
  cat("Three-way SE: original",round(b0["se"],3),
      "; expanded",round(b1["se"],3),
      "; inflation",round(b1["se"]/b0["se"],2),"x\n")
  cat("Weighted partial correlation (Post x M vs triple):",
      round(corr_adjusted,3),"; VIF proxy:",round(vif_proxy,2),"\n")
  rm(D,S,base,full)
  gc(verbose=FALSE)
}
main <- rbindlist(all_main)
fwrite(main,file.path(OUT,"01_PRECISION_AND_COLLINEARITY.csv"))
fwrite(rbindlist(all_overlap),file.path(OUT,"02_OVERLAP_DETAILS.csv"))
fwrite(rbindlist(all_ranges),file.path(OUT,"03_ACTUAL_ELI_MODERATOR_RANGES.csv"))
cat("\n=== RESULTS (all from verified FINAL30) ===\n")
print(main[,.(window,moderator,original_beta,original_se,
              expanded_beta,expanded_se,expanded_p,se_inflation,
              weighted_partial_corr,vif_proxy_from_partial_corr)],digits=4)
cat("\nSaved to:",OUT,"\n")
cat("Interpretation: large absolute partial correlation / VIF suggests\n",
    "insufficient independent variation to separate the two interactions.\n",
    "Small partial correlation with large SE increase suggests looking beyond\n",
    "simple overlap (e.g. sparse outcome or weak identifying variation).\n",
    "Neither result establishes whether H2's mechanism is true or false.\n",sep="")
