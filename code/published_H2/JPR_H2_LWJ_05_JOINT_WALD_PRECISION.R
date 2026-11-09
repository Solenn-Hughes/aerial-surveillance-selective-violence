# 05: joint Wald tests of all three differential ELI moderation coefficients.
# Tests are descriptive, not a solution to causal identification.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
OUT5 <- file.path(OUT,'05_joint_wald_precision')
dir.create(OUT5,recursive=TRUE,showWarnings=FALSE)
R <- list(); E <- list()
for (w in names(DIAG_WINDOWS)) {
  f <- diag_fit(DIAG_WINDOWS[[w]])$fit
  b <- coef(f)
  V <- vcov(f,vcov=~cell_id+month)
  terms <- DIAG_TRIPLES
  if(!all(terms %in% names(b))) stop('Joint moderation terms missing')
  vv <- V[terms,terms,drop=FALSE]
  bb <- b[terms]
  eig <- eigen(vv,symmetric=TRUE,only.values=TRUE)$values
  if(any(!is.finite(eig)) || min(eig)<=0) {
    warning('Joint covariance not positive definite in ',w,'; Wald test not calculated')
    p <- NA_real_;stat <- NA_real_
  } else {
    stat <- as.numeric(t(bb)%*%solve(vv,bb))
    p <- pchisq(stat,df=length(terms),lower.tail=FALSE)
  }
  R[[length(R)+1L]] <- data.table(window=w, test='All status-specific differential ELI moderation terms are zero',
                                 statistic=stat,df=length(terms),p=p, n=nobs(f),
                                 min_cov_eigenvalue=min(eig),max_cov_eigenvalue=max(eig))
  C <- cov2cor(vv)
  C <- as.data.table(C,keep.rownames='term')
  C[,window:=w]
  E[[length(E)+1L]] <- C
}
print(rbindlist(R))
fwrite(rbindlist(R),file.path(OUT5,'01_joint_wald.csv'))
fwrite(rbindlist(E,fill=TRUE),file.path(OUT5,'02_cluster_covariance_correlations.csv'))
cat('Joint Wald files: ',OUT5,'\n')
