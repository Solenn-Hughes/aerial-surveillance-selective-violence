# ============================================================
# R09y_pretrend_diagnostics_m0boot_v2.R
# Pretrend diagnostics for m0_boot ES (fixes Wald NAs)
# - ES baseline ref=-1
# - ES alt ref=-6:-2
# - ES alt ref=-6:-2 + donut(-1,0,1)
#
# Computes:
# - Joint Wald test on leads k in [-24,-2] excluding ref months
#   (manual Wald with PSD-stabilization)
# - Linear slope test on lead coefficients (weighted by 1/SE^2)
#
# Outputs audits to outputs/paper1/audit/
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

# -------------------------
# 0) Config + helpers
# -------------------------
cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found.")
source(cfg_path)

if (!exists("write_audit_dt")) {
  write_audit_dt <- function(dt, stem) {
    out_dir <- file.path("outputs","paper1","audit")
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    f <- file.path(out_dir, paste0(stem, "_", ts, ".csv"))
    fwrite(as.data.table(dt), f)
    message("Wrote audit: ", f)
    invisible(f)
  }
}

assert_file <- function(p) if (!file.exists(p)) stop("INTEGRITY FAIL: Missing file: ", p)
assert_cols2 <- function(dt, cols, nm="dt") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("INTEGRITY FAIL: ", nm, " missing cols: ", paste(miss, collapse=", "))
}

# -------------------------
# 1) Load panel
# -------------------------
P_PANEL <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds")
assert_file(P_PANEL)
dt <- readRDS(P_PANEL); setDT(dt)

if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
if (exists("as_month_idate")) dt[, month := as_month_idate(month)] else dt[, month := as.IDate(month)]
dt[, month_int := year(month)*12L + month(month)]
dt[, theatre_id := as.character(theatre_id)]

req <- c("pg_id","theatre_id","month_int","treated","event_time","eli_feas","events_type3")
assert_cols2(dt, req, "ALLCW_posteli_v2")

# m0_boot controls
CTRL_M0 <- c("nlights_calib_mean","nlights_mean","nlights_min","prec_gpcp","temp","water_ih")
missing_ctrl <- setdiff(CTRL_M0, names(dt))
if (length(missing_ctrl)) stop("INTEGRITY FAIL: missing m0_boot controls: ", paste(missing_ctrl, collapse=", "))

# -------------------------
# 2) Prepare ES variables (never-treated retained)
# -------------------------
window <- 60L
d <- copy(dt)

d[, ylog := log1p(events_type3)]
d[, event_time_es := as.integer(event_time)]
d[treated == 0L, event_time_es := 0L]
d[treated == 1L & (event_time_es < -24L | event_time_es > window), event_time_es := NA_integer_]
d[, eli_feas_t := fifelse(treated == 1L, eli_feas, 0)]

# -------------------------
# 3) Helpers: extract ES coefs and compute tests
# -------------------------
extract_es <- function(m) {
  ct <- as.data.table(coeftable(m), keep.rownames="term")
  es <- ct[grepl("^event_time_es::\\-?\\d+:eli_feas_t$", term)]
  if (!nrow(es)) stop("INTEGRITY FAIL: no ES terms matched event_time_es::<k>:eli_feas_t")
  es[, k := as.integer(sub("^event_time_es::(-?\\d+):eli_feas_t$", "\\1", term))]
  setorder(es, k)
  es
}

# PSD-stabilized inverse solve for Wald
wald_leads_manual <- function(m, ref_vec, lead_min=-24L, lead_max=-2L, ridge=1e-10) {
  es <- extract_es(m)
  leads <- es[k >= lead_min & k <= lead_max & !(k %in% ref_vec)]
  terms <- leads$term
  
  if (length(terms) < 5) {
    return(data.table(n_leads=length(terms), wald_stat=NA_real_, wald_df=length(terms), wald_p=NA_real_, vcov_min_eig=NA_real_, wald_note="Too few lead terms"))
  }
  
  b <- coef(m)[terms]
  V <- vcov(m)[terms, terms, drop=FALSE]
  
  # Symmetrize
  V <- (V + t(V))/2
  
  # Eigen diagnostics + ridge to ensure PSD
  eig <- tryCatch(eigen(V, symmetric=TRUE, only.values=TRUE)$values, error=function(e) NA_real_)
  min_eig <- if (all(is.na(eig))) NA_real_ else min(eig)
  
  V_adj <- V
  if (!is.na(min_eig) && min_eig < ridge) {
    V_adj <- V + diag(abs(min_eig) + ridge, nrow(V))
  }
  
  # Wald statistic
  stat <- tryCatch(as.numeric(t(b) %*% solve(V_adj) %*% b), error=function(e) NA_real_)
  df <- length(b)
  p  <- if (is.na(stat)) NA_real_ else 1 - pchisq(stat, df=df)
  
  data.table(
    n_leads = df,
    wald_stat = stat,
    wald_df = df,
    wald_p = p,
    vcov_min_eig = min_eig,
    wald_note = if (!is.na(min_eig) && min_eig < ridge) "VCOV ridge-adjusted" else ""
  )
}

slope_test_leads <- function(es, ref_vec, lead_min=-24L, lead_max=-2L) {
  leads <- es[k >= lead_min & k <= lead_max & !(k %in% ref_vec)]
  leads <- leads[!is.na(`Std. Error`) & `Std. Error` > 0]
  if (nrow(leads) < 5) return(data.table(slope=NA_real_, slope_se=NA_real_, slope_p=NA_real_, n=nrow(leads)))
  
  w <- 1/(leads$`Std. Error`^2)
  fit <- lm(Estimate ~ k, data=leads, weights=w)
  cf <- summary(fit)$coefficients
  
  data.table(
    slope = as.numeric(cf["k","Estimate"]),
    slope_se = as.numeric(cf["k","Std. Error"]),
    slope_p = as.numeric(cf["k","Pr(>|t|)"]),
    n = nrow(leads)
  )
}

fit_es <- function(ref_vec, donut=FALSE, label="") {
  x <- copy(d)
  if (donut) x[treated==1L & event_time_es %in% c(-1L, 0L, 1L), event_time_es := NA_integer_]
  
  rhs <- paste(
    c(sprintf("i(event_time_es, eli_feas_t, ref = c(%s))", paste(ref_vec, collapse=",")),
      CTRL_M0),
    collapse=" + "
  )
  
  fml <- as.formula(paste0("ylog ~ ", rhs, " | pg_id + theatre_id^month_int"))
  m <- feols(fml, data=x, vcov=~theatre_id)
  
  es <- extract_es(m)
  w  <- wald_leads_manual(m, ref_vec=ref_vec)
  sl <- slope_test_leads(es, ref_vec=ref_vec)
  
  diag <- cbind(
    data.table(spec=label, ref=paste(ref_vec, collapse=","), donut=donut),
    w, sl,
    data.table(nobs=nobs(m))
  )
  
  list(model=m, es=es, diag=diag)
}

# -------------------------
# 4) Run three specs
# -------------------------
res1 <- fit_es(ref_vec=c(-1L), donut=FALSE, label="ES baseline ref=-1")
res2 <- fit_es(ref_vec=as.integer(-6:-2), donut=FALSE, label="ES alt ref=-6:-2")
res3 <- fit_es(ref_vec=as.integer(-6:-2), donut=TRUE, label="ES alt ref=-6:-2 + donut(-1,0,1)")

diag <- rbindlist(list(res1$diag, res2$diag, res3$diag), fill=TRUE)
print(diag)
write_audit_dt(diag, "R09y_pretrend_diagnostics_m0boot_v2")

es_all <- rbindlist(list(
  cbind(res1$es, spec="ES baseline ref=-1"),
  cbind(res2$es, spec="ES alt ref=-6:-2"),
  cbind(res3$es, spec="ES alt ref=-6:-2 + donut")
), fill=TRUE)
write_audit_dt(es_all, "R09y_pretrend_diagnostics_m0boot_v2_es_coefs")

cat("\n✅ Done. The printed diag table now includes valid Wald p-values.\n")