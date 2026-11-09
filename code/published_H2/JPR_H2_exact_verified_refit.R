# JPR H2: exact FINAL30 replication, then lower-order interaction sensitivity.
# Does not overwrite your published source scripts, models, figures, or tables.
# Input: final_analysis_panel.rds from FINAL30, with spatial ADM2 H2 moderators.
# Tested provenance: FINAL30/H2 IMPROVEMENT.R run_h2_model() and supplied results CSVs.

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

THESIS <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
PANEL <- file.path(THESIS, "outputs/final_airfield_screen_20260809_FINAL30/11_final_analysis_panel.rds")
OUT <- file.path(THESIS, "outputs/JPR_H2_EXACT_REPLICATION_AND_HIERARCHY_CHECK")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(PANEL)) stop("Final panel not found: ", PANEL)

D <- as.data.table(readRDS(PANEL))
D[, month := as.IDate(month)]
D[, cell_id := as.factor(cell_id)]

moderators <- c("pre_contested_share_spatial", "pre_taliban_activity_share_spatial", "pre_contested_intensity_spatial")
mandatory <- c("events_type3", "post_closure", "cell_id", "month", moderators)
missing <- setdiff(mandatory, names(D))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse=", "))

# Both ELI variables occur in the FINAL30 RDS. Only the value reproducing ALL
# published spatial H2 coefficients is permitted for the corrected refits.
eli_candidates <- intersect(c("eli_feas", "eli_final"), names(D))
if (!length(eli_candidates)) stop("No ELI column found")
controls <- intersect(c("prec_gpcp", "temp", "nlights_calib_mean"), names(D))
expected <- data.table(
  window = rep(c("full", "drop_anticipation"), each=3),
  moderator = rep(moderators, 2),
  expected_estimate = c(-16.1023802013192, -17.6667522797993, -3.12067660696125,
                        -17.4082078077671, -19.0707077670960, -3.38898087940946),
  expected_se = c(4.27319350902624, 4.6031800283220, 0.880878277343997,
                  4.86739390345531, 5.26815958334537, 0.984560054320593),
  expected_n = c(rep(54312L, 3), rep(49956L, 3))
)

extract_term <- function(m, term) {
  ct <- as.data.table(fixest::coeftable(m, vcov=~cell_id+month), keep.rownames="term")
  if (!term %in% ct$term) return(data.table(estimate=NA_real_, se=NA_real_, p=NA_real_))
  row <- ct[ct[["term"]] == term]
  data.table(estimate=as.numeric(row[["Estimate"]][1]),
             se=as.numeric(row[["Std. Error"]][1]),
             p=as.numeric(row[names(row)[grepl("^Pr\\(", names(row))][1]][1]))
}

fit_once <- function(dt, eli, mod, add_lower_order=FALSE) {
  # Exact original specification in H2 IMPROVEMENT.R (run_h2_model):
  # events_type3 ~ post_closure:ELI + post_closure:ELI:M + controls | cell_id + month
  # Fully specified sensitivity: additionally includes post_closure:M.
  rhs <- c(paste0("post_closure:",eli),
           paste0("post_closure:",eli,":",mod),
           if(add_lower_order) paste0("post_closure:",mod),
           controls)
  fml <- as.formula(paste0("events_type3 ~ ",paste(rhs,collapse=" + ")," | cell_id + month"))
  fit <- fixest::fepois(fml, data=dt, notes=FALSE)
  key <- paste0("post_closure:",eli,":",mod)
  est <- extract_term(fit, key)
  lower <- if(add_lower_order) extract_term(fit, paste0("post_closure:",mod)) else data.table(estimate=NA_real_,se=NA_real_,p=NA_real_)
  data.table(estimate=est$estimate, se=est$se, p=est$p,
             lower_order_estimate=lower$estimate, lower_order_se=lower$se,
             lower_order_p=lower$p, nobs=as.integer(nobs(fit)),
             collinear_terms=paste(fit$collin.var, collapse="; "))
}

windows <- list(
  full=D,
  drop_anticipation=D[month < as.IDate("2020-03-01") | month >= as.IDate("2021-07-01")]
)

cat("FINAL30 panel rows:",nrow(D),"\n")
cat("Controls:",paste(controls, collapse=", "),"\n")
cat("Checking all 6 published H2 models before running corrected variants...\n")
checks <- list(); k <- 0L
for(eli in eli_candidates) for(w in names(windows)) for(mod in moderators) {
  k <- k+1L
  a <- tryCatch(fit_once(windows[[w]],eli,mod,FALSE),error=function(e) {
    data.table(estimate=NA_real_,se=NA_real_,p=NA_real_,lower_order_estimate=NA_real_,
               lower_order_se=NA_real_,lower_order_p=NA_real_,nobs=NA_integer_,
               collinear_terms=paste("ERROR",conditionMessage(e)))
  })
  checks[[k]] <- cbind(data.table(eli=eli,window=w,moderator=mod),a)
}
check_dt <- rbindlist(checks,fill=TRUE)
check_dt <- merge(check_dt,expected,by=c("window","moderator"),all.x=TRUE,sort=FALSE)
check_dt[, reproduces := !is.na(estimate) & !is.na(se) & !is.na(nobs) &
            abs(estimate-expected_estimate)<0.002 & abs(se-expected_se)<0.002 & nobs==expected_n]
fwrite(check_dt,file.path(OUT,"01_REPRODUCTION_OF_PUBLISHED_H2.csv"))
print(check_dt[,.(eli,window,moderator,estimate,se,nobs,reproduces)],row.names=FALSE)

eli_ok <- check_dt[,.(matches=sum(reproduces),total=.N),by=eli][matches==6 & total==6,eli]
if(!length(eli_ok)) {
  cat("\nSTOP: Could not reproduce all 6 published H2 models.\n",
      "No interpretation of the corrected specification is valid yet.\n",
      "Inspect the reproduction CSV and the exact source data lineage.\n",sep="")
  stop("Original H2 replication not verified. Corrected models were NOT run.")
}
eli <- if ("eli_feas" %in% eli_ok) "eli_feas" else eli_ok[1]
if (length(eli_ok)>1L) cat("Multiple ELI columns reproduce the results; using the original H2 source variable: ",eli,"\n",sep="")
cat("\nAll six original H2 models VERIFIED with ",eli,".\n",sep="")

comparisons <- list(); j <- 0L
for(w in names(windows)) for(mod in moderators) {
  original <- check_dt[check_dt[["eli"]] == eli & window == w & moderator == mod][1]
  corrected <- fit_once(windows[[w]],eli,mod,TRUE)
  j <- j+1L
  comparisons[[j]] <- data.table(
    window=w, moderator=mod, eli=eli,
    original_estimate=original$estimate, original_se=original$se,
    original_p=original$p, original_nobs=original$nobs,
    corrected_triple_estimate=corrected$estimate,
    corrected_triple_se=corrected$se,
    corrected_triple_p=corrected$p,
    corrected_lower_order_estimate=corrected$lower_order_estimate,
    corrected_lower_order_se=corrected$lower_order_se,
    corrected_lower_order_p=corrected$lower_order_p,
    corrected_nobs=corrected$nobs,
    same_sample_n=(corrected$nobs==original$nobs),
    collinear_terms=corrected$collinear_terms
  )
}
results <- rbindlist(comparisons, fill=TRUE)
fwrite(results,file.path(OUT,"02_EXACT_PAIRED_H2_RESULTS.csv"))
print(results[,.(window,moderator,original_estimate,original_p,
                 corrected_triple_estimate,corrected_triple_se,
                 corrected_triple_p,corrected_lower_order_p,
                 original_nobs,corrected_nobs)],row.names=FALSE)
cat("\nDone. Original results preserved.\nResults folder: ",OUT,"\n",sep="")
