# JPR H2 forensic repair: reproduce FINAL30 six coefficients, then test Post x moderator
# Does NOT modify any original scripts, data or manuscript.
# Based on actual final results (07_H2_anticipation_window_original_spec.csv):
# all H2 terms use post_closure:eli_final:<spatial_moderator>.
# Errors are written verbatim (never silently converted to NA).

suppressPackageStartupMessages({ library(data.table); library(fixest) })
ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
BASE <- file.path(ROOT, "outputs/final_airfield_screen_20260809_FINAL30")
OUT  <- file.path(ROOT, "outputs/JPR_H2_FINAL30_REPAIRED_REPLICATION")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

panel_options <- c(
  file.path(BASE, "11_final_analysis_panel.rds"),
  file.path(ROOT, "outputs/paper1_geo_afghanistan_executable_acled/diagnostics_h2_spatial_adm2/panel_core_with_spatial_h2_diagnostic.rds")
)
mods <- c("pre_contested_share_spatial", "pre_taliban_activity_share_spatial", "pre_contested_intensity_spatial")
expected <- data.table(
  window=rep(c("full", "drop_anticipation"),each=3),
  moderator=rep(mods,2),
  expected_estimate=c(-16.1023802013192,-17.6667522797993,-3.12067660696125,-17.4082078077671,-19.070707767096,-3.38898087940946),
  expected_se=c(4.27319350902624,4.6031800283220,0.880878277343997,4.86739390345531,5.26815958334537,0.984560054320593),
  expected_n=c(rep(54312L,3),rep(49956L,3))
)

extract_coef <- function(m, term) {
  ct <- data.table::as.data.table(fixest::coeftable(m, vcov=~cell_id+month), keep.rownames="term")
  if (!(term %in% ct$term)) stop("Term not found: ",term, "; available terms: ",paste(ct$term,collapse="; "))
  row <- ct[ct$term==term,]
  pcol <- grep("^Pr\\(", names(row), value=TRUE)
  if (length(pcol)!=1L) stop("Could not identify p-value column: ",paste(names(row),collapse=", "))
  list(estimate=as.numeric(row[["Estimate"]][1]),se=as.numeric(row[["Std. Error"]][1]),p=as.numeric(row[[pcol]][1]))
}

fit_model <- function(dt, eli, mod, corrected=FALSE) {
  controls <- intersect(c("prec_gpcp","temp","nlights_calib_mean"),names(dt))
  rhs <- c(paste0("post_closure:",eli),
           paste0("post_closure:",eli,":",mod),
           if(corrected) paste0("post_closure:",mod),controls)
  form <- stats::as.formula(paste("events_type3 ~",paste(rhs,collapse=" + "),"| cell_id + month"))
  # The ORIGINAL H2 code calls fepois(fml, data=dt), with no 'notes' parameter.
  model <- fixest::fepois(form, data=dt)
  triple <- extract_coef(model, paste0("post_closure:",eli,":",mod))
  lower <- if (corrected) extract_coef(model,paste0("post_closure:",mod)) else list(estimate=NA_real_,se=NA_real_,p=NA_real_)
  data.table(estimate=triple$estimate,se=triple$se,p=triple$p,
             lower_estimate=lower$estimate,lower_se=lower$se,lower_p=lower$p,
             nobs=as.integer(stats::nobs(model)),formula=paste(deparse(form),collapse=""))
}

attempts <- list(); paired <- list(); ka <- 0L; kp <- 0L
for (path in panel_options) {
  if (!file.exists(path)) {cat("Skipping absent panel: ",path,"\n",sep=""); next}
  cat("\nPANEL: ",path,"\n",sep="")
  D <- data.table::as.data.table(readRDS(path))
  D[,month:=data.table::as.IDate(month)]
  D[,cell_id:=as.factor(cell_id)]
  cat("Rows: ",nrow(D),"; columns: ",ncol(D),"\n",sep="")
  needed <- c("events_type3","post_closure","cell_id","month",mods)
  miss <- setdiff(needed,names(D))
  if(length(miss)){cat("Missing columns: ",paste(miss,collapse=", "),"\n",sep="");next}
  # FINAL30 published output explicitly calls this 'eli_final', unlike the earlier diagnostic code.
  eli_candidates <- intersect(c("eli_final","eli_feas"),names(D))
  if(!length(eli_candidates)){cat("Missing ELI\n");next}
  windows <- list(full=D,drop_anticipation=D[month < data.table::as.IDate("2020-03-01") | month >= data.table::as.IDate("2021-07-01")])
  for(eli in eli_candidates){
    local <- list(); k <- 0L
    for(w in names(windows)) for(mod in mods){
      k <- k+1L
      a <- tryCatch(fit_model(windows[[w]],eli,mod,FALSE),error=function(e) {
        msg <- conditionMessage(e)
        cat("ERROR [",basename(path),", ",eli,", ",w,", ",mod,"]: ",msg,"\n",sep="")
        data.table(estimate=NA_real_,se=NA_real_,p=NA_real_,lower_estimate=NA_real_,lower_se=NA_real_,lower_p=NA_real_,nobs=NA_integer_,formula="")
      })
      local[[k]] <- cbind(data.table(panel=path,eli=eli,window=w,moderator=mod),a)
    }
    x <- data.table::rbindlist(local,fill=TRUE)
    x <- merge(x,expected,by=c("window","moderator"),all.x=TRUE,sort=FALSE)
    x[,reproduces:=!is.na(estimate)&!is.na(se)&!is.na(nobs)&
      abs(estimate-expected_estimate)<0.002 & abs(se-expected_se)<0.002 & nobs==expected_n]
    ka<-ka+1L; attempts[[ka]]<-x
    cat("Replicated ",sum(x$reproduces)," / 6 using ",eli,"\n",sep="")
    if(all(x$reproduces)){
      cat("\nVERIFIED: six FINAL30 H2 models replicated. Running fully specified sensitivity...\n")
      for(w in names(windows)) for(mod in mods){
        original<-x[x$window==w & x$moderator==mod,][1]
        correction<-tryCatch(fit_model(windows[[w]],eli,mod,TRUE),error=function(e){
          cat("CORRECTED FIT ERROR [",w,",",mod,"]: ",conditionMessage(e),"\n",sep="")
          NULL
        })
        if(is.null(correction)) next
        kp<-kp+1L
        paired[[kp]]<-data.table(panel=path,eli=eli,window=w,moderator=mod,
          original_triple=original$estimate,original_se=original$se,original_p=original$p,original_n=original$nobs,
          corrected_triple=correction$estimate,corrected_se=correction$se,corrected_p=correction$p,
          corrected_post_m=correction$lower_estimate,corrected_post_m_se=correction$lower_se,
          corrected_post_m_p=correction$lower_p,corrected_n=correction$nobs,
          same_n=(original$nobs==correction$nobs))
      }
      break
    }
  }
  if(length(paired)) break
}

if(!length(attempts)) stop("No compatible panel was available; see printed paths.")
all_checks <- data.table::rbindlist(attempts,fill=TRUE)
data.table::fwrite(all_checks,file.path(OUT,"01_ALL_REPRODUCTION_ATTEMPTS.csv"))
if(length(paired)){
  ans<-data.table::rbindlist(paired,fill=TRUE)
  data.table::fwrite(ans,file.path(OUT,"02_CORRECTED_VS_ORIGINAL_H2.csv"))
  print(ans[,.(window,moderator,original_triple,corrected_triple,corrected_se,corrected_p,corrected_post_m_p,original_n,corrected_n,same_n)],row.names=FALSE)
  cat("\nFINISHED. Upload 02_CORRECTED_VS_ORIGINAL_H2.csv from: ",OUT,"\n",sep="")
} else {
  cat("\nNo panel/ELI combination reproduced all six FINAL30 estimates.\n")
  cat("The first ERROR message above is the actual failure to diagnose; no coefficients have been changed.\n")
  cat("Upload 01_ALL_REPRODUCTION_ATTEMPTS.csv and copy the first ERROR line.\n")
}
