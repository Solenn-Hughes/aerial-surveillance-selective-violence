# JPR H2: direct hierarchical sensitivity using the original spatial H2 panel.
# Does not overwrite original data, code, or the manuscript.
# The exact ORIGINAL source is H2 IMPROVEMENT.R / run_h2_model().
# Its reported H2 estimates are -16.102380, -17.666752, -3.120677.

suppressPackageStartupMessages({library(data.table); library(fixest)})
root <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
folder <- file.path(root,"outputs/paper1_geo_afghanistan_executable_acled/diagnostics_h2_spatial_adm2")
input <- file.path(folder,"panel_core_with_spatial_h2_diagnostic.rds")
output <- file.path(root,"outputs/JPR_H2_DIRECT_SENSITIVITY")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
if (!file.exists(input)) stop("Original spatial H2 panel not found at: ",input,
  "\nUse the source panel from H2 IMPROVEMENT.R; do not substitute a FINAL30 panel.")
d <- as.data.table(readRDS(input))
if (!inherits(d$month,"IDate")) d[,month:=as.IDate(month)]
d[,cell_id:=as.factor(cell_id)]
mods <- c("pre_contested_share_spatial", "pre_taliban_activity_share_spatial", "pre_contested_intensity_spatial")
need <- c("events_type3","post_closure","eli_feas","cell_id","month",mods)
missing <- setdiff(need,names(d))
if (length(missing)) stop("Missing original model variable(s): ",paste(missing,collapse=", "))
ctrl <- intersect(c("prec_gpcp","temp","nlights_calib_mean"),names(d))
cat("Original spatial panel:", input,"\n", nrow(d),"rows; controls:",paste(ctrl,collapse=", "),"\n")

# Coefficients and inference are extracted independently of data.table's auto-renaming rules.
get_term <- function(m, term) {
  ct <- fixest::coeftable(m,vcov=~cell_id+month)
  if (!(term %in% rownames(ct))) return(list(est=NA_real_,se=NA_real_,p=NA_real_,error=paste("Term absent",term,
    "Available:",paste(rownames(ct),collapse="; "))))
  i <- which(rownames(ct)==term)[1]
  list(est=unname(ct[i,1]),se=unname(ct[i,2]),p=unname(ct[i,4]),error="")
}
run_one <- function(data,window,moderator,variant){
  eli <- "eli_feas"; post <- "post_closure"
  triple <- paste(post,eli,moderator,sep=":")
  lower <- paste(post,moderator,sep=":")
  rhs <- c(paste(post,eli,sep=":"), if (variant=="full_hierarchy") lower,
           triple,ctrl)
  fm <- as.formula(paste("events_type3 ~",paste(rhs,collapse=" + "),"| cell_id + month"))
  head <- data.table(window=window,moderator=moderator,specification=variant,
                     formula=paste(deparse(fm),collapse=""))
  m <- tryCatch(fixest::fepois(fm,data=data),error=function(e)e)
  if(inherits(m,"error")) {
    cat("MODEL ERROR",window,moderator,variant,":",conditionMessage(m),"\n")
    return(cbind(head,data.table(triple=NA_real_,se=NA_real_,p=NA_real_,
                   post_m=NA_real_,post_m_se=NA_real_,post_m_p=NA_real_,n=NA_integer_,
                   error=conditionMessage(m))))
  }
  a <- get_term(m,triple)
  b <- if(variant=="full_hierarchy") get_term(m,lower) else list(est=NA_real_,se=NA_real_,p=NA_real_,error="")
  cbind(head,data.table(triple=a$est,se=a$se,p=a$p,post_m=b$est,
                         post_m_se=b$se,post_m_p=b$p,n=as.integer(nobs(m)),
                         error=paste(c(a$error,b$error)[nzchar(c(a$error,b$error))],collapse="; ")))
}

windows <- list(full=d, drop_anticipation=d[month<as.IDate("2020-03-01") | month>=as.IDate("2021-07-01")])
results <- list();k <- 0L
for (w in names(windows)) for (mod in mods) for (variant in c("original","full_hierarchy")) {
 k<-k+1L;results[[k]]<-run_one(windows[[w]],w,mod,variant)
}
r <- rbindlist(results,fill=TRUE)
original_expected <- data.table(window=rep(c("full","drop_anticipation"),each=3),
 moderator=rep(mods,2),
 published=c(-16.1023802013192,-17.6667522797993,-3.12067660696125,
             -17.4082078077671,-19.070707767096,-3.38898087940946),
 published_n=c(rep(54312L,3),rep(49956L,3)))
r<-merge(r,original_expected,by=c("window","moderator"),all.x=TRUE,sort=FALSE)
r[,original_reproduced:=ifelse(specification=="original",!is.na(triple)&!is.na(n)&
 abs(triple-published)<0.001 & n==published_n,NA)]
r[,sort_window:=match(window,c("full","drop_anticipation"))]
r[,sort_mod:=match(moderator,mods)]
r[,sort_variant:=match(specification,c("original","full_hierarchy"))]
setorder(r,sort_window,sort_mod,sort_variant)
r[,c("sort_window","sort_mod","sort_variant"):=NULL]
fwrite(r,file.path(output,"H2_original_vs_full_hierarchy.csv"))
print(r[,.(window,moderator,specification,triple,se,p,post_m,post_m_p,n,original_reproduced,error)],row.names=FALSE)
cat("\nSaved: ",file.path(output,"H2_original_vs_full_hierarchy.csv"),"\n",sep="")
if(!all(r[specification=="original",original_reproduced])) {
 cat("WARNING: Original coefficients did not reproduce exactly. Do not replace manuscript numbers.\n")
} else cat("All six original estimates reproduced. Hierarchical sensitivity results above are directly comparable.\n")
