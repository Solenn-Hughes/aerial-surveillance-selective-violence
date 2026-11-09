# ============================================================
# R09c_postcheck_models.R — Post-run sanity checks for Step 9 v2
# - Loads latest Step 9 v2 models RDS
# - Recomputes lead-joint Wald tests robustly
# - Reports ES sample sizes + cluster counts
# - Writes a clean checkpoint audit CSV
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))
suppressMessages(library(fixest))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found: R00_config.R or R/00_config.R")
source(cfg_path)

latest_file <- function(dir, pattern) {
  ff <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(ff)) stop("No files matched in: ", dir, " pattern=", pattern)
  ff[which.max(file.info(ff)$mtime)]
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- robustly extract (stat, df, p) from fixest::wald() across versions
as_wald_fields <- function(w) {
  # If it's a matrix/data.frame
  if (is.matrix(w) || is.data.frame(w)) {
    wd <- as.data.frame(w)
    nms <- tolower(names(wd))
    pick <- function(cands) {
      id <- which(nms %in% cands)
      if (!length(id)) return(NULL)
      wd[[id[1]]][1]
    }
    stat <- pick(c("stat","statistic","f","chi2","chisq","wald"))
    pval <- pick(c("p","p.value","pval","p_value","pvalue","pr(>f)","pr(>chisq)"))
    df1  <- pick(c("df","df1","df_num","dof","parameter"))
    return(list(stat=as.numeric(stat), df=as.integer(df1), p=as.numeric(pval)))
  }
  
  # If it's a list (common)
  if (is.list(w)) {
    stat <- w$stat %||% w$statistic %||% w$Fstat %||% w$chi2 %||% w$chisq
    pval <- w$p.value %||% w$p %||% w$pval %||% w$`Pr(>F)` %||% w$`Pr(>Chisq)`
    df1  <- w$df %||% w$df1 %||% w$parameter
    return(list(stat=as.numeric(stat), df=as.integer(df1), p=as.numeric(pval)))
  }
  
  # Fallback: unknown structure
  list(stat=NA_real_, df=NA_integer_, p=NA_real_)
}

wald_leads_joint <- function(m, lead_min=-24L, lead_max=-2L) {
  cn <- names(coef(m))
  if (is.null(cn) || !length(cn)) return(list(n_leads=0L, stat=NA_real_, df=NA_integer_, p=NA_real_, note="no coef"))
  et_terms <- cn[grepl("^event_time_es::", cn)]
  if (!length(et_terms)) return(list(n_leads=0L, stat=NA_real_, df=NA_integer_, p=NA_real_, note="no event_time terms"))
  
  k <- suppressWarnings(as.integer(sub("^event_time_es::\\(?(-?\\d+)\\)?.*$", "\\1", et_terms)))
  keep <- et_terms[!is.na(k) & k >= lead_min & k <= lead_max]
  keep <- intersect(keep, cn)
  if (!length(keep)) return(list(n_leads=0L, stat=NA_real_, df=NA_integer_, p=NA_real_, note="no leads kept"))
  
  w <- tryCatch(wald(m, keep), error = function(e) e)
  if (inherits(w, "error")) return(list(n_leads=length(keep), stat=NA_real_, df=NA_integer_, p=NA_real_, note=paste0("wald error: ", w$message)))
  
  wf <- as_wald_fields(w)
  list(n_leads=length(keep), stat=wf$stat, df=wf$df, p=wf$p, note="ok")
}

parse_es_key <- function(key) {
  # expected: outcome_eventstudy_feols_log1p__FE__WINDOW
  # e.g. events_ged_eventstudy_feols_log1p__global_month__w24
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  out <- list(outcome=NA_character_, fe=NA_character_, window=NA_character_)
  out$window <- tail(parts, 1)
  out$fe     <- tail(parts, 2)[1]
  out$outcome <- sub("_eventstudy_feols_log1p.*$", "", parts[1])
  out
}

# ---- load latest models RDS produced by Step 9 v2
mod_path <- latest_file(DIRS$out_paper1, "^R09_models_event_study_v2_models_.*\\.rds$")
all_results <- readRDS(mod_path)
cat("Loaded models RDS:\n", mod_path, "\n")

# ---- build ES sample stats by reloading panels (fast)
pick_existing <- function(paths, label){
  for (p in paths) if (file.exists(p)) return(p)
  stop("No existing path found for ", label, ". Tried:\n", paste(paths, collapse = "\n"))
}

P_ALLCW <- pick_existing(
  c(file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli_v2.rds"),
    file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1_posteli.rds")),
  "ALLCW posteli"
)
P_STRICT <- pick_existing(
  c(file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli_v2.rds"),
    file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1_posteli.rds")),
  "STRICT posteli"
)

load_panel <- function(path) {
  dt <- readRDS(path); setDT(dt)
  if (!("pg_id" %in% names(dt)) && "gid" %in% names(dt)) dt[, pg_id := gid]
  dt[, pg_id := as.integer(pg_id)]
  dt[, month := as_month_idate(month)]
  if (!("theatre_id" %in% names(dt))) {
    if ("country" %in% names(dt)) dt[, theatre_id := as.character(country)]
    else if ("country_id" %in% names(dt)) dt[, theatre_id := paste0("country_id_", as.integer(country_id))]
    else stop("No theatre_id/country/country_id.")
  }
  dt[, month_int := (year(month) * 12L + month(month))]
  dt
}

panels <- list(
  ALLCW  = load_panel(P_ALLCW),
  STRICT = load_panel(P_STRICT)
)

get_es_sample_stats <- function(dt, window) {
  dt[, event_time_es := as.integer(event_time)]
  dt[event_time_es < window[1] | event_time_es > window[2], event_time_es := NA_integer_]
  dt_es <- dt[!is.na(event_time_es)]
  data.table(
    es_rows = nrow(dt_es),
    es_pg_clusters = uniqueN(dt_es$pg_id),
    es_month_clusters = uniqueN(dt_es$month_int),
    es_theatres = uniqueN(dt_es$theatre_id),
    share_zero_type3 = if ("events_type3" %in% names(dt_es)) mean(dt_es$events_type3 == 0) else NA_real_,
    share_zero_ged   = if ("events_ged" %in% names(dt_es)) mean(dt_es$events_ged == 0) else NA_real_
  )
}

win_map <- list(w24=c(-24L,24L), w60=c(-24L,60L))

rows <- list()

for (panel_name in names(all_results)) {
  es_list <- all_results[[panel_name]]$models_es
  if (!length(es_list)) next
  
  for (k in names(es_list)) {
    m <- es_list[[k]]
    meta <- parse_es_key(k)
    
    w <- win_map[[meta$window]]
    if (is.null(w)) next
    
    wl <- wald_leads_joint(m, lead_min = w[1], lead_max = -2L)
    
    ss <- get_es_sample_stats(panels[[panel_name]], w)
    
    rows[[length(rows)+1]] <- cbind(
      data.table(
        panel = panel_name,
        outcome = meta$outcome,
        fe_mode = meta$fe,
        window = meta$window,
        model_key = k,
        nobs = nobs(m),
        n_leads = wl$n_leads,
        wald_stat = wl$stat,
        wald_df = wl$df,
        wald_p = wl$p,
        pass_pretrend = !is.na(wl$p) & wl$p > 0.05,
        note = wl$note
      ),
      ss
    )
  }
}

chk <- rbindlist(rows, fill = TRUE)
setorder(chk, panel, outcome, fe_mode, window)

print(chk[, .(panel,outcome,fe_mode,window,nobs,es_rows,es_pg_clusters,es_month_clusters,es_theatres,n_leads,wald_stat,wald_df,wald_p,pass_pretrend)])

OUT <- file.path(DIRS$audit, paste0("R09c_postcheck_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
fwrite(chk, OUT)
cat("\nWrote postcheck:\n", OUT, "\n")
