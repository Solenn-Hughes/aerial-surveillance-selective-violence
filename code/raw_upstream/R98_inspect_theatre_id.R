# ============================================================
# R98_inspect_theatre_id.R — inspect what IDs exist for "theatre"
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found: R00_config.R or R/00_config.R")
source(cfg_path)

P1_ALLCW  <- file.path(DIRS$panels, "panel_pg_month_ALLCW_paper1.rds")
P1_STRICT <- file.path(DIRS$panels, "panel_pg_month_ALLCW_ACD_STRICT_paper1.rds")
GED_CLEAN <- file.path(DIRS$base,   "ged_clean.rds")

for (pp in c(P1_ALLCW, P1_STRICT, GED_CLEAN)) {
  if (!file.exists(pp)) stop("Missing required file: ", pp)
}

load_dt <- function(path, label){
  dt <- readRDS(path); setDT(dt)
  data.table(dataset = label,
             n = nrow(dt),
             k = ncol(dt),
             cols = paste(names(dt), collapse = " | "))
}

# Candidate ID-like columns
pick_id_cols <- function(nms){
  nms[grepl("theatre|conflict|episode|dyad|side|actor|country|iso|gwno|^id$|_id$",
            nms, ignore.case = TRUE)]
}

p_allcw  <- readRDS(P1_ALLCW);  setDT(p_allcw)
p_strict <- readRDS(P1_STRICT); setDT(p_strict)
ged      <- readRDS(GED_CLEAN); setDT(ged)

# Standard key check (should pass)
if (!("pg_id" %in% names(p_allcw)) && "gid" %in% names(p_allcw)) p_allcw[, pg_id := gid]
if (!("pg_id" %in% names(p_strict)) && "gid" %in% names(p_strict)) p_strict[, pg_id := gid]

assert_cols(p_allcw,  c("pg_id","month"), "ALLCW_paper1")
assert_cols(p_strict, c("pg_id","month"), "STRICT_paper1")

# Summaries
out <- rbindlist(list(
  data.table(dataset="ALLCW_paper1",  id_cols = paste(pick_id_cols(names(p_allcw)),  collapse=", ")),
  data.table(dataset="STRICT_paper1", id_cols = paste(pick_id_cols(names(p_strict)), collapse=", ")),
  data.table(dataset="GED_clean",     id_cols = paste(pick_id_cols(names(ged)),      collapse=", "))
), fill=TRUE)

# For the most important candidates, print distinct counts
count_candidates <- function(dt, label){
  cand <- intersect(c("country_id","country","iso3","conflict_id","conflict_new_id","episode_id","dyad_id"),
                    names(dt))
  if (!length(cand)) return(data.table(dataset=label, var=character(), uniqueN=integer()))
  rbindlist(lapply(cand, function(v) data.table(dataset=label, var=v, uniqueN=uniqueN(dt[[v]]))))
}

counts <- rbindlist(list(
  count_candidates(p_allcw,  "ALLCW_paper1"),
  count_candidates(p_strict, "STRICT_paper1"),
  count_candidates(ged,      "GED_clean")
), fill=TRUE)

audit <- merge(out, counts, by="dataset", all=TRUE)

write_audit_dt(audit, "R98_inspect_theatre_id")
print(audit)
cat("\nIf you paste the printed table above, I can tell you exactly what theatre_id we can use.\n")
