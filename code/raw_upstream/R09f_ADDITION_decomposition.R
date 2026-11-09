
# ============================================================
# ADDITION TO R09f_robustness_decay_radius.R
# Paste this at the END of R09f (after the decay sweep).
# 
# Runs eli_feas DECOMPOSED into components:
#   (a) Distance decay only (eli_feas_range)
#   (b) Terrain openness only (eli_feas_terrain)
#   (c) Full composite (eli_feas = range × terrain)
#
# This confirms β is not artefactual to index construction.
# ============================================================

cat("\n====== ELI_FEAS COMPONENT DECOMPOSITION ======\n")

# Ensure components are in panel
if (!("eli_feas_terrain" %in% names(dt))) {
  eli <- readRDS(file.path(DIRS$eli, "eli_feasibility_pg.rds")); setDT(eli)
  eli2 <- eli[, .(pg_id = as.integer(gid), eli_feas_terrain)]
  if ("eli_feas_terrain" %in% names(dt)) dt[, eli_feas_terrain := NULL]
  dt <- merge(dt, eli2, by = "pg_id", all.x = TRUE)
}

if (!("eli_feas_range" %in% names(dt))) {
  # Recompute at baseline R_KM = 400
  dt[, eli_feas_range := exp(-min_km_airbase / 400)]
}

ylog <- paste0("log1p_", OUTCOME)
if (!(ylog %in% names(dt))) dt[, (ylog) := log1p(get(OUTCOME))]

decomp_results <- list()

# (a) Distance only
f_dist <- as.formula(paste0(ylog, " ~ post_eli:eli_feas_range + ", rhs_controls, " | ", FE_SPEC))
m_dist <- tryCatch(feols(f_dist, data = dt, vcov = vcov_arg), error = function(e) NULL)

# (b) Terrain only
f_terr <- as.formula(paste0(ylog, " ~ post_eli:eli_feas_terrain + ", rhs_controls, " | ", FE_SPEC))
m_terr <- tryCatch(feols(f_terr, data = dt, vcov = vcov_arg), error = function(e) NULL)

# (c) Full composite (already run as baseline, but repeat for clean comparison)
f_full <- as.formula(paste0(ylog, " ~ post_eli:eli_feas + ", rhs_controls, " | ", FE_SPEC))
m_full <- tryCatch(feols(f_full, data = dt, vcov = vcov_arg), error = function(e) NULL)

extract_coef <- function(m, label) {
  if (is.null(m)) return(data.table(component = label, coef = NA, se = NA, p = NA, nobs = NA))
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  hit <- ct[grepl("eli_feas|post_eli", term)][1]
  data.table(component = label, coef = round(hit$Estimate, 5),
             se = round(hit$`Std. Error`, 5), p = round(hit$`Pr(>|t|)`, 4), nobs = nobs(m))
}

decomp_out <- rbindlist(list(
  extract_coef(m_dist, "Distance decay only (eli_feas_range)"),
  extract_coef(m_terr, "Terrain openness only (eli_feas_terrain)"),
  extract_coef(m_full, "Full composite (eli_feas = range × terrain)")
))

cat("\nComponent Decomposition Results:\n")
print(decomp_out)
cat("\nIf all three are positive and significant: both components contribute.\n")
cat("If only one is significant: that component drives the result.\n")

write_audit_dt(decomp_out, "R09f_component_decomposition")
