# ============================================================================
# JPR H2 RE-TEST — SOURCE-VERIFIED, NO AUTOMATIC MANUSCRIPT CHANGES
# 2026-11-06
#
# Source basis: FINAL30/staged_original_code/H2 IMPROVEMENT.R, run_h2_model().
# Main H2 published results: FINAL30/08_FINAL30_CANONICAL_H2_spatial_results.csv
# Plus anticipation models: FINAL30/07_H2_anticipation_window_original_spec.csv
#
# Goal: for each of 3 spatial moderators, on SAME panel and same window,
# estimate original restricted PPML and its one-term-expanded sensitivity.
# NO assertion that the expanded specification must be preferred.
# NO coefficient replacement unless the original reproduces published values.
# NO overwriting original panels, manuscript, or model outputs.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
OUT <- file.path(ROOT, "outputs", "JPR_H2_PROPER_RETEST")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PANEL_FINAL30 <- file.path(ROOT,
  "outputs/final_airfield_screen_20260809_FINAL30/11_final_analysis_panel.rds")
PANEL_SPATIAL_SOURCE <- file.path(ROOT,
  "outputs/paper1_geo_afghanistan_executable_acled/diagnostics_h2_spatial_adm2/panel_core_with_spatial_h2_diagnostic.rds")

# Compare BOTH real panel versions. Do not mix columns from separate panels.
panels <- c(FINAL30 = PANEL_FINAL30, SPATIAL_ORIGINAL = PANEL_SPATIAL_SOURCE)
panels <- panels[file.exists(panels)]
if (!length(panels)) stop("Both known H2 panels are missing. Check ROOT/path settings.")

MODS <- c("pre_contested_share_spatial",
          "pre_taliban_activity_share_spatial",
          "pre_contested_intensity_spatial")
CONTROL_CANDIDATES <- c("prec_gpcp", "temp", "nlights_calib_mean")

# Exact saved FINAL30 estimates from included provenance files. These are
# comparison benchmarks, not pre-filled re-estimation output.
expected <- data.table(
  window = rep(c("full", "drop_anticipation"), each = 3L),
  moderator = rep(MODS, 2L),
  expected_estimate = c(-16.1023802013192, -17.6667522797993,
                        -3.12067660696125, -17.4082078077671,
                        -19.0707077670960, -3.38898087940946),
  expected_se = c(4.27319350902624, 4.60318002832200,
                  0.880878277343997, 4.86739390345531,
                  5.26815958334537, 0.984560054320593),
  expected_n = rep(c(54312L, 49956L), each = 3L)
)

# Avoid selecting the wrong ELI silently. The published anticipation output
# refers to 'eli_final'; the spatial diagnostic source uses 'eli_feas'.
eli_candidates <- function(panel_name, data_names) {
  preference <- if (panel_name == "FINAL30") c("eli_final", "eli_feas")
                else c("eli_feas", "eli_final")
  intersect(preference, data_names)
}

pull_coef <- function(model, term) {
  tab <- fixest::coeftable(model, vcov = ~ cell_id + month)
  if (!term %in% rownames(tab)) return(list(est = NA_real_, se = NA_real_,
    p = NA_real_, message = paste0("TERM NOT FOUND: ", term,
                                 " [available: ", paste(rownames(tab), collapse = " | "), "]")))
  x <- as.numeric(tab[term, 1:4])
  list(est = x[1], se = x[2], p = x[4], message = "")
}

safe_obs <- function(model) {
  tryCatch(as.integer(fixest::obs(model)), error = function(e) integer(0))
}

fit_pair <- function(d, panel_name, eli, window, moderator, controls) {
  post <- "post_closure"
  triple <- paste(post, eli, moderator, sep = ":")
  lower <- paste(post, moderator, sep = ":")
  treatment <- paste(post, eli, sep = ":")
  # This is the ORIGINAL exact functional form, not a generic factorial model:
  original_rhs <- c(treatment, triple, controls)
  expanded_rhs <- c(treatment, lower, triple, controls)
  rhs_versions <- list(original_restricted = original_rhs,
                       sensitivity_add_post_m = expanded_rhs)
  fitted <- list()
  results <- list()
  for (variant in names(rhs_versions)) {
    fml <- as.formula(paste("events_type3 ~",
                 paste(rhs_versions[[variant]], collapse = " + "),
                 "| cell_id + month"))
    cat("  ",panel_name," / ",eli," / ",window," / ",moderator," / ",variant,"\n",sep="")
    ans <- tryCatch(
      withCallingHandlers(fixest::fepois(fml, data = d),
        warning = function(w) {
          cat("    WARNING: ",conditionMessage(w),"\n",sep="")
          invokeRestart("muffleWarning")
        }),
      error = function(e) e)
    if (inherits(ans, "error")) {
      err <- conditionMessage(ans)
      cat("    MODEL ERROR: ",err,"\n",sep="")
      results[[variant]] <- data.table(
        panel = panel_name, eli = eli, window = window,
        moderator = moderator, specification = variant,
        term = triple, triple_est = NA_real_, triple_se = NA_real_,
        triple_p = NA_real_, post_m_est = NA_real_, post_m_se = NA_real_,
        post_m_p = NA_real_, nobs = NA_integer_,
        formula = paste(deparse(fml), collapse = ""),
        error = err)
      next
    }
    fitted[[variant]] <- ans
    tri <- pull_coef(ans, triple)
    low <- if (variant == "sensitivity_add_post_m") pull_coef(ans, lower)
           else list(est = NA_real_, se = NA_real_, p = NA_real_, message = "")
    errs <- c(tri$message,low$message)
    errs <- errs[nzchar(errs)]
    results[[variant]] <- data.table(
      panel = panel_name, eli = eli, window = window,
      moderator = moderator, specification = variant,
      term = triple, triple_est = tri$est, triple_se = tri$se,
      triple_p = tri$p, post_m_est = low$est, post_m_se = low$se,
      post_m_p = low$p, nobs = as.integer(nobs(ans)),
      formula = paste(deparse(fml), collapse = ""),
      error = paste(errs, collapse = "; "))
  }
  rows <- rbindlist(results, fill = TRUE)
  if (all(names(rhs_versions) %in% names(fitted))) {
    oa <- safe_obs(fitted$original_restricted)
    ob <- safe_obs(fitted$sensitivity_add_post_m)
    # obs indices are panel-relative for models fit using the same data object.
    same_sample <- length(oa) > 0L && length(ob) > 0L &&
      identical(sort(oa), sort(ob))
    rows[, same_estimation_rows := same_sample]
    cat("    Both models fit. Same estimation rows: ",same_sample,"\n",sep="")
  } else {
    rows[, same_estimation_rows := NA]
  }
  rows
}

out <- list(); kk <- 0L
inventory <- list()
for (panel_name in names(panels)) {
  p <- panels[[panel_name]]
  D <- as.data.table(readRDS(p))
  needed <- c("events_type3", "post_closure", "cell_id", "month", MODS)
  miss <- setdiff(needed, names(D))
  if (length(miss)) {
    cat("SKIPPING ",panel_name," — missing: ",paste(miss, collapse=", "),"\n",sep="")
    next
  }
  D[, month := as.IDate(month)]
  D[, cell_id := as.factor(cell_id)]
  controls <- intersect(CONTROL_CANDIDATES, names(D))
  usable_eli <- eli_candidates(panel_name, names(D))
  if (!length(usable_eli)) {
    cat("SKIPPING ",panel_name," — no ELI column available\n",sep="")
    next
  }
  cat("\nPANEL ",panel_name,": ",nrow(D)," rows, ",ncol(D)," columns; ELI = ",
      paste(usable_eli,collapse=", "),"; controls = ",
      paste(controls,collapse=", "),"\n",sep="")
  inventory[[panel_name]] <- data.table(
    panel=panel_name, path=p, total_rows=nrow(D),
    total_cells=uniqueN(D$cell_id),
    total_months=uniqueN(D$month),
    eli_columns=paste(usable_eli,collapse=";"),
    controls=paste(controls,collapse=";"),
    post_values=paste(sort(unique(D$post_closure)), collapse=";"),
    treatment_on_month=paste(sort(unique(as.character(D[post_closure == 1, month])))[1],collapse=""))
  windows <- list(full = D,
    drop_anticipation = D[month < as.IDate("2020-03-01") |
                            month >= as.IDate("2021-07-01")])
  for (eli in usable_eli) for (win in names(windows)) for (mod in MODS) {
    kk <- kk + 1L
    out[[kk]] <- fit_pair(windows[[win]], panel_name, eli, win, mod, controls)
  }
}
if (!length(out)) stop("No compatible panel/model combination could be fitted.")
results <- rbindlist(out, fill=TRUE)
results <- merge(results, expected, by=c("window","moderator"), all.x=TRUE,sort=FALSE)
results[, reproduces_published := fifelse(specification == "original_restricted",
  !is.na(triple_est) & !is.na(triple_se) & !is.na(nobs) &
  abs(triple_est - expected_estimate) < 0.001 &
  abs(triple_se - expected_se) < 0.001 & nobs == expected_n,
  NA)]
setcolorder(results, c("panel","eli","window","moderator","specification",
  "triple_est","triple_se","triple_p","post_m_est","post_m_se","post_m_p",
  "nobs","same_estimation_rows","expected_estimate","expected_se","expected_n",
  "reproduces_published","formula","error","term"))
fwrite(results, file.path(OUT, "01_ALL_H2_PAIRED_RESULTS.csv"))
if (length(inventory)) fwrite(rbindlist(inventory,fill=TRUE),
  file.path(OUT,"02_PANEL_PROVENANCE.csv"))

# Compact readable comparison. Keep the published value next to each pair.
comparison <- copy(results[, .(panel,eli,window,moderator,specification,
  triple_est,triple_se,triple_p,post_m_est,post_m_p,nobs,same_estimation_rows,
  expected_estimate,reproduces_published,error)])
fwrite(comparison,file.path(OUT,"03_SHORT_H2_COMPARISON.csv"))
cat("\n================ REPRODUCTION CHECK ================\n")
print(results[specification == "original_restricted",
  .(panel,eli,window,moderator,estimate=round(triple_est,4),
    published=round(expected_estimate,4),n=nobs,match=reproduces_published,error)],
  row.names=FALSE)
cat("\n================ PAIRED SENSITIVITY CHECK ================\n")
print(results[,.(panel,eli,window,moderator,specification,
  triple=round(triple_est,4),se=round(triple_se,4),
  p=signif(triple_p,3),post_m=round(post_m_est,4),
  n=nobs,same_estimation_rows)], row.names=FALSE)
cat("\nFiles saved to: ",OUT,"\n",sep="")
cat("Only a panel/ELI where the original model reproduces published estimates",
    " supports a claim about the *published* H2 sensitivity.\n")
cat("All other pairs are diagnostic only, even if their observation counts match.\n")
