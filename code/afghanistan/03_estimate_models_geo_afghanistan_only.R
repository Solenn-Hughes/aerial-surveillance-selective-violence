# ============================================================
# 03_estimate_models_geo_afghanistan_only_EXECUTABLE_ACLED.R
# Afghanistan executable branch:
# - H2 from ACLED district contestedness
# - H3 local-cluster headline
# ============================================================

source("C:/Users/USERNAME/OneDrive/Desktop/THESIS/Proper code/super proper/00_config_geo_afghanistan_only.R")
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

log_msg("03_estimate_models_geo_afghanistan_only: start")

get_model_nobs <- function(model) {
  out <- tryCatch(as.integer(stats::nobs(model)), error = function(e) NA_integer_)
  if (length(out) == 0L || is.na(out)) {
    out <- tryCatch(as.integer(model$nobs), error = function(e) NA_integer_)
  }
  out
}

fit_fixest <- function(fun_name, fml, data) {
  warns <- character(); err <- NULL
  fit <- withCallingHandlers(
    tryCatch({
      switch(fun_name,
             ppml = fepois(fml, data = data),
             ols_level = feols(fml, data = data),
             stop("Unsupported fun_name: ", fun_name))
    }, error = function(e) { err <<- e; NULL }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(model = fit, error = err, warnings = unique(warns))
}

safe_coeftable <- function(model, vcov_spec) {
  warns <- character(); err <- NULL
  ct <- withCallingHandlers(
    tryCatch(fixest::coeftable(model, vcov = vcov_spec), error = function(e) { err <<- e; NULL }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = ct, error = err, warnings = unique(warns))
}

available_controls <- function(dt) intersect(c("prec_gpcp", "temp", "nlights_calib_mean"), names(dt))
should_retry_without_controls <- function(msg) grepl("collinear|singular", msg %||% "", ignore.case = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

build_formula <- function(fun_name, y, rhs_terms, fe) {
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & nzchar(rhs_terms)]
  rhs <- if (length(rhs_terms) == 0L) "1" else paste(rhs_terms, collapse = " + ")
  as.formula(paste0(y, " ~ ", rhs, " | ", fe))
}

estimate_model <- function(fun_name, dt, y, treatment_terms, fe, primary_cluster, fallback_cluster = NULL,
                           controls = available_controls(dt), force_no_controls = FALSE) {
  use_dt <- dt[!is.na(get(y))]
  if (nrow(use_dt) == 0L) return(list(model = NULL, coeftable = NULL, status = "error", reason = "No usable rows", warnings = character(), cluster_used = NA_character_, vcov_fallback = NA, controls_dropped = NA))
  use_controls <- if (isTRUE(force_no_controls)) character() else controls
  rhs_terms_full <- unique(c(treatment_terms, use_controls))
  fml_full <- build_formula(fun_name, y, rhs_terms_full, fe)
  fit1 <- fit_fixest(fun_name, fml_full, use_dt)
  controls_dropped <- FALSE
  if (!is.null(fit1$error) && length(use_controls) > 0L && should_retry_without_controls(fit1$error$message)) {
    fml_retry <- build_formula(fun_name, y, treatment_terms, fe)
    fit2 <- fit_fixest(fun_name, fml_retry, use_dt)
    if (is.null(fit2$error)) { fit1 <- fit2; controls_dropped <- TRUE }
  }
  if (!is.null(fit1$error)) return(list(model = NULL, coeftable = NULL, status = "error", reason = fit1$error$message, warnings = fit1$warnings, cluster_used = NA_character_, vcov_fallback = NA, controls_dropped = controls_dropped))
  ct1 <- safe_coeftable(fit1$model, primary_cluster)
  primary_cluster_txt <- paste(deparse(primary_cluster), collapse = "")
  if (!is.null(ct1$error)) return(list(model = NULL, coeftable = NULL, status = "error", reason = ct1$error$message, warnings = unique(c(fit1$warnings, ct1$warnings)), cluster_used = NA_character_, vcov_fallback = NA, controls_dropped = controls_dropped))
  list(model = fit1$model, coeftable = ct1$value, status = "ok", reason = NA_character_, warnings = unique(c(fit1$warnings, ct1$warnings)), cluster_used = primary_cluster_txt, vcov_fallback = FALSE, controls_dropped = controls_dropped)
}

extract_results <- function(fit, family, panel, outcome, spec_name, tier, result_block) {
  ct <- as.data.table(fit$coeftable, keep.rownames = "term")
  if ("Estimate" %in% names(ct)) setnames(ct, "Estimate", "estimate")
  if ("Std. Error" %in% names(ct)) setnames(ct, "Std. Error", "std_error")
  if ("t value" %in% names(ct)) setnames(ct, "t value", "stat")
  if ("z value" %in% names(ct)) setnames(ct, "z value", "stat")
  if ("Pr(>|t|)" %in% names(ct)) setnames(ct, "Pr(>|t|)", "p_value")
  if ("Pr(>|z|)" %in% names(ct)) setnames(ct, "Pr(>|z|)", "p_value")
  ct[, `:=`(result_block = result_block, theatre = TARGET_THEATRE, panel = panel, outcome = outcome, family = family, spec = spec_name, tier = tier, nobs = get_model_nobs(fit$model), vcov_fallback = fit$vcov_fallback)]
  ct[, .(result_block, theatre, panel, outcome, family, spec, tier, term, estimate, std_error, stat, p_value, nobs, vcov_fallback)]
}

results <- list(); model_audit <- list()
register_fit <- function(result_block, fit, data_used, family, panel, outcome, spec_name, tier, id_var) {
  model_audit[[result_block]] <<- data.table(result_block = result_block, theatre = TARGET_THEATRE, panel = panel, outcome = outcome, family = family, spec = spec_name, tier = tier, status = fit$status, reason = ifelse(is.null(fit$reason), NA_character_, fit$reason), rows_in_data = nrow(data_used), nobs_estimation = if (!is.null(fit$model)) get_model_nobs(fit$model) else NA_integer_, nonzero_outcomes = if (outcome %in% names(data_used)) sum(data_used[[outcome]] > 0, na.rm = TRUE) else NA_integer_, n_units = if (id_var %in% names(data_used)) uniqueN(data_used[[id_var]]) else NA_integer_, n_months = if ("month" %in% names(data_used)) uniqueN(data_used$month) else NA_integer_, cluster = fit$cluster_used, vcov_fallback = fit$vcov_fallback, controls_dropped = fit$controls_dropped, collinear_vars = if (!is.null(fit$model) && length(fit$model$collin.var) > 0L) paste(fit$model$collin.var, collapse = ", ") else NA_character_, warning_text = if (length(fit$warnings) > 0L) paste(unique(fit$warnings), collapse = " | ") else NA_character_)
  if (identical(fit$status, "ok") && !is.null(fit$coeftable)) results[[result_block]] <<- extract_results(fit, family, panel, outcome, spec_name, tier, result_block)
}

core <- as.data.table(readRDS(file.path(DIRS$panels, "panel_core_cell_month.rds")))
core[, month := as.IDate(month)]
core[, cell_id := as.factor(cell_id)]

fit <- estimate_model("ppml", core, "events_type3", c("post_closure:eli_feas"), "cell_id + month", ~ cell_id + month, ~ cell_id)
register_fit("afghanistan_h1_core_ppml", fit, core, "PPML", "core", "events_type3", "H1: post_closure x eli_feas", "main", "cell_id")

# ACLED-based H2 main
if ("pre_contested_share" %in% names(core)) {
  fit <- estimate_model("ppml", core, "events_type3",
                        c("post_closure:eli_feas", "post_closure:eli_feas:pre_contested_share"),
                        "cell_id + month", ~ cell_id + month, ~ cell_id)
  register_fit("afghanistan_h2_contested_share_ppml", fit, core, "PPML", "core", "events_type3", "H2 main: post_closure x eli_feas x pre_contested_share", "main", "cell_id")
}
if ("pre_taliban_activity_share" %in% names(core)) {
  fit <- estimate_model("ppml", core, "events_type3",
                        c("post_closure:eli_feas", "post_closure:eli_feas:pre_taliban_activity_share"),
                        "cell_id + month", ~ cell_id + month, ~ cell_id)
  register_fit("afghanistan_h2_taliban_activity_share_ppml", fit, core, "PPML", "core", "events_type3", "H2 appendix: post_closure x eli_feas x pre_taliban_activity_share", "appendix", "cell_id")
}
if ("pre_contested_intensity" %in% names(core)) {
  fit <- estimate_model("ppml", core, "events_type3",
                        c("post_closure:eli_feas", "post_closure:eli_feas:pre_contested_intensity"),
                        "cell_id + month", ~ cell_id + month, ~ cell_id)
  register_fit("afghanistan_h2_contested_intensity_ppml", fit, core, "PPML", "core", "events_type3", "H2 appendix: post_closure x eli_feas x pre_contested_intensity", "appendix", "cell_id")
}
if ("pre_weak_control_or_contested_legacy" %in% names(core)) {
  fit <- estimate_model("ppml", core, "events_type3",
                        c("post_closure:eli_feas", "post_closure:eli_feas:pre_weak_control_or_contested_legacy"),
                        "cell_id + month", ~ cell_id + month, ~ cell_id)
  register_fit("afghanistan_h2_legacy_proxy_ppml", fit, core, "PPML", "core", "events_type3", "H2 legacy appendix", "appendix", "cell_id")
}

for (yy in c("events_protests", "events_riots", "events_type1")) {
  if (yy %in% names(core)) {
    fit <- estimate_model("ppml", core, yy, c("post_closure:eli_feas"), "cell_id + month", ~ cell_id + month, ~ cell_id)
    block <- switch(yy, events_protests = "afghanistan_nc_protests_ppml", events_riots = "afghanistan_nc_riots_ppml", events_type1 = "afghanistan_contrast_type1_ppml")
    spec <- switch(yy, events_protests = "negative control: protests", events_riots = "negative control: riots", events_type1 = "contrast outcome: type1 battles")
    register_fit(block, fit, core, "PPML", "core", yy, spec, "support", "cell_id")
  }
}

if ("low_visibility_q25" %in% names(core)) {
  fit <- estimate_model("ppml", core, "events_type3", c("post_closure:eli_feas", "post_closure:eli_feas:low_visibility_q25"), "cell_id + month", ~ cell_id + month, ~ cell_id)
  register_fit("afghanistan_nc_low_visibility_ppml", fit, core, "PPML", "core", "events_type3", "negative-control geography", "support", "cell_id")
}



# Event-study / pretrend diagnostic for H1
parse_rel_month_term <- function(term, rel_min = EVENT_STUDY_MIN, rel_max = EVENT_STUDY_MAX) {
  nums <- regmatches(term, gregexpr("-?[0-9]+", term, perl = TRUE))[[1]]
  if (length(nums) == 0L) return(NA_integer_)
  vals <- suppressWarnings(as.integer(nums))
  vals <- vals[!is.na(vals)]
  vals <- vals[vals >= rel_min & vals <= rel_max]
  if (length(vals) == 0L) return(NA_integer_)
  vals[length(vals)]
}

safe_vcov_matrix <- function(model, vcov_spec = ~ cell_id + month) {
  V <- tryCatch({
    sm <- summary(model, vcov = vcov_spec)
    out <- sm$cov.scaled
    if (is.null(out)) out <- tryCatch(stats::vcov(sm), error = function(e) NULL)
    if (is.null(out)) return(NULL)
    as.matrix(out)
  }, error = function(e) NULL)
  V
}

compute_joint_pretrend_test <- function(model, rel_map, vcov_spec = ~ cell_id + month) {
  coef_names <- names(stats::coef(model))
  use_idx <- which(rel_map$rel_month < -1L & rel_map$term %in% coef_names)
  if (length(use_idx) == 0L) {
    return(data.table(
      test = "joint_pretrend",
      n_leads = 0L,
      wald_stat = NA_real_,
      df = 0L,
      p_value = NA_real_,
      vcov_type = paste(deparse(vcov_spec), collapse = "")
    ))
  }
  keep_terms <- rel_map$term[use_idx]
  b <- stats::coef(model)[keep_terms]
  V_full <- safe_vcov_matrix(model, vcov_spec = vcov_spec)
  if (is.null(V_full)) {
    return(data.table(
      test = "joint_pretrend",
      n_leads = length(b),
      wald_stat = NA_real_,
      df = length(b),
      p_value = NA_real_,
      vcov_type = paste(deparse(vcov_spec), collapse = "")
    ))
  }
  V <- V_full[keep_terms, keep_terms, drop = FALSE]
  stat <- pval <- NA_real_
  ok <- TRUE
  if (any(is.na(b)) || any(!is.finite(b)) || any(!is.finite(V))) ok <- FALSE
  if (ok) {
    stat <- tryCatch(as.numeric(t(b) %*% qr.solve(V, b)), error = function(e) NA_real_)
    pval <- if (is.finite(stat)) stats::pchisq(stat, df = length(b), lower.tail = FALSE) else NA_real_
  }
  data.table(
    test = "joint_pretrend",
    n_leads = length(b),
    wald_stat = stat,
    df = length(b),
    p_value = pval,
    vcov_type = paste(deparse(vcov_spec), collapse = "")
  )
}

build_event_study_outputs <- function(model, out_stub, support_dt = NULL,
                                      rel_min = EVENT_STUDY_MIN, rel_max = EVENT_STUDY_MAX) {
  ct <- as.data.table(fixest::coeftable(model, vcov = ~ cell_id + month), keep.rownames = "term")
  if ("Estimate" %in% names(ct)) setnames(ct, "Estimate", "estimate")
  if ("Std. Error" %in% names(ct)) setnames(ct, "Std. Error", "std_error")
  if ("Pr(>|t|)" %in% names(ct)) setnames(ct, "Pr(>|t|)", "p_value")
  if ("Pr(>|z|)" %in% names(ct)) setnames(ct, "Pr(>|z|)", "p_value")

  es <- copy(ct)
  es[, rel_month := vapply(term, parse_rel_month_term, integer(1), rel_min = rel_min, rel_max = rel_max)]
  es <- es[!is.na(rel_month)]
  if (nrow(es) == 0L) {
    write_audit_text(c("Event study estimated, but no relative-month terms were parsed."), paste0(out_stub, "_audit"))
    return(invisible(NULL))
  }

  es[, `:=`(
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error,
    is_pre = rel_month < -1L,
    is_post = rel_month >= 0L
  )]
  es <- es[order(rel_month)]
  es_plot <- rbind(
    es[, .(term, rel_month, estimate, std_error, p_value, ci_low, ci_high)],
    data.table(term = "(reference -1)", rel_month = -1L, estimate = 0, std_error = NA_real_,
               p_value = NA_real_, ci_low = 0, ci_high = 0),
    fill = TRUE
  )
  es_plot <- unique(es_plot, by = c("term", "rel_month"))
  setorder(es_plot, rel_month)

  fwrite(es_plot, file.path(DIRS$tables, paste0(out_stub, "_coefficients.csv")))
  if (!is.null(support_dt)) fwrite(support_dt, file.path(DIRS$audit, paste0(out_stub, "_support.csv")))

  pre <- es[rel_month < -1L]
  joint_pre <- compute_joint_pretrend_test(model, es[, .(term, rel_month)], vcov_spec = ~ cell_id + month)
  fwrite(joint_pre, file.path(DIRS$audit, paste0(out_stub, "_joint_pretrend_test.csv")))

  pretrend_lines <- c(
    paste0("Event-study diagnostic: ", out_stub),
    paste0("Window: [", rel_min, ", ", rel_max, "] with reference month -1."),
    paste0("Estimated bins parsed: ", paste(es_plot$rel_month, collapse = ", ")),
    paste0("Number of pre-period lead coefficients (< -1): ", nrow(pre)),
    paste0("Pre-period leads with p < 0.10: ", sum(pre$p_value < 0.10, na.rm = TRUE)),
    paste0("Pre-period leads with p < 0.05: ", sum(pre$p_value < 0.05, na.rm = TRUE)),
    paste0("Largest absolute pre-period coefficient: ", if (nrow(pre) > 0L) sprintf("%.4f", max(abs(pre$estimate), na.rm = TRUE)) else "NA"),
    paste0("Joint pretrend Wald chi2: ", if (is.finite(joint_pre$wald_stat)) sprintf("%.4f", joint_pre$wald_stat) else "NA"),
    paste0("Joint pretrend df: ", joint_pre$df),
    paste0("Joint pretrend p-value: ", if (is.finite(joint_pre$p_value)) sprintf("%.4f", joint_pre$p_value) else "NA")
  )
  write_audit_text(pretrend_lines, paste0(out_stub, "_audit"))

  png(file.path(DIRS$tables, paste0(out_stub, "_plot.png")), width = 1100, height = 700, res = 120)
  plot(es_plot$rel_month, es_plot$estimate,
       ylim = range(c(es_plot$ci_low, es_plot$ci_high), na.rm = TRUE),
       xlab = paste0("Months relative to ", format(SHOCK_DATE, "%B %Y")),
       ylab = "Coefficient on high-ELI cell × relative month",
       main = paste0("Afghanistan / Bagram event-study diagnostic (", rel_min, " to ", rel_max, ")"),
       pch = 19)
  abline(h = 0, lty = 2)
  abline(v = -1, lty = 3)
  segments(es_plot$rel_month, es_plot$ci_low, es_plot$rel_month, es_plot$ci_high)
  dev.off()
  invisible(es_plot)
}

run_event_study_block <- function(dt, rel_min, rel_max, out_stub, result_block, spec_label) {
  es_dt <- copy(dt)
  es_dt[, month_num := month_int_from_idate(month)]
  shock_num <- month_int_from_idate(SHOCK_DATE)
  es_dt[, rel_month := month_num - shock_num]
  es_dt <- es_dt[rel_month >= rel_min & rel_month <= rel_max]
  if (!"eli_high_theatre" %in% names(es_dt) || all(is.na(es_dt$eli_high_theatre))) {
    med_eli <- median(es_dt$eli_feas, na.rm = TRUE)
    es_dt[, eli_high_theatre := as.integer(eli_feas >= med_eli)]
  }
  es_support <- es_dt[, .(
    treated_positive = sum(events_type3 > 0 & eli_high_theatre == 1L, na.rm = TRUE),
    treated_events = sum(events_type3[eli_high_theatre == 1L], na.rm = TRUE),
    treated_cells_nonzero = uniqueN(cell_id[events_type3 > 0 & eli_high_theatre == 1L]),
    total_events = sum(events_type3, na.rm = TRUE)
  ), by = rel_month]
  keep_rel <- es_support[treated_positive >= 5 | rel_month == -1L, rel_month]
  es_dt <- es_dt[rel_month %in% keep_rel]

  if (nrow(es_dt) == 0L) {
    write_audit_text(c("Event study skipped.", "No rows remained after support filtering."), paste0(out_stub, "_audit"))
    return(invisible(NULL))
  }

  fit <- estimate_model("ppml", es_dt, "events_type3", c("i(rel_month, eli_high_theatre, ref = -1)"),
                        "cell_id + month", ~ cell_id + month, ~ cell_id,
                        controls = character(0), force_no_controls = TRUE)
  register_fit(result_block, fit, es_dt, "PPML_event_study", "core", "events_type3", spec_label, "support", "cell_id")
  if (identical(fit$status, "ok") && !is.null(fit$model)) {
    build_event_study_outputs(fit$model, out_stub, es_support, rel_min = rel_min, rel_max = rel_max)
  } else {
    write_audit_text(c("Event study failed.", paste0("Reason: ", fit$reason %||% "unknown")), paste0(out_stub, "_audit"))
  }
  invisible(NULL)
}

# Event-study / pretrend diagnostics for H1
run_event_study_block(
  dt = core,
  rel_min = EVENT_STUDY_MIN,
  rel_max = EVENT_STUDY_MAX,
  out_stub = "afghanistan_h1_event_study",
  result_block = "afghanistan_h1_event_study_ppml",
  spec_label = "H1 event study: rel_month x eli_high_theatre"
)

run_event_study_block(
  dt = core,
  rel_min = -6L,
  rel_max = 6L,
  out_stub = "afghanistan_h1_event_study_short",
  result_block = "afghanistan_h1_event_study_short_ppml",
  spec_label = "H1 short-window event study: rel_month x eli_high_theatre"
)

# H1 robustness excluding the very large +6 month from the event-study plot
core_drop_p6 <- copy(core)
core_drop_p6[, rel_month_h1 := month_int_from_idate(month) - month_int_from_idate(SHOCK_DATE)]
core_drop_p6 <- core_drop_p6[rel_month_h1 != 6L]
fit <- estimate_model("ppml", core_drop_p6, "events_type3", c("post_closure:eli_feas"), "cell_id + month", ~ cell_id + month, ~ cell_id)
register_fit("afghanistan_h1_core_ppml_drop_rel_month_6", fit, core_drop_p6, "PPML", "core", "events_type3", "H1 robustness: drop relative month +6", "appendix", "cell_id")

placebo <- copy(core); placebo[, placebo_post := as.integer(month >= PLACEBO_DATE)]
fit <- estimate_model("ppml", placebo, "events_type3", c("placebo_post:eli_feas"), "cell_id + month", ~ cell_id + month, ~ cell_id)
register_fit("afghanistan_placebo_time_ppml", fit, placebo, "PPML", "core", "events_type3", "time placebo", "support", "cell_id")

h3_local <- copy(core[low_eli_cell == 1L & h3_local_sample == 1L & (adjacent_0_50 == 1L | fringe_50_100 == 1L)])
if (nrow(h3_local) > 0L && "source_cluster_id" %in% names(h3_local)) {
  h3_local[, source_cluster_id := as.factor(source_cluster_id)]
  h3_local[, source_cluster_month := interaction(source_cluster_id, month, drop = TRUE)]
  fit <- estimate_model("ppml", h3_local, "events_type3", c("post_closure:adjacent_0_50"), "cell_id + source_cluster_month", ~ cell_id + source_cluster_id, ~ cell_id)
  register_fit("afghanistan_h3_local_cluster_ppml", fit, h3_local, "PPML", "core", "events_type3", "H3 main local source-cluster", "main", "cell_id")
}

if (all(c("refuge_ring_0_50", "refuge_ring_50_100", "refuge_ring_100_150") %in% names(core))) {
  h3_global <- copy(core[low_eli_cell == 1L])
  fit <- estimate_model("ppml", h3_global, "events_type3", c("post_closure:refuge_ring_0_50", "post_closure:refuge_ring_50_100", "post_closure:refuge_ring_100_150"), "cell_id + month", ~ cell_id + month, ~ cell_id)
  register_fit("afghanistan_h3_global_rings_ppml", fit, h3_global, "PPML", "core", "events_type3", "H3 appendix global rings", "appendix", "cell_id")
}

if (length(results) == 0L) stop("No models were successfully estimated.")
all_results <- rbindlist(results, fill = TRUE); audit_dt <- rbindlist(model_audit, fill = TRUE)
fwrite(all_results, file.path(DIRS$tables, "all_model_results.csv"))
fwrite(audit_dt, file.path(DIRS$audit, "model_estimation_audit.csv"))
write_audit_text(c("Estimation summary (executable ACLED branch)", "H1 main retained.", "H1 event-study / pretrend diagnostics added and exported in both full and short windows.", "A joint pretrend Wald test is written for each event-study specification.", "H1 robustness excluding relative month +6 is also estimated.", "H2 main now uses ACLED district contestedness.", "H3 main is local source-cluster spillover.", "Negative controls: protests, riots, low-visibility geography, time placebo.", "Type1 retained as contrast outcome."), "estimation_summary_executable_acled")
log_msg("03_estimate_models_geo_afghanistan_only_EXECUTABLE_ACLED: done")
