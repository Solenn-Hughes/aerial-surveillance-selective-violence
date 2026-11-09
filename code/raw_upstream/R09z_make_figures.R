# ============================================================
# R09z_make_figures_v2.R
# Same as R09z_make_figures.R but fixes axes/labels (esp Fig 4)
# Outputs:
#  outputs/paper1/figures/Figure1_EventStudy_Baseline.png
#  outputs/paper1/figures/Figure2_EventStudy_AltRef_Donut.png
#  outputs/paper1/figures/Figure3_DecayRadius_Robustness.png
#  outputs/paper1/figures/Figure4_PretrendControlled_Robustness.png
#  outputs/paper1/audit/R09z_figures_manifest_*.csv
# ============================================================

rm(list = ls())
suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
if (!file.exists(cfg_path)) stop("Config not found at R00_config.R or R/00_config.R")
source(cfg_path)

AUDIT_DIR <- if (exists("DIRS") && "audit" %in% names(DIRS)) DIRS$audit else file.path("outputs","paper1","audit")
FIG_DIR   <- file.path("outputs","paper1","figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

latest_file <- function(dir, pattern) {
  ff <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(ff)) return(NA_character_)
  ff[which.max(file.info(ff)$mtime)]
}
assert_file <- function(p) if (is.na(p) || !nzchar(p) || !file.exists(p)) stop("INTEGRITY FAIL: Missing file: ", p)

write_manifest <- function(dt, stem="R09z_figures_manifest") {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out <- file.path(AUDIT_DIR, paste0(stem, "_", ts, ".csv"))
  fwrite(dt, out)
  message("Wrote manifest: ", out)
  invisible(out)
}

ES_COEFS <- latest_file(AUDIT_DIR, "^R09y_pretrend_diagnostics_m0boot_v2_es_coefs_.*\\.csv$")
DECAY    <- latest_file(AUDIT_DIR, "^R09f_m0boot_decay_radius_.*\\.csv$")
TAB5     <- latest_file(AUDIT_DIR, "^R09h_table5_full_export_.*\\.csv$")
if (is.na(TAB5)) TAB5 <- latest_file(AUDIT_DIR, "^R09h_m0boot_pretrend_controlled_.*\\.csv$")

assert_file(ES_COEFS); assert_file(DECAY); assert_file(TAB5)

es_all <- fread(ES_COEFS)
decay  <- fread(DECAY)
tab5   <- fread(TAB5)

# ---- Standardize ES columns
if (!("k" %in% names(es_all))) {
  if (!("term" %in% names(es_all))) stop("INTEGRITY FAIL: ES coef file missing both k and term.")
  es_all[, k := as.integer(sub("^event_time_es::(-?\\d+):eli_feas_t$", "\\1", term))]
}
if (!("Estimate" %in% names(es_all)) && "estimate" %in% names(es_all)) setnames(es_all, "estimate", "Estimate")
if (!("Std. Error" %in% names(es_all))) {
  if ("std.error" %in% names(es_all)) setnames(es_all, "std.error", "Std. Error")
  if ("se" %in% names(es_all)) setnames(es_all, "se", "Std. Error")
}
if (!("spec" %in% names(es_all))) stop("INTEGRITY FAIL: ES coef file missing 'spec'.")
if (!all(c("Estimate","Std. Error","k") %in% names(es_all))) stop("INTEGRITY FAIL: ES coef file missing Estimate/Std. Error/k.")

# ---- Plot helpers
plot_es <- function(es_dt, file_out, main, subtitle=NULL) {
  es_dt <- copy(es_dt)
  setorder(es_dt, k)
  es_dt[, ci_lo := Estimate - 1.96 * `Std. Error`]
  es_dt[, ci_hi := Estimate + 1.96 * `Std. Error`]
  
  png(file_out, width=1800, height=1100, res=200)
  op <- par(mar=c(5,5,4,2) + 0.1)
  
  ylim <- range(c(es_dt$ci_lo, es_dt$ci_hi), finite=TRUE)
  plot(es_dt$k, es_dt$Estimate,
       type="n",
       xlim=c(-24, 60),
       ylim=ylim,
       xlab="Event time (months relative to onset)",
       ylab="Event-time coefficient on ELI feasibility (relative to ref)",
       main=main)
  
  abline(h=0, lty=2)
  abline(v=0, lty=3)
  
  segments(es_dt$k, es_dt$ci_lo, es_dt$k, es_dt$ci_hi)
  points(es_dt$k, es_dt$Estimate, pch=16)
  
  if (!is.null(subtitle)) mtext(subtitle, side=3, line=0.5, cex=0.9)
  
  par(op)
  dev.off()
}

plot_points_ci_discrete <- function(labels, est, se, file_out, main, ylab) {
  x <- seq_along(est)
  ci_lo <- est - 1.96*se
  ci_hi <- est + 1.96*se
  
  png(file_out, width=1800, height=1100, res=200)
  op <- par(mar=c(7,5,4,2) + 0.1)  # bigger bottom margin for labels
  
  ylim <- range(c(ci_lo, ci_hi), finite=TRUE)
  plot(x, est, type="n", xaxt="n",
       xlab="", ylab=ylab, ylim=ylim, main=main)
  axis(1, at=x, labels=labels, las=2, cex.axis=0.95)  # vertical labels
  abline(h=0, lty=2)
  
  segments(x, ci_lo, x, ci_hi)
  points(x, est, pch=16)
  
  par(op)
  dev.off()
}

# ---- Figure 1 & 2 specs
spec_baseline <- unique(es_all[grepl("baseline", spec, ignore.case=TRUE)]$spec)
if (!length(spec_baseline)) stop("Could not find baseline ES spec in ES coef file.")
spec_baseline <- spec_baseline[1]

spec_donut <- unique(es_all[grepl("donut", spec, ignore.case=TRUE)]$spec)
if (!length(spec_donut)) stop("Could not find donut ES spec in ES coef file.")
spec_donut <- spec_donut[1]

es1 <- es_all[spec == spec_baseline]
es2 <- es_all[spec == spec_donut]

F1 <- file.path(FIG_DIR, "Figure1_EventStudy_Baseline.png")
F2 <- file.path(FIG_DIR, "Figure2_EventStudy_AltRef_Donut.png")
plot_es(es1, F1, "Figure 1. Event study (baseline reference)", paste0("Spec: ", spec_baseline))
plot_es(es2, F2, "Figure 2. Event study (alt reference + donut)", paste0("Spec: ", spec_donut))

# ---- Figure 3 (Decay radius)
if (!all(c("R_km","coef","se") %in% names(decay))) stop("Decay file missing R_km/coef/se.")
setorder(decay, R_km)

F3 <- file.path(FIG_DIR, "Figure3_DecayRadius_Robustness.png")
plot_points_ci_discrete(
  labels = paste0(decay$R_km, " km"),
  est = decay$coef,
  se  = decay$se,
  file_out = F3,
  main = "Figure 3. Decay radius robustness (m0_boot spec)",
  ylab = "Coefficient on post_eli × eli_feas_R"
)

# ---- Figure 4 (Pretrend-controlled robustness)
if ("row" %in% names(tab5)) {
  t5 <- tab5[row == "post_eli × eli_feas"]
  ord <- c("M1 baseline","M2 + lin trend","M3 + lin+quad")
  t5[, spec := factor(spec, levels=ord)]
  setorder(t5, spec)
  labels <- as.character(t5$spec)
  est <- t5$coef
  se  <- t5$se
} else {
  # fallback simple file
  t5 <- copy(tab5)
  t5[, spec := as.character(spec)]
  t5[, ord_i := fifelse(grepl("^M1", spec), 1L,
                        fifelse(grepl("^M2", spec), 2L,
                                fifelse(grepl("^M3", spec), 3L, 99L)))]
  setorder(t5, ord_i)
  labels <- t5$spec
  est <- t5$coef
  se  <- t5$se
}

F4 <- file.path(FIG_DIR, "Figure4_PretrendControlled_Robustness.png")
plot_points_ci_discrete(
  labels = labels,
  est = est,
  se  = se,
  file_out = F4,
  main = "Figure 4. Pretrend-controlled robustness (m0_boot spec)",
  ylab = "Coefficient on post_eli × eli_feas"
)

manifest <- data.table(
  figure = c("Figure 1","Figure 2","Figure 3","Figure 4"),
  file = c(F1, F2, F3, F4),
  source_data = c(ES_COEFS, ES_COEFS, DECAY, TAB5),
  spec_used = c(spec_baseline, spec_donut, "R=200/400/800", "M1/M2/M3")
)
write_manifest(manifest)

cat("\n✅ Figures written to: ", FIG_DIR, "\n", sep="")
print(manifest[, .(figure, file, spec_used)])