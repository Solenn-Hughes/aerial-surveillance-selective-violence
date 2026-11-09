# ============================================================
# R11_check_paper1_results.R  (run anytime)
# - Reads latest Step 9/9b outputs from disk
# - Summarises pretrends (wald), main DID coef, and dispersion sanity
# ============================================================

suppressMessages(library(data.table))

cfg_path <- if (file.exists("R00_config.R")) "R00_config.R" else "R/00_config.R"
source(cfg_path)

latest_file <- function(dir, pattern) {
  ff <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (!length(ff)) stop("No files match pattern in ", dir, ": ", pattern)
  ff[which.max(file.info(ff)$mtime)]
}

cat("\n--- Latest Step 9 checkpoint (if exists) ---\n")
ckpt_pat <- "^R09c_checkpoint_fixed_ES_controls_.*\\.csv$"
ckpt_files <- list.files(DIRS$audit, pattern = ckpt_pat, full.names = TRUE)
if (length(ckpt_files)) {
  ckpt <- fread(latest_file(DIRS$audit, ckpt_pat))
  setorder(ckpt, panel, outcome, fe_mode, window)
  print(ckpt[, .(panel, outcome, fe_mode, window, n_leads, wald_stat, wald_df, wald_p)])
  cat("\nRule of thumb: wald_p > 0.10 => pretrends not rejected.\n")
} else {
  cat("No checkpoint file found yet. (Make R09 write one; see next section.)\n")
}

cat("\n--- Dispersion sanity (latest R09b summaries) ---\n")
disp_pat <- "^R09b_dispersion_summary_ALLCW_v2_.*\\.csv$"
disp_files <- list.files(DIRS$audit, pattern = disp_pat, full.names = TRUE)
if (length(disp_files)) {
  # print both type3 and ged summaries if they exist
  for (f in disp_files) {
    d <- fread(f)
    cat("\nFile:", f, "\n")
    print(d)
  }
} else {
  cat("No dispersion summaries found.\n")
}

cat("\n--- Dispersion conditional on sum_y > 0 (recommended check) ---\n")
for (y in c("events_type3", "events_ged")) {
  rds <- file.path(DIRS$panels, paste0("dispersion_theatre_month_ALLCW_v2_", y, ".rds"))
  if (!file.exists(rds)) next
  d <- readRDS(rds); setDT(d)
  
  d_pos <- d[sum_y > 0]
  if (!nrow(d_pos)) next
  
  cat("\nOutcome:", y, "\n")
  cat("Share sum_y==0:", mean(d$sum_y == 0), "\n")
  cat("Mean dispersion treated pre  (sum_y>0):", d_pos[treated==1 & post_eli==0, mean(dispersion)], "\n")
  cat("Mean dispersion treated post (sum_y>0):", d_pos[treated==1 & post_eli==1, mean(dispersion)], "\n")
}

cat("\n✅ Done.\n")
