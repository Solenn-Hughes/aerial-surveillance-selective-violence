# ============================================================
# H2 ANTICIPATION-WINDOW ROBUSTNESS — ORIGINAL DISSERTATION SPEC
# ============================================================
# Purpose:
#   1) Replicate the three FINAL H2 models on the full sample.
#   2) Re-estimate the SAME H2 models after dropping the entire
#      anticipation window: March 2020 through June 2021.
#
# IMPORTANT:
#   This deliberately preserves the dissertation H2 specification:
#       post_closure:eli_final
#       post_closure:eli_final:moderator
#   It DOES NOT add post_closure:moderator.
#
# Expected full-sample replication targets:
#   pre-contested share       = -16.1023802
#   Taliban activity share    = -17.6667523
#   contested intensity       =  -3.1206766
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

try(fixest::setFixest_notes(FALSE), silent = TRUE)

# ------------------------------------------------------------
# 0) PATHS — same final panel used by the canonical FINAL30 pipeline
# ------------------------------------------------------------
ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
FINAL_DIR <- file.path(
  ROOT,
  "outputs",
  "final_airfield_screen_20260809_FINAL30"
)
PANEL_PATH <- file.path(FINAL_DIR, "11_final_analysis_panel.rds")
OUT_PATH <- file.path(FINAL_DIR, "H2_anticipation_window_original_spec.csv")

if (!file.exists(PANEL_PATH)) {
  stop(
    "Final panel not found at:\n", PANEL_PATH,
    "\n\nChange only ROOT at the top of this script if your THESIS folder is on another drive."
  )
}

# ------------------------------------------------------------
# 1) LOAD AND NORMALISE FINAL PANEL
# ------------------------------------------------------------
dt <- as.data.table(readRDS(PANEL_PATH))

if (!"cell_id_num" %in% names(dt)) {
  if (!"cell_id" %in% names(dt)) stop("Panel has neither cell_id_num nor cell_id.")
  dt[, cell_id_num := suppressWarnings(as.integer(as.character(cell_id)))]
}

if (!"month" %in% names(dt)) stop("Panel has no month variable.")
dt[, month := as.IDate(month)]

if (!"post_closure" %in% names(dt)) {
  dt[, post_closure := as.integer(month >= as.IDate("2021-07-01"))]
} else {
  dt[, post_closure := as.integer(post_closure)]
}

required <- c(
  "events_type3",
  "eli_final",
  "pre_contested_share_spatial",
  "pre_taliban_activity_share_spatial",
  "pre_contested_intensity_spatial"
)
missing <- setdiff(required, names(dt))
if (length(missing)) stop("Missing required variables: ", paste(missing, collapse = ", "))

controls <- c("prec_gpcp", "temp", "nlights_calib_mean")
controls <- controls[controls %in% names(dt)]

mods <- c(
  "pre_contested_share_spatial",
  "pre_taliban_activity_share_spatial",
  "pre_contested_intensity_spatial"
)

mod_labels <- c(
  "Pre-contested share",
  "Pre-shock Taliban activity share",
  "Pre-contested intensity"
)

# ------------------------------------------------------------
# 2) ESTIMATION HELPER
# ------------------------------------------------------------
fit_one_h2 <- function(data, mod, sample_label) {

  rhs <- c(
    "post_closure:eli_final",
    paste0("post_closure:eli_final:", mod),
    controls
  )

  fml <- as.formula(
    paste0(
      "events_type3 ~ ",
      paste(rhs, collapse = " + "),
      " | cell_id_num + month"
    )
  )

  m <- fepois(fml, data = data, notes = FALSE)

  # Match the dissertation/canonical inference: cell + month clustering.
  # Fall back to cell clustering only if the two-way VCOV cannot be used.
  warnings_tw <- character()
  ct <- withCallingHandlers(
    tryCatch(
      coeftable(m, vcov = ~ cell_id_num + month),
      error = function(e) NULL
    ),
    warning = function(w) {
      warnings_tw <<- c(warnings_tw, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  cluster_used <- "cell + month"

  if (is.null(ct) || any(grepl("not positive definite", warnings_tw, ignore.case = TRUE))) {
    ct <- coeftable(m, vcov = ~ cell_id_num)
    cluster_used <- "cell fallback"
  }

  ct <- as.data.table(ct, keep.rownames = "term")
  pcol <- grep("^Pr\\(", names(ct), value = TRUE)[1]

  target <- paste0("post_closure:eli_final:", mod)
  row <- ct[term == target]

  if (nrow(row) != 1L) {
    stop(
      "Could not uniquely find H2 term: ", target,
      "\nAvailable terms:\n", paste(ct$term, collapse = "\n")
    )
  }

  data.table(
    sample = sample_label,
    moderator = mod,
    moderator_label = mod_labels[match(mod, mods)],
    term = target,
    estimate = as.numeric(row[["Estimate"]][1]),
    std_error = as.numeric(row[["Std. Error"]][1]),
    p_value = as.numeric(row[[pcol]][1]),
    nobs = as.integer(nobs(m)),
    cluster = cluster_used
  )
}

# ------------------------------------------------------------
# 3) FULL-SAMPLE REPLICATION OF ORIGINAL H2
# ------------------------------------------------------------
full_results <- rbindlist(
  lapply(mods, function(mod) fit_one_h2(dt, mod, "Full sample")),
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------
# 4) THOMAS TIMING ROBUSTNESS:
#    keep pre-Doha OR post-closure; drop Mar 2020-Jun 2021
# ------------------------------------------------------------
gap <- dt[
  month < as.IDate("2020-03-01") |
    month >= as.IDate("2021-07-01")
]

gap_results <- rbindlist(
  lapply(mods, function(mod) fit_one_h2(gap, mod, "Drop Mar 2020-Jun 2021")),
  use.names = TRUE,
  fill = TRUE
)

results <- rbindlist(
  list(full_results, gap_results),
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------
# 5) REPLICATION CHECK
# ------------------------------------------------------------
expected <- data.table(
  moderator = mods,
  expected_main = c(
    -16.1023802013192,
    -17.6667522797993,
    -3.12067660696125
  )
)

rep_check <- merge(
  full_results[, .(moderator, estimate)],
  expected,
  by = "moderator",
  all.x = TRUE
)
rep_check[, abs_difference := abs(estimate - expected_main)]

cat("\n============================================================\n")
cat("ORIGINAL H2 — FULL-SAMPLE REPLICATION CHECK\n")
cat("============================================================\n")
print(rep_check)

cat("\n============================================================\n")
cat("ORIGINAL H2 — ANTICIPATION-WINDOW RESULTS\n")
cat("Drop March 2020 through June 2021\n")
cat("============================================================\n")
print(gap_results)

fwrite(results, OUT_PATH)
cat("\nSaved results to:\n", OUT_PATH, "\n")
