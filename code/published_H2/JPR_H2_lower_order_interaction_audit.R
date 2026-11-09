# ============================================================================
# JPR H2 SPECIFICATION FORENSICS + PAIRED PPML RE-ESTIMATION
# Purpose: determine if the FINAL H2 specifications include Post x moderator;
#          compare the original model with a full hierarchical interaction model.
# No existing analysis files are overwritten. No original R scripts are sourced.
# ============================================================================

ROOT <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
PANEL_FILE <- ""             # Optional: enter an exact .rds/.csv/.parquet panel path.
RUN_REFITS <- TRUE            # Set FALSE to perform only the source-code audit.

# Optional: fill these only if the script cannot unambiguously find the columns.
OUTCOME_COL <- ""            # Usually events_type3
POST_COL <- ""               # Usually shock_bagram or post_bagram
ELI_COL <- ""                # Usually eli_feas
CELL_COL <- ""               # Usually cell_id
MONTH_COL <- ""              # Usually month
MODERATOR_COLS <- c(
  pre_contested_share = "",
  pre_taliban_activity_share = "",
  pre_contested_intensity_log = ""
)
CONTROL_COLS <- NULL          # NULL = use available canonical three controls.

if (!dir.exists(ROOT)) stop("THESIS root not found. Update ROOT at the top of the script.")
OUT <- file.path(ROOT, "outputs", "jpr_h2_lower_order_audit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("\nJPR H2 audit output folder:\n", OUT, "\n", sep = "")

safe_lines <- function(path) tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"),
                                       error = function(e) character())
write_out <- function(obj, filename) write.csv(obj, file.path(OUT, filename),
                                                row.names = FALSE, na = "")

# ---- 1. Audit the actual source code on your drive -------------------------
source_files <- list.files(ROOT, pattern = "\\.(r|rmd|qmd)$", recursive = TRUE,
                           full.names = TRUE, ignore.case = TRUE)
source_files <- source_files[!grepl("[/\\\\](renv|packrat|node_modules|\\.git)[/\\\\]",
                                   source_files, ignore.case = TRUE)]

signature <- "(16[.]102|17[.]667|3[.]121|17[.]408|19[.]071|3[.]389|pre.?contest|taliban.{0,25}(share|activ)|post.{0,45}eli.{0,45}(contest|taliban)|feglm\\s*\\(|fepois\\s*\\(|feols\\s*\\(|H2.{0,25}(model|spec|moderator))"
script_rows <- list(); hit_rows <- list()
for (f in source_files) {
  x <- safe_lines(f)
  if (!length(x)) next
  hits <- grep(signature, x, ignore.case = TRUE, perl = TRUE)
  bn <- basename(f)
  rank <- 3L * length(hits) +
    50L * grepl("H2.*IMPROVEMENT|ANTICIPATION.*WINDOW|H2.*ORIGINAL", bn, ignore.case = TRUE) +
    12L * grepl("H2|FINAL30|ROBUST|SUPERPROPER", f, ignore.case = TRUE) +
    20L * any(grepl("16[.]102|17[.]667|3[.]121", x))
  if (!length(hits) && rank == 0L) next
  script_rows[[length(script_rows) + 1L]] <- data.frame(
    file = f, score = rank, matches = length(hits),
    modified = as.character(file.info(f)$mtime), stringsAsFactors = FALSE
  )
  if (length(hits)) {
    for (i in head(hits, 120L)) {
      hit_rows[[length(hit_rows) + 1L]] <- data.frame(
        file = f, line = i, code = x[i],
        context_before = if (i <= 1L) "" else paste(x[max(1L, i - 2L):(i - 1L)], collapse = " / "),
        context_after = if (i >= length(x)) "" else paste(x[(i + 1L):min(length(x), i + 3L)], collapse = " / "),
        stringsAsFactors = FALSE
      )
    }
  }
}
scripts <- if (length(script_rows)) do.call(rbind, script_rows) else data.frame(
  file = character(), score = integer(), matches = integer(), modified = character())
if (nrow(scripts)) scripts <- scripts[order(-scripts$score, -scripts$matches), ]
write_out(scripts, "01_ranked_H2_source_files.csv")

hits <- if (length(hit_rows)) do.call(rbind, hit_rows) else data.frame(
  file = character(), line = integer(), code = character(),
  context_before = character(), context_after = character())
write_out(hits, "02_H2_source_lines_with_context.csv")
cat("\nTop likely H2 scripts (read these before interpreting the estimates):\n")
print(utils::head(scripts[, c("file", "score", "matches"), drop = FALSE], 12L),
      row.names = FALSE, right = FALSE)
cat("\nSource-code audit written. Do NOT assume a matching file is the final model\n",
    "until its path and output coefficients match the current manuscript.\n", sep = "")

# ---- 2. Inventory candidate saved panels ----------------------------------
data_files <- list.files(ROOT, pattern = "\\.(rds|csv|parquet|fst)$",
                         recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(data_files)) {
  data_files <- data_files[!grepl("[/\\\\](renv|packrat|node_modules|\\.git)[/\\\\]",
                                 data_files, ignore.case = TRUE)]
  info <- file.info(data_files)
  inventory <- data.frame(
    path = data_files,
    size_mb = round(info$size / 1024^2, 2),
    modified = as.character(info$mtime),
    stringsAsFactors = FALSE
  )
  inventory$score <- 10L * grepl("panel", inventory$path, ignore.case = TRUE) +
    10L * grepl("final30|final_airfield|superproper", inventory$path, ignore.case = TRUE) +
    5L * grepl("afghan|bagram|h2", inventory$path, ignore.case = TRUE) -
    10L * grepl("audit|model|figure|bootstrap", inventory$path, ignore.case = TRUE)
  inventory <- inventory[order(-inventory$score, -as.POSIXct(inventory$modified)), ]
} else {
  inventory <- data.frame(path = character(), size_mb = numeric(),
                          modified = character(), score = integer())
}
write_out(inventory, "03_candidate_data_files.csv")
cat("\nTop candidate panel files:\n")
print(utils::head(inventory[, c("path", "size_mb", "score"), drop = FALSE], 18L),
      row.names = FALSE, right = FALSE)

if (!RUN_REFITS) {
  cat("\nRUN_REFITS=FALSE: source inspection completed.\n")
} else {
  if (!nzchar(PANEL_FILE) && interactive() && nrow(inventory)) {
    cat("\nSelect the exact final Afghanistan H2 panel, NOT a model-results RDS.\n")
    top_n <- min(20L, nrow(inventory))
    for (j in seq_len(top_n)) cat(sprintf("%2d: %s\n", j, inventory$path[j]))
    answer <- trimws(readline("Panel number (1-20), F to browse, or Enter to skip: "))
    if (grepl("^[0-9]+$", answer) && as.integer(answer) >= 1L &&
        as.integer(answer) <= top_n) {
      PANEL_FILE <- inventory$path[as.integer(answer)]
    } else if (toupper(answer) == "F") {
      PANEL_FILE <- tryCatch(file.choose(), error = function(e) "")
    }
  }

  if (!nzchar(PANEL_FILE) || !file.exists(PANEL_FILE)) {
    cat("\nNo panel selected. Code and panel inventory are complete.\n",
        "To re-estimate, choose the final panel or set PANEL_FILE and source again.\n", sep = "")
  } else {
    if (!requireNamespace("fixest", quietly = TRUE)) stop("Install fixest: install.packages('fixest')")
    ext <- tolower(tools::file_ext(PANEL_FILE))
    panel <- switch(ext,
      rds = readRDS(PANEL_FILE),
      csv = {
        if (requireNamespace("data.table", quietly = TRUE)) data.table::fread(PANEL_FILE)
        else read.csv(PANEL_FILE, stringsAsFactors = FALSE)
      },
      parquet = {
        if (!requireNamespace("arrow", quietly = TRUE)) stop("Install arrow to read parquet.")
        arrow::read_parquet(PANEL_FILE)
      },
      fst = {
        if (!requireNamespace("fst", quietly = TRUE)) stop("Install fst to read fst.")
        fst::read_fst(PANEL_FILE)
      },
      stop("Unsupported panel extension: ", ext)
    )
    if (!is.data.frame(panel)) stop("Selected file is not a data-frame panel: ", PANEL_FILE)
    nms <- names(panel)
    write_out(data.frame(column = nms, class = vapply(panel, function(z)
      paste(class(z), collapse = "/"), character(1))), "04_selected_panel_columns.csv")
    cat("\nPanel loaded:", PANEL_FILE, "\nRows:", nrow(panel), " Columns:", length(nms), "\n")

    choose_one <- function(manual, alternatives, label) {
      if (nzchar(manual)) {
        if (!manual %in% nms) stop(label, " column not found: ", manual)
        return(manual)
      }
      exact <- alternatives[alternatives %in% nms]
      if (length(exact)) return(exact[1L])
      stop("Could not identify ", label, ". See 04_selected_panel_columns.csv.\n",
           "Enter the correct column in the configuration section and rerun.")
    }
    outcome <- choose_one(OUTCOME_COL, c("events_type3", "ucdp_type3_events", "type3_events"), "outcome")
    post <- choose_one(POST_COL, c("shock_bagram", "post_bagram", "post_closure", "post"), "post")
    eli <- choose_one(ELI_COL, c("eli_feas", "eli_feasibility", "ELI", "eli"), "ELI")
    cell <- choose_one(CELL_COL, c("cell_id", "pg_id", "priogrid_gid"), "cell ID")
    month <- choose_one(MONTH_COL, c("month", "month_date", "month_id"), "month")

    if ("theatre_id" %in% names(panel)) {
      af <- panel[!is.na(panel$theatre_id) &
                  tolower(as.character(panel$theatre_id)) == "afghanistan", , drop = FALSE]
      if (nrow(af)) panel <- af
    }

    # Avoid silently choosing between similarly named moderators: exact aliases
    # have priority. Unclear matches are reported and skipped.
    aliases <- list(
      pre_contested_share = c("pre_contested_share", "pre_contested_month_share", "pre_contestedness_share", "h2_pre_contested_share"),
      pre_taliban_activity_share = c("pre_taliban_activity_share", "pre_taliban_active_share", "pre_taliban_share", "h2_pre_taliban_activity_share"),
      pre_contested_intensity_log = c("pre_contested_intensity_log", "pre_contested_log_intensity", "log_pre_contested_intensity", "pre_contested_intensity_logged")
    )
    patterns <- c(pre_contested_share = "contest.*share|share.*contest",
                  pre_taliban_activity_share = "taliban.*(share|activ)|(share|activ).*taliban",
                  pre_contested_intensity_log = "contest.*(intens|log)|(intens|log).*contest")
    mods <- MODERATOR_COLS
    for (j in names(mods)) {
      if (nzchar(mods[j])) next
      exact <- aliases[[j]][aliases[[j]] %in% nms]
      if (length(exact)) {
        mods[j] <- exact[1L]
      } else {
        maybe <- grep(patterns[j], nms, value = TRUE, ignore.case = TRUE, perl = TRUE)
        if (length(maybe) == 1L) mods[j] <- maybe
        else cat("\nModerator", j, "is ambiguous/missing. Candidates:",
                 paste(maybe, collapse = ", "), "\n")
      }
    }
    mods <- mods[nzchar(mods)]
    write_out(data.frame(moderator_type = names(MODERATOR_COLS),
                         selected_column = unname(vapply(names(MODERATOR_COLS),
                           function(k) if (k %in% names(mods)) mods[[k]] else "", character(1)))),
              "05_H2_variable_mapping.csv")
    cat("\nSelected H2 moderators:\n"); print(mods)
    if (!length(mods)) stop("No moderator identified. Fill MODERATOR_COLS using the panel column inventory.")

    canonical_controls <- c("prec_gpcp", "temp", "nlights_calib_mean")
    controls <- if (is.null(CONTROL_COLS)) intersect(canonical_controls, names(panel)) else CONTROL_COLS
    if (!all(controls %in% names(panel))) stop("One or more CONTROL_COLS not found.")
    if (length(controls) != 3L) warning("Fewer than 3 canonical controls found: check exact FINAL model controls.")

    # Quoting protects column names that contain dots/spaces; fixest interaction
    # notation is used only for safe syntactic names.
    needed <- unique(c(outcome, post, eli, cell, month, controls, unname(mods)))
    if (any(make.names(needed) != needed)) warning("Unusual column names: inspect formulas.")
    panel <- as.data.frame(panel)
    if (inherits(panel[[month]], "factor")) panel[[month]] <- as.character(panel[[month]])
    month_values <- suppressWarnings(as.Date(as.character(panel[[month]])))
    if (all(is.na(month_values))) {
      month_values <- tryCatch(suppressWarnings(as.Date(panel[[month]])), error = function(e) rep(as.Date(NA), nrow(panel)))
    }
    if (all(is.na(month_values))) warning("Could not parse month as date. Anticipation-window check skipped.")

    model_rows <- list()
    for (j in names(mods)) {
      mod <- mods[[j]]
      for (window in c("full", "exclude_Mar2020_Jun2021")) {
        if (window != "full" && all(is.na(month_values))) next
        keep <- rep(TRUE, nrow(panel))
        if (window != "full") keep <- !is.na(month_values) &
          (month_values < as.Date("2020-03-01") | month_values >= as.Date("2021-07-01"))
        cols <- unique(c(outcome, post, eli, cell, month, mod, controls))
        keep <- keep & stats::complete.cases(panel[, cols, drop = FALSE])
        dat <- panel[keep, , drop = FALSE]
        if (!nrow(dat)) next
        if (length(unique(dat[[post]])) < 2) {
          warning("Only one post value in ", j, ", ", window, ". Skipping.")
          next
        }
        triple <- paste(post, eli, mod, sep = ":")
        eli_term <- paste(post, eli, sep = ":")
        post_mod <- paste(post, mod, sep = ":")
        rhs_old <- c(eli_term, triple, controls)
        rhs_new <- c(eli_term, post_mod, triple, controls)
        form_old <- stats::as.formula(paste(outcome, "~", paste(rhs_old, collapse = " + "),
                                             "|", cell, "+", month))
        form_new <- stats::as.formula(paste(outcome, "~", paste(rhs_new, collapse = " + "),
                                             "|", cell, "+", month))
        two_way <- stats::as.formula(paste("~", cell, "+", month))
        for (kind in c("reported_structure", "add_post_x_moderator")) {
          fm <- if (kind == "reported_structure") form_old else form_new
          err <- ""
          fit <- tryCatch(fixest::feglm(fm, data = dat, family = "poisson"),
                          error = function(e) {err <<- conditionMessage(e); NULL})
          if (is.null(fit)) {
            model_rows[[length(model_rows) + 1L]] <- data.frame(
              moderator = j, moderator_column = mod, window = window, specification = kind,
              term = triple, estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
              nobs = NA_integer_, candidate_rows = nrow(dat),
              dropped_by_collinearity = "", error = err, formula = paste(deparse(fm),collapse=""))
            next
          }
          ct <- tryCatch(fixest::coeftable(fit, vcov = two_way),
                         error = function(e) {err <<- conditionMessage(e); NULL})
          if (is.null(ct)) {
            model_rows[[length(model_rows) + 1L]] <- data.frame(
              moderator = j, moderator_column = mod, window = window, specification = kind,
              term = triple, estimate = NA_real_, std_error = NA_real_, p_value = NA_real_,
              nobs = stats::nobs(fit), candidate_rows = nrow(dat),
              dropped_by_collinearity = paste(fit$collin.var, collapse = ";"),
              error = err, formula = paste(deparse(fm),collapse=""))
            next
          }
          for (term in c(triple, post_mod)) {
            cn <- rownames(ct)
            tokens <- strsplit(cn, ":", fixed = TRUE)
            wanted <- strsplit(term, ":", fixed = TRUE)[[1L]]
            match_ix <- which(vapply(tokens,
              function(v) length(v) == length(wanted) && setequal(v, wanted), logical(1)))
            estimate <- se <- pvalue <- NA_real_
            if (length(match_ix) == 1L) {
              estimate <- unname(ct[match_ix, 1]); se <- unname(ct[match_ix, 2]);
              pvalue <- unname(ct[match_ix, ncol(ct)])
            }
            model_rows[[length(model_rows) + 1L]] <- data.frame(
              moderator = j, moderator_column = mod, window = window, specification = kind,
              term = term, estimate = estimate, std_error = se, p_value = pvalue,
              nobs = stats::nobs(fit), candidate_rows = nrow(dat),
              dropped_by_collinearity = paste(fit$collin.var, collapse = ";"),
              error = err, formula = paste(deparse(fm), collapse = ""))
          }
        }
      }
    }
    model_compare <- if (length(model_rows)) do.call(rbind, model_rows) else data.frame()
    write_out(model_compare, "06_paired_H2_PPML_results.csv")
    expected <- data.frame(
      moderator = rep(c("pre_contested_share", "pre_taliban_activity_share", "pre_contested_intensity_log"), 2),
      window = rep(c("full", "exclude_Mar2020_Jun2021"), each = 3),
      manuscript_estimate = c(-16.102, -17.667, -3.121, -17.408, -19.071, -3.389),
      manuscript_N = rep(c(54312L, 49956L), each = 3),
      stringsAsFactors = FALSE
    )
    if (nrow(model_compare)) {
      baselines <- model_compare[model_compare$specification == "reported_structure" &
                                 grepl(":.*:", model_compare$term),
                                 c("moderator", "window", "estimate", "nobs"), drop = FALSE]
      validation <- merge(expected, baselines, by = c("moderator", "window"), all.x = TRUE)
      validation$matches_manuscript <- !is.na(validation$estimate) &
        !is.na(validation$nobs) & (abs(validation$estimate - validation$manuscript_estimate) <= 0.005) &
        (validation$nobs == validation$manuscript_N)
      write_out(validation, "07_baseline_reproduction_check.csv")
      cat("\nDoes the baseline replicate the manuscript?\n")
      print(validation[, c("moderator", "window", "matches_manuscript"), drop = FALSE], row.names = FALSE)
    }
    cat("\nPAIRED H2 RESULTS (inspect N and baseline coefficients first):\n")
    print(model_compare[, intersect(c("moderator", "window", "specification", "term",
                                     "estimate", "std_error", "p_value", "nobs"),
                                   names(model_compare)), drop = FALSE], row.names = FALSE)
    cat("\nImportant: if the reported-structure coefficients and N do NOT reproduce\n",
        "Table 2 / Table 4, this is not the canonical final panel/specification.\n",
        "Find the exact final H2 R script and adjust control/sample definitions first.\n", sep = "")
  }
}
cat("\nDone. Send 01_ranked_H2_source_files.csv, 02_H2_source_lines_with_context.csv,\n",
    "05_H2_variable_mapping.csv, 06_paired_H2_PPML_results.csv and\n",
    "07_baseline_reproduction_check.csv (if created).\n")
