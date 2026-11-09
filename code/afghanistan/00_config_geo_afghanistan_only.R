# ============================================================
# 00_config_geo_afghanistan_only_EXECUTABLE_ACLED.R
# SELF-CONTAINED config + helpers for the executable Afghanistan branch.
# Does NOT depend on any previous thesis config script.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -------------------------
# 0) ROOT PATHS
# -------------------------
BASE_DIR <- "C:/Users/USERNAME/OneDrive/Desktop/THESIS"

DIRS <- list(
  data_raw = file.path(BASE_DIR, "data", "raw"),
  data_processed = file.path(BASE_DIR, "data", "processed"),
  base = file.path(BASE_DIR, "outputs", "paper1_geo_afghanistan_executable_acled"),
  panels = file.path(BASE_DIR, "outputs", "paper1_geo_afghanistan_executable_acled", "panels"),
  tables = file.path(BASE_DIR, "outputs", "paper1_geo_afghanistan_executable_acled", "tables"),
  audit = file.path(BASE_DIR, "outputs", "paper1_geo_afghanistan_executable_acled", "audit")
)
for (p in DIRS) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# 1) GLOBAL SETTINGS
# -------------------------
TARGET_THEATRE <- "afghanistan"
SHOCK_DATE <- as.IDate("2021-07-01")
SHOCK_VAR <- "shock_bagram"
SHOCK_LABEL <- "Bagram closure"
SHOCK_SITE_LAT <- 34.9461
SHOCK_SITE_LON <- 69.2647
DONUT_RADIUS_KM <- 50
EVENT_STUDY_MIN <- -12L
EVENT_STUDY_MAX <- 12L
PLACEBO_DATE <- as.IDate("2020-07-01")
H3_ADJ_KM <- 50
H3_LOCAL_MAX_KM <- 100
H3_LOCAL_MIN_FRINGE_KM <- 50
LOW_VIS_MAIN_Q <- 0.25
LOW_VIS_STRICT_Q <- 0.10

# Backbone panel used as the Afghanistan source panel
FOUNDATION_PATH <- file.path(BASE_DIR, "data", "processed", "panels", "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds")
if (!file.exists(FOUNDATION_PATH)) {
  FOUNDATION_PATH <- file.path(BASE_DIR, "outputs", "paper1", "panels", "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds")
}

# Raw-data expectations
AFG_DISTRICT_BOUNDARIES_PATH <- file.path(BASE_DIR, "data", "raw", "afg_admin", "afg_admin2.shp")
AFG_DISTRICT_NAME_CROSSWALK_PATH <- file.path(BASE_DIR, "data", "raw", "afg_admin", "afg_district_name_crosswalk.csv")
RAW_ACLED_PATTERNS <- c("acled.*\\.csv$")
RAW_GED_PATTERNS <- c("GEDEvent.*\\.csv$", "ged.*\\.csv$")

DISTRICT_ID_CANDIDATES <- c("adm2_pcode", "adm2_code", "dist_id", "district_id", "adm2_code_en")
DISTRICT_NAME_CANDIDATES <- c("adm2_en", "adm2_name", "district", "district_name", "name_2")
PROVINCE_NAME_CANDIDATES <- c("adm1_en", "adm1_name", "province", "province_name", "name_1")

ACLED_DISTRICT_CANDIDATES <- c("admin2", "district", "district_name", "adm2_en", "adm2_name")
ACLED_PROVINCE_CANDIDATES <- c("admin1", "province", "province_name", "adm1_en", "adm1_name")
ACLED_DATE_CANDIDATES <- c("event_date", "eventdate", "date", "eventdate_date")
ACLED_EVENT_CANDIDATES <- c("event_type")
ACLED_SUBEVENT_CANDIDATES <- c("sub_event_type")
ACLED_ACTOR_CANDIDATES <- c("actor1", "actor2", "assoc_actor_1", "assoc_actor_2")

# -------------------------
# 2) HELPERS
# -------------------------
`%notin%` <- function(x, table) !(x %in% table)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

safe_name <- function(x) {
  x <- trimws(as.character(x))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

normalize_text <- function(x) {
  x <- trimws(as.character(x))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

pick_col <- function(nms, candidates, required = TRUE) {
  sn <- safe_name(nms)
  sc <- safe_name(candidates)
  m <- match(sc, sn)
  m <- m[!is.na(m)]
  if (!length(m)) {
    if (required) stop("Could not find required column among candidates: ", paste(candidates, collapse = ", "))
    return(NA_character_)
  }
  nms[m[1]]
}

clean_theatre <- function(x) normalize_text(x)

coerce_month <- function(x) {
  if (inherits(x, "IDate")) return(as.IDate(format(x, "%Y-%m-01")))
  if (inherits(x, "Date")) return(as.IDate(format(x, "%Y-%m-01")))
  if (is.numeric(x)) {
    # support YYYYMM integers
    xs <- as.character(as.integer(x))
    ok <- grepl("^[0-9]{6}$", xs)
    out <- rep(as.IDate(NA), length(xs))
    out[ok] <- as.IDate(paste0(substr(xs[ok],1,4), "-", substr(xs[ok],5,6), "-01"))
    return(out)
  }
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  # try YYYY-MM-DD / YYYY-MM / MM/DD/YYYY
  out <- suppressWarnings(as.IDate(x))
  bad <- is.na(out) & !is.na(x)
  if (any(bad)) {
    # YYYY-MM
    ok2 <- grepl("^[0-9]{4}-[0-9]{2}$", x[bad])
    out[bad][ok2] <- as.IDate(paste0(x[bad][ok2], "-01"))
    bad2 <- is.na(out) & !is.na(x)
    if (any(bad2)) {
      parsed <- suppressWarnings(as.Date(x[bad2], format = "%m/%d/%Y"))
      out[bad2] <- as.IDate(parsed)
    }
  }
  as.IDate(format(out, "%Y-%m-01"))
}

month_int_from_idate <- function(x) {
  x <- as.IDate(x)
  as.integer(format(x, "%Y")) * 12L + as.integer(format(x, "%m"))
}

assert_unique_keys <- function(dt, keys, name = "data") {
  if (!is.data.table(dt)) dt <- as.data.table(dt)
  dup <- dt[, .N, by = keys][N > 1L]
  if (nrow(dup)) stop("Duplicate keys in ", name, ": ", paste(keys, collapse = ", "))
  invisible(TRUE)
}

assert_cols <- function(dt, cols, name = "data") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Missing columns in ", name, ": ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

read_any_dt <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(as.data.table(readRDS(path)))
  if (ext %in% c("csv", "txt", "tsv")) return(as.data.table(fread(path, showProgress = FALSE)))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read: ", path)
    return(as.data.table(readxl::read_excel(path)))
  }
  stop("Unsupported file type: ", path)
}

find_first_recursive <- function(root, patterns) {
  if (!dir.exists(root)) return(NA_character_)
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  for (p in patterns) {
    hit <- files[grepl(p, basename(files), ignore.case = TRUE)]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371.0
  to_rad <- pi / 180
  lat1 <- lat1 * to_rad; lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad; lon2 <- lon2 * to_rad
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

compute_pg_id <- function(lat, lon) {
  lat <- as.numeric(lat)
  lon <- as.numeric(lon)
  out <- rep(NA_integer_, length(lat))
  ok <- !is.na(lat) & !is.na(lon) & lat >= -90 & lat < 90 & lon >= -180 & lon < 180
  if (any(ok)) {
    lat_i <- floor((lat[ok] + 90) / 0.5)
    lon_i <- floor((lon[ok] + 180) / 0.5)
    out[ok] <- as.integer(lat_i * 720L + lon_i + 1L)
  }
  out
}

write_audit_text <- function(lines, filename_stub) {
  out <- file.path(DIRS$audit, paste0(filename_stub, ".txt"))
  writeLines(as.character(lines), out)
  invisible(out)
}

safe_merge <- function(x, y, by, all.x = TRUE, sort = FALSE, name_x = "x", name_y = "y") {
  out <- merge(x, y, by = by, all.x = all.x, sort = sort)
  assert_unique_keys(out, by, paste0("merged_", name_x, "_", name_y))
  out
}

CAUSAL_UPGRADE_TAG <- "Afghanistan executable ACLED-based H2 + local-cluster H3 + protests negative controls (self-contained)"
