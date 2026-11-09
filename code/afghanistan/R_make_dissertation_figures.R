# ============================================================
# R_make_dissertation_figures_FINAL.R
# Final polished dissertation figures
# Bagram / Afghanistan ISR design
# ============================================================

# -------------------------
# 0) SETUP
# -------------------------

THESIS_DIR <- Sys.getenv(
  "THESIS_DIR",
  unset = "C:/Users/USERNAME/OneDrive/Desktop/THESIS"
)

OUT_DIR <- file.path(THESIS_DIR, "outputs", "dissertation_figures_FINAL")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

needed <- c(
  "data.table", "ggplot2", "sf", "rnaturalearth",
  "rnaturalearthdata", "scales"
)

to_install <- needed[!needed %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(scales)
})

theme_set(theme_minimal(base_size = 12))

save_fig <- function(plot, filename, width = 8.5, height = 5.5) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

pick_var <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit)) hit[1] else NA_character_
}

first_non_na <- function(x) {
  idx <- which(!is.na(x))
  if (!length(idx)) {
    if (is.character(x)) return(NA_character_)
    if (is.integer(x)) return(NA_integer_)
    if (is.numeric(x)) return(NA_real_)
    if (is.logical(x)) return(NA)
    return(NA)
  }
  x[idx[1]]
}

cat("\nWriting final figures to:\n", OUT_DIR, "\n")

# -------------------------
# 1) LOAD PANEL
# -------------------------

panel_candidates <- c(
  file.path(THESIS_DIR, "outputs", "paper1_geo_afghanistan_executable_acled", "panels", "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds"),
  file.path(THESIS_DIR, "outputs", "paper1_geo_afghanistan_executable_acled", "panels", "panel_core_cell_month.rds"),
  file.path(THESIS_DIR, "data", "processed", "panels", "panel_pg_month_ALLCW_paper1_posteli_v2_H2H3.rds"),
  file.path(THESIS_DIR, "data", "processed", "panels", "panel_core_cell_month.rds")
)

panel_path <- first_existing(panel_candidates)

if (is.na(panel_path)) {
  stop(
    "Could not find the main panel. Edit panel_candidates in this script.\nLooked in:\n",
    paste(panel_candidates, collapse = "\n")
  )
}

cat("\nUsing panel:\n", panel_path, "\n")

dt <- as.data.table(readRDS(panel_path))

if (!"month" %in% names(dt)) stop("Panel must contain month.")
dt[, month := as.Date(month)]

required_core <- c("cell_id", "month", "lat", "lon", "eli_feas")
missing_core <- setdiff(required_core, names(dt))
if (length(missing_core)) {
  stop("Panel is missing required columns: ", paste(missing_core, collapse = ", "))
}

cell_cols <- c(
  "cell_id", "lat", "lon", "eli_feas",
  "pre_contested_share", "pre_contested_share_adm2",
  "adm2_pre_contested_share", "pre_taliban_activity_share",
  "pre_contested_intensity"
)
cell_cols <- intersect(cell_cols, names(dt))

cells <- dt[
  ,
  lapply(.SD, first_non_na),
  by = cell_id,
  .SDcols = setdiff(cell_cols, "cell_id")
]

cells <- cells[!is.na(lat) & !is.na(lon) & !is.na(eli_feas)]

# -------------------------
# 2) BASEMAP AND AIRBASES
# -------------------------

afg <- rnaturalearth::ne_countries(
  country = "Afghanistan",
  scale = "medium",
  returnclass = "sf"
)

airbases <- data.table(
  name = c("Bagram", "Kandahar", "Shindand"),
  lon = c(69.264, 65.847, 62.261),
  lat = c(34.946, 31.506, 33.392)
)

# -------------------------
# 3) PALETTE
# -------------------------

COL_HIGH   <- "#2F4B7C"  # muted navy
COL_LOW    <- "#E8E1D5"  # light beige
COL_ACCENT <- "#A23E3E"  # muted red
COL_MID    <- "#7EAAC8"  # muted blue
COL_LIGHT  <- "#DDEAF6"  # pale blue
COL_DARK   <- "#25364A"  # dark navy

# ============================================================
# FIGURE 1: TIMELINE
# ============================================================

timeline <- data.table(
  date = as.Date(c("2020-02-29", "2021-04-14", "2021-07-02", "2021-08-15")),
  event = c(
    "Doha Agreement",
    "Biden announcement",
    "Bagram closure",
    "Taliban takeover"
  ),
  note = c(
    "anticipation benchmark",
    "withdrawal confirmed",
    "main treatment date",
    "state collapse"
  ),
  type = c("anticipation", "context", "treatment", "context"),
  y = c(1.20, 0.80, 1.20, 0.80)
)

fig1 <- ggplot(timeline, aes(x = date)) +
  geom_hline(yintercept = 1, linewidth = 0.6, color = "grey25") +
  geom_segment(aes(xend = date, y = 1, yend = y), linewidth = 0.45, color = "grey35") +
  geom_point(aes(y = 1, color = type), size = 3.5) +
  geom_label(
    aes(y = y, label = paste0(event, "\n", note)),
    size = 3.2,
    label.size = 0.22,
    label.padding = grid::unit(0.18, "lines"),
    fill = "white",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "anticipation" = COL_MID,
      "treatment" = COL_ACCENT,
      "context" = "grey35"
    ),
    guide = "none"
  ) +
  scale_x_date(
    date_breaks = "3 months",
    date_labels = "%b\n%Y",
    limits = as.Date(c("2020-01-15", "2021-09-15"))
  ) +
  coord_cartesian(ylim = c(0.60, 1.40), clip = "off") +
  labs(
    title = "Withdrawal timeline and treatment timing",
    subtitle = "Bagram closure is the main treatment date; Doha is the anticipation benchmark.",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 25, 10, 10)
  )

save_fig(fig1, "Figure_1_FINAL_withdrawal_timeline.png", width = 9, height = 4.5)

# ============================================================
# FIGURE 2: HIGH/LOW ELI MAP
# ============================================================

eli_med <- median(cells$eli_feas, na.rm = TRUE)

cells[, eli_group := fifelse(
  eli_feas >= eli_med,
  "Higher ISR feasibility",
  "Lower ISR feasibility"
)]

cells[, eli_group := factor(
  eli_group,
  levels = c("Lower ISR feasibility", "Higher ISR feasibility")
)]

fig2_audit <- cells[
  ,
  .(
    cells = .N,
    min_eli = min(eli_feas, na.rm = TRUE),
    median_eli = median(eli_feas, na.rm = TRUE),
    max_eli = max(eli_feas, na.rm = TRUE)
  ),
  by = eli_group
]

cat("\n================ FIGURE 2 ELI MAP AUDIT ================\n")
print(fig2_audit)

bad_high <- cells[eli_group == "Higher ISR feasibility" & eli_feas < eli_med]
bad_low  <- cells[eli_group == "Lower ISR feasibility" & eli_feas >= eli_med]

if (nrow(bad_high) > 0 | nrow(bad_low) > 0) {
  stop("Figure 2 classification error: some cells are assigned to the wrong ELI group.")
} else {
  cat("\nFigure 2 audit passed: darker cells are exactly cells with eli_feas >= theatre median.\n")
}

fwrite(fig2_audit, file.path(OUT_DIR, "Figure_2_ELI_map_audit.csv"))

fig2 <- ggplot() +
  geom_sf(data = afg, fill = "grey97", color = "grey50", linewidth = 0.4) +
  geom_tile(
    data = cells,
    aes(x = lon, y = lat, fill = eli_group),
    width = 0.48,
    height = 0.48,
    color = "white",
    linewidth = 0.10,
    alpha = 0.92
  ) +
  geom_point(
    data = airbases,
    aes(x = lon, y = lat),
    size = 2.8,
    shape = 21,
    fill = "white",
    color = COL_ACCENT,
    stroke = 1
  ) +
  geom_label(
    data = airbases,
    aes(x = lon, y = lat, label = name),
    nudge_y = 0.35,
    size = 3,
    label.size = 0.15,
    fill = "white",
    label.padding = grid::unit(0.12, "lines")
  ) +
  scale_fill_manual(
    values = c(
      "Lower ISR feasibility" = COL_LOW,
      "Higher ISR feasibility" = COL_HIGH
    ),
    name = NULL
  ) +
  coord_sf(
    xlim = range(cells$lon, na.rm = TRUE) + c(-0.7, 0.7),
    ylim = range(cells$lat, na.rm = TRUE) + c(-0.7, 0.7),
    expand = FALSE
  ) +
  labs(
    title = "High and low predetermined ISR feasibility across Afghanistan",
    subtitle = "Darker cells are above the theatre median of ELI feasibility.",
    x = NULL,
    y = NULL,
    caption = "ELI combines terrain observability and distance-decayed proximity to relevant U.S. or coalition airpower hubs. White areas are outside the analysed PRIO-GRID cell set."
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    plot.caption = element_text(size = 8.8, color = "grey35")
  )

save_fig(fig2, "Figure_2_FINAL_high_low_eli_map.png", width = 8.8, height = 6.5)

# ============================================================
# FIGURE 3: H2 MODERATOR DISTRIBUTION
# ============================================================

h2_var <- pick_var(
  names(cells),
  c(
    "pre_contested_share",
    "pre_contested_share_adm2",
    "adm2_pre_contested_share"
  )
)

if (is.na(h2_var)) {
  warning("No H2 pre-contested share variable found. Figure 3 skipped.")
} else {
  
  h2_data <- cells[!is.na(get(h2_var))]
  h2_data[, h2_share := get(h2_var)]
  
  fig3_audit <- h2_data[
    ,
    .(
      cells = .N,
      min = min(h2_share, na.rm = TRUE),
      p25 = quantile(h2_share, 0.25, na.rm = TRUE),
      median = median(h2_share, na.rm = TRUE),
      mean = mean(h2_share, na.rm = TRUE),
      p75 = quantile(h2_share, 0.75, na.rm = TRUE),
      max = max(h2_share, na.rm = TRUE),
      zero_cells = sum(h2_share == 0, na.rm = TRUE),
      zero_share = mean(h2_share == 0, na.rm = TRUE)
    )
  ]
  
  cat("\n================ FIGURE 3 H2 MODERATOR AUDIT ================\n")
  cat("H2 variable used:", h2_var, "\n")
  print(fig3_audit)
  
  fwrite(fig3_audit, file.path(OUT_DIR, "Figure_3_H2_moderator_audit.csv"))
  
  fig3 <- ggplot(h2_data, aes(x = h2_share)) +
    geom_histogram(
      bins = 25,
      fill = COL_MID,
      color = "white",
      linewidth = 0.25
    ) +
    geom_vline(
      xintercept = median(h2_data$h2_share, na.rm = TRUE),
      linetype = "dashed",
      linewidth = 0.55,
      color = COL_ACCENT
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      breaks = pretty_breaks(n = 6)
    ) +
    labs(
      title = "Pre-shock contestedness used in H2",
      subtitle = "Distribution of the ADM2 pre-shock contested-share moderator.",
      x = "Pre-shock contested share",
      y = "Number of PRIO-GRID cells",
      caption = "Dashed line marks the median cell value. The outcome remains UCDP-GED type-3; ACLED is used only to construct the pre-shock moderator."
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 8.8, color = "grey35")
    )
  
  save_fig(fig3, "Figure_3_FINAL_h2_moderator_distribution.png", width = 8.5, height = 5.2)
}

# ============================================================
# FIGURE 4: COEFFICIENT SUMMARY
# ============================================================

coef_dt <- data.table(
  group = c(
    "H1 and timing",
    "H1 and timing",
    "H1 and timing",
    "H2 mechanism",
    "H2 mechanism",
    "H2 mechanism",
    "Checks",
    "Checks",
    "Checks"
  ),
  model = c(
    "H1 main",
    "H1 drop month +6",
    "Doha timing H1",
    "H2 pre-contested share",
    "H2 Taliban-active share",
    "H2 contested intensity",
    "ACLED protests",
    "UCDP type-1 battles",
    "Low-visibility terrain"
  ),
  estimate = c(
    -3.055,
    -3.132,
    -2.082,
    -15.670,
    -17.062,
    -3.058,
    0.835,
    -0.744,
    1.548
  ),
  se = c(
    1.245,
    1.254,
    0.957,
    4.218,
    4.501,
    0.883,
    1.108,
    1.113,
    0.659
  )
)

coef_dt[, `:=`(
  ci_low = estimate - 1.96 * se,
  ci_high = estimate + 1.96 * se
)]

coef_dt[, group := factor(
  group,
  levels = c("H1 and timing", "H2 mechanism", "Checks")
)]

coef_dt[, model := factor(model, levels = rev(model))]

cat("\n================ FIGURE 4 COEFFICIENT AUDIT ================\n")
print(coef_dt)
fwrite(coef_dt, file.path(OUT_DIR, "Figure_4_coefficient_summary_data.csv"))

fig4 <- ggplot(coef_dt, aes(x = estimate, y = model, color = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, color = "grey45") +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.16,
    linewidth = 0.5
  ) +
  geom_point(size = 2.4) +
  facet_wrap(~ group, ncol = 1, scales = "free") +
  scale_color_manual(
    values = c(
      "H1 and timing" = COL_HIGH,
      "H2 mechanism" = COL_ACCENT,
      "Checks" = "grey35"
    ),
    guide = "none"
  ) +
  labs(
    title = "Main estimates and diagnostic checks",
    subtitle = "Points show PPML estimates; bars show approximate 95% confidence intervals.",
    x = "Coefficient estimate",
    y = NULL,
    caption = "H2 coefficients are triple interactions and use a different scale from H1 and diagnostic checks."
  ) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.caption = element_text(size = 8.8, color = "grey35")
  )

save_fig(fig4, "Figure_4_FINAL_coefficient_summary.png", width = 8.8, height = 7.0)

# ============================================================
# FIGURE 5: H3 SCHEMATIC, FIXED LABELS
# ============================================================

grid <- CJ(x = 1:7, y = 1:7)
grid[, zone := "Outside H3 sample"]

# 2x2 high-ELI source cluster
grid[x %in% 4:5 & y %in% 4:5, zone := "High-ELI source cluster"]

# Adjacent ring around source cluster
grid[
  x %in% 3:6 & y %in% 3:6 &
    !(x %in% 4:5 & y %in% 4:5),
  zone := "Adjacent low-ELI cells, 0–50 km"
]

# Farther fringe ring
grid[
  x %in% 2:7 & y %in% 2:7 &
    !(x %in% 3:6 & y %in% 3:6),
  zone := "Fringe low-ELI cells, 50–100 km"
]

plot_grid <- grid[zone != "Outside H3 sample"]

plot_grid[, zone := factor(
  zone,
  levels = c(
    "Fringe low-ELI cells, 50–100 km",
    "Adjacent low-ELI cells, 0–50 km",
    "High-ELI source cluster"
  )
)]

fig5 <- ggplot(plot_grid, aes(x = x, y = y, fill = zone)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.98, height = 0.98) +
  
  # Centre label sits fully inside the dark source cluster
  annotate(
    "label",
    x = 4.5,
    y = 4.5,
    label = "High-ELI\nsource cluster\n(anchor)",
    size = 3.4,
    fontface = "bold",
    fill = "white",
    label.size = 0.2
  ) +
  
  # Callout labels placed outside the coloured cells
  annotate(
    "label",
    x = 8.15,
    y = 5.75,
    label = "Fringe low-ELI cells\n50–100 km\ncomparison group",
    size = 3.2,
    fill = "white",
    label.size = 0.2
  ) +
  annotate(
    "segment",
    x = 7.45,
    xend = 6.55,
    y = 5.75,
    yend = 6.15,
    arrow = arrow(length = grid::unit(0.13, "inches")),
    linewidth = 0.45
  ) +
  
  annotate(
    "label",
    x = 8.15,
    y = 3.15,
    label = "Adjacent low-ELI cells\n0–50 km\nH3 group",
    size = 3.2,
    fill = "white",
    label.size = 0.2
  ) +
  annotate(
    "segment",
    x = 7.45,
    xend = 6.05,
    y = 3.15,
    yend = 3.45,
    arrow = arrow(length = grid::unit(0.13, "inches")),
    linewidth = 0.45
  ) +
  
  annotate(
    "text",
    x = 4.5,
    y = 1.45,
    label = "H3 compares the adjacent low-ELI ring with the farther low-ELI fringe,\nwithin the same source-cluster month.",
    size = 3.2
  ) +
  
  scale_fill_manual(
    values = c(
      "Fringe low-ELI cells, 50–100 km" = COL_LIGHT,
      "Adjacent low-ELI cells, 0–50 km" = COL_MID,
      "High-ELI source cluster" = COL_DARK
    ),
    name = NULL
  ) +
  coord_equal(xlim = c(1.65, 9.25), ylim = c(1.25, 7.35), expand = FALSE, clip = "off") +
  labs(
    title = "Illustrative schematic of the H3 local design",
    subtitle = "H3 compares nearby and farther low-ELI cells around the same high-ELI cluster.",
    x = NULL,
    y = NULL,
    caption = "Schematic only. It explains the comparison structure and is not drawn to geographic scale."
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(hjust = 0, margin = margin(b = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 9.5),
    plot.caption = element_text(size = 8.8, color = "grey35"),
    plot.margin = margin(12, 16, 12, 16)
  )

save_fig(fig5, "Figure_5_FINAL_h3_schematic.png", width = 9.5, height = 6.2)

# -------------------------
# DONE
# -------------------------

cat("\nDone. Final figures saved here:\n")
cat(OUT_DIR, "\n\n")
print(list.files(OUT_DIR, pattern = "\\.png$", full.names = TRUE))