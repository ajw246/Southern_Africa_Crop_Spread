suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(terra)
})

# Set warning handling and validate required inputs

without_vectorised_note <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Vectorized input to `element_text", conditionMessage(w),
              fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}
FILE_SITES <- file.path("data", "PaperA_fig_sites.csv")
OUT_PNG <- file.path("figures", "figure_S1_climate_zone_time_slices.png")
OUT_PDF <- file.path("figures", "figure_S1_climate_zone_time_slices.pdf")
if (!file.exists(FILE_SITES)) {
  stop("Cannot find ", FILE_SITES, ". Run 02_figure_datasets.R first.",
       call. = FALSE)
}

# Set climate periods, zone labels and plot colours

PERIODS <- c("1901_1930", "1931_1960", "1961_1990", "1991_2020")
REF_PERIOD <- "1991_2020"
kg_path <- function(period) {
  file.path("data", "koppen_geiger", period, "koppen_geiger_0p00833333.tif")
}
for (p in PERIODS) {
  if (!file.exists(kg_path(p))) {
    stop("Cannot find the ", p, " climate map at ", kg_path(p), ".\n",
         "All four periods are needed for this figure.  They are not committed\n",
         "to this repository - they are 50 MB of third-party raster from\n",
         "Beck et al. (2023), Scientific Data 10:724,\n",
         "https://doi.org/10.1038/s41597-023-02549-6.",
         "\n\nRun this once to fetch them, then run this script again:\n",
         '  source("00_fetch_koppen_geiger.R")', call. = FALSE)
  }
}
KG_ZONE <- c(
  "1" = "Af",  "2" = "Am",  "3" = "Aw",  "4" = "BWh", "5" = "BWk", "6" = "BSh",
  "7" = "BSk", "8" = "Csa", "9" = "Csb", "10" = "Csc", "11" = "Cwa",
  "12" = "Cwb", "13" = "Cwc", "14" = "Cfa", "15" = "Cfb", "16" = "Cfc",
  "17" = "Dsa", "18" = "Dsb", "19" = "Dsc", "20" = "Dsd", "21" = "Dwa",
  "22" = "Dwb", "23" = "Dwc", "24" = "Dwd", "25" = "Dfa", "26" = "Dfb",
  "27" = "Dfc", "28" = "Dfd", "29" = "ET",  "30" = "EF")
major_class <- function(code) substr(code, 1, 1)
KG_COLS <- c(
  Aw  = rgb( 70, 170, 250, maxColorValue = 255),
  BWh = rgb(255,   0,   0, maxColorValue = 255),
  BSh = rgb(245, 165,   0, maxColorValue = 255),
  BSk = rgb(255, 220, 100, maxColorValue = 255),
  Cwa = rgb(150, 255, 150, maxColorValue = 255),
  Cwb = rgb(100, 200, 100, maxColorValue = 255),
  Cfa = rgb(200, 255,  80, maxColorValue = 255),
  Cfb = rgb(100, 255,  80, maxColorValue = 255))
PERIOD_LABELS <- c(
  "1901_1930" = "1901–1930",
  "1931_1960" = "1931–1960",
  "1961_1990" = "1961–1990",
  "1991_2020" = "1991–2020\n(main analysis)")

# Read site records and extract zones from each climate map

message("Reading the site table ...")
sites <- read_csv(FILE_SITES, show_col_types = FALSE) |>
  select(site_name, country, latitude, longitude, kg_zone) |>
  arrange(site_name)
stopifnot(!any(is.na(sites$latitude)), !any(is.na(sites$longitude)),
          !any(duplicated(sites$site_name)))
points <- vect(as.data.frame(sites[, c("longitude", "latitude")]),
               geom = c("longitude", "latitude"), crs = "EPSG:4326")
read_period <- function(period) {
  grid   <- rast(kg_path(period))
  values <- as.integer(terra::extract(grid, points, method = "simple")[[2]])
  tibble(site_name = sites$site_name,
         period    = period,
         kg_zone   = unname(KG_ZONE[as.character(values)]))
}
message("Looking the 35 sites up in each of the four climate maps ...")
zones_long <- bind_rows(lapply(PERIODS, function(p) {
  message("  ", p)
  read_period(p)
}))
stopifnot(!any(is.na(zones_long$kg_zone)))

# Compare period zones and summarize site stability

ref_here <- zones_long |>
  filter(period == REF_PERIOD) |>
  arrange(match(site_name, sites$site_name)) |>
  pull(kg_zone)
stopifnot(identical(ref_here, sites$kg_zone))
message("  the 1991-2020 column matches the zone recorded in ", FILE_SITES)
zones_wide <- zones_long |>
  pivot_wider(names_from = period, values_from = kg_zone)
site_stability <- zones_wide |>
  rowwise() |>
  mutate(
    zones_seen   = paste(sort(unique(c_across(all_of(PERIODS)))),
                         collapse = " | "),
    n_zones      = n_distinct(c_across(all_of(PERIODS))),
    major_change = n_distinct(major_class(c_across(all_of(PERIODS)))) > 1) |>
  ungroup() |>
  left_join(select(sites, site_name, country), by = "site_name")
n_sites   <- nrow(site_stability)
n_moving  <- sum(site_stability$n_zones > 1)
n_major   <- sum(site_stability$major_change)
message(sprintf("  %d of %d sites change zone; %d cross a major boundary",
                n_moving, n_sites, n_major))
stopifnot(n_major <= n_moving)
row_order <- site_stability |>
  mutate(rank = (n_zones > 1) + major_change) |>
  arrange(rank, n_zones, site_name)
site_levels <- row_order$site_name
moving_sites <- site_stability$site_name[site_stability$n_zones > 1]
major_sites  <- site_stability$site_name[site_stability$major_change]
y_face   <- ifelse(site_levels %in% major_sites,  "bold",  "plain")
y_colour <- ifelse(site_levels %in% moving_sites, "black", "grey55")

# Prepare the climate stability plot

plot_df <- zones_long |>
  mutate(site_name = factor(site_name, levels = site_levels),
         period    = factor(period,    levels = PERIODS))
message("Drawing ...")
subtitle_text <- sprintf(
  paste("Beck et al. (2023) 1-km maps. %d/%d sites change zone across the",
        "four historical slices; %d cross a major A/B/C boundary (bold)."),
  n_moving, n_sites, n_major)
figure_S1 <- without_vectorised_note(ggplot(plot_df, aes(x = period, y = site_name, fill = kg_zone)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = kg_zone), size = 2.5, colour = "grey10") +
  scale_fill_manual(values = KG_COLS, name = "K-G zone (Beck 2023)") +
  scale_x_discrete(labels = PERIOD_LABELS, position = "top") +
  labs(
    title    = paste("Köppen-Geiger time-slice sensitivity for the",
                     "35 dated sites"),
    subtitle = subtitle_text,
    x = NULL, y = NULL,
    caption  = paste("Reference = 1991–2020 (used in main analysis).",
                     "All four slices are 20th-century observational climate,",
                     "not Iron-Age.")) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y   = element_text(size = 8, face = y_face, colour = y_colour),
    axis.text.x   = element_text(size = 8.5, face = "bold"),
    panel.grid    = element_blank(),
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, colour = "grey25"),
    plot.caption  = element_text(size = 7, colour = "grey40"),
    legend.position = "right"))

# Save Figure S1 and report zone changes

message("Saving ...")
ggsave(OUT_PNG, figure_S1, width = 8.5, height = 8, units = "in",
       dpi = 300, bg = "white")
ggsave(OUT_PDF, figure_S1, width = 8.5, height = 8, units = "in",
       bg = "white", device = grDevices::pdf, encoding = "WinAnsi.enc")
message("")
message(sprintf("Figure S1: %d sites x %d climate maps.", n_sites, length(PERIODS)))
message(sprintf("  %d sites keep the same zone in all four; %d move; %d of those",
                n_sites - n_moving, n_moving, n_major))
message("  change the first letter of the code.")
message("")
message("Sites that do NOT keep the same zone across the four maps:")
print(as.data.frame(
  site_stability |>
    filter(n_zones > 1) |>
    arrange(desc(major_change), site_name) |>
    select(site_name, country, zones_seen, major_change)),
  row.names = FALSE)
message("")
message("Stability grouped by the zone the 1991-2020 map gives:")
print(as.data.frame(
  site_stability |>
    group_by(reference_zone = .data[[REF_PERIOD]]) |>
    summarise(sites = n(), of_which_move = sum(n_zones > 1), .groups = "drop") |>
    arrange(desc(sites))),
  row.names = FALSE)
stopifnot(
  n_sites  == 35,
  n_moving == 13,
  n_major  == 5,
  all(unique(zones_long$kg_zone) %in% names(KG_COLS)))
message("")
message("Saved ", OUT_PNG, " and ", OUT_PDF, ".")
