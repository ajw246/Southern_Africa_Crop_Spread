# Load packages and validate data and climate grid

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(terra)
  library(tidyterra)
  library(ggnewscale)
  library(rnaturalearth)
})
for (pkg in c("ragg", "rnaturalearthdata")) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop("The '", pkg, "' package is needed. Install it with:\n",
         '  install.packages("', pkg, '")', call. = FALSE)
  }
}
FILE_RECORDS <- file.path("data", "PaperA_fig_plant_records.csv")
FILE_BINS    <- file.path("data", "PaperA_fig_sample_bins.csv")
FILE_SITES   <- file.path("data", "PaperA_fig_sites.csv")
KG_TIF <- file.path("data", "koppen_geiger", "1991_2020",
                    "koppen_geiger_0p00833333.tif")
OUT_PNG <- file.path("figures", "figure_5_crop_maps_by_period.png")
OUT_PDF <- file.path("figures", "figure_5_crop_maps_by_period.pdf")
for (f in c(FILE_RECORDS, FILE_BINS, FILE_SITES)) {
  if (!file.exists(f)) {
    stop("Cannot find ", f, ". Run 02_figure_datasets.R first.", call. = FALSE)
  }
}
if (!file.exists(KG_TIF)) {
  stop("Cannot find the climate grid at ", KG_TIF, ".\n",
       "It is not committed to this repository - it is 50 MB of third-party\n",
       "raster from Beck et al. (2023), Scientific Data 10:724,\n",
       "https://doi.org/10.1038/s41597-023-02549-6.",
       "\n\nRun this once to fetch it, then run this script again:\n",
       '  source("00_fetch_koppen_geiger.R")', call. = FALSE)
}
suppressMessages(sf_use_s2(FALSE))

# Set crop selection, map bounds and climate legend

quietly <- function(expr) suppressWarnings(suppressMessages(expr))
BASE_FAMILY <- "Helvetica"
CROP_SPECIES <- c("Sorghum bicolor",
                  "Cenchrus americanus",
                  "Eleusine coracana",
                  "Vigna subterranea",
                  "Vigna unguiculata")
BIN_ORDER <- c("1 - 500 CE", "501 - 1000 CE", "1001 - 1500 CE")
LON_MIN <-  11
LON_MAX <-  41
LAT_MIN <- -35.5
LAT_MAX <-  -6
LAT_MAX_GRID <- -5.8
ALBERS <- paste("+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=-20 +lon_0=25",
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
beck_legend <- tibble::tribble(
  ~code, ~zone,  ~description,                              ~r,  ~g,  ~b,
   1,    "Af",  "Tropical, rainforest",                      0,   0, 255,
   2,    "Am",  "Tropical, monsoon",                         0, 120, 255,
   3,    "Aw",  "Tropical, savannah",                       70, 170, 250,
   4,    "BWh", "Arid, desert, hot",                       255,   0,   0,
   5,    "BWk", "Arid, desert, cold",                      255, 150, 150,
   6,    "BSh", "Arid, steppe, hot",                       245, 165,   0,
   7,    "BSk", "Arid, steppe, cold",                      255, 220, 100,
   8,    "Csa", "Temperate, dry summer, hot summer",       255, 255,   0,
   9,    "Csb", "Temperate, dry summer, warm summer",      200, 200,   0,
  11,    "Cwa", "Temperate, dry winter, hot summer",       150, 255, 150,
  12,    "Cwb", "Temperate, dry winter, warm summer",      100, 200, 100,
  14,    "Cfa", "Temperate, no dry season, hot summer",    200, 255,  80,
  15,    "Cfb", "Temperate, no dry season, warm summer",   100, 255,  80,
  29,    "ET",  "Polar, tundra",                           178, 178, 178,
  30,    "EF",  "Polar, frost",                            102, 102, 102
) |>
  mutate(hex = rgb(r, g, b, maxColorValue = 255))
beck_cols <- setNames(beck_legend$hex, beck_legend$zone)
KG_LABEL_WRAP <- 24
wrap_label <- function(x, width = KG_LABEL_WRAP) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}
zone_label <- setNames(
  wrap_label(paste0(beck_legend$zone, ": ", beck_legend$description)),
  beck_legend$zone)
KG_LEGEND_DROP <- c("ET", "EF", "Af", "Am")
KG_LEGEND_ORDER <- c("Aw", "BWh", "BWk", "BSh", "BSk",
                     "Csa", "Csb", "Cwa", "Cwb", "Cfa", "Cfb")
SUIT_LEVELS <- c("Suitable", "Not suitable")
SUIT_COLS   <- c("Suitable" = "#4A8B8B", "Not suitable" = "#FFFFFF8C")

# Read figure tables and prepare crop occurrences

message("Reading the three tables ...")
plant_records <- read_csv(FILE_RECORDS, show_col_types = FALSE)
sample_bins   <- read_csv(FILE_BINS,    show_col_types = FALSE)
sites         <- read_csv(FILE_SITES,   show_col_types = FALSE)
crop_records <- plant_records |>
  filter(taxon_merged_figs45 %in% CROP_SPECIES,
         !is_post_1500_taxon,
         !is_genus_exclude)
occurrences <- crop_records |>
  select(species = taxon_merged_figs45, common_name, site_name, sample_number) |>
  inner_join(select(sample_bins, site_name, sample_number, bin_label, bin_index),
             by = c("site_name", "sample_number"),
             relationship = "many-to-many") |>
  filter(bin_label %in% BIN_ORDER) |>
  inner_join(select(sites, site_name, latitude, longitude, kg_zone),
             by = "site_name") |>
  distinct(species, common_name, site_name, bin_label, bin_index,
           latitude, longitude, kg_zone)
stopifnot(!any(is.na(occurrences$latitude)),
          !any(is.na(occurrences$kg_zone)))
ecocrop_codes <- plant_records |>
  filter(taxon_merged_figs45 %in% CROP_SPECIES, !is.na(ecocrop_kg_codes)) |>
  distinct(species = taxon_merged_figs45, codes = ecocrop_kg_codes)
stopifnot(!any(duplicated(ecocrop_codes$species)),
          setequal(ecocrop_codes$species, CROP_SPECIES))
split_codes <- function(codes) {
  z <- trimws(unlist(strsplit(codes, "[;,]")))
  unique(z[nzchar(z)])
}
suit_zones <- setNames(lapply(ecocrop_codes$codes, split_codes),
                       ecocrop_codes$species)
occurrences <- occurrences |>
  mutate(suitability = ifelse(
           mapply(function(sp, z) z %in% suit_zones[[sp]], species, kg_zone),
           "Suitable", "Not suitable"))
crop_labels <- occurrences |>
  distinct(species, common_name)
crop_labels <- crop_labels[match(CROP_SPECIES, crop_labels$species), ]
occurrences <- occurrences |>
  mutate(common_name = factor(common_name, levels = crop_labels$common_name),
         bin_label   = factor(bin_label,   levels = BIN_ORDER),
         suitability = factor(suitability, levels = SUIT_LEVELS))
occurrences <- occurrences |>
  arrange(species, site_name, bin_index)
message(sprintf("  %d dots: %d sites, %d crops, %d windows",
                nrow(occurrences), n_distinct(occurrences$site_name),
                n_distinct(occurrences$common_name),
                n_distinct(occurrences$bin_label)))
message("  site climate zones present: ",
        paste(sort(unique(occurrences$kg_zone)), collapse = ", "))

# Project crop sites and prepare the basemap

projected <- occurrences |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(ALBERS) |>
  st_coordinates()
occurrences$x_albers <- projected[, "X"]
occurrences$y_albers <- projected[, "Y"]
message("Building the basemap ...")
n <- 160
ring <- rbind(
  cbind(seq(LON_MIN, LON_MAX, length.out = n), LAT_MAX),
  cbind(LON_MAX, seq(LAT_MAX, LAT_MIN, length.out = n)),
  cbind(seq(LON_MAX, LON_MIN, length.out = n), LAT_MIN),
  cbind(LON_MIN, seq(LAT_MIN, LAT_MAX, length.out = n)))
ring <- rbind(ring, ring[1, , drop = FALSE])
window_albers <- st_sfc(st_polygon(list(unname(ring))), crs = 4326) |>
  st_transform(ALBERS)
page <- st_bbox(window_albers)
X_RANGE <- c(page[["xmin"]], page[["xmax"]])
Y_RANGE <- c(page[["ymin"]], page[["ymax"]])
countries <- quietly(
  ne_countries(scale = "medium", returnclass = "sf") |>
    st_make_valid() |>
    st_crop(st_bbox(c(xmin = LON_MIN, xmax = LON_MAX,
                      ymin = LAT_MIN, ymax = LAT_MAX), crs = st_crs(4326))) |>
    st_transform(ALBERS))
kg_grid <- rast(KG_TIF) |>
  crop(ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX_GRID))
coltab(kg_grid) <- NULL
kg_albers <- project(kg_grid, ALBERS, method = "near")
codes_present <- sort(unique(values(kg_albers, mat = FALSE)))
codes_present <- codes_present[!is.na(codes_present) &
                                 codes_present %in% beck_legend$code]
levels(kg_albers) <- as.data.frame(
  beck_legend |> filter(code %in% codes_present) |> select(value = code, zone))
names(kg_albers) <- "zone"
kg_visible <- terra::mask(kg_albers,
                          terra::vect(st_sf(geometry = window_albers)))
visible_codes <- sort(unique(values(kg_visible, mat = FALSE)))
visible_codes <- visible_codes[!is.na(visible_codes) &
                                 visible_codes %in% beck_legend$code]
zones_present <- beck_legend |>
  filter(code %in% visible_codes, !zone %in% KG_LEGEND_DROP) |>
  pull(zone)
zones_present <- KG_LEGEND_ORDER[KG_LEGEND_ORDER %in% zones_present]
message("  climate zones in the window: ", paste(zones_present, collapse = ", "))

# Draw the crop maps for each time period

message("Drawing ...")
land_layer <- suppressMessages(
  geom_spatraster(data = kg_albers, aes(fill = zone),
                  alpha = 0.42, maxcell = 5e5))
figure_5 <- ggplot() +
  land_layer +
  geom_sf(data = countries, fill = NA, colour = "grey55", linewidth = 0.25) +
  scale_fill_manual(values = beck_cols, breaks = zones_present,
                    labels = zone_label[zones_present], drop = FALSE,
                    na.translate = FALSE, name = "Köppen–Geiger zones:",
                    guide = guide_legend(order = 1, ncol = 3, byrow = TRUE,
                                         title.position = "left",
                                         title.vjust = 0.5,
                                         theme = theme(
                                           legend.text = element_text(
                                             size = 7.5, lineheight = 0.95),
                                           legend.key.spacing.x = unit(8, "pt"),
                                           legend.key.spacing.y = unit(2, "pt")),
                                         override.aes = list(alpha = 1,
                                                             colour = "grey40"))) +
  new_scale_fill() +
  geom_point(data = occurrences,
             aes(x = x_albers, y = y_albers, fill = suitability),
             shape = 21, colour = "black", size = 2.4, stroke = 0.6,
             alpha = 0.75) +
  scale_fill_manual(values = SUIT_COLS, breaks = SUIT_LEVELS, drop = FALSE,
                    na.translate = FALSE,
                    name = "Crop suitability for K-G zone (ECOCROP):",
                    guide = guide_legend(order = 2, nrow = 1,
                                         override.aes = list(alpha = 1,
                                                             colour = "black"))) +
  coord_sf(crs = ALBERS, xlim = X_RANGE, ylim = Y_RANGE, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  facet_grid(common_name ~ bin_label, switch = "y", drop = FALSE) +
  theme_minimal(base_size = 9, base_family = BASE_FAMILY) +
  theme(
    panel.background  = element_rect(fill = "#EAF2F8", colour = NA),
    panel.border      = element_rect(colour = "grey55", fill = NA,
                                     linewidth = 0.4),
    panel.grid        = element_blank(),
    panel.spacing     = unit(2, "pt"),
    axis.text         = element_blank(),
    axis.ticks        = element_blank(),
    strip.text.x      = element_text(size = 8, face = "bold"),
    strip.text.y.left = element_text(size = 8, face = "bold", angle = 0),
    strip.placement   = "outside",
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.location   = "plot",
    legend.title      = element_text(size = 8.5, face = "bold"),
    legend.text       = element_text(size = 8),
    plot.title        = element_text(size = 11, face = "bold"),
    plot.subtitle     = element_text(size = 8.5),
    plot.margin       = margin(6, 8, 6, 6))

# Save Figure 5 as PNG and PDF

message("Saving ...")
ggsave(OUT_PNG, figure_5, device = ragg::agg_png,
       width = 9.5, height = 9.7, units = "in", dpi = 300, bg = "white")
