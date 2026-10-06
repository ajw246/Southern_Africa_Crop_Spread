# Load packages and validate figure inputs

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
for (pkg in c("ragg", "rnaturalearthdata", "geodata", "scales")) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop("The '", pkg, "' package is needed. Install it with:\n",
         '  install.packages("', pkg, '")', call. = FALSE)
  }
}

# Set output paths and map parameters

FILE_RECORDS <- file.path("data", "PaperA_fig_plant_records.csv")
FILE_BINS    <- file.path("data", "PaperA_fig_sample_bins.csv")
FILE_SITES   <- file.path("data", "PaperA_fig_sites.csv")
CACHE      <- file.path("data", "map_data_cache")
ELEV_CACHE <- file.path(CACHE, "elevation_region.tif")
OUT_PNG <- file.path("figures", "figure_S2_crop_maps_by_elevation.png")
OUT_PDF <- file.path("figures", "figure_S2_crop_maps_by_elevation.pdf")
for (f in c(FILE_RECORDS, FILE_BINS, FILE_SITES)) {
  if (!file.exists(f)) {
    stop("Cannot find ", f, ". Run 02_figure_datasets.R first.", call. = FALSE)
  }
}
suppressMessages(sf_use_s2(FALSE))
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
LAT_MAX_GRID <- -4
ALBERS <- paste("+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=-20 +lon_0=25",
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
ELEV_LOW  <- "#EFE9E1"
ELEV_HIGH <- "#3B3024"
SUIT_LEVELS <- c("Suitable", "Not suitable")
SUIT_COLS   <- c("Suitable" = "#4A8B8B", "Not suitable" = "#FFFFFF8C")

# Read crop records and calculate site elevations

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
  inner_join(select(sites, site_name, latitude, longitude),
             by = "site_name") |>
  distinct(species, common_name, site_name, bin_label, bin_index,
           latitude, longitude)
stopifnot(!any(is.na(occurrences$latitude)),
          !any(is.na(occurrences$longitude)))
message(sprintf("  %d dots: %d sites, %d crops, %d windows",
                nrow(occurrences), n_distinct(occurrences$site_name),
                n_distinct(occurrences$species),
                n_distinct(occurrences$bin_label)))
if (!file.exists(ELEV_CACHE)) {
  dir.create(CACHE, showWarnings = FALSE)
  message("Downloading the WorldClim elevation grid (once, about 5 MB) ...")
  world <- quietly(geodata::elevation_global(res = 5, path = CACHE))
  writeRaster(crop(world, ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX_GRID)),
              ELEV_CACHE, overwrite = TRUE)
  message("  cached as ", ELEV_CACHE)
}
elev_lonlat <- rast(ELEV_CACHE)
names(elev_lonlat) <- "elevation"
site_xy <- as.data.frame(occurrences[, c("longitude", "latitude")])
elev_m <- terra::extract(elev_lonlat, site_xy, method = "bilinear")[, 2]
missing <- which(is.na(elev_m))
if (length(missing) > 0) {
  elev_m[missing] <- terra::extract(elev_lonlat,
                                    site_xy[missing, , drop = FALSE],
                                    method = "simple")[, 2]
}
elev_m[is.na(elev_m)] <- 0
occurrences$elev_m <- pmax(elev_m, 0)
message(sprintf("  sites sit between %.0f and %.0f m above sea level",
                min(occurrences$elev_m), max(occurrences$elev_m)))

# Compare site elevations with crop growing ranges

alt_ranges <- plant_records |>
  filter(taxon_merged_figs45 %in% CROP_SPECIES,
         !is.na(ecocrop_alt_min_m), !is.na(ecocrop_alt_max_m)) |>
  distinct(species = taxon_merged_figs45,
           alt_min = ecocrop_alt_min_m,
           alt_max = ecocrop_alt_max_m)
stopifnot(!any(duplicated(alt_ranges$species)),
          setequal(alt_ranges$species, CROP_SPECIES),
          all(alt_ranges$alt_min < alt_ranges$alt_max))
occurrences <- occurrences |>
  left_join(alt_ranges, by = "species") |>
  mutate(suitability = ifelse(elev_m >= alt_min & elev_m <= alt_max,
                              "Suitable", "Not suitable"))
stopifnot(!any(is.na(occurrences$suitability)))
crop_labels <- occurrences |>
  distinct(species, common_name)
crop_labels <- crop_labels[match(CROP_SPECIES, crop_labels$species), ]
occurrences <- occurrences |>
  mutate(common_name = factor(common_name, levels = crop_labels$common_name),
         bin_label   = factor(bin_label,   levels = BIN_ORDER),
         suitability = factor(suitability, levels = SUIT_LEVELS))
occurrences <- occurrences |>
  arrange(species, site_name, bin_index)

# Prepare the projected elevation basemap

message("Building the basemap ...")
projected <- occurrences |>
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
  st_transform(ALBERS) |>
  st_coordinates()
occurrences$x_albers <- projected[, "X"]
occurrences$y_albers <- projected[, "Y"]
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
elev_albers <- project(elev_lonlat, ALBERS, method = "bilinear")
elev_albers <- clamp(elev_albers, lower = 0, values = TRUE)
elev_albers <- terra::mask(elev_albers, terra::vect(countries))
elev_range <- minmax(elev_albers)[, 1]
message(sprintf("  land in the window runs from %.0f to %.0f m",
                elev_range[1], elev_range[2]))
ELEV_LIMIT <- ceiling(elev_range[2] / 100) * 100

# Draw crop maps by period and elevation suitability

message("Drawing ...")
land_layer <- suppressMessages(
  geom_spatraster(data = elev_albers, aes(fill = elevation), maxcell = 5e5))
figure_S2 <- ggplot() +
  land_layer +
  geom_sf(data = countries, fill = NA, colour = "grey55", linewidth = 0.25) +
  scale_fill_gradient(
    low = ELEV_LOW, high = ELEV_HIGH, na.value = "transparent",
    name = "Elevation (m)",
    limits = c(0, ELEV_LIMIT),
    breaks = scales::breaks_width(750),
    guide = guide_colourbar(order = 1, barwidth = unit(200, "pt"),
                            barheight = unit(7, "pt"), title.vjust = 1)) +
  ggnewscale::new_scale_fill() +
  geom_point(data = occurrences,
             aes(x = x_albers, y = y_albers, fill = suitability),
             shape = 21, colour = "black", size = 2.4, stroke = 0.6,
             alpha = 0.75) +
  scale_fill_manual(values = SUIT_COLS, breaks = SUIT_LEVELS, drop = FALSE,
                    na.translate = FALSE,
                    name = "Site within crop's FAO ECOCROP altitude range:",
                    guide = guide_legend(order = 2, nrow = 1,
                                         override.aes = list(
                                           alpha = 1, colour = "black",
                                           fill = unname(SUIT_COLS[SUIT_LEVELS])))) +
  coord_sf(crs = ALBERS, xlim = X_RANGE, ylim = Y_RANGE, expand = FALSE) +
  labs(x = NULL, y = NULL,
       title = paste("Introduced crops: presence by site and 500-yr bin",
                     "(sample-level chronology)")) +
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
    legend.title      = element_text(size = 8.5, face = "bold"),
    legend.text       = element_text(size = 8),
    plot.title        = element_text(size = 11, face = "bold"),
    plot.subtitle     = element_text(size = 8.5),
    plot.margin       = margin(6, 8, 6, 6))

# Save Figure S2 as PNG and PDF

message("Saving ...")
ggsave(OUT_PNG, figure_S2, device = ragg::agg_png,
       width = 9.5, height = 9.2, units = "in", dpi = 300, bg = "white")
crop_sites <- occurrences |> distinct(common_name, site_name, suitability)
message("")
message(sprintf("Figure S2: %d dots - %d sites, %d crop-site pairs.",
                nrow(occurrences), n_distinct(occurrences$site_name),
                nrow(crop_sites)))
message("")
message("Dots per crop per window:")
print(table(occurrences$common_name, occurrences$bin_label))
message("")
message("Altitude suitability, counted by crop-site pair:")
print(table(crop_sites$common_name, crop_sites$suitability))
message("")
message("Site altitudes, by crop:")
print(as.data.frame(
  occurrences |>
    group_by(Crop = common_name) |>
    summarise(sites = n_distinct(site_name),
              lowest_m  = round(min(elev_m)),
              highest_m = round(max(elev_m)),
              fao_range = paste0(first(alt_min), "-", first(alt_max), " m"),
              .groups = "drop")),
  row.names = FALSE)
stopifnot(
  nrow(occurrences) == 57,
  n_distinct(occurrences$site_name) == 27,
  nrow(crop_sites) == 52,
  all(occurrences$suitability == "Suitable"))
message("")
message("Saved ", OUT_PNG, " and ", OUT_PDF, ".")
