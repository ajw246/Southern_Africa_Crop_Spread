# Load packages and validate site and climate grid inputs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(sf)
  library(terra)
  library(tidyterra)
  library(rnaturalearth)
  library(ggspatial)
  library(ggrepel)
})
for (pkg in c("ragg", "rnaturalearthdata", "marmap")) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop("The '", pkg, "' package is needed. Install it with:\n",
         '  install.packages("', pkg, '")', call. = FALSE)
  }
}

# Set message handling and map parameters

without_scale_note <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("Scale on map varies", conditionMessage(m), fixed = TRUE)) {
      invokeRestart("muffleMessage")
    }
  })
}
FILE_SITES <- file.path("data", "PaperA_fig_sites.csv")
KG_TIF <- file.path("data", "koppen_geiger", "1991_2020",
                    "koppen_geiger_0p00833333.tif")
CACHE <- file.path("data", "map_data_cache")
dir.create(CACHE, showWarnings = FALSE)
OUT_PNG <- file.path("figures", "figure_S3_site_map_climate_zones.png")
OUT_PDF <- file.path("figures", "figure_S3_site_map_climate_zones.pdf")
if (!file.exists(FILE_SITES)) {
  stop("Cannot find ", FILE_SITES, ". Run 02_figure_datasets.R first.",
       call. = FALSE)
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
quietly <- function(expr) suppressWarnings(suppressMessages(expr))
LON_MIN <-  11
LON_MAX <-  41
LAT_MIN <- -35.5
LAT_MAX <-  -5.8
window_ext <- ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX)
window_box <- st_bbox(c(xmin = LON_MIN, xmax = LON_MAX,
                        ymin = LAT_MIN, ymax = LAT_MAX),
                      crs = st_crs(4326))
ZONE_AREA_MIN_PCT <- 0.10
KG_LEGEND_ORDER <- c("Aw", "BWh", "BWk", "BSh", "BSk",
                     "Csa", "Csb", "Cwa", "Cwb", "Cfa", "Cfb")
beck_legend <- tibble::tribble(
  ~code, ~zone,  ~description,                                ~r,  ~g,  ~b,
   1,  "Af",  "Tropical, rainforest",                        0,   0, 255,
   2,  "Am",  "Tropical, monsoon",                           0, 120, 255,
   3,  "Aw",  "Tropical, savannah",                         70, 170, 250,
   4,  "BWh", "Arid, desert, hot",                         255,   0,   0,
   5,  "BWk", "Arid, desert, cold",                        255, 150, 150,
   6,  "BSh", "Arid, steppe, hot",                         245, 165,   0,
   7,  "BSk", "Arid, steppe, cold",                        255, 220, 100,
   8,  "Csa", "Temperate, dry summer, hot summer",         255, 255,   0,
   9,  "Csb", "Temperate, dry summer, warm summer",        200, 200,   0,
  10,  "Csc", "Temperate, dry summer, cold summer",        150, 150,   0,
  11,  "Cwa", "Temperate, dry winter, hot summer",         150, 255, 150,
  12,  "Cwb", "Temperate, dry winter, warm summer",        100, 200, 100,
  13,  "Cwc", "Temperate, dry winter, cold summer",         50, 150,  50,
  14,  "Cfa", "Temperate, no dry season, hot summer",      200, 255,  80,
  15,  "Cfb", "Temperate, no dry season, warm summer",     100, 255,  80,
  16,  "Cfc", "Temperate, no dry season, cold summer",      50, 200,   0,
  17,  "Dsa", "Cold, dry summer, hot summer",              255,   0, 255,
  18,  "Dsb", "Cold, dry summer, warm summer",             200,   0, 200,
  19,  "Dsc", "Cold, dry summer, cold summer",             150,  50, 150,
  20,  "Dsd", "Cold, dry summer, very cold winter",        150, 100, 150,
  21,  "Dwa", "Cold, dry winter, hot summer",              170, 175, 255,
  22,  "Dwb", "Cold, dry winter, warm summer",              90, 120, 220,
  23,  "Dwc", "Cold, dry winter, cold summer",              75,  80, 180,
  24,  "Dwd", "Cold, dry winter, very cold winter",         50,   0, 135,
  25,  "Dfa", "Cold, no dry season, hot summer",             0, 255, 255,
  26,  "Dfb", "Cold, no dry season, warm summer",           55, 200, 255,
  27,  "Dfc", "Cold, no dry season, cold summer",            0, 125, 125,
  28,  "Dfd", "Cold, no dry season, very cold winter",       0,  70,  95,
  29,  "ET",  "Polar, tundra",                             178, 178, 178,
  30,  "EF",  "Polar, frost",                              102, 102, 102
) |>
  mutate(hex = rgb(r, g, b, maxColorValue = 255))
beck_cols  <- setNames(beck_legend$hex, beck_legend$zone)
code_to_zone <- setNames(beck_legend$zone, beck_legend$code)
zone_label <- setNames(paste0(beck_legend$zone, ": ", beck_legend$description),
                       beck_legend$zone)
LAKE_FILL    <- "#1C86D6"
LAKE_OUTLINE <- "#0D4C8B"
RIVER_COLOUR <- "#3A6BA5"

# Read sites and verify their recorded climate zones

message("Reading the site table ...")
sites <- read_csv(FILE_SITES, show_col_types = FALSE)
stopifnot(!any(is.na(sites$latitude)), !any(is.na(sites$longitude)),
          !any(duplicated(sites$site_name)))
sites <- sites |>
  mutate(kg_zone = factor(kg_zone, levels = beck_legend$zone))
sites_sf <- st_as_sf(sites, coords = c("longitude", "latitude"),
                     crs = 4326, remove = FALSE)
message(sprintf("  %d sites in %d countries",
                nrow(sites), n_distinct(sites$country)))
message("Reading the climate grid ...")
kg_grid <- rast(KG_TIF) |> crop(window_ext)
coltab(kg_grid) <- NULL
site_codes <- terra::extract(kg_grid, as.matrix(st_coordinates(sites_sf)))[, 1]
site_zone_from_grid <- unname(code_to_zone[as.character(site_codes)])
n_matching <- sum(site_zone_from_grid == as.character(sites$kg_zone))
message(sprintf("  %d of %d sites sit on the zone the site table records",
                n_matching, nrow(sites)))
stopifnot(n_matching == nrow(sites))
zone_area <- freq(kg_grid) |>
  as_tibble() |>
  mutate(pct = 100 * count / sum(count))
qualifying_codes <- zone_area$value[zone_area$pct >= ZONE_AREA_MIN_PCT]
site_codes_in_key <- beck_legend$code[beck_legend$zone %in% site_zone_from_grid]
key_zones <- unname(code_to_zone[as.character(
  sort(unique(c(qualifying_codes, site_codes_in_key))))])
key_zones <- key_zones[!is.na(key_zones)]
present_codes <- sort(unique(values(kg_grid, mat = FALSE)))
present_codes <- present_codes[!is.na(present_codes) &
                                 present_codes %in% beck_legend$code]
levels(kg_grid) <- as.data.frame(
  beck_legend |> filter(code %in% present_codes) |> select(value = code, zone))
names(kg_grid) <- "zone"
missing_from_order <- setdiff(key_zones, KG_LEGEND_ORDER)
if (length(missing_from_order) > 0) {
  stop("These zones qualify for the key but are not in KG_LEGEND_ORDER: ",
       paste(missing_from_order, collapse = ", "),
       ".\nAdd them there, in the position they should read.", call. = FALSE)
}
key_zones <- KG_LEGEND_ORDER[KG_LEGEND_ORDER %in% key_zones]
message(sprintf("  all land coloured; %d zones in the key: %s",
                length(key_zones), paste(key_zones, collapse = ", ")))

# Prepare climate zone, land and ocean map layers

OCEAN_CACHE <- file.path(CACHE, "ocean_bathy_rgb.tif")
if (file.exists(OCEAN_CACHE)) {
  ocean_rgb <- rast(OCEAN_CACHE)
} else {
  message("Fetching sea-floor depths from NOAA (once) ...")
  old_wd <- setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)
  bathy <- try(marmap::getNOAA.bathy(lon1 = LON_MIN, lon2 = LON_MAX,
                                     lat1 = LAT_MIN, lat2 = LAT_MAX,
                                     resolution = 4,
                                     keep = FALSE), silent = TRUE)
  setwd(old_wd)
  if (inherits(bathy, "try-error")) {
    stop("Could not download the sea-floor depths from NOAA.\n",
         "Check that this computer is online and that ",
         "https://www.ncei.noaa.gov is reachable, then run the script again.\n",
         "The error was:\n", attr(bathy, "condition")$message, call. = FALSE)
  }
  relief <- crop(rast(marmap::as.raster(bathy)), window_ext)
  depth  <- ifel(relief < 0, relief, NA)
  deepest <- as.numeric(global(depth, "min", na.rm = TRUE))
  scaled  <- (depth - deepest) / (0 - deepest)
  ramp <- colorRampPalette(c("#16466F", "#21618E", "#3585B5",
                             "#5DA2CE", "#8FC0E2", "#BFDDF0"))(101)
  step <- app(scaled, function(v) ifelse(is.na(v), NA_integer_,
                                         as.integer(round(v * 100) + 1)))
  ramp_rgb <- col2rgb(ramp)
  base_r <- classify(step, cbind(1:101, ramp_rgb[1, ]))
  base_g <- classify(step, cbind(1:101, ramp_rgb[2, ]))
  base_b <- classify(step, cbind(1:101, ramp_rgb[3, ]))
  slope  <- terrain(depth, "slope",  unit = "radians")
  aspect <- terrain(depth, "aspect", unit = "radians")
  light  <- shade(slope, aspect, angle = 45, direction = 315)
  factor <- app(light, function(v) ifelse(is.na(v), 1, 0.90 + 0.20 * v))
  shade_channel <- function(channel) {
    app(c(channel, factor), function(v) {
      x <- v[1] * v[2]
      if (is.na(x)) NA else min(max(x, 0), 255)
    })
  }
  ocean_rgb <- c(shade_channel(base_r),
                 shade_channel(base_g),
                 shade_channel(base_b))
  names(ocean_rgb) <- c("r", "g", "b")
  ocean_rgb <- mask(ocean_rgb, depth)
  writeRaster(ocean_rgb, OCEAN_CACHE, overwrite = TRUE)
}

# Build the basemap and place country labels

message("Assembling the basemap ...")
lakes_cache <- file.path(CACHE, "lakes_region.gpkg")
if (file.exists(lakes_cache)) {
  lakes <- st_read(lakes_cache, quiet = TRUE)
} else {
  message("  downloading Natural Earth lakes (once) ...")
  lakes <- quietly(
    ne_download(scale = 10, type = "lakes", category = "physical",
                returnclass = "sf") |>
      st_make_valid() |>
      st_crop(window_box))
  st_write(lakes, lakes_cache, quiet = TRUE, delete_dsn = TRUE)
}
rivers_cache <- file.path(CACHE, "rivers_region.gpkg")
if (file.exists(rivers_cache)) {
  rivers <- st_read(rivers_cache, quiet = TRUE)
} else {
  message("  downloading Natural Earth rivers (once) ...")
  rivers <- quietly(
    ne_download(scale = 10, type = "rivers_lake_centerlines",
                category = "physical", returnclass = "sf") |>
      st_make_valid() |>
      st_crop(window_box))
  st_write(rivers, rivers_cache, quiet = TRUE, delete_dsn = TRUE)
}
if ("scalerank" %in% names(rivers)) {
  rivers <- rivers |>
    filter(suppressWarnings(as.numeric(scalerank)) <= 8)
}
countries_world <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_make_valid()
countries <- quietly(st_crop(countries_world, window_box))
NORTHERN_STUDY_COUNTRIES <- c("Angola", "Zambia", "Malawi", "Mozambique")
NORTHERN_NEIGHBOURS      <- c("Democratic Republic of the Congo",
                              "United Republic of Tanzania")
pieces <- quietly(st_cast(st_union(
  countries_world |> filter(admin %in% NORTHERN_STUDY_COUNTRIES)), "POLYGON"))
piece_area <- vapply(seq_along(pieces), function(i) {
  b <- st_bbox(pieces[i]); as.numeric((b[3] - b[1]) * (b[4] - b[2]))
}, numeric(1))
mainland <- pieces[which.max(piece_area)]
neighbours <- quietly(
  st_union(countries_world |> filter(admin %in% NORTHERN_NEIGHBOURS)))
study_area_north <- quietly(
  st_boundary(mainland) |>
    st_intersection(st_buffer(neighbours, 0.2)) |>
    st_collection_extract("LINESTRING") |>
    st_union() |>
    st_cast("MULTILINESTRING") |>
    st_line_merge())
country_labels <- tibble::tribble(
  ~label,          ~lon,   ~lat,
  "ANGOLA",        16.5,  -11.0,
  "ZAMBIA",        25.3,  -14.3,
  "MALAWI",        34.8,  -13.5,
  "NAMIBIA",       17.8,  -25.6,
  "BOTSWANA",      23.4,  -23.3,
  "ZIMBABWE",      27.1,  -18.3,
  "MOZAMBIQUE",    35.3,  -19.4,
  "SOUTH AFRICA",  22.6,  -30.4
)
avoid_anchors <- country_labels |>
  mutate(half = 0.32 + 0.085 * nchar(label)) |>
  rowwise() |>
  reframe(longitude = lon + c(-half, 0, half),
          latitude  = lat,
          site_name = "")
repel_data <- bind_rows(
  sites |> select(longitude, latitude, site_name),
  avoid_anchors)

# Draw and save Figure S3

message("Drawing ...")
zone_layer <- suppressMessages(
  geom_spatraster(data = kg_grid, aes(fill = zone), alpha = 0.5, maxcell = 3e6))
figure_S3 <- ggplot() +
  geom_spatraster_rgb(data = ocean_rgb, maxcell = 6e6) +
  geom_sf(data = countries, fill = "white", colour = NA) +
  zone_layer +
  geom_sf(data = lakes, fill = LAKE_FILL, colour = LAKE_OUTLINE,
          linewidth = 0.35) +
  geom_sf(data = rivers, colour = RIVER_COLOUR, linewidth = 0.35,
          alpha = 0.85) +
  geom_sf(data = countries, fill = NA, colour = "white", linewidth = 0.4) +
  geom_sf(data = study_area_north, colour = "white",
          linewidth = 1.7, linetype = "dashed", lineend = "round") +
  geom_text(data = country_labels, aes(x = lon, y = lat, label = label),
            colour = "grey40", size = 3.6, fontface = "bold", alpha = 0.55) +
  geom_sf(data = sites_sf, aes(fill = kg_zone),
          shape = 21, colour = "black", size = 2.8, stroke = 0.7) +
  geom_text_repel(
    data = repel_data, aes(x = longitude, y = latitude, label = site_name),
    size = 4.0, colour = "grey10", fontface = "bold",
    bg.color = "white", bg.r = 0.16,
    min.segment.length = 0, segment.size = 0.4, segment.colour = "grey30",
    box.padding = 0.3, point.padding = 0.5,
    force = 2.0, force_pull = 1.0, max.overlaps = Inf,
    max.iter = 200000, max.time = Inf, seed = 1
  ) +
  scale_fill_manual(
    values      = beck_cols,
    breaks      = key_zones,
    labels      = zone_label[key_zones],
    drop        = TRUE,
    na.translate = FALSE,
    name = "Köppen-Geiger climate zone\n(Beck et al. 2023, 1991-2020)",
    guide = guide_legend(order = 1, ncol = 3, byrow = TRUE,
                         title.position = "top",
                         override.aes = list(shape = 22, size = 5, alpha = 1,
                                             colour = "grey40"))) +
  coord_sf(xlim = c(LON_MIN, LON_MAX), ylim = c(LAT_MIN, LAT_MAX),
           expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.22, text_cex = 0.7,
                   line_width = 0.5, bar_cols = c("grey30", "white")) +
  annotation_north_arrow(location = "br", which_north = "true",
                         height = unit(0.95, "cm"), width = unit(0.95, "cm"),
                         pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         style = north_arrow_fancy_orienteering()) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "#16466F", colour = NA),
    panel.grid       = element_line(colour = "grey92", linewidth = 0.15),
    panel.border     = element_rect(colour = "grey35", fill = NA,
                                    linewidth = 0.6),
    legend.position  = "bottom",
    legend.key       = element_rect(fill = "white", colour = NA),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 10.5),
    axis.text        = element_text(size = 10, colour = "grey40"))

# Write the map files and summarize site zones

message("Saving ...")
without_scale_note(
  ggsave(OUT_PNG, figure_S3, width = 11, height = 12.8, units = "in",
         dpi = 300, bg = "white", device = ragg::agg_png))
message("")
message(sprintf("Figure S3: %d sites, %d zones in the key.",
                nrow(sites), length(key_zones)))
message("")
message("Sites per climate zone:")
print(as.data.frame(
  sites |>
    count(zone = kg_zone, name = "sites", .drop = TRUE) |>
    filter(sites > 0) |>
    mutate(description = zone_label[as.character(zone)]) |>
    select(zone, sites, description)),
  row.names = FALSE)
message("")
message("Share of the land in the window covered by each zone in the key:")
print(as.data.frame(
  zone_area |>
    transmute(zone = unname(code_to_zone[as.character(value)]),
              land_pct = round(pct, 2)) |>
    filter(zone %in% key_zones) |>
    arrange(match(zone, KG_LEGEND_ORDER))),
  row.names = FALSE)
stopifnot(
  nrow(sites)   == 35,
  n_matching    == 35,
  identical(key_zones, c("Aw", "BWh", "BWk", "BSh", "BSk",
                         "Csa", "Csb", "Cwa", "Cwb", "Cfa", "Cfb")))
message("")
message("Saved ", OUT_PNG, " and ", OUT_PDF, ".")
