# Load packages and validate inputs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(terra)
  library(tidyterra)
  library(rnaturalearth)
  library(ggspatial)
  library(ggrepel)
})
for (pkg in c("marmap", "ragg", "rnaturalearthdata", "geodata")) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop("The '", pkg, "' package is needed. Install it with:\n",
         '  install.packages("', pkg, '")', call. = FALSE)
  }
}

# Set output paths, map bounds and projection

FILE_SITES <- file.path("data", "PaperA_fig_sites.csv")
OUT_PNG <- file.path("figures", "figure_2_site_map.png")
OUT_PDF <- file.path("figures", "figure_2_site_map.pdf")
if (!file.exists(FILE_SITES)) {
  stop("Cannot find ", FILE_SITES, ". Run 02_figure_datasets.R first.",
       call. = FALSE)
}
CACHE <- file.path("data", "map_data_cache")
dir.create(CACHE, showWarnings = FALSE)
suppressMessages(sf_use_s2(FALSE))
terraOptions(progress = 0)
quietly <- function(expr) suppressWarnings(suppressMessages(expr))
LON_MIN <-   9
LON_MAX <-  44.1
LAT_MIN <- -35.5
LAT_MAX <-  -5.4
ALBERS <- paste("+proj=aea +lat_1=-18 +lat_2=-32 +lat_0=-20 +lon_0=25",
                "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
project_edge <- function(lon, lat) {
  pts <- sf::st_sfc(sf::st_multipoint(cbind(lon, lat)), crs = 4326)
  sf::st_coordinates(sf::st_transform(pts, ALBERS))
}
n_pts    <- 400
edge_top <- project_edge(seq(LON_MIN, LON_MAX, length.out = n_pts), LAT_MAX)
edge_bot <- project_edge(seq(LON_MIN, LON_MAX, length.out = n_pts), LAT_MIN)
edge_lft <- project_edge(LON_MIN, seq(LAT_MIN, LAT_MAX, length.out = n_pts))
edge_rgt <- project_edge(LON_MAX, seq(LAT_MIN, LAT_MAX, length.out = n_pts))
X_RANGE <- c(max(edge_lft[, 1]), min(edge_rgt[, 1]))
Y_RANGE <- c(max(edge_bot[, 2]), min(edge_top[, 2]))
message(sprintf("Map panel: %.0f x %.0f km (height / width = %.3f)",
                diff(X_RANGE) / 1000, diff(Y_RANGE) / 1000,
                diff(Y_RANGE) / diff(X_RANGE)))
RELIEF_VEIL <- 0.28
HS_STRENGTH <- 0.28
Z_EXAG <- 12
HS_CLIP <- c(0.02, 0.98)
RELIEF_BREAKS <- c(-6000,    -4000,    -2000,     -500,     -100,       -1,
                       0,      200,      500,      800,     1100,     1400,
                    1700,     2100,     2600,     3500)
RELIEF_COLS   <- c("#5B7F9B","#6B8DA7","#7C9DB4","#96B2C4","#B3C8D6","#CFDFE9",
                   "#C2CFB2","#CBD4B6","#D6D9BB","#DEDCC0","#DED4B4","#D6C6A6",
                   "#CBB699","#BDA68B","#AE957E","#A08872")

# Project archaeological sites and check map coverage

sites <- read_csv(FILE_SITES, show_col_types = FALSE) |>
  select(site_name, country, latitude, longitude)
sites_sf <- st_as_sf(sites, coords = c("longitude", "latitude"),
                     crs = 4326, remove = FALSE) |>
  st_transform(ALBERS)
site_xy <- st_coordinates(sites_sf)
sites_labels <- sites |>
  mutate(x = site_xy[, 1], y = site_xy[, 2])
outside <- site_xy[, 1] < X_RANGE[1] | site_xy[, 1] > X_RANGE[2] |
           site_xy[, 2] < Y_RANGE[1] | site_xy[, 2] > Y_RANGE[2]
if (any(outside)) {
  stop("These sites fall outside the map window: ",
       paste(sites$site_name[outside], collapse = ", "),
       ".\nWiden LON_MIN / LON_MAX / LAT_MIN / LAT_MAX in Section 2.",
       call. = FALSE)
}
message(sprintf("Sites: %d, all inside the map window.", nrow(sites)))

# Load or fetch terrain relief and prepare the shaded basemap

ETOPO_CACHE <- file.path(
  CACHE, sprintf("etopo2022_30s_%g_%g_%g_%g.tif",
                 LON_MIN, LON_MAX, LAT_MIN, LAT_MAX))
if (file.exists(ETOPO_CACHE)) {
  etopo <- rast(ETOPO_CACHE)
} else {
  message("Fetching the ETOPO 2022 relief grid from NOAA. ",
          "This happens once and takes about a minute ...")
  old_wd <- setwd(tempdir())
  on.exit(setwd(old_wd), add = TRUE)
  fetch_tile <- function(lon1, lon2, lat1, lat2) {
    rast(marmap::as.raster(
      marmap::getNOAA.bathy(lon1 = lon1, lon2 = lon2,
                            lat1 = lat1, lat2 = lat2,
                            resolution = 0.5,
                            keep = FALSE)))
  }
  etopo <- try(fetch_tile(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX), silent = TRUE)
  if (inherits(etopo, "try-error")) {
    message("  The single request was refused. Retrying as four quarters ...")
    lon_mid <- (LON_MIN + LON_MAX) / 2
    lat_mid <- (LAT_MIN + LAT_MAX) / 2
    etopo <- try(terra::merge(
      fetch_tile(LON_MIN, lon_mid, LAT_MIN, lat_mid),
      fetch_tile(lon_mid, LON_MAX, LAT_MIN, lat_mid),
      fetch_tile(LON_MIN, lon_mid, lat_mid, LAT_MAX),
      fetch_tile(lon_mid, LON_MAX, lat_mid, LAT_MAX)), silent = TRUE)
  }
  if (inherits(etopo, "try-error")) {
    stop("Could not download the ETOPO 2022 relief grid from NOAA.\n",
         "Check that this computer is online and that ",
         "https://www.ncei.noaa.gov is reachable, then run the script again.\n",
         "The error was:\n", attr(etopo, "condition")$message, call. = FALSE)
  }
  setwd(old_wd)
  crs(etopo)   <- "EPSG:4326"
  names(etopo) <- "elevation"
  writeRaster(round(etopo), ETOPO_CACHE, overwrite = TRUE, datatype = "INT2S",
              gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"))
  etopo <- rast(ETOPO_CACHE)
}
names(etopo) <- "elevation"
message(sprintf("Relief grid: %d x %d cells, %.0f m to %.0f m.",
                ncol(etopo), nrow(etopo),
                minmax(etopo)[1, 1], minmax(etopo)[2, 1]))
RELIEF_RGB_CACHE <- file.path(CACHE, "relief_rgb_albers.tif")
SETTINGS_FILE    <- paste0(RELIEF_RGB_CACHE, ".settings")
settings_now <- paste(Z_EXAG, HS_STRENGTH, RELIEF_VEIL,
                      paste(HS_CLIP,       collapse = ","),
                      paste(RELIEF_BREAKS, collapse = ","),
                      paste(RELIEF_COLS,   collapse = ","),
                      LON_MIN, LON_MAX, LAT_MIN, LAT_MAX, sep = "|")
settings_match <- file.exists(SETTINGS_FILE) &&
  identical(readLines(SETTINGS_FILE, warn = FALSE), settings_now)
if (file.exists(RELIEF_RGB_CACHE) && settings_match) {
  relief_rgb <- rast(RELIEF_RGB_CACHE)
} else {
  if (file.exists(RELIEF_RGB_CACHE)) {
    message("The terrain settings have changed - rebuilding the basemap image.")
  }
  message("Building the shaded relief basemap. ",
          "This happens once and takes about a minute ...")
  relief <- project(etopo, ALBERS, method = "bilinear")
  exaggerated <- relief * Z_EXAG
  slope       <- terrain(exaggerated, "slope",  unit = "radians")
  aspect      <- terrain(exaggerated, "aspect", unit = "radians")
  hillshade   <- Reduce(`+`, lapply(
    c(270, 315, 0, 45),
    function(direction) shade(slope, aspect, angle = 45, direction = direction)
  )) / 4
  sampled   <- spatSample(hillshade, 2e5, na.rm = TRUE, method = "regular")[, 1]
  cut_offs  <- quantile(sampled, HS_CLIP, names = FALSE, na.rm = TRUE)
  hillshade <- clamp((hillshade - cut_offs[1]) / (cut_offs[2] - cut_offs[1]),
                     0, 1, values = TRUE)
  hillshade <- ifel(is.na(hillshade), 0.5, hillshade)
  N_STEPS <- 512
  heights <- seq(min(RELIEF_BREAKS), max(RELIEF_BREAKS), length.out = N_STEPS)
  ramp    <- col2rgb(RELIEF_COLS)
  lookup  <- lapply(1:3, function(channel)
    approx(RELIEF_BREAKS, ramp[channel, ], xout = heights, rule = 2)$y)
  capped <- clamp(relief, min(RELIEF_BREAKS), max(RELIEF_BREAKS), values = TRUE)
  step   <- round((capped - min(RELIEF_BREAKS)) /
                  (max(RELIEF_BREAKS) - min(RELIEF_BREAKS)) * (N_STEPS - 1)) + 1
  shading <- 1 + HS_STRENGTH * (hillshade - 0.5)
  channels <- lapply(1:3, function(channel) {
    flat   <- classify(step, cbind(seq_len(N_STEPS), lookup[[channel]]))
    shaded <- flat * shading
    clamp(shaded * (1 - RELIEF_VEIL) + 255 * RELIEF_VEIL, 0, 254, values = TRUE)
  })
  relief_rgb <- rast(channels)
  names(relief_rgb) <- c("red", "green", "blue")
  writeRaster(relief_rgb, RELIEF_RGB_CACHE, overwrite = TRUE, datatype = "INT1U",
              gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  writeLines(settings_now, SETTINGS_FILE)
  relief_rgb <- rast(RELIEF_RGB_CACHE)
}
coastline <- st_as_sf(as.contour(etopo, levels = 0)) |>
  st_set_crs(4326) |>
  st_transform(ALBERS)
window_box <- st_bbox(c(xmin = LON_MIN, xmax = LON_MAX,
                        ymin = LAT_MIN, ymax = LAT_MAX), crs = st_crs(4326))
WINDOW_TAG <- sprintf("%g_%g_%g_%g", LON_MIN, LON_MAX, LAT_MIN, LAT_MAX)
lakes_cache <- file.path(CACHE, paste0("lakes_region_", WINDOW_TAG, ".gpkg"))
if (file.exists(lakes_cache)) {
  lakes <- st_read(lakes_cache, quiet = TRUE)
} else {
  message("Downloading Natural Earth lakes (once) ...")
  lakes <- quietly(
    ne_download(scale = 10, type = "lakes", category = "physical",
                returnclass = "sf") |>
      st_make_valid() |>
      st_crop(window_box))
  st_write(lakes, lakes_cache, quiet = TRUE, delete_dsn = TRUE)
}
rivers_cache <- file.path(CACHE, paste0("rivers_region_", WINDOW_TAG, ".gpkg"))
if (file.exists(rivers_cache)) {
  rivers <- st_read(rivers_cache, quiet = TRUE)
} else {
  message("Downloading Natural Earth rivers (once) ...")
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
countries <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_make_valid()

# Build country outlines and define the study region

NORTHERN_STUDY_COUNTRIES <- c("Angola", "Zambia", "Malawi", "Mozambique")
NORTHERN_NEIGHBOURS      <- c("Democratic Republic of the Congo",
                              "United Republic of Tanzania")
pieces <- quietly(st_cast(st_union(
  countries |> filter(admin %in% NORTHERN_STUDY_COUNTRIES)), "POLYGON"))
piece_area <- vapply(seq_along(pieces), function(i) {
  b <- st_bbox(pieces[i]); as.numeric((b[3] - b[1]) * (b[4] - b[2]))
}, numeric(1))
mainland <- pieces[which.max(piece_area)]
neighbours <- quietly(
  st_union(countries |> filter(admin %in% NORTHERN_NEIGHBOURS)))
study_area_north <- quietly(
  st_boundary(mainland) |>
    st_intersection(st_buffer(neighbours, 0.2)) |>
    st_collection_extract("LINESTRING") |>
    st_union() |>
    st_cast("MULTILINESTRING") |>
    st_line_merge())
BORDER_INSET <- 0.02
africa <- countries |> filter(continent == "Africa")
borders_internal <- quietly(
  st_boundary(africa) |>
    st_intersection(st_buffer(st_union(africa), -BORDER_INSET)) |>
    st_collection_extract("LINESTRING") |>
    st_union() |>
    st_crop(window_box))
WC_RES          <- 5
WRZ_THRESHOLD   <- 66
SRZ_THRESHOLD   <- 33
MIN_ANNUAL_MM   <- 20
SMOOTH_WINDOW   <- 5
MIN_LINE_KM     <- 100

# Calculate rainfall zones and isohyets

message("Working out the rainfall seasons (WorldClim 2.1) ...")
prec <- crop(
  geodata::worldclim_global(var = "prec", res = WC_RES, path = CACHE),
  ext(LON_MIN, LON_MAX, LAT_MIN, LAT_MAX))
annual_mm    <- sum(prec)
winter_share <- 100 * sum(prec[[4:9]]) / annual_mm
winter_share <- ifel(annual_mm < MIN_ANNUAL_MM, NA, winter_share)
winter_smooth <- mask(
  focal(winter_share, w = SMOOTH_WINDOW, fun = mean, na.rm = TRUE),
  winter_share)
rain_zone <- classify(winter_smooth,
                      rbind(c(-1, SRZ_THRESHOLD,   1),
                            c(SRZ_THRESHOLD, WRZ_THRESHOLD, 2),
                            c(WRZ_THRESHOLD, 101,  3)))
blobs     <- patches(ifel(rain_zone >= 2, 1, NA), directions = 8)
wrz_blobs <- unique(values(mask(blobs, ifel(rain_zone == 3, 1, NA)), na.rm = TRUE))
connected <- ifel(blobs %in% wrz_blobs, 1, NA)
dropped   <- as.numeric(global(ifel(rain_zone >= 2, 1, 0), "sum", na.rm = TRUE) -
                        global(connected, "sum", na.rm = TRUE))
rain_zone <- ifel(is.na(connected), 1, rain_zone)
guard <- quietly(st_buffer(st_union(st_as_sf(as.polygons(connected))), 0.2))
rain_lines_proj <- lapply(c(SRZ_THRESHOLD, WRZ_THRESHOLD), function(level) {
  parts <- quietly(
    st_as_sf(as.contour(winter_smooth, levels = level)) |>
      st_set_crs(4326) |>
      st_geometry() |>
      st_cast("LINESTRING"))
  parts <- parts[lengths(quietly(st_intersects(parts, guard))) > 0]
  parts <- st_transform(parts, ALBERS)
  parts[as.numeric(st_length(parts)) >= MIN_LINE_KM * 1000]
})
rain_lines_proj <- do.call(c, rain_lines_proj)
zone_share       <- as.data.frame(freq(rain_zone))
zone_share$zone  <- c("summer", "year-round", "winter")[zone_share$value]
zone_share$pct   <- 100 * zone_share$count / sum(zone_share$count)
message(sprintf("  land in each rainfall zone: %s",
                paste(sprintf("%s %.1f%%", zone_share$zone, zone_share$pct),
                      collapse = ", ")))
message(sprintf("  %.0f cells of unconnected northern wet-season land returned to the summer zone",
                dropped))
message(sprintf("  %d isohyet segments drawn, %.0f km in total",
                length(rain_lines_proj),
                sum(as.numeric(st_length(rain_lines_proj))) / 1000))
lakes_proj      <- st_transform(lakes,            ALBERS)
rivers_proj     <- st_transform(rivers,           ALBERS)
north_edge_proj <- st_transform(study_area_north, ALBERS)
borders_proj    <- st_transform(borders_internal, ALBERS)

# Style the map and assemble the figure

apply_veil <- function(cols, amount = RELIEF_VEIL) {
  m <- col2rgb(cols) * (1 - amount) + 255 * amount
  rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}
RELIEF_COLS_VEILED <- apply_veil(RELIEF_COLS)
DEEP_SEA_VEILED    <- RELIEF_COLS_VEILED[1]
LAKE_FILL    <- "#7FAECF"
LAKE_OUTLINE <- "#5589A9"
RIVER_BLUE   <- "#6A9CBF"
COAST_GREY   <- "grey20"
STUDY_GREY   <- "grey45"
BORDER_GREY  <- "grey55"
BORDER_WIDTH <- 0.25
RAIN_COL        <- STUDY_GREY
RAIN_LINE_WIDTH <- 1.5
RAIN_LABEL_SIZE <- 5.3
rain_labels <- tibble::tribble(
  ~zone,        ~label,  ~lon,   ~lat,
  "winter",     "WRZ",    19.00, -33.00,
  "year-round", "YRZ",    21.45, -33.10,
  "summer",     "SRZ",    22.85, -31.25
)
rain_label_xy <- st_coordinates(st_transform(
  st_as_sf(rain_labels, coords = c("lon", "lat"), crs = 4326, remove = FALSE),
  ALBERS))
rain_labels <- rain_labels |>
  mutate(x = rain_label_xy[, 1], y = rain_label_xy[, 2])
landed <- c("summer", "year-round", "winter")[
  terra::extract(rain_zone, as.matrix(rain_labels[, c("lon", "lat")]))[, 1]]
message(sprintf("  rainfall names: %s",
                paste(sprintf("%s -> %s", rain_labels$zone,
                              ifelse(is.na(landed), "off the grid", landed)),
                      collapse = ", ")))
if (!identical(landed, rain_labels$zone)) {
  message("  NOTE: a rainfall name is not sitting in its own zone - ",
          "move it in `rain_labels` above.")
}
LABEL_HALF_M_PER_CHAR <- 28000
rain_anchors <- rain_labels |>
  mutate(half = LABEL_HALF_M_PER_CHAR *
           vapply(strsplit(label, "\n"), function(l) max(nchar(l)), numeric(1))) |>
  rowwise() |>
  reframe(x = x + c(-half, 0, half), y = y, site_name = "")
sites_labels <- bind_rows(sites_labels, rain_anchors)
rain_gap_km <- (as.numeric(st_distance(
  st_as_sf(rain_labels, coords = c("x", "y"), crs = ALBERS),
  st_union(rain_lines_proj))) -
  LABEL_HALF_M_PER_CHAR * nchar(rain_labels$label)) / 1000
message(sprintf("  rainfall names clear of the nearest isohyet: %s",
                paste(sprintf("%s %.0f km", rain_labels$label, rain_gap_km),
                      collapse = ", ")))
if (any(rain_gap_km < 0)) {
  message("  NOTE: an isohyet runs through a rainfall name - ",
          "move it in `rain_labels` above.")
}
apart_km <- function(a, b) sqrt((rain_labels$x[a] - rain_labels$x[b])^2 +
                                (rain_labels$y[a] - rain_labels$y[b])^2) / 1000
message(sprintf("  rainfall names apart: WRZ-YRZ %.0f km, YRZ-SRZ %.0f km",
                apart_km(1, 2), apart_km(2, 3)))
figure_2 <- ggplot() +
  geom_spatraster_rgb(data = relief_rgb, maxcell = 1.5e7) +
  geom_sf(data = coastline, colour = COAST_GREY, linewidth = 0.4) +
  geom_sf(data = borders_proj, colour = BORDER_GREY,
          linewidth = BORDER_WIDTH) +
  geom_sf(data = lakes_proj, fill = LAKE_FILL, colour = LAKE_OUTLINE,
          linewidth = 0.3) +
  geom_sf(data = rivers_proj, colour = RIVER_BLUE, linewidth = 0.5) +
  geom_sf(data = rain_lines_proj, colour = RAIN_COL,
          linewidth = RAIN_LINE_WIDTH, linetype = "dotted") +
  geom_text(data = rain_labels, aes(x = x, y = y, label = label),
            colour = RAIN_COL, fontface = "bold", size = RAIN_LABEL_SIZE) +
  geom_sf(data = north_edge_proj, colour = STUDY_GREY,
          linewidth = 1.5, linetype = "dashed", lineend = "round") +
  geom_sf(data = sites_sf, shape = 21, fill = "black", colour = "white",
          size = 3.9, stroke = 1.0) +
  geom_text_repel(
    data = sites_labels, aes(x = x, y = y, label = site_name),
    size = 4.2, colour = "grey10", fontface = "bold",
    bg.color = "white", bg.r = 0.16,
    min.segment.length = 0, segment.size = 0.4, segment.colour = "grey30",
    box.padding = 0.3, point.padding = 0.5,
    force = 2.0, force_pull = 1.0, max.overlaps = Inf,
    max.iter = 200000, max.time = Inf, seed = 1
  ) +
  geom_point(data = data.frame(x = X_RANGE[1], y = Y_RANGE[1], height = 0),
             aes(x = x, y = y, fill = height),
             shape = 21, size = 0, stroke = 0, alpha = 0) +
  scale_fill_gradientn(
    colours = RELIEF_COLS_VEILED,
    values  = scales::rescale(RELIEF_BREAKS),
    limits  = range(RELIEF_BREAKS),
    breaks  = c(-4000, -2000, 0, 1500, 3000),
    labels  = function(v) format(v, big.mark = ",", trim = TRUE),
    name    = "Elevation / sea-floor depth (m)\n(ETOPO 2022, NOAA NCEI)",
    guide   = guide_colourbar(
      frame.colour = "grey35", ticks.colour = "grey35",
      theme = theme(legend.title       = element_text(hjust = 0.5),
                    legend.key.width   = unit(330, "pt"),
                    legend.key.height  = unit(10,  "pt")))
  ) +
  coord_sf(crs = ALBERS, xlim = X_RANGE, ylim = Y_RANGE,
           expand = FALSE, datum = 4326) +
  annotation_scale(location = "bl", width_hint = 0.22, text_cex = 0.7,
                   line_width = 0.5, bar_cols = c("grey30", "white")) +
  annotation_north_arrow(location = "br", which_north = "true",
                         height = unit(0.95, "cm"), width = unit(0.95, "cm"),
                         pad_x  = unit(0.25, "cm"), pad_y = unit(0.25, "cm"),
                         style  = north_arrow_fancy_orienteering()) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background      = element_rect(fill = DEEP_SEA_VEILED, colour = NA),
    panel.grid            = element_line(colour = "#55555540", linewidth = 0.25),
    panel.border          = element_rect(colour = "grey35", fill = NA,
                                         linewidth = 0.6),
    legend.position       = "bottom",
    legend.title.position = "top",
    legend.title          = element_text(size = 11, face = "bold"),
    legend.text           = element_text(size = 10),
    axis.text             = element_text(size = 10, colour = "grey40")
  )
FIG_WIDTH  <- 10
FIG_HEIGHT <- 10.8

# Save the site map as PNG and PDF

ggsave(OUT_PNG, figure_2, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 300,
       bg = "white", device = ragg::agg_png)
pdf_tmp <- tempfile(fileext = ".pdf")
ggsave(pdf_tmp, figure_2, width = FIG_WIDTH, height = FIG_HEIGHT, bg = "white")
if (!file.copy(pdf_tmp, OUT_PDF, overwrite = TRUE)) {
  stop("Could not write ", OUT_PDF, ". Is it open in another program?",
       call. = FALSE)
}
unlink(pdf_tmp)
message("Saved ", OUT_PNG)
message("Saved ", OUT_PDF)
message("Done.")
