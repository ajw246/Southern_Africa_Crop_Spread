# Load packages and validate figure inputs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(ggnewscale)
})

# Set output paths and climate zone order

FILE_RECORDS <- file.path("data", "PaperA_fig_plant_records.csv")
FILE_BINS    <- file.path("data", "PaperA_fig_sample_bins.csv")
FILE_SITES   <- file.path("data", "PaperA_fig_sites.csv")
OUT_PNG <- file.path("figures", "figure_4_crops_and_climate_zones.png")
OUT_PDF <- file.path("figures", "figure_4_crops_and_climate_zones.pdf")
for (f in c(FILE_RECORDS, FILE_BINS, FILE_SITES)) {
  if (!file.exists(f)) {
    stop("Cannot find ", f, ". Run 02_figure_datasets.R first.", call. = FALSE)
  }
}

# Prepare crop records and climate zone summaries

plant_records <- read_csv(FILE_RECORDS, show_col_types = FALSE)
sample_bins   <- read_csv(FILE_BINS,    show_col_types = FALSE)
sites         <- read_csv(FILE_SITES,   show_col_types = FALSE)
PLANT_GROUP_ORDER <- c("Introduced from wider Africa",
                       "Introduced from outside Africa")
KG_ZONE_ORDER <- c("BWh", "BWk", "BSh", "BSk",
                   "Cwa", "Cwb", "Cwc",
                   "Cfa", "Cfb", "Cfc",
                   "Csa", "Csb", "Csc")
ENV_OF_ZONE <- c("BWh" = "Arid",
                 "BSh" = "Steppe",   "BSk" = "Steppe",
                 "Cwa" = "Temperate", "Cwb" = "Temperate", "Cfa" = "Temperate")
SORGHUM_MERGED_LABEL <- "Sorghum*"
crop_records <- plant_records |>
  filter(plant_group %in% PLANT_GROUP_ORDER,
         !is_post_1500_taxon,
         !is_genus_exclude)
taxon_labels <- crop_records |>
  group_by(taxon = taxon_merged_figs45) |>
  summarise(plant_group = first(plant_group),
            taxon_label = first(common_name),
            n_labels    = n_distinct(common_name),
            n_groups    = n_distinct(plant_group),
            .groups     = "drop") |>
  mutate(taxon_label = if_else(taxon == "Sorghum bicolor",
                               SORGHUM_MERGED_LABEL, taxon_label))
stopifnot(all(taxon_labels$n_labels == 1), all(taxon_labels$n_groups == 1))
taxon_labels <- select(taxon_labels, taxon, plant_group, taxon_label)
occurrences <- crop_records |>
  select(taxon = taxon_merged_figs45, site_name, sample_number) |>
  inner_join(select(sample_bins, site_name, sample_number, bin_label, bin_index),
             by = c("site_name", "sample_number"),
             relationship = "many-to-many") |>
  inner_join(select(sites, site_name, kg_zone), by = "site_name") |>
  distinct(taxon, site_name, kg_zone, bin_label, bin_index)
stopifnot(!any(is.na(occurrences$kg_zone)))
sites_per_cell <- occurrences |>
  group_by(taxon, kg_zone) |>
  summarise(n_sites = n_distinct(site_name), .groups = "drop")
sites_per_slice <- occurrences |>
  group_by(taxon, kg_zone, bin_label, bin_index) |>
  summarise(n_slice = n_distinct(site_name), .groups = "drop")
zone_totals <- count(sites, kg_zone, name = "n_total_sites")
zones_with_finds <- unique(sites_per_cell$kg_zone)
zones_dropped    <- setdiff(zone_totals$kg_zone, zones_with_finds)
if (length(zones_dropped) > 0) {
  for (z in zones_dropped) {
    message("Zone dropped - has dated sites but no introduced crop find: ", z,
            " (", zone_totals$n_total_sites[zone_totals$kg_zone == z],
            " site(s): ",
            paste(sites$site_name[sites$kg_zone == z], collapse = ", "), ")")
  }
}
active_zones <- KG_ZONE_ORDER[KG_ZONE_ORDER %in% zones_with_finds]
taxon_order <- sites_per_cell |>
  group_by(taxon) |>
  summarise(total_sites = sum(n_sites), .groups = "drop") |>
  left_join(taxon_labels, by = "taxon") |>
  mutate(plant_group = factor(plant_group, levels = PLANT_GROUP_ORDER)) |>
  arrange(plant_group, desc(total_sites), taxon_label)
n_taxa_total   <- nrow(taxon_order)
n_outside_taxa <- sum(taxon_order$plant_group == "Introduced from outside Africa")
separator_y    <- n_outside_taxa + 0.5
ecocrop_codes <- crop_records |>
  distinct(taxon = taxon_merged_figs45, codes = ecocrop_kg_codes) |>
  filter(!is.na(codes))
stopifnot(!any(duplicated(ecocrop_codes$taxon)))
zone_is_listed <- function(codes, zone) {
  if (is.na(codes) || startsWith(codes, "(")) return(NA)
  listed <- trimws(unlist(strsplit(codes, "[;,]")))
  zone %in% listed[nzchar(listed)]
}
tile_draw_order <- taxon_order$taxon[
  order(match(taxon_order$plant_group, rev(PLANT_GROUP_ORDER)))
]
heatmap_grid <- expand.grid(kg_zone = active_zones,
                            taxon   = tile_draw_order,
                            stringsAsFactors = FALSE) |>
  as_tibble() |>
  select(taxon, kg_zone) |>
  left_join(taxon_order,    by = "taxon") |>
  left_join(sites_per_cell, by = c("taxon", "kg_zone")) |>
  left_join(ecocrop_codes,  by = "taxon") |>
  mutate(
    n_sites = coalesce(n_sites, 0L),
    listed  = mapply(zone_is_listed, codes, kg_zone),
    ecocrop_status = case_when(
      is.na(listed) ~ "No ECOCROP data",
      listed        ~ "compatible",
      TRUE          ~ "incompatible"
    )
  )

# Build the heatmap, climate suitability and environmental bands

R_MIN <- 0.18
R_MAX <- 0.40
n_max <- max(heatmap_grid$n_sites)
cell_totals <- heatmap_grid |>
  filter(n_sites > 0) |>
  mutate(r = R_MIN + (R_MAX - R_MIN) * (n_sites - 1) / max(1, n_max - 1))
pie_slices <- cell_totals |>
  select(taxon, kg_zone, r) |>
  inner_join(sites_per_slice, by = c("taxon", "kg_zone")) |>
  arrange(taxon, kg_zone, bin_index) |>
  group_by(taxon, kg_zone) |>
  mutate(prop      = n_slice / sum(n_slice),
         arc_end   = cumsum(prop) * 2 * pi,
         arc_start = lag(arc_end, default = 0)) |>
  ungroup()
stopifnot(nrow(pie_slices) > 0)
ECOCROP_COLS <- c("compatible"      = "#B8D8D8",
                  "incompatible"    = "white",
                  "No ECOCROP data" = "#E8E8E8")
BIN_COLS <- c("500 BCE - 0 CE" = "grey85",
              "1 - 500 CE"     = "grey70",
              "501 - 1000 CE"  = "grey40",
              "1001 - 1500 CE" = "grey20")
ENV_COLS <- c("Arid"      = "#E74C3C",
              "Steppe"    = "#E29A3C",
              "Temperate" = "#3C9A5F",
              "Other"     = "#BDBDBD")
ENV_ALPHA <- 0.30
bin_levels <- sample_bins |>
  distinct(bin_label, bin_index) |>
  arrange(bin_index) |>
  pull(bin_label)
stopifnot(setequal(bin_levels, names(BIN_COLS)))
as_plot_factors <- function(d) {
  d |> mutate(
    taxon_label = factor(taxon_label, levels = rev(taxon_order$taxon_label)),
    kg_zone     = factor(kg_zone,     levels = active_zones)
  )
}
heatmap_plot <- heatmap_grid |>
  as_plot_factors() |>
  mutate(ecocrop_status = factor(ecocrop_status, levels = names(ECOCROP_COLS)))
cell_totals <- as_plot_factors(cell_totals)
pie_slices <- pie_slices |>
  left_join(select(taxon_order, taxon, taxon_label), by = "taxon") |>
  as_plot_factors() |>
  mutate(bin_label = factor(bin_label, levels = bin_levels))
BAND_YMIN <- -0.65
BAND_YMAX <-  0.45
TOP_YMIN  <- n_taxa_total + 0.62
TOP_YMAX  <- n_taxa_total + 0.86
env_band <- tibble(
  kg_zone = active_zones,
  x       = seq_along(active_zones),
  env     = unname(ENV_OF_ZONE[active_zones]),
  label   = sprintf("%s\n(n=%d)", active_zones,
                    zone_totals$n_total_sites[match(active_zones,
                                                    zone_totals$kg_zone)])
) |>
  mutate(env = coalesce(env, "Other"))
env_blocks <- env_band |>
  arrange(x) |>
  mutate(block_id = cumsum(env != lag(env, default = ""))) |>
  group_by(block_id, env) |>
  summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, .groups = "drop")
status_present <- levels(droplevels(heatmap_plot$ecocrop_status))
status_outline <- unname(c("compatible"      = "white",
                           "incompatible"    = "grey50",
                           "No ECOCROP data" = "grey50")[status_present])
figure_4 <- ggplot(heatmap_plot, aes(x = kg_zone, y = taxon_label)) +
  annotate("rect", xmin = 0.5, xmax = length(active_zones) + 0.5,
           ymin = 0.5, ymax = separator_y, fill = "#F0F0F0") +
  annotate("rect", xmin = 0.5, xmax = length(active_zones) + 0.5,
           ymin = separator_y, ymax = n_taxa_total + 0.5, fill = "#FAFAFA") +
  geom_tile(aes(fill = ecocrop_status), colour = "white", linewidth = 0.7) +
  scale_fill_manual(
    values = ECOCROP_COLS,
    name   = "FAO ECOCROP\nclimate suitability",
    guide  = guide_legend(override.aes = list(colour = status_outline),
                          order = 1)
  ) +
  new_scale_fill() +
  geom_arc_bar(
    data = pie_slices,
    aes(x0 = kg_zone, y0 = taxon_label, r0 = 0, r = r,
        start = arc_start, end = arc_end, fill = bin_label),
    inherit.aes = FALSE, colour = "black", linewidth = 0.3
  ) +
  scale_fill_manual(
    values = BIN_COLS,
    name   = "500-year bin",
    guide  = guide_legend(override.aes = list(colour = "white"), order = 2)
  ) +
  geom_text(data = cell_totals, aes(label = n_sites),
            size = 5, colour = "white", fontface = "bold") +
  annotate("segment", x = 0.5, xend = length(active_zones) + 0.5,
           y = separator_y, yend = separator_y,
           colour = "grey45", linewidth = 0.5) +
  new_scale_fill() +
  geom_rect(data = env_blocks,
            aes(xmin = xmin, xmax = xmax,
                ymin = BAND_YMIN, ymax = BAND_YMAX, fill = env),
            inherit.aes = FALSE, alpha = ENV_ALPHA,
            colour = "white", linewidth = 0.7) +
  geom_rect(data = env_blocks,
            aes(xmin = xmin, xmax = xmax,
                ymin = TOP_YMIN, ymax = TOP_YMAX, fill = env),
            inherit.aes = FALSE, alpha = ENV_ALPHA,
            colour = "white", linewidth = 0.7) +
  scale_fill_manual(
    values = ENV_COLS,
    name   = "Environment",
    breaks = c("Arid", "Steppe", "Temperate"),
    guide  = guide_legend(override.aes = list(alpha = ENV_ALPHA,
                                              colour = "white"), order = 3)
  ) +
  geom_text(data = env_band,
            aes(x = x, y = (BAND_YMIN + BAND_YMAX) / 2, label = label),
            inherit.aes = FALSE, size = 3.8, lineheight = 0.9,
            colour = "grey15") +
  coord_fixed(ratio = 1, clip = "off") +
  scale_x_discrete(labels = NULL, drop = FALSE,
                   expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(expand = expansion(add = c(1.6, 1.0))) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title.y          = element_blank(),
    axis.text.x           = element_blank(),
    axis.ticks.x          = element_blank(),
    axis.text.y           = element_text(size = 11.5, hjust = 1),
    axis.title.x          = element_text(size = 11, margin = margin(t = 6)),
    panel.grid            = element_blank(),
    panel.border          = element_rect(colour = "white", fill = NA,
                                         linewidth = 0),
    legend.position       = "right",
    legend.key.size       = unit(0.9, "lines"),
    legend.title          = element_text(size = 10.5),
    legend.text           = element_text(size = 9.5),
    plot.title            = element_text(size = 14, face = "bold"),
    plot.subtitle         = element_text(size = 11, colour = "grey30",
                                         lineheight = 1.3),
    plot.margin           = margin(t = 8, r = 8, b = 10, l = 8)
  ) +
  labs(
    x        = "Koppen-Geiger climate zone (Beck et al. 2023)  |  n = total sites",
    title    = "Introduced taxa across climate zones",
    subtitle = "Southern African archaeobotanical sites, 500 BCE - 1500 CE "
  )

# Save Figure 4 as PNG and PDF

ggsave(OUT_PNG, figure_4, width = 7, height = 10, units = "in",
       dpi = 300, device = "png")
