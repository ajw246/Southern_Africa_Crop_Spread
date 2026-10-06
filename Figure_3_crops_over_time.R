# Load packages and validate figure inputs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(scales)
})

# Set output paths and figure parameters

FILE_RECORDS <- file.path("data", "PaperA_fig_plant_records.csv")
FILE_BINS    <- file.path("data", "PaperA_fig_sample_bins.csv")
OUT_PNG <- file.path("figures", "figure_3_crops_over_time.png")
OUT_PDF <- file.path("figures", "figure_3_crops_over_time.pdf")
for (f in c(FILE_RECORDS, FILE_BINS)) {
  if (!file.exists(f)) {
    stop("Cannot find ", f, ". Run 02_figure_datasets.R first.", call. = FALSE)
  }
}

# Prepare crop counts across archaeological periods

plant_records <- read_csv(FILE_RECORDS, show_col_types = FALSE)
sample_bins   <- read_csv(FILE_BINS,    show_col_types = FALSE)
TARGET_BINS <- c("500_BC_to_AD_0", "AD_1_to_500", "AD_501_to_1000", "AD_1001_to_1500")
FIG3_PHASES <- c("Later Stone Age", "Early Iron Age",
                 "Middle Iron Age", "Late Iron Age", "Iron Age")
NATIVE_CROPS_KEEP <- c("Citrullus lanatus", "Lagenaria siceraria")
INTRODUCED <- c("Introduced from wider Africa", "Introduced from outside Africa")
PHASE_COLOURS <- c("Later Stone Age" = "#661100",
                   "Iron Age"        = "#DDCC77")
fig3_records <- plant_records |>
  filter(
    !is_post_1500_taxon,
    plant_group %in% INTRODUCED | taxon_plot_id_fig3 %in% NATIVE_CROPS_KEEP
  ) |>
  inner_join(
    sample_bins,
    by = c("site_name", "sample_number"),
    relationship = "many-to-many"
  ) |>
  filter(
    phase_final %in% FIG3_PHASES,
    chron_cat   %in% TARGET_BINS
  ) |>
  mutate(
    phase_plot = if_else(phase_final == "Later Stone Age",
                         "Later Stone Age", "Iron Age"),
    family = if_else(is.na(family), "Unknown", family)
  )
bar_heights <- fig3_records |>
  group_by(chron_cat, phase_plot, taxon_plot_id_fig3, taxon_plot_name_fig3) |>
  summarise(n_sites = n_distinct(site_name), .groups = "drop")
taxon_order <- fig3_records |>
  group_by(taxon_plot_id_fig3) |>
  summarise(
    taxon_name  = first(taxon_plot_name_fig3),
    plant_group = first(plant_group),
    family      = min(family),
    first_bin   = min(match(chron_cat, TARGET_BINS)),
    total_sites = n_distinct(site_name),
    .groups     = "drop"
  ) |>
  arrange(first_bin, desc(total_sites), family, taxon_name)
axis_labels <- taxon_order |>
  mutate(
    colour = if_else(plant_group %in% INTRODUCED, "black", "grey60"),
    label  = sprintf("<span style='color:%s'>%s</span>", colour, taxon_name)
  )
axis_labels <- setNames(axis_labels$label, axis_labels$taxon_plot_id_fig3)
bar_heights <- bar_heights |>
  mutate(
    taxon_plot_id_fig3 = factor(taxon_plot_id_fig3,
                                levels = taxon_order$taxon_plot_id_fig3),
    chron_cat  = factor(chron_cat, levels = rev(TARGET_BINS)),
    phase_plot = factor(phase_plot, levels = c("Later Stone Age", "Iron Age"))
  )
figure_3 <- ggplot(
  bar_heights,
  aes(x = taxon_plot_id_fig3, y = n_sites, fill = phase_plot, group = phase_plot)
) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.5) +
  geom_col(
    colour = NA,
    width = 0.75,
    position = position_dodge2(preserve = "single", padding = 0)
  ) +
  scale_fill_manual(values = PHASE_COLOURS, name = "") +
  scale_x_discrete(labels = axis_labels) +
  facet_grid(
    chron_cat ~ .,
    labeller = labeller(chron_cat = c(
      `500_BC_to_AD_0` = "500 BCE–0 CE",
      AD_1_to_500      = "1–500 CE",
      AD_501_to_1000   = "501–1000 CE",
      AD_1001_to_1500  = "1001–1500 CE"
    ))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.16)),
    breaks = scales::breaks_pretty()
  ) +
  labs(x = NULL, y = "Number of sites at which taxa present", title = "") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_classic(base_size = 22) +
  theme(
    axis.text.x        = element_markdown(angle = 90, hjust = 1, vjust = 0.5,
                                          size = 12),
    axis.text.y        = element_text(size = 10),
    axis.title.y       = element_text(size = 12),
    text               = element_text(size = 12),
    panel.grid.major.y = element_line(colour = "grey98", linewidth = 0.3),
    panel.spacing      = unit(0.3, "lines"),
    legend.position    = "bottom",
    legend.box         = "vertical",
    legend.text        = element_text(size = 12),
    legend.title       = element_text(size = 12),
    legend.title.position = "top",
    legend.key.size    = unit(3, "mm"),
    legend.key.spacing.y = unit(0, "pt"),
    strip.text.y       = element_text(angle = 0, size = 12),
    strip.background   = element_rect(colour = NA, fill = NA)
  )

# Save Figure 3 as PNG and PDF

ggsave(OUT_PNG, figure_3, width = 200, height = 200, units = "mm",
       dpi = 300, device = "png")
ggsave(OUT_PDF, figure_3, width = 200, height = 200, units = "mm",
       device = grDevices::pdf, encoding = "WinAnsi.enc")
