# Load packages and locate the analysis dataset

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})
FILE_ANALYSIS <- file.path("data", "PaperA_analysis_dataset.csv")
OUT_SITES   <- file.path("data", "PaperA_fig_sites.csv")
OUT_BINS    <- file.path("data", "PaperA_fig_sample_bins.csv")
OUT_RECORDS <- file.path("data", "PaperA_fig_plant_records.csv")
if (!file.exists(FILE_ANALYSIS)) {
  stop("Cannot find ", FILE_ANALYSIS, ". Run 01_data_wrangling.R first.")
}

# Select dated records and build the figure input tables

analysis <- read_csv(FILE_ANALYSIS, guess_max = Inf, show_col_types = FALSE)
TARGET_BINS <- c("500_BC_to_AD_0", "AD_1_to_500", "AD_501_to_1000", "AD_1001_to_1500")
FIGURE_PLANT_GROUPS <- c(
  "Introduced from wider Africa",
  "Introduced from outside Africa",
  "Native subsistence plant"
)
in_scope <- analysis |>
  filter(status == "binned", chron_cat %in% TARGET_BINS)
stopifnot(all(in_scope$in_dated_sites_cohort))
fig_sites <- in_scope |>
  distinct(site_name, country, latitude, longitude, kg_zone) |>
  arrange(country, site_name)
stopifnot(!any(duplicated(fig_sites$site_name)))
fig_sample_bins <- in_scope |>
  distinct(site_name, sample_number, chron_cat, bin_label, bin_index, phase_final) |>
  arrange(site_name, sample_number, bin_index)
stopifnot(!any(duplicated(fig_sample_bins[, c("sample_number", "chron_cat")])))
fig_plant_records <- in_scope |>
  filter(
    plant_group %in% FIGURE_PLANT_GROUPS,
    !is.na(accepted_name_short)
  ) |>
  distinct(
    site_name, sample_number,
    accepted_name_short,
    common_name = top_common_name,
    family, plant_group,
    taxon_merged_figs45, taxon_plot_id_fig3, taxon_plot_name_fig3,
    is_post_1500_taxon, is_genus_exclude,
    ecocrop_kg_codes, ecocrop_alt_min_m, ecocrop_alt_max_m
  ) |>
  arrange(site_name, sample_number, accepted_name_short)
stopifnot(
  all(fig_plant_records$sample_number %in% fig_sample_bins$sample_number),
  all(fig_plant_records$site_name %in% fig_sites$site_name),
  all(fig_sample_bins$site_name %in% fig_sites$site_name)
)
write_csv(fig_sites,         OUT_SITES)
write_csv(fig_sample_bins,   OUT_BINS)
write_csv(fig_plant_records, OUT_RECORDS)

# Derive crop occurrences and Figure 3 records

INTRODUCED <- c("Introduced from wider Africa", "Introduced from outside Africa")
crop_occurrences <- fig_plant_records |>
  filter(plant_group %in% INTRODUCED,
         !is_post_1500_taxon, !is_genus_exclude) |>
  inner_join(fig_sample_bins, by = c("site_name", "sample_number"),
             relationship = "many-to-many") |>
  distinct(taxon_merged_figs45, site_name, chron_cat)
NATIVE_CROPS_KEEP <- c("Citrullus lanatus", "Lagenaria siceraria")
FIG3_PHASES <- c("Later Stone Age", "Early Iron Age",
                 "Middle Iron Age", "Late Iron Age", "Iron Age")
fig3_records <- fig_plant_records |>
  filter(!is_post_1500_taxon,
         plant_group %in% INTRODUCED | taxon_plot_id_fig3 %in% NATIVE_CROPS_KEEP) |>
  inner_join(fig_sample_bins, by = c("site_name", "sample_number"),
             relationship = "many-to-many") |>
  filter(phase_final %in% FIG3_PHASES)

# Print dataset summaries and check expected counts

line <- function() cat(strrep("=", 74), "\n")
cat("\n"); line()
cat("PAPER A FIGURE DATASETS\n"); line()
cat(sprintf("%-30s %5d rows x %2d cols  -> %s\n", "Sites",
            nrow(fig_sites), ncol(fig_sites), OUT_SITES))
cat(sprintf("%-30s %5d rows x %2d cols  -> %s\n", "Samples x time bins",
            nrow(fig_sample_bins), ncol(fig_sample_bins), OUT_BINS))
cat(sprintf("%-30s %5d rows x %2d cols  -> %s\n", "Plant finds",
            nrow(fig_plant_records), ncol(fig_plant_records), OUT_RECORDS))
cat("\nSites per time bin (Figure 4's denominators):\n")
fig_sample_bins |>
  distinct(site_name, bin_label, bin_index) |>
  count(bin_index, bin_label, name = "sites") |>
  arrange(bin_index) |> select(bin_label, sites) |>
  as.data.frame() |> print(row.names = FALSE)
cat(sprintf("\nFigures 4/5/S2: %d crop-at-a-site-in-a-bin occurrences, %d crops\n",
            nrow(crop_occurrences), n_distinct(crop_occurrences$taxon_merged_figs45)))
cat(sprintf("Figure 3:       %d plants across %d sites\n",
            n_distinct(fig3_records$taxon_plot_id_fig3),
            n_distinct(fig3_records$site_name)))
cat("Figure 3 plants: ",
    paste(sort(unique(fig3_records$taxon_plot_name_fig3)), collapse = ", "), "\n")
expected <- list(sites = 35L, sample_bins = 121L, plant_records = 230L,
                 crop_occurrences = 70L, crop_taxa = 10L, fig3_taxa = 12L,
                 fig3_sites = 30L)
observed <- list(sites = nrow(fig_sites), sample_bins = nrow(fig_sample_bins),
                 plant_records = nrow(fig_plant_records),
                 crop_occurrences = nrow(crop_occurrences),
                 crop_taxa = n_distinct(crop_occurrences$taxon_merged_figs45),
                 fig3_taxa = n_distinct(fig3_records$taxon_plot_id_fig3),
                 fig3_sites = n_distinct(fig3_records$site_name))
ok <- TRUE
for (nm in names(expected)) {
  if (observed[[nm]] != expected[[nm]]) {
    ok <- FALSE
    warning(sprintf("'%s' is %d but %d was expected for this database version. ",
                    nm, observed[[nm]], expected[[nm]]),
            "If the database has been updated this is fine - but check the figures.")
  }
}