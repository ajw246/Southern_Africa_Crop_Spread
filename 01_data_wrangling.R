# Load packages, define input files and validate required inputs

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(terra)
})
FILE_SOURCE    <- file.path("data", "PaperA_source_data.csv")
FILE_C14       <- file.path("data", "PaperA_radiocarbon_calibrated.csv")
FILE_EXCLUDED  <- file.path("data", "PaperA_excluded_sites.csv")
FILE_KG_RASTER <- file.path("data", "koppen_geiger", "1991_2020",
                            "koppen_geiger_0p00833333.tif")
FILE_OUTPUT    <- file.path("data", "PaperA_analysis_dataset.csv")
for (f in c(FILE_SOURCE, FILE_C14, FILE_EXCLUDED)) {
  if (!file.exists(f)) {
    stop("Required input file not found: ", f,
         "\nMake sure R's working directory is the folder containing this script",
         "\n(opening CropSpread_FigureCode.Rproj does that for you).")
  }
}
if (!file.exists(FILE_KG_RASTER)) {
  stop("The Koppen-Geiger climate map was not found at:\n  ", FILE_KG_RASTER,
       "\n\nRun this once to fetch it, then run this script again:\n",
       '  source("00_fetch_koppen_geiger.R")')
}

# Set chronology bins and crop classification keys

chron_bins_500 <- tibble(
  chron_cat = c("pre_500_BC", "500_BC_to_AD_0", "AD_1_to_500",
                "AD_501_to_1000", "AD_1001_to_1500", "AD_1501_to_2000"),
  bin_lower = c(-100000, -500,   1,  501, 1001, 1501),
  bin_upper = c(   -501,    0, 500, 1000, 1500, 2000)
)
TARGET_BINS <- c("500_BC_to_AD_0", "AD_1_to_500", "AD_501_to_1000", "AD_1001_to_1500")
BIN_LABELS  <- c(
  "500_BC_to_AD_0"  = "500 BCE - 0 CE",
  "AD_1_to_500"     = "1 - 500 CE",
  "AD_501_to_1000"  = "501 - 1000 CE",
  "AD_1001_to_1500" = "1001 - 1500 CE"
)
BIN_IDX <- setNames(seq_along(TARGET_BINS), TARGET_BINS)
figure3_taxa <- tribble(
  ~accepted_name_short,  ~common_name_expected, ~plant_group_expected,
  "Citrullus lanatus",   "Tsamma melon",        "Native subsistence plant",
  "Lagenaria siceraria", "Bottle gourd",        "Native subsistence plant",
  "Ricinus communis",    "Castor-oil plant",    "Introduced from wider Africa",
  "Sorghum bicolor",     "Sorghum",             "Introduced from wider Africa",
  "Cenchrus americanus", "Pearl millet",        "Introduced from wider Africa",
  "Eleusine coracana",   "Finger millet",       "Introduced from wider Africa",
  "Vigna unguiculata",   "Cowpea",              "Introduced from wider Africa",
  "Vigna subterranea",   "Bambara groundnut",   "Introduced from wider Africa",
  "Adansonia digitata",  "Baobab",              "Introduced from wider Africa",
  "Gossypium herbaceum", "Cotton",              "Introduced from wider Africa",
  "Vigna radiata",       "Mung bean",           "Introduced from outside Africa",
  "Panicum miliaceum",   "Broomcorn millet",    "Introduced from outside Africa"
)
POST_1500_TAXA <- c(
  "Zea mays",
  "Arachis hypogaea",
  "Phaseolus vulgaris",
  "Cucurbita maxima",
  "Helianthus annuus"
)
GENUS_EXCLUDE <- c("Cannabis")
GENUS_MERGE <- c("Sorghum" = "Sorghum bicolor")
taxon_collapse_key <- tribble(
  ~accepted_name_short, ~taxon_plot_id,     ~taxon_plot_name_override,
  "Sorghum",            "Sorghum combined", "Sorghum",
  "Sorghum bicolor",    "Sorghum combined", "Sorghum",
  "Ficus",              "Ficus combined",   "Fig",
  "Ficus ingens",       "Ficus combined",   "Fig"
)
ecocrop_kg_codes_key <- tribble(
  ~accepted_name_short,             ~ecocrop_kg_codes,
  "Adansonia digitata",             "Aw; BSh, BSk",
  "Amaranthus cruentus",            "Aw; Af, Am; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Arachis hypogaea",               "Aw; BSh, BSk; Cfa, Cfb",
  "Berchemia discolor",             "Aw; BSh, BSk; Csa, Csb",
  "Cenchrus americanus",            "Aw; Af, Am; BWh, BWk; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Citrullus lanatus",              "Aw; Af, Am; Cfa, Cfb; Csa, Csb; Cfb, Cfc",
  "Colophospermum mopane",          "Aw; BWh, BWk; BSh, BSk",
  "Eleusine coracana",              "Aw; BWh, BWk; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Gossypium herbaceum",            "Aw; Af, Am; BWh, BWk; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Helianthus annuus",              "Aw; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb; Cfb, Cfc; Dfa, Dfb, Dsa, Dsb; Dfa, Dfb, Dfc; Dwa, Dwb, Dwc",
  "Lagenaria siceraria",            "Aw; Af, Am; BWh, BWk; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Panicum miliaceum",              "Aw; BWh, BWk; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Parinari curatellifolia",        "Aw; Af, Am; BWh, BWk; BSh, BSk; Cfa, Cfb; Csa, Csb",
  "Ricinus communis",               "Aw; BSh, BSk; Cfa, Cfb; Csa, Csb",
  "Sclerocarya birrea",             "Aw; Af, Am; BSh, BSk; Csa, Csb",
  "Sclerocarya birrea subsp. afra", "Aw; Af, Am; BSh, BSk; Csa, Csb",
  "Sorghum bicolor",                "Aw; BSh, BSk; Csa, Csb",
  "Tylosema esculentum",            "Aw; Af, Am; BWh, BWk; BSh, BSk",
  "Uapaca kirkiana",                "Aw",
  "Vigna radiata",                  "Aw; Cfa, Cfb",
  "Vigna subterranea",              "Aw; Csa, Csb",
  "Vigna unguiculata",              "Aw; Af, Am; BSh, BSk; Cfa, Cfb; Csa, Csb; Cwa, Cwb",
  "Ziziphus mucronata",             "Aw; BSh, BSk"
)
ecocrop_altitude <- tribble(
  ~accepted_name_short,   ~ecocrop_id, ~ecocrop_alt_min_m, ~ecocrop_alt_max_m,
  "Sorghum bicolor",           48747L,                  0,               2500,
  "Cenchrus americanus",        8418L,                  0,               1800,
  "Eleusine coracana",          5657L,                  0,               2500,
  "Vigna unguiculata",          2153L,                  0,               2000,
  "Vigna subterranea",         10830L,                200,               1400
)
NON_IDENTIFICATIONS <- c("Adiagnostic")
BIN_OVERLAP_THRESHOLD_YEARS <- 100
MAX_SPAN_YEARS              <- 1100
MAX_BINS_PER_SAMPLE         <- 2
EIA_LIA_BOUNDARY_CE <- 1000
POST_WINDOW_PHASES <- c("1800-present")
NORMALISE_STRATA <- FALSE

# Read and inspect the raw sample and plant records

message("1/9  Reading the raw data ...")
raw <- read_csv(FILE_SOURCE, guess_max = Inf, show_col_types = FALSE)
metadata_raw <- raw |>
  select(
    site_name, sample_number, country, latitude, longitude,
    described_chronology, stratigraphy, age_estimate,
    archbot_material, kg_zone_recorded
  ) |>
  distinct()
if (any(duplicated(metadata_raw$sample_number))) {
  stop("These samples have inconsistent information across their rows: ",
       paste(metadata_raw$sample_number[duplicated(metadata_raw$sample_number)],
             collapse = ", "))
}
plant_records <- raw |>
  filter(!is.na(species_archaeobotany)) |>
  select(sample_number, species_archaeobotany, presence_absence) |>
  distinct()
taxa_recorded <- raw |>
  filter(!is.na(species_archaeobotany)) |>
  select(
    species_archaeobotany, accepted_name_short, family, plant_group,
    top_common_name
  ) |>
  distinct()
if (any(duplicated(taxa_recorded$species_archaeobotany))) {
  stop("These recorded plant names carry more than one botanical identity: ",
       paste(taxa_recorded$species_archaeobotany[
         duplicated(taxa_recorded$species_archaeobotany)], collapse = ", "))
}
message(sprintf("     %d samples, %d plant records, %d distinct plant names",
                nrow(metadata_raw), nrow(plant_records), nrow(taxa_recorded)))

# Remove charcoal samples and tidy sample metadata

message("2/9  Checking no charcoal samples are present ...")
n_charcoal <- sum(metadata_raw$archbot_material == "Charcoal", na.rm = TRUE)
if (n_charcoal > 0) {
  stop(n_charcoal, " charcoal samples are present in ", FILE_SOURCE,
       ". They must be removed before this script is run.")
}
message(sprintf("     none found; %d samples remain, of these %d record their material",
                nrow(metadata_raw), sum(!is.na(metadata_raw$archbot_material))))
message("3/9  Tidying sample information ...")
clean_number <- function(x) suppressWarnings(as.numeric(str_trim(as.character(x))))
samples <- metadata_raw |>
  mutate(
    site_name_clean = str_squish(str_to_lower(site_name)),
    latitude            = clean_number(latitude),
    longitude           = clean_number(longitude),
    described_chronology_raw = described_chronology,
    described_chronology = case_when(
      described_chronology == "1800s"           ~ "1800-present",
      described_chronology == "Modern"          ~ "1800-present",
      described_chronology == "Early Stone Age" ~ "Earlier Stone Age",
      TRUE                                      ~ described_chronology
    )
  )
stopifnot(!any(duplicated(samples$sample_number)))

# Extract climate zones for archaeological sites

message("4/9  Reading climate zones from the Koppen-Geiger map ...")
KG_CODE_LOOKUP <- c(
  "Af","Am","Aw","BWh","BWk","BSh","BSk","Csa","Csb","Csc",
  "Cwa","Cwb","Cwc","Cfa","Cfb","Cfc","Dsa","Dsb","Dsc","Dsd",
  "Dwa","Dwb","Dwc","Dwd","Dfa","Dfb","Dfc","Dfd","ET","EF"
)
kg_raster <- terra::rast(FILE_KG_RASTER)
site_points <- samples |>
  distinct(site_name, site_name_clean, latitude, longitude) |>
  filter(!is.na(latitude), !is.na(longitude))
stopifnot(nrow(site_points) == n_distinct(samples$site_name))
kg_values <- terra::extract(
  kg_raster,
  as.matrix(site_points[, c("longitude", "latitude")])
)[[1]]
site_climate <- site_points |>
  mutate(
    kg_code = as.integer(kg_values),
    kg_zone = KG_CODE_LOOKUP[kg_code]
  ) |>
  select(site_name_clean, kg_zone)
stopifnot(!any(is.na(site_climate$kg_zone)))
samples <- samples |> left_join(site_climate, by = "site_name_clean")
n_checked <- sum(!is.na(samples$kg_zone_recorded))
n_diff    <- sum(samples$kg_zone != samples$kg_zone_recorded, na.rm = TRUE)
message(sprintf("     %d sites checked against the recorded zone; %d disagree",
                n_distinct(samples$site_name[!is.na(samples$kg_zone_recorded)]),
                n_distinct(samples$site_name[
                  !is.na(samples$kg_zone_recorded) &
                  samples$kg_zone != samples$kg_zone_recorded])))
if (n_diff > 0) {
  print(samples |> filter(kg_zone != kg_zone_recorded) |>
          distinct(site_name, kg_zone, kg_zone_recorded))
}

# Classify taxa and attach growing limits

message("5/9  Adding taxon groupings, flags and growing limits ...")
taxa_recorded <- taxa_recorded |>
  mutate(
    accepted_name_short = str_squish(str_replace_all(accepted_name_short,
                                                     "\u00a0", " "))
  )
ecocrop_key <- ecocrop_kg_codes_key |>
  full_join(ecocrop_altitude, by = "accepted_name_short") |>
  mutate(ecocrop_url = if_else(
    is.na(ecocrop_id), NA_character_,
    paste0("https://ecocrop.apps.fao.org/ecocrop/srv/en/dataSheet?id=", ecocrop_id)))
ecocrop_unmatched <- setdiff(ecocrop_key$accepted_name_short,
                             taxa_recorded$accepted_name_short)
if (length(ecocrop_unmatched) > 0) {
  stop("These ECOCROP entries in Section 2.5 match no taxon in the database: ",
       paste(ecocrop_unmatched, collapse = ", "))
}
taxa <- taxa_recorded |>
  left_join(ecocrop_key,        by = "accepted_name_short") |>
  left_join(taxon_collapse_key, by = "accepted_name_short") |>
  mutate(
    plant_group_filled = if_else(is.na(plant_group) | plant_group == "",
                                 "Other (mixed flora)", plant_group),
    is_figure3_taxon         = accepted_name_short %in% figure3_taxa$accepted_name_short,
    is_post_1500_taxon       = accepted_name_short %in% POST_1500_TAXA,
    is_genus_exclude         = accepted_name_short %in% GENUS_EXCLUDE,
    is_indet_label           = species_archaeobotany %in% NON_IDENTIFICATIONS |
                               is.na(accepted_name_short),
    taxon_merged_figs45 = recode(accepted_name_short, !!!GENUS_MERGE),
    taxon_plot_id_fig3   = coalesce(taxon_plot_id, accepted_name_short),
    taxon_plot_name_fig3 = coalesce(taxon_plot_name_override,
                                    top_common_name,
                                    accepted_name_short)
  ) |>
  select(-taxon_plot_id, -taxon_plot_name_override)
figure3_check <- figure3_taxa |>
  left_join(taxa |> distinct(accepted_name_short, plant_group),
            by = "accepted_name_short")
if (any(is.na(figure3_check$plant_group)) ||
    !all(figure3_check$plant_group_expected == figure3_check$plant_group)) {
  print(figure3_check |> filter(is.na(plant_group) |
                                plant_group_expected != plant_group))
  stop("The Figure 3 taxon classification written into Section 2.2 disagrees ",
       "with the plant_group column in ", FILE_SOURCE, ". Resolve before continuing.")
}
message("     Figure 3 taxon classification agrees with the raw data.")
message("6/9  Checking the plant records ...")
stopifnot(
  all(plant_records$presence_absence == 1),
  all(plant_records$species_archaeobotany %in% taxa$species_archaeobotany),
  all(plant_records$sample_number %in% samples$sample_number)
)
message(sprintf("     %d plant records across %d samples",
                nrow(plant_records), n_distinct(plant_records$sample_number)))

# Read radiocarbon dates and resolve sample ages

message("7/9  Reading radiocarbon dates ...")
radiocarbon <- read_csv(FILE_C14, guess_max = Inf, show_col_types = FALSE) |>
  mutate(
    sample_id     = str_remove(name, "^R_Date\\s+"),
    date_suffix   = str_extract(sample_id, "[ab]$"),
    sample_number = str_remove(sample_id, "[ab]$"),
    date_suffix   = if_else(is.na(date_suffix), "a", date_suffix)
  ) |>
  select(sample_number, date_suffix, lower_cal, upper_cal) |>
  pivot_wider(names_from = date_suffix, values_from = c(lower_cal, upper_cal),
              names_sep = "_", values_fn = first)
normalise_stratum <- function(s) {
  if (NORMALISE_STRATA) str_replace(str_trim(s), "(?<=^\\d{1,2})[a-zA-Z].*$", "")
  else str_trim(s)
}
sample_dates <- samples |>
  select(site_name, sample_number, age_estimate, stratigraphy,
         described_chronology) |>
  mutate(stratum = if_else(is.na(stratigraphy), NA_character_,
                           normalise_stratum(as.character(stratigraphy)))) |>
  left_join(radiocarbon, by = "sample_number") |>
  rowwise() |>
  mutate(
    c14_lo = { v <- c_across(any_of(c("lower_cal_a", "lower_cal_b")))
               v <- v[!is.na(v)]; if (!length(v)) NA_real_ else min(v) },
    c14_hi = { v <- c_across(any_of(c("upper_cal_a", "upper_cal_b")))
               v <- v[!is.na(v)]; if (!length(v)) NA_real_ else max(v) }
  ) |>
  ungroup() |>
  mutate(
    has_c14 = !is.na(c14_lo) | !is.na(c14_hi),
    c14_lo_eff = coalesce(c14_lo, c14_hi),
    c14_hi_eff = coalesce(c14_hi, c14_lo),
    c14_single_bound     = has_c14 & (is.na(c14_lo) | is.na(c14_hi)),
    c14_n_determinations = (!is.na(coalesce(lower_cal_a, upper_cal_a))) +
                           (if ("lower_cal_b" %in% names(radiocarbon))
                              (!is.na(coalesce(lower_cal_b, upper_cal_b))) else 0L)
  )
parse_age_estimate <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(c(NA_real_, NA_real_))
  x0 <- str_to_upper(str_squish(x))
  x0 <- str_replace_all(x0, "±", "+-")
  x0 <- str_replace_all(x0, "–|-", "-")
  x0 <- str_replace_all(x0, "\\bTO\\b", "-")
  x0 <- str_replace_all(x0, "\\bCIRCA\\b", "")
  x0 <- str_replace_all(x0, "\\bCA\\.?\\s*", "")
  x0 <- str_replace_all(x0, "\\bC\\.?\\s+", "")
  x0 <- str_squish(x0)
  bp_to_year <- function(bp) 1950 - bp
  if (str_detect(x0, "\\d+\\s*\\+-\\s*\\d+")) {
    nums <- as.numeric(str_extract_all(x0, "\\d+")[[1]])
    if (length(nums) >= 2) {
      centre <- nums[1]; err <- nums[2]
      if (str_detect(x0, "\\bBP\\b")) {
        y1 <- bp_to_year(centre + err); y2 <- bp_to_year(centre - err)
      } else if (str_detect(x0, "\\bBC\\b|\\bBCE\\b")) {
        y1 <- -(centre + err); y2 <- -(centre - err)
      } else {
        y1 <- centre - err; y2 <- centre + err
      }
      return(c(min(y1, y2), max(y1, y2)))
    }
  }
  if (str_detect(x0, "\\d+\\s*(BC|BCE).*(AD|CE)\\s*\\d+")) {
    nums <- as.numeric(str_extract_all(x0, "\\d+")[[1]])
    if (length(nums) >= 2) {
      return(c(min(-abs(nums[1]), abs(nums[2])), max(-abs(nums[1]), abs(nums[2]))))
    }
  }
  nums <- as.numeric(str_extract_all(x0, "\\d+")[[1]])
  if (length(nums) == 0) return(c(NA_real_, NA_real_))
  if (str_detect(x0, "\\bBP\\b")) {
    if (length(nums) == 1) { y <- bp_to_year(nums[1]); return(c(y, y)) }
    return(c(min(bp_to_year(nums[1]), bp_to_year(nums[2])),
             max(bp_to_year(nums[1]), bp_to_year(nums[2]))))
  }
  if (str_detect(x0, "\\bBC\\b|\\bBCE\\b")) {
    if (length(nums) == 1) { y <- -abs(nums[1]); return(c(y, y)) }
    return(c(min(-abs(nums[1]), -abs(nums[2])), max(-abs(nums[1]), -abs(nums[2]))))
  }
  if (str_detect(x0, "\\bAD\\b|\\bCE\\b")) {
    if (length(nums) == 1) { y <- abs(nums[1]); return(c(y, y)) }
    return(c(min(abs(nums[1]), abs(nums[2])), max(abs(nums[1]), abs(nums[2]))))
  }
  if (length(nums) == 1) return(c(nums[1], nums[1]))
  c(min(nums[1], nums[2]), max(nums[1], nums[2]))
}
estimate_matrix <- t(vapply(sample_dates$age_estimate, parse_age_estimate, numeric(2)))
sample_dates <- sample_dates |>
  mutate(est_lo = estimate_matrix[, 1], est_hi = estimate_matrix[, 2])
classify_range <- function(lower, upper,
                           threshold = BIN_OVERLAP_THRESHOLD_YEARS,
                           span_max  = MAX_SPAN_YEARS,
                           max_bins  = MAX_BINS_PER_SAMPLE) {
  b <- chron_bins_500
  out <- function(status, bins = character(0), detail = "")
    list(status = status, bins = bins, detail = detail)
  if (is.na(lower) && is.na(upper)) return(out("no_date", detail = "no usable date"))
  if (is.na(lower)) lower <- upper
  if (is.na(upper)) upper <- lower
  if (lower > upper) { tmp <- lower; lower <- upper; upper <- tmp }
  if (lower == upper) {
    idx <- which(lower >= b$bin_lower & lower <= b$bin_upper)
    if (length(idx) == 0) return(out("out_of_window", detail = "point date outside all bins"))
    cc <- b$chron_cat[idx[1]]
    if (!cc %in% TARGET_BINS) return(out("out_of_window", detail = paste("point in", cc)))
    return(out("kept", cc, "point"))
  }
  overlaps <- pmax(0, pmin(upper, b$bin_upper) - pmax(lower, b$bin_lower))
  if (sum(overlaps) == 0) return(out("out_of_window", detail = "range outside all bins"))
  occupies <- overlaps > threshold
  if (sum(occupies) > max_bins)
    return(out("imprecise",
               detail = sprintf("overlaps %d bins (>%dyr each)", sum(occupies), threshold)))
  if (sum(occupies) > 1 && (upper - lower) >= span_max)
    return(out("imprecise",
               detail = sprintf("span %d yr >= %d cap", upper - lower, span_max)))
  keep <- which(occupies)
  if (length(keep) == 0) keep <- which.max(overlaps)
  keep <- keep[order(overlaps[keep], decreasing = TRUE)][seq_len(min(max_bins, length(keep)))]
  cc <- b$chron_cat[keep]
  cc_in_window <- cc[cc %in% TARGET_BINS]
  if (length(cc_in_window) == 0)
    return(out("out_of_window", detail = paste("bins:", paste(cc, collapse = ","))))
  cc_in_window <- cc_in_window[order(BIN_IDX[cc_in_window])]
  out("kept", cc_in_window, "range")
}
adjacent_ok <- function(bins) {
  idx <- sort(unique(BIN_IDX[bins]))
  idx <- idx[!is.na(idx)]
  length(idx) >= 1 && length(idx) <= 2 && (length(idx) == 1 || diff(range(idx)) == 1)
}

# Assign samples to chronological bins

message("8/9  Assigning samples to chronological bins ...")
sample_dates <- sample_dates |>
  rowwise() |>
  mutate(
    direct = list(if (has_c14) classify_range(c14_lo_eff, c14_hi_eff)
                  else list(status = "no_c14", bins = character(0), detail = "no C14")),
    estim  = list(classify_range(est_lo, est_hi))
  ) |>
  ungroup() |>
  mutate(direct_status = map_chr(direct, "status"),
         estim_status  = map_chr(estim,  "status"))
sample_dates <- sample_dates |>
  mutate(own_precise_bins = pmap(
    list(has_c14, direct_status, map(direct, "bins"), estim_status, map(estim, "bins")),
    function(hc, d_status, d_bins, e_status, e_bins) {
      if (isTRUE(hc) && d_status == "kept") d_bins
      else if (e_status == "kept")          e_bins
      else                                  character(0)
    }))
site_span <- sample_dates |>
  group_by(site_name) |>
  summarise(
    site_bins = list(sort(unique(unlist(own_precise_bins)))),
    site_has_date_outside_window = any(has_c14 & direct_status == "out_of_window"),
    .groups = "drop"
  ) |>
  mutate(
    site_n_bins     = map_int(site_bins, length),
    site_single_bin = site_n_bins == 1 & !site_has_date_outside_window,
    site_single_cc  = map2_chr(site_bins, site_single_bin,
                               function(b, single) if (isTRUE(single)) b[[1]] else NA_character_)
  )
sample_dates <- sample_dates |>
  left_join(site_span |> select(site_name, site_single_bin, site_single_cc),
            by = "site_name")
donors <- sample_dates |>
  filter(has_c14, direct_status == "kept", !is.na(stratum)) |>
  mutate(direct_bins = map(direct, "bins")) |>
  group_by(site_name, stratum) |>
  summarise(
    donor_bins    = list(sort(unique(unlist(direct_bins)))),
    donor_lo      = min(c14_lo_eff, na.rm = TRUE),
    donor_hi      = max(c14_hi_eff, na.rm = TRUE),
    donor_samples = paste(sort(unique(sample_number)), collapse = ","),
    .groups = "drop"
  )
sample_dates <- sample_dates |>
  left_join(donors, by = c("site_name", "stratum"))
sample_dates <- sample_dates |>
  mutate(
    .phase_midpoint = case_when(
      has_c14 & !is.na(c14_lo_eff) ~ (c14_lo_eff + c14_hi_eff) / 2,
      !is.na(donor_lo)             ~ (donor_lo + donor_hi) / 2,
      TRUE                         ~ NA_real_
    ),
    .phase_source = case_when(
      has_c14 & !is.na(c14_lo_eff) ~ "own_C14",
      !is.na(donor_lo)             ~ "stratum_C14",
      TRUE                         ~ NA_character_
    ),
    phase_basis = case_when(
      is.na(described_chronology) |
        described_chronology != "Iron Age" ~ "database_label",
      is.na(.phase_midpoint)               ~ "unresolved",
      TRUE                                 ~ .phase_source
    ),
    phase_final = case_when(
      phase_basis == "database_label"      ~ described_chronology,
      phase_basis == "unresolved"          ~ "Iron Age",
      .phase_midpoint < EIA_LIA_BOUNDARY_CE ~ "Early Iron Age",
      TRUE                                 ~ "Late Iron Age"
    )
  ) |>
  select(-.phase_midpoint, -.phase_source)
decide_bin <- function(has_c14, direct, estim, est_lo, est_hi,
                       donor_bins, donor_lo, donor_hi, donor_samples,
                       site_single_bin, site_single_cc, described_chronology) {
  excluded <- function(reason)
    tibble(tier = "excluded", confidence = NA_character_, basis = "none",
           bins = list(character(0)), donor = NA_character_, reason = reason)
  if (has_c14) {
    if (direct$status == "kept") {
      n <- length(direct$bins)
      return(tibble(tier = if (n == 1) "T1_direct_single" else "T2_direct_double",
                    confidence = "high", basis = "C14",
                    bins = list(direct$bins), donor = NA_character_,
                    reason = "direct C14 date"))
    }
    if (direct$status == "out_of_window")
      return(excluded(paste0("directly dated but outside 500 BCE-1500 CE (",
                             direct$detail, ")")))
  }
  if (!isTRUE(has_c14) && estim$status == "kept" && length(estim$bins) == 1)
    return(tibble(tier = "T1_estimate_single", confidence = "medium",
                  basis = "age_estimate", bins = list(estim$bins),
                  donor = NA_character_, reason = "single-bin age estimate"))
  if (!is.na(described_chronology) && described_chronology %in% POST_WINDOW_PHASES)
    return(excluded(paste0("recorded phase '", described_chronology,
                           "' lies outside 500 BCE-1500 CE")))
  if (isTRUE(site_single_bin) && !is.na(site_single_cc))
    return(tibble(tier = "T1b_site_single_bin", confidence = "medium",
                  basis = "site_single_bin", bins = list(site_single_cc),
                  donor = NA_character_, reason = "whole site confined to one bin"))
  if (!is.null(donor_bins) && length(donor_bins) >= 1 && all(!is.na(donor_bins))) {
    if (adjacent_ok(donor_bins)) {
      consistent <- is.na(est_lo) || is.na(est_hi) ||
                    (min(est_hi, donor_hi) >= max(est_lo, donor_lo))
      if (consistent)
        return(tibble(tier = "T3_stratum_inherit", confidence = "medium",
                      basis = "stratigraphy", bins = list(donor_bins),
                      donor = donor_samples,
                      reason = paste0("inherited from stratum donor(s) ", donor_samples)))
      else
        return(excluded(paste0("stratum donor(s) ", donor_samples,
                               " inconsistent with own estimate")))
    }
  }
  if (estim$status == "kept" && length(estim$bins) %in% 1:2)
    return(tibble(tier = "T4_estimate_range", confidence = "low",
                  basis = "age_estimate", bins = list(estim$bins),
                  donor = NA_character_, reason = "age estimate (<=2 bins)"))
  excluded(case_when(
    has_c14 && direct$status == "imprecise" ~
      paste0("C14 date too imprecise (", direct$detail, ")"),
    estim$status == "imprecise" ~
      paste0("age estimate too imprecise (", estim$detail, ")"),
    estim$status == "out_of_window" ~
      paste0("age estimate outside window (", estim$detail, ")"),
    TRUE ~ "no usable date"
  ))
}
sample_bins <- sample_dates |>
  mutate(decision = pmap(list(has_c14, direct, estim, est_lo, est_hi,
                              donor_bins, donor_lo, donor_hi, donor_samples,
                              site_single_bin, site_single_cc,
                              described_chronology), decide_bin)) |>
  unnest(decision)
sample_bins <- sample_bins |>
  rowwise() |>
  mutate(
    placement_outside_window = {
      lo <- if (has_c14) c14_lo_eff else est_lo
      hi <- if (has_c14) c14_hi_eff else est_hi
      if (tier != "excluded" || is.na(lo)) NA_character_
      else {
        mid <- (lo + hi) / 2
        idx <- which(mid >= chron_bins_500$bin_lower & mid <= chron_bins_500$bin_upper)
        cc  <- if (length(idx)) chron_bins_500$chron_cat[idx[1]] else NA_character_
        if (!is.na(cc) && !cc %in% TARGET_BINS) cc else NA_character_
      }
    },
    range_used_source = case_when(
      basis == "C14"             ~ "C14",
      basis == "age_estimate"    ~ "age_estimate",
      basis == "stratigraphy"    ~ "donor_stratum",
      basis == "site_single_bin" ~ "site_single_bin",
      TRUE                       ~ "none"
    ),
    range_used_lo = case_when(
      range_used_source == "C14"           ~ c14_lo_eff,
      range_used_source == "age_estimate"  ~ est_lo,
      range_used_source == "donor_stratum" ~ donor_lo,
      TRUE                                 ~ NA_real_
    ),
    range_used_hi = case_when(
      range_used_source == "C14"           ~ c14_hi_eff,
      range_used_source == "age_estimate"  ~ est_hi,
      range_used_source == "donor_stratum" ~ donor_hi,
      TRUE                                 ~ NA_real_
    )
  ) |>
  ungroup()
sample_bins_long <- sample_bins |>
  mutate(n_bins = map_int(bins, length),
         bins   = map(bins, function(x) if (length(x)) x else NA_character_)) |>
  unnest(bins) |>
  rename(chron_cat = bins) |>
  mutate(
    bin_label = if_else(is.na(chron_cat), NA_character_, BIN_LABELS[chron_cat]),
    bin_index = match(chron_cat, TARGET_BINS),
    status    = if_else(tier == "excluded", "excluded", "binned")
  ) |>
  select(site_name, sample_number, stratum,
         status, tier, confidence, basis, chron_cat, bin_label, bin_index,
         n_bins, donor, reason, direct_status, estim_status,
         site_single_bin, site_single_cc,
         phase_final, phase_basis, placement_outside_window,
         c14_lo, c14_hi, c14_lo_eff, c14_hi_eff, c14_single_bound,
         c14_n_determinations, est_lo, est_hi, donor_lo, donor_hi,
         range_used_lo, range_used_hi, range_used_source,
         any_of(c("lower_cal_a", "upper_cal_a", "lower_cal_b", "upper_cal_b")))
dated_sites <- sample_bins_long |>
  filter(status == "binned", chron_cat %in% TARGET_BINS) |>
  group_by(site_name) |>
  summarise(
    site_n_dated_samples = n_distinct(sample_number),
    site_bins_in_range   = paste(sort(unique(chron_cat)), collapse = "; "),
    .groups = "drop"
  ) |>
  mutate(in_dated_sites_cohort = TRUE)
message(sprintf("     %d sites have at least one sample dated to 500 BCE - 1500 CE",
                nrow(dated_sites)))

# Assemble, validate and save the analysis dataset

message("9/9  Assembling the final table ...")
plant_records_full <- plant_records |>
  left_join(taxa, by = "species_archaeobotany")
samples_without_plants <- setdiff(samples$sample_number,
                                  unique(plant_records_full$sample_number))
blank_rows <- tibble(sample_number = samples_without_plants,
                     species_archaeobotany = NA_character_,
                     presence_absence = NA_integer_)
analysis_dataset <- sample_bins_long |>
  left_join(samples |> select(-site_name), by = "sample_number") |>
  left_join(dated_sites, by = "site_name") |>
  left_join(bind_rows(plant_records_full, blank_rows),
            by = "sample_number", relationship = "many-to-many") |>
  mutate(
    in_dated_sites_cohort = coalesce(in_dated_sites_cohort, FALSE),
    in_chron_cohort       = status == "binned",
    in_fig3_cohort = status == "binned" &
                     !is.na(accepted_name_short) &
                     !is_post_1500_taxon &
                     phase_final %in% c("Later Stone Age", "Early Iron Age",
                                        "Middle Iron Age", "Late Iron Age",
                                        "Iron Age")
  ) |>
  select(
    site_name, site_name_clean, sample_number,
    country, latitude, longitude,
    kg_zone, kg_zone_recorded,
    described_chronology_raw, described_chronology,
    stratigraphy, stratum, age_estimate,
    any_of(c("lower_cal_a", "upper_cal_a", "lower_cal_b", "upper_cal_b")),
    c14_n_determinations, c14_lo, c14_hi, c14_lo_eff, c14_hi_eff,
    c14_single_bound, est_lo, est_hi, donor_lo, donor_hi,
    range_used_lo, range_used_hi, range_used_source,
    status, tier, confidence, basis, chron_cat, bin_label, bin_index,
    n_bins, donor, reason, direct_status, estim_status,
    site_single_bin, site_single_cc, placement_outside_window,
    phase_final, phase_basis,
    archbot_material,
    presence_absence,
    species_archaeobotany, accepted_name_short, family,
    plant_group, plant_group_filled, top_common_name,
    is_figure3_taxon, is_post_1500_taxon, is_genus_exclude,
    is_indet_label, taxon_merged_figs45, taxon_plot_id_fig3, taxon_plot_name_fig3,
    ecocrop_id, ecocrop_url, ecocrop_kg_codes,
    ecocrop_alt_min_m, ecocrop_alt_max_m,
    in_dated_sites_cohort, site_n_dated_samples, site_bins_in_range,
    in_chron_cohort, in_fig3_cohort
  ) |>
  arrange(country, site_name, sample_number, bin_index, species_archaeobotany)
n_samples  <- n_distinct(analysis_dataset$sample_number)
n_binned   <- n_distinct(analysis_dataset$sample_number[analysis_dataset$status == "binned"])
n_sites    <- n_distinct(analysis_dataset$site_name)
n_cohort   <- n_distinct(analysis_dataset$site_name[analysis_dataset$in_dated_sites_cohort])
n_records  <- sum(!is.na(analysis_dataset$species_archaeobotany))
stopifnot(
  !any(analysis_dataset$archbot_material == "Charcoal", na.rm = TRUE),
  n_samples == nrow(samples),
  !any(duplicated(analysis_dataset[, c("sample_number", "chron_cat",
                                       "species_archaeobotany")])),
  all(analysis_dataset$chron_cat[analysis_dataset$status == "binned"] %in% TARGET_BINS),
  !any(is.na(analysis_dataset$latitude)),
  !any(is.na(analysis_dataset$kg_zone))
)
expected <- list(samples = 131, binned_samples = 107, cohort_sites = 35, rows = 661)
observed <- list(samples = n_samples, binned_samples = n_binned,
                 cohort_sites = n_cohort, rows = nrow(analysis_dataset))
for (nm in names(expected)) {
  if (observed[[nm]] != expected[[nm]]) {
    warning(sprintf(
      "'%s' is %d but %d was expected for this database version. ",
      nm, observed[[nm]], expected[[nm]]),
      "If the database has been updated this is fine - but check the figures.")
  }
}
write_csv(analysis_dataset, FILE_OUTPUT)
line <- function() cat(strrep("=", 74), "\n")
cat("\n"); line()
cat("PAPER A ANALYSIS DATASET\n"); line()
cat(sprintf("Written to            : %s\n", FILE_OUTPUT))
cat(sprintf("Rows x columns        : %d x %d\n",
            nrow(analysis_dataset), ncol(analysis_dataset)))
cat(sprintf("Samples               : %d  (charcoal excluded upstream)\n",
            n_samples))
cat(sprintf("  ... dated to a bin  : %d\n", n_binned))
cat(sprintf("Sites                 : %d  (%d in the dated cohort)\n", n_sites, n_cohort))
cat(sprintf("Plant records         : %d\n", n_records))
excluded <- read_csv(FILE_EXCLUDED, show_col_types = FALSE)
cat(sprintf(
  "\nSelection             : %d sites analysed + %d screened out = %d in the\n",
  n_cohort, nrow(excluded), n_cohort + nrow(excluded)))
cat(sprintf(
  "                        source database (%d samples here, %d screened out).\n",
  n_samples, sum(excluded$n_samples)))
cat("                        Reasons per excluded site: ", FILE_EXCLUDED, "\n", sep = "")
cat("\nSamples per chronological bin (a sample may count in two):\n")
analysis_dataset |>
  filter(status == "binned") |>
  distinct(sample_number, bin_label, bin_index) |>
  count(bin_index, bin_label, name = "samples") |>
  arrange(bin_index) |>
  select(bin_label, samples) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\nSites per chronological bin:\n")
analysis_dataset |>
  filter(status == "binned") |>
  distinct(site_name, bin_label, bin_index) |>
  count(bin_index, bin_label, name = "sites") |>
  arrange(bin_index) |>
  select(bin_label, sites) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\nHow each dated sample got its date:\n")
analysis_dataset |>
  distinct(sample_number, tier, confidence) |>
  count(tier, confidence, name = "samples") |>
  arrange(desc(samples)) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\nChronological phase (note: Middle Iron Age is kept as recorded):\n")
analysis_dataset |>
  distinct(sample_number, phase_final) |>
  count(phase_final, name = "samples") |>
  arrange(desc(samples)) |>
  as.data.frame() |>
  print(row.names = FALSE)
line()
line()
check <- suppressWarnings(read_csv(FILE_OUTPUT, guess_max = Inf, show_col_types = FALSE))
stopifnot(nrow(check) == nrow(analysis_dataset),
          identical(names(check), names(analysis_dataset)),
          is.character(check$donor))
message("Output file verified: reads back with all columns intact.")
