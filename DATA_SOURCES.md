# Data sources and provenance

This file records the sources, versions, uses and licence conditions for data used in the figures. Third-party datasets are fetched from their publishers at run time and are not committed to the repository, except where stated.

## Environmental data and suitability records

| Source | Use and version | Access, licence and attribution |
|---|---|---|
| **Köppen–Geiger maps** — Beck et al. (2023), *Scientific Data* 10, 724. [Article](https://doi.org/10.1038/s41597-023-02549-6); [data](https://doi.org/10.6084/m9.figshare.21789074). | Global 30-arc-second grids for 1901–1930, 1931–1960, 1961–1990 and 1991–2020. Used for site climate classifications and Figures 5, S1 and S3. Figshare article version 2 (corrected release, 14 January 2026). | CC BY 4.0. Retrieved 5 June 2026. `00_fetch_koppen_geiger.R` downloads the archive and verifies its MD5 (`7fc2f5a15d4f5fe0ce59c9a9b502aa09`) before extracting the four periods and legend. The grids are not redistributed. Cite Beck et al. and the data DOI when using them. |
| **ETOPO 2022** — NOAA National Centers for Environmental Information. [Product](https://doi.org/10.25921/fd45-gt74). | Bedrock relief: 30 arc-seconds for Figure 2 and 60 arc-seconds for Figure S3. | Public domain. Retrieved on first run through `marmap`; the requested map areas are extracted by the service, so outputs are not pinned to a fixed file. Cached locally and not redistributed. |
| **WorldClim 2.1** — Fick & Hijmans (2017), *International Journal of Climatology* 37, 4302–4315. [Article](https://doi.org/10.1002/joc.5086); [data](https://www.worldclim.org/data/worldclim21.html). | Elevation and monthly precipitation at 5 arc-minutes. Elevation is used in Figure S2; precipitation underlies the rainfall-zone boundaries in Figure 2. | Academic and other non-commercial use; redistribution and commercial use require prior permission. Elevation retrieved 29 August 2026; precipitation retrieved 14 September 2026. Downloaded on first run, cached locally and not redistributed. |
| **Natural Earth.** [Source](https://www.naturalearthdata.com/). | Rivers and lakes at 1:10 million; country outlines at 1:50 million, used in Figures 2 and S3. River and lake releases are not pinned; bundled outlines use `rnaturalearthdata` 1.0.0 and `rnaturalearth` 1.2.0. | Public domain. Retrieved on first run and cached locally; not redistributed. |
| **Rainfall seasonality definition** — Chase & Meadows (2007), *Earth-Science Reviews* 84, 103–138. [Article](https://doi.org/10.1016/j.earscirev.2007.06.002). | The study's definition is applied to WorldClim monthly precipitation: summer rainfall exceeds 66% of annual precipitation from October to March; winter rainfall exceeds 66% from April to September; the intervening belt is year-round rainfall. Figure 2 implements these thresholds. | This is a definition, not a source dataset; no map is copied or redistributed. The source paper does not identify the precipitation dataset behind its own map, so the boundaries here are a new application of its definition to WorldClim 2.1. |
| **FAO ECOCROP.** [Live database](https://gaez.fao.org/pages/ecocrop-find-plant). | Climate-zone and altitude tolerances for the crops shown in Figures 4, 5, and S2. Ten crop records were retrieved 10 June 2026 and rechecked 29 August 2026. The database has no fixed version. | FAO terms allow attributed, non-commercial use; see [FAO terms](https://www.fao.org/contact-us/terms/en/). Values are transcribed into lookup tables in `01_data_wrangling.R`; no database files are included. Records used: sorghum 48747; pearl millet 8418; finger millet 5657; cowpea 2153; Bambara groundnut 10830; baobab 2807; cotton 6493; castor-oil plant 1866; mung bean 2150; broomcorn millet 8280. |

ECOCROP zone groups are expanded in the code to match the Köppen–Geiger sub-zones used in the figures. Blank minimum altitude values are treated as 0 m. Altitude ranges are transcribed only for the five crops shown in Figure S2; other crops carry no altitude value in the analysis tables. Suitability denotes whether a site falls within a database range, not whether cultivation was possible there.

## Archaeobotanical and radiocarbon data

We compiled the archaeobotanical records from published excavation reports and theses. The paper and Supplementary Materials provide the relevant references and dating methods. 

| File | Content and provenance |
|---|---|
| `data/PaperA_source_data.csv` | Published plant records and associated site and sample information, standardised for analysis (574 rows). |
| `data/PaperA_radiocarbon_calibrated.csv` | Calibrated radiocarbon determinations (47 rows); calibration was undertaken in OxCal. |
| `data/PaperA_excluded_sites.csv` | The 37 screened-out sites and reasons for exclusion. |
| `data/PaperA_analysis_dataset.csv` | Analysis table produced by `01_data_wrangling.R`, including chronological assignments and exclusions. |
| `data/PaperA_fig_sites.csv`, `data/PaperA_fig_sample_bins.csv`, `data/PaperA_fig_plant_records.csv` | Figure tables produced by `02_figure_datasets.R`. They are also committed so the figures can be generated without rebuilding the pipeline. |

Site coordinates retain the precision reported in their source publications. The climate-zone assignments for the 35 included sites are derived from the 1991–2020 Köppen–Geiger grid.