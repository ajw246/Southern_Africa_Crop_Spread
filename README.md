# The spread of crops into southern Africa was a mosaic

Code, data and figure files accompanying the paper *“The spread of crops into southern Africa was a mosaic”*. The study synthesises published archaeobotanical evidence from 35 sites dated 500 BCE–1500 CE and examines the timing and distribution of introduced crops in relation to climate and crop-growing conditions.

**Authors:** Alice J. Williams, Mica B. Jones, Russell Kapumha, Michael Charles, Amy Bogaard, Shadreck Chirikure.

## Abstract

The spread of domesticated crops remains understudied across large parts of the world. In Southern Africa, crop cultivation began relatively recently and coexisted with foraging over the last two thousand years, but direct archaeobotanical investigation has been limited. We synthesized published plant remains from 35 sites dated 500 BCE -1500 CE distributed across southern Africa and compared occurrences of introduced crops with modern climate zones and crop-suitability ranges. Introduced cereals appear in the first half of the first millennium CE, and legume crops in its second half. Introduced crops were recorded only at Iron Age sites and not at Later Stone Age sites, despite periods of contemporaneous occupation. Crop occurrences were geographically dispersed, extending beyond conducive zones, indicating that factors other than climatic suitability played a role in their sequence and distribution. Our results exclude the arrival of a uniform crop package and instead demonstrate a mosaic of temporally staggered and locally variable crop dispersals, alongside established foraging lifeways. 

## Repository contents

- `data/` — the compiled archaeobotanical records, calibrated radiocarbon dates, analysis table, and figure-specific tables.
- `figures/` — PNG and PDF versions of Figures 2–5 and Supplementary Figures S1–S3.
- `00_fetch_koppen_geiger.R` — downloads and verifies the Köppen–Geiger grids used in the analysis.
- `01_data_wrangling.R` and `02_figure_datasets.R` — assign samples to chronological bins and prepare the tables required to generate the figures.
- `Figure_*.R` — scripts to regenerate the main and supplementary figures.

## Reproducing the analysis

Open `CropSpread_FigureCode.Rproj` in RStudio and run the scripts from the repository root. To rebuild the analysis tables, run `00_fetch_koppen_geiger.R`, followed by `01_data_wrangling.R` and `02_figure_datasets.R`, then run the relevant figure scripts. The scripts require R 4.5.2 and the packages listed in their opening sections. Figure scripts may download third-party map layers on first use. The Köppen–Geiger grids are not included in the repository and downloaded in `00_fetch_koppen_geiger.R`.

The analysis uses four 500-year bins from 500 BCE to 1500 CE. Chronological assignment and uncertainty rules are described in Supplementary Materials, SM 2, and recorded in the analysis table. Climate comparisons use modern Köppen–Geiger classifications and FAO ECOCROP suitability ranges.

## Data, citation and reuse

Data provenance, source versions, retrieval details and third-party licences are documented in [`DATA_SOURCES.md`](DATA_SOURCES.md).

Please cite:
- Alice J. Williams, Mica B. Jones, Russell Kapumha, Michael Charles, Amy Bogaard, Shadreck Chirikure. *The spread of crops into southern Africa was a mosaic*. Repository citation. 
- Paper URL and DOI to be added in due course.

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Contact Details
[Alice.Williams@arch.ox.ac.uk](mailto:Alice.Williams@arch.ox.ac.uk)
